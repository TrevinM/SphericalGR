//================================================
//
// File that contains routines for class ScalarField
//
//================================================
#include "Grid.h"
#include "Matter.h"
#include <ctime>
#include "tensors.h"
//
//================================================
// Initialize scalar fields
//================================================
// #define FLAT

void ScalarField::Initialize(state *s, curvature *c, diagnostics *d) {
  //
  // set up initial data 
  //
  // NOTE: we'll assume that time derivative and shift vanish initially
  // 
  indata->Initialize_Scalar_Field(last->sf, last->pi);
  //
  // also: set matter state inter to last so that correct
  // values are used in initial evaluation of ADM sources
  //
  inter->equals(last);
  //
  // a2_center_new = Compute_A2();
  // rho_center_new = rho_center;
  // rho_center_old = rho_center_new;
  // rho_center_previous = rho_center_new;
  sf_center_new = last->sf(0.0,N_g,N_g);
  sf_center_old = sf_center_new;
  sf_center_previous = sf_center_new;
};
//===============================================
// Compute RHS sides for matter equations
//===============================================
void ScalarField::Compute_RHS(state *s, curvature *c, double time) {
  //
  // compute time derivative of vector potential A
  // 
  dot_sf(inter, s, time);
  //
  // compute time derivative of electric field E
  // 
  dot_pi(inter, s, c, time);
};
//
//===============================================
// derivative of scalar field itself
//===============================================
//
void ScalarField::dot_sf(scalar_state *m, state *s, double time) {
  for (int i = N_g; i < N_r - N_g; i++ ) {
    const double rl = grid->r(i);
    for (int j = N_g; j < N_t - N_g; j++ ) {
      const double sinthetal = grid->sintheta(j);
      const double costhetal = grid->costheta(j);
      //      const double cotthetal = costhetal / sinthetal;
      const double rst = rl*sinthetal;
      for (int k = N_g; k < N_p - N_g; k++ ) {
	//
	// first Lie derivative of sf along shift
	//
	const double br = s->shift_r(i,j,k);
	const double bt = s->shift_t(i,j,k) / rl;
	const double bp = s->shift_p(i,j,k) / rst;
   	derivs->sf[i][j][k] = br*m->sf.dr(i,j,k,br) +  
	  bt*m->sf.dtheta(i,j,k,bt) +  
	  bp*m->sf.dphi(i,j,k,bp);
	//
	// and compute RHS...
	//
	derivs->sf[i][j][k] -= s->lapse(i,j,k) * m->pi(i,j,k);
      }
    }
  }
};
//===============================================
// derivative of time derivative pi
//===============================================
void ScalarField::dot_pi(scalar_state *m, state *s, curvature *c, double time) {
  const double a_friedmann = cosmology->a(time);
  //
  // just to make sure...
  //    
  c->det.fill_ghosts();
  //
  // and now...
  //
  for (int i = N_g; i < N_r - N_g; i++ ) {
    const double rl = grid->r(i);
    for (int j = N_g; j < N_t-N_g; j++ ) {
      const double sinthetal = grid->sintheta(j);
      const double costhetal = grid->costheta(j);
      const double rst = rl*sinthetal;
      const double sin2theta = sinthetal * sinthetal;
      //      const double cotthetal = costhetal / sinthetal;
      for (int k = N_g; k < N_p-N_g; k++ ) {  
	//
	// Compute Lie-derivative of Pi along shift and add alpha K pi term
	// 
	const double br = s->shift_r(i,j,k);
	const double bt = s->shift_t(i,j,k) / rl;
	const double bp = s->shift_p(i,j,k) / rst;
	//
	derivs->pi[i][j][k] = br * m->pi.dr(i,j,k,br) + 
	  bt * m->pi.dtheta(i,j,k,bt) + 
	  bp * m->pi.dphi(i,j,k,bp) +
	  s->lapse(i,j,k) * s->K(i,j,k) * m->pi(i,j,k) + s->lapse(i,j,k) * potential->dVdsf(m->sf[i][j][k]);
#ifdef FLAT
	const double e4p = 1.0;
	tensor gup(1.0, 0.0, 0.0, 1.0/(rl*rl), 0.0, 1.0/(rst*rst));
	vect dphi(0.0, 0.0, 0.0);
	vect Gam( -2.0 / rl,
		  - costhetal/(rl * rl * sinthetal),
		  0.0 );
#else
	const double e4p = exp(4.0 * s->phi(i,j,k)) * a_friedmann*a_friedmann;
	// NOTE: gup is *not* rescaled, so in contractions also want *unrescaled* quantities...
	tensor gup(c->gup_rr(i,j,k), c->gup_rt(i,j,k), c->gup_rp(i,j,k),
		   c->gup_tt(i,j,k), c->gup_tp(i,j,k), c->gup_pp(i,j,k));
	vect dphi(s->phi.dr(i,j,k),s->phi.dtheta(i,j,k),s->phi.dphi(i,j,k));
	// conformal \bar \Gamma^i, computed from \Lambda^i and \hat \Gamma^i
	vect Gam( ( s->lam_r(i,j,k) - rl * c->gup_tt(i,j,k) - rl*sin2theta * c->gup_pp(i,j,k) ) / e4p,
		  ( s->lam_t(i,j,k)/rl + 2.0 * c->gup_rt(i,j,k) / rl - sinthetal*costhetal * c->gup_pp(i,j,k) ) / e4p,
		  ( s->lam_p(i,j,k)/rst + 2.0 * c->gup_rp(i,j,k) / rl + 2.0 * costhetal/sinthetal * c->gup_tp(i,j,k) ) / e4p );
#endif
	vect psi(m->sf.dr(i,j,k), m->sf.dtheta(i,j,k), m->sf.dphi(i,j,k));
	vect dlapse(s->lapse.dr(i,j,k),s->lapse.dtheta(i,j,k),s->lapse.dphi(i,j,k));
	tensor dpsi(m->sf.ddr(i,j,k), 
		    m->sf.drdtheta(i,j,k), 
		    m->sf.drdphi(i,j,k),
		    m->sf.drdtheta(i,j,k),
		    m->sf.ddtheta(i,j,k),
		    m->sf.dthetadphi(i,j,k),
		    m->sf.drdphi(i,j,k),
		    m->sf.dthetadphi(i,j,k),
		    m->sf.ddphi(i,j,k));
	// now conformal factor terms to compute physical \Gamma^i 
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++)
	    Gam[a] -= 2.0*gup[a][b] / e4p * dphi[b];
	//
	for (int a = 0; a < 3; a++) { 
	  derivs->pi[i][j][k] += s->lapse(i,j,k) * Gam[a] * psi[a];
	  for (int b = 0; b < 3; b++)
	    derivs->pi[i][j][k] -= gup[a][b] / e4p * ( s->lapse(i,j,k) * dpsi[a][b] + psi[a] * dlapse[b]);
	}
	// finally add Kreiss-Oliger...
	derivs->pi[i][j][k] += eta_KO * m->pi.KO(i,j,k);
      }
    }
  }
};
//
//================================================
// 
// Compute ADM sources 
//
// NOTE: will compute unrescaled sources, as expected in field equations
//================================================
//
void ScalarField::ADM_Sources(state * s, curvature * c)
{
  for (int i = N_g; i < N_r - N_g; i++) {   
    //    const double rl = grid->r(i);
    for (int j = N_g; j < N_t - N_g; j++) {
      //      const double sinthetal = grid->sintheta(j);
      //      const double costhetal = grid->costheta(j);
      //      const double rst = rl * sinthetal; 
      for (int k = N_g; k < N_p - N_g; k++) {     
	const double e4p = exp(4.0 * s->phi(i,j,k));
	tensor gup(c->gup_rr(i,j,k), c->gup_rt(i,j,k), c->gup_rp(i,j,k),
		   c->gup_tt(i,j,k), c->gup_tp(i,j,k), c->gup_pp(i,j,k));
	tensor g_conf = gup.inverse();
	//
	// Compute auxiliary quantity \gamma^{ij} \psi_i \psi_j
	//
	// NOTE: gup is *not* rescaled, so in contractions also want *unrescaled* quantities...
	//	vect psi((*psi_r)(i,j,k), rl * (*psi_t)(i,j,k), rst * (*psi_p)(i,j,k));
	vect psi(inter->sf.dr(i,j,k), inter->sf.dtheta(i,j,k), inter->sf.dphi(i,j,k));
	double psi2 = 0.0;
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++)
	    psi2 += gup[a][b] / e4p * psi[a] * psi[b];
	double pi2 = inter->pi(i,j,k) * inter->pi(i,j,k);
  	//
  	// Now compute ADM density...
  	// 
	adm_sources->rho_ADM[i][j][k] = 0.5 * ( pi2 + psi2 ) + potential->V(inter->sf[i][j][k]);
	//
  	// ... fluxes ...
  	//
	adm_sources->S_r[i][j][k] = inter->pi(i,j,k) * psi[0] ;
  	adm_sources->S_t[i][j][k] = inter->pi(i,j,k) * psi[1] ;
  	adm_sources->S_p[i][j][k] = inter->pi(i,j,k) * psi[2] ;
  	//
  	// ... stresses ...
  	//
  	adm_sources->S_rr[i][j][k] = psi[0] * psi[0] + 0.5 * e4p * g_conf[0][0] * ( pi2 - psi2 ) - e4p * g_conf[0][0] * potential->V(inter->sf[i][j][k]);
	adm_sources->S_rt[i][j][k] = psi[0] * psi[1] + 0.5 * e4p * g_conf[0][1] * ( pi2 - psi2 ) - e4p * g_conf[0][1] * potential->V(inter->sf[i][j][k]);
  	adm_sources->S_rp[i][j][k] = psi[0] * psi[2] + 0.5 * e4p * g_conf[0][2] * ( pi2 - psi2 ) - e4p * g_conf[0][2] * potential->V(inter->sf[i][j][k]);
  	adm_sources->S_tt[i][j][k] = psi[1] * psi[1] + 0.5 * e4p * g_conf[1][1] * ( pi2 - psi2 ) - e4p * g_conf[1][1] * potential->V(inter->sf[i][j][k]);
  	adm_sources->S_tp[i][j][k] = psi[1] * psi[2] + 0.5 * e4p * g_conf[1][2] * ( pi2 - psi2 ) - e4p * g_conf[1][2] * potential->V(inter->sf[i][j][k]);
  	adm_sources->S_pp[i][j][k] = psi[2] * psi[2] + 0.5 * e4p * g_conf[2][2] * ( pi2 - psi2 ) - e4p * g_conf[2][2] * potential->V(inter->sf[i][j][k]);
  	//
  	// ... and trace of stress: 
  	//
  	adm_sources->trace_S[i][j][k] = 1.5 * pi2 - 0.5 * psi2 - 3.0 * potential->V(inter->sf[i][j][k]);
      } 
    }
  }
  //
  adm_sources->rho_ADM.fill_ghosts();
  rho_center = adm_sources->rho_ADM(0.0,N_g,N_g);
  if (rho_center > rho_c_max) rho_c_max = rho_center;
  drhoddr = adm_sources->rho_ADM.ddr(N_g,N_g,N_g);
};
