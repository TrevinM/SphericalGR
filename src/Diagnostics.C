//
//==================================================================================
// Routines that belong to class Diagnostics - compute constraint violations
//==================================================================================
//
#include "tensors.h"
#include "Manager.h"

//=========================================================================================
// Compute violations of Hamiltonian constraint
//
// provide current state, curvature, excise (decides whether Hamiltonian is evaluated in interior
// of a black hole), and the fraction of the grid to which the constraint violation is evaluated (in
// terms of coordinate radius).  NOTE: time needed for cosmology term only
// 
//=========================================================================================
double diagnostics::Hamiltonian(state *s, curvature *curve, auxiliary *aux,
				ADM_Source_Terms * adm,
				bool excise, double frac_grid, double time) {
  //
  // make sure X is updated
  // 
  const double a_friedmann = (*cosmology).a(time);
  aux->Compute_X(s, a_friedmann);
  //
  //   const double small_number = 1.e-12;
#pragma omp parallel for collapse(3)
  for (int i = N_g; i < N_r-N_g; i++) {
    for (int j = N_g; j < N_t-N_g; j++) {
      for (int k = N_g; k < N_p-N_g; k++) {
	const double rl = grid->r(i);
	const double r2 = rl * rl;
	const double stl = grid->sintheta(j);
	const double st2 = stl * stl;
	const double ctl = grid->costheta(j);
	bool inside_horizon = horizonfinder->InsideHorizon(rl,j,k);
	if ( excise && inside_horizon ) Ham[i][j][k] = 0.0;
	else {
	  //
	  // compute trace of R and D^i phi D_i phi, and Laplace operator of phi
	  //	
	  tensor gup(curve->gup_rr[i][j][k],curve->gup_rt[i][j][k],curve->gup_rp[i][j][k],
		     curve->gup_tt[i][j][k],curve->gup_tp[i][j][k],curve->gup_pp[i][j][k]);
	  vect D_phi;
	  tensor DD_phi;
	  vect D_X(aux->X.dr(i,j,k), aux->X.dtheta(i,j,k), aux->X.dphi(i,j,k));
	  tensor DD_X(aux->X.ddr(i,j,k), aux->X.drdtheta(i,j,k), aux->X.drdphi(i,j,k),
		      aux->X.ddtheta(i,j,k), aux->X.dthetadphi(i,j,k), aux->X.ddphi(i,j,k));
	  const double factor1 = a_friedmann * exp(2.0*s->phi(i,j,k))/2.0;
	  const double factor2 = a_friedmann * a_friedmann * exp(4.0*s->phi(i,j,k))/2.0;
	  //
	  D_phi[0] = - factor1 * D_X[0];
	  D_phi[1] = - factor1 * D_X[1];
	  D_phi[2] = - factor1 * D_X[2];
	  //
	  DD_phi[0][0] = - factor1 * DD_X[0][0] + factor2 * D_X[0] * D_X[0]; 
	  DD_phi[0][1] = - factor1 * DD_X[0][1] + factor2 * D_X[0] * D_X[1]; 
	  DD_phi[1][0] = - factor1 * DD_X[1][0] + factor2 * D_X[1] * D_X[0]; 
	  DD_phi[0][2] = - factor1 * DD_X[0][2] + factor2 * D_X[0] * D_X[2]; 
	  DD_phi[2][0] = - factor1 * DD_X[2][0] + factor2 * D_X[2] * D_X[0]; 
	  DD_phi[1][1] = - factor1 * DD_X[1][1] + factor2 * D_X[1] * D_X[1]; 
	  DD_phi[1][2] = - factor1 * DD_X[1][2] + factor2 * D_X[1] * D_X[2]; 
	  DD_phi[2][1] = - factor1 * DD_X[2][1] + factor2 * D_X[2] * D_X[1]; 
	  DD_phi[2][2] = - factor1 * DD_X[2][2] + factor2 * D_X[2] * D_X[2]; 
	  //
	  // store Delta Gamma^i_{jk} in Delta_Gam
	  rank3tens Delta_Gam(curve->DG_r_rr(i,j,k),curve->DG_r_rt(i,j,k),curve->DG_r_rp(i,j,k),
			      curve->DG_r_tt(i,j,k),curve->DG_r_tp(i,j,k),curve->DG_r_pp(i,j,k),
			      curve->DG_t_rr(i,j,k),curve->DG_t_rt(i,j,k),curve->DG_t_rp(i,j,k),
			      curve->DG_t_tt(i,j,k),curve->DG_t_tp(i,j,k),curve->DG_t_pp(i,j,k),
			      curve->DG_p_rr(i,j,k),curve->DG_p_rt(i,j,k),curve->DG_p_rp(i,j,k),
			      curve->DG_p_tt(i,j,k),curve->DG_p_tp(i,j,k),curve->DG_p_pp(i,j,k));
	  //
	  // add Delta connection terms to second derivatives
	  //
	  for (int a = 0; a < 3; a++)
	    for (int b = 0; b < 3; b++) 
	      for (int c = 0; c < 3; c++) 
		DD_phi[a][b]   -= D_phi[c] * Delta_Gam[c][a][b];
	  //
	  // finally add reference connection terms
	  //
	  DD_phi[0][1] -= D_phi[1] / rl;
	  DD_phi[1][0] -= D_phi[1] / rl;
	  DD_phi[0][2] -= D_phi[2] / rl;
	  DD_phi[2][0] -= D_phi[2] / rl;
	  DD_phi[1][1] += D_phi[0] * rl;
	  DD_phi[1][2] -= D_phi[2] * ctl / stl;
	  DD_phi[2][1] -= D_phi[2] * ctl / stl;
	  DD_phi[2][2] += D_phi[0] * rl * st2 + D_phi[1] * stl * ctl;	
	  //
	  // now compute Laplace operator of phi, and DphiDphi
	  //
	  double DphiDphi = 0.0;
	  double Lap_phi = 0.0;
	  //
	  // For better comparison with old code, use derivatives of phi computed
	  // from phi directly, rather than from X, in computation of DphiDphi
	  //
	  D_phi[0] = s->phi.dr(i,j,k);
	  D_phi[1] = s->phi.dtheta(i,j,k);
	  D_phi[2] = s->phi.dphi(i,j,k);
	  for (int a = 0; a < 3; a++)
	    for (int b = 0; b < 3; b++) { 
	      Lap_phi += gup[a][b] * DD_phi[a][b];
	      DphiDphi += gup[a][b] * D_phi[a] * D_phi[b];
	    }
	  //
	  // store extrinsic curvature (lower indices) in A (not rescaled)
	  //
	  tensor A(s->a_rr(i,j,k), rl*s->a_rt(i,j,k), rl*stl*s->a_rp(i,j,k),
		   r2*s->a_tt(i,j,k), r2*stl*s->a_tp(i,j,k), r2*st2*s->a_pp(i,j,k));
	  //
	  // storage for A_{ac} A^c_{~b}
	  //
	  tensor AdotA;
	  //
	  // trace of AdotA
	  //
	  double A2_trace = 0.0;
	  // 
	  // now loop over free indices of A2
	  //
	  for (int a = 0; a < 3; a++) 
	    for (int b = 0; b < 3; b++) { 
	      AdotA[a][b] = 0.0;
	      for (int c = 0; c < 3; c++) 
		for (int d = 0; d < 3; d++)
		  AdotA[a][b] += gup[c][d] * A[a][c] * A[d][b];
	      A2_trace += gup[a][b] * AdotA[a][b];
	    }
	  //
	  // see eq. (16) in B09
	  //
	  double KplusA = s->K(i,j,k) + curve->trace_A(i,j,k) + 2.0 * s->Theta(i,j,k);
	  const double psim4 = exp(-4.0*s->phi(i,j,k)) / (a_friedmann*a_friedmann); 
	  const double K_parts = 2.0/3.0*KplusA*KplusA + curve->trace_A(i,j,k)*curve->trace_A(i,j,k)/3.0 - A2_trace;
	  const double R_parts = psim4 * ( curve->trace_R(i,j,k) - 8.0 * DphiDphi - 8.0 * Lap_phi );
	  // note cosmological term
	  const double rho_parts = - 16.0 * PI * adm->rho_ADM(i,j,k) - 2.0 * (*cosmology).Lambda() ;
	  // const double R_scaling = psim4 * ( abs(trace_R(i,j,k)) + 
	  // 				     8.0 * abs(DphiDphi)  + 
	  // 				     8.0 * abs(D2_phi(i,j,k)) );
	  // const double scaling = abs(K_parts) + R_scaling + abs(rho_parts) + small_number;
	  //	  Ham[i][j][k] = ( K_parts + R_parts + rho_parts ) / scaling;
	  Ham[i][j][k] = ( K_parts + R_parts + rho_parts );
	  // Ham[i][j][k] = 2.0/3.0*KplusA*KplusA + trace_A(i,j,k)*trace_A(i,j,k)/3.0 - A2(i,j,k)
	  //   + exp(-4.0*phi_o(i,j,k)) * ( trace_R(i,j,k) - 8.0 * DphiDphi - 8.0 * D2_phi(i,j,k) )
	  //   - 16.0 * PI * rho_ADM(i,j,k);
	}
	// done...
      }
    }
  }
  return Ham.L2_norm(frac_grid);
};
//================================================
//
// Momentum constraint
//
//================================================
double diagnostics::MomentumConstraint(state * s, curvature * curve, state * t, 
				       ADM_Source_Terms * adm, 
				       double & Mom_r_norm, double & Mom_t_norm, double & Mom_p_norm,
				       bool excise, double frac_grid, double time) {
  const double a_friedmann = (*cosmology).a(time);
  //
  // and make sure det has its ghosts filled
  // 
  curve->det.fill_ghosts();
  //
  // first loop over all gridpoints to compute \bar A^{ij} 
  //
#pragma omp parallel for collapse(3)  
  for (int i = N_g; i < N_r-N_g; i++) {   
    for (int j = N_g; j < N_t-N_g; j++) {
      for (int k = N_g; k < N_p-N_g; k++) {
	const double rl = grid->r(i);
	const double r2 = rl*rl;
	const double sintheta = grid->sintheta(j);
	const double sin2theta = sintheta*sintheta;
  	//
  	// store rescaled \bar A_{ij} in A_down, and (also rescaled!) \bar g^{ij} in g_up 
  	//
  	tensor A_down(s->a_rr(i,j,k), s->a_rt(i,j,k), s->a_rp(i,j,k),
  		      s->a_tt(i,j,k), s->a_tp(i,j,k), s->a_pp(i,j,k));
  	tensor gup(curve->gup_rr(i,j,k),    curve->gup_rt(i,j,k)*rl,          curve->gup_rp(i,j,k)*rl*sintheta, 
  		   curve->gup_tt(i,j,k)*r2, curve->gup_tp(i,j,k)*r2*sintheta, curve->gup_pp(i,j,k)*r2*sin2theta);
  	//
  	// raise both indices of \bar A_{ij} and store in A_up (result is rescaled!)
  	//
  	tensor A_up;
  	for (int a = 0; a < 3; a++)
  	  for (int b = a; b < 3; b++) {
  	    A_up[a][b] = 0.0;
  	    for (int c = 0; c < 3; c++) 
  	      for (int d = 0; d < 3; d++) 
  		A_up[a][b] += gup[a][c] * gup[b][d] * A_down[c][d];
  	  }	
  	//
  	// store result in temporary storage...
  	//
  	t->a_rr[i][j][k] = A_up[0][0];
  	t->a_rt[i][j][k] = A_up[0][1];
  	t->a_rp[i][j][k] = A_up[0][2];
  	t->a_tt[i][j][k] = A_up[1][1];
  	t->a_tp[i][j][k] = A_up[1][2];
  	t->a_pp[i][j][k] = A_up[2][2];
      }
    }
  }
  t->a_rr.fill_ghosts();
  t->a_rt.fill_ghosts();
  t->a_rp.fill_ghosts();
  t->a_tt.fill_ghosts();
  t->a_tp.fill_ghosts();
  t->a_pp.fill_ghosts();
  //
  // Now we can take derivatives and compute momentum constraint...
  //
#pragma omp parallel for collapse(3)  
  for (int i = N_g; i < N_r-2*N_g; i++) {   
    for (int j = N_g; j < N_t-N_g; j++) {
      for (int k = N_g; k < N_p-N_g; k++) {
	const double rl = grid->r(i);
	const double r2 = rl*rl;
	const double sintheta = grid->sintheta(j);
	const double sin2theta = sintheta*sintheta;
	const double costheta = grid->costheta(j);
	const double cottheta = costheta/sintheta;
	bool inside_horizon = horizonfinder->InsideHorizon(rl,j,k);
	if (excise && inside_horizon ) {
	  Mom_r[i][j][k] = 0.0;
	  Mom_t[i][j][k] = 0.0;
	  Mom_p[i][j][k] = 0.0;
	} else {
	  //
	  // compute *rescaled* divergence of \hat D_j \bar A^{ij} (stored in a_ij_t)
	  //
	  vect Mom;
	  Mom[0] = 
	    t->a_rr.dr(i,j,k) + t->a_rt.dtheta(i,j,k)/rl + t->a_rp.dphi(i,j,k)/(rl*sintheta)
	    + ( 2.0 * t->a_rr(i,j,k) - t->a_tt(i,j,k) - t->a_pp(i,j,k) + cottheta * t->a_rt(i,j,k) ) / rl;
	  Mom[1] =
	    t->a_rt.dr(i,j,k) + t->a_tt.dtheta(i,j,k)/rl + t->a_tp.dphi(i,j,k)/(rl*sintheta)
	    + ( 3.0 * t->a_rt(i,j,k) + cottheta * ( t->a_tt(i,j,k) - t->a_pp(i,j,k) ) ) / rl;
	  Mom[2] =
	    t->a_rp.dr(i,j,k) + t->a_tp.dtheta(i,j,k)/rl + t->a_pp.dphi(i,j,k)/(rl*sintheta)
	    + ( 3.0 * t->a_rp(i,j,k) + 2.0 * cottheta * t->a_tp(i,j,k) ) / rl;
	  //	  if (i == 100 && j == 2 && k == 2 && fabs(Mom[1]) > 1.e-13 ) cout << " Mom_t = " << endl;
	  //
	  // compute *rescaled* derivatives...
	  //
	  vect dphi(s->phi.dr(i,j,k), s->phi.dtheta(i,j,k)/rl, s->phi.dphi(i,j,k)/(rl*sintheta));
	  vect dK(  ( s->K.dr(i,j,k)     + 2.0 * s->Theta.dr(i,j,k) ),
		    ( s->K.dtheta(i,j,k) + 2.0 * s->Theta.dtheta(i,j,k) ) / rl,
		    ( s->K.dphi(i,j,k)   + 2.0 * s->Theta.dphi(i,j,k) ) /(rl*sintheta));
	  vect ddet(curve->det.dr(i,j,k),   curve->det.dtheta(i,j,k)/rl,   curve->det.dphi(i,j,k)/(rl*sintheta));
	  //
	  // store rescaled vectors... (note: S_i not rescaled in rest of code, so will rescale here)
	  //
	  vect S_low(adm->S_r(i,j,k),   adm->S_t(i,j,k)/rl, adm->S_p(i,j,k)/(rl*sintheta));
	  //
	  // ... tensors... (again, want rescaled versions...)
	  //
	  tensor A_up(t->a_rr(i,j,k), t->a_rt(i,j,k), t->a_rp(i,j,k),
		      t->a_tt(i,j,k), t->a_tp(i,j,k), t->a_pp(i,j,k));
	  tensor gup(curve->gup_rr(i,j,k),    curve->gup_rt(i,j,k)*rl,          curve->gup_rp(i,j,k)*rl*sintheta, 
		     curve->gup_tt(i,j,k)*r2, curve->gup_tp(i,j,k)*r2*sintheta, curve->gup_pp(i,j,k)*r2*sin2theta);
	  //
	  // raise indices of S (with physical metric)
	  //
	  vect S_up;
	  const double psim4 = exp(-4.0*s->phi(i,j,k)) / (a_friedmann*a_friedmann); 
	  for (int a = 0; a < 3; a++) {
	    S_up[a] = 0.0;
	    for (int b = 0; b < 3; b++) {
	      S_up[a] += psim4 * gup[a][b] * S_low[b];
	    }
	  }
	  //
	  // ... and connection coefficients (need to rescale these...
	  //
	  rank3tens Delta_Gam(curve->DG_r_rr(i,j,k),             curve->DG_r_rt(i,j,k)/rl,             curve->DG_r_rp(i,j,k)/(rl*sintheta),
			      curve->DG_r_tt(i,j,k)/r2,          curve->DG_r_tp(i,j,k)/(r2*sintheta),  curve->DG_r_pp(i,j,k)/(r2*sin2theta),
			      curve->DG_t_rr(i,j,k)*rl,          curve->DG_t_rt(i,j,k),                curve->DG_t_rp(i,j,k)/sintheta,
			      curve->DG_t_tt(i,j,k)/rl,          curve->DG_t_tp(i,j,k)/(rl*sintheta),  curve->DG_t_pp(i,j,k)/(rl*sin2theta),
			      curve->DG_p_rr(i,j,k)*rl*sintheta, curve->DG_p_rt(i,j,k)*sintheta,       curve->DG_p_rp(i,j,k),
			      curve->DG_p_tt(i,j,k)/rl*sintheta, curve->DG_p_tp(i,j,k)/rl,             curve->DG_p_pp(i,j,k)/(rl*sintheta));
	  //
	  // now add other terms:
	  // 
	  for (int a = 0; a < 3; a++) {
	    for (int b = 0; b < 3; b++) {
	      Mom[a] += A_up[a][b] * ( 6.0 * dphi[b] + 0.5 * ddet[b]/curve->det(i,j,k) )
		- 2.0/3.0 * gup[a][b] * dK[b];
	      for (int c = 0; c < 3; c++)
		Mom[a] += A_up[b][c] * Delta_Gam[a][b][c];
	    }
	  }
	  Mom_r[i][j][k] = psim4 * Mom[0] - 8.0*PI*S_up[0]; 
	  Mom_t[i][j][k] = psim4 * Mom[1] - 8.0*PI*S_up[1]; 
	  Mom_p[i][j][k] = psim4 * Mom[2] - 8.0*PI*S_up[2]; 
	}
      }
    }
  }
  //
  // now compute norms
  //
  Mom_r_norm = Mom_r.L2_norm(frac_grid);
  Mom_t_norm = Mom_t.L2_norm(frac_grid);
  Mom_p_norm = Mom_p.L2_norm(frac_grid);
  return sqrt(Mom_r_norm*Mom_r_norm + Mom_t_norm*Mom_t_norm + Mom_p_norm*Mom_p_norm);
}
//================================================
//
// Connection function constraint
//
//================================================
double diagnostics::ConnectionFunctionConstraint(state * s, curvature *c,
						 double & CFC_r_norm, double & CFC_t_norm, 
						 double & CFC_p_norm) {
  for (int i = N_g; i < N_r-N_g; i++) {
    const double rl = grid->r(i);
    for (int j = N_g; j < N_t-N_g; j++) {
      const double stl = grid->sintheta(j);
      for (int k = N_g; k < N_p-N_g; k++) {
	// tensor gup(c->gup_rr[i][j][k],c->gup_rt[i][j][k],c->gup_rp[i][j][k],
	// 	   c->gup_tt[i][j][k],c->gup_tp[i][j][k],c->gup_pp[i][j][k]);
  	// rank3tens Delta_Gam(c->DG_r_rr(i,j,k),c->DG_r_rt(i,j,k),c->DG_r_rp(i,j,k),
	// 		    c->DG_r_tt(i,j,k),c->DG_r_tp(i,j,k),c->DG_r_pp(i,j,k),
	// 		    c->DG_t_rr(i,j,k),c->DG_t_rt(i,j,k),c->DG_t_rp(i,j,k),
	// 		    c->DG_t_tt(i,j,k),c->DG_t_tp(i,j,k),c->DG_t_pp(i,j,k),
	// 		    c->DG_p_rr(i,j,k),c->DG_p_rt(i,j,k),c->DG_p_rp(i,j,k),
	// 		    c->DG_p_tt(i,j,k),c->DG_p_tp(i,j,k),c->DG_p_pp(i,j,k));
	// vect Lam(s->lam_r(i,j,k), s->lam_t(i,j,k)/rl, s->lam_p(i,j,k)/(rl * stl));
	// vect CFC;
	// for (int a = 0; a < 3; a++) {
	//   CFC[a] = Lam[a];
	//   for (int b = 0; b < 3; b++) 
	//     for (int c = 0; c < 3; c++) {
	//       CFC[a] -= gup[b][c] * Delta_Gam[a][b][c];  
	//     }	  
	// }
	// CFC_r[i][j][k] = CFC[0];
	// CFC_t[i][j][k] = CFC[1];
	// CFC_p[i][j][k] = CFC[2];
	CFC_r[i][j][k] = s->lam_r(i,j,k) - c->DG_r(i,j,k);
	CFC_t[i][j][k] = s->lam_t(i,j,k) - c->DG_t(i,j,k);
	CFC_p[i][j][k] = s->lam_p(i,j,k) - c->DG_p(i,j,k);
      }
    }
  }
  CFC_r_norm = CFC_r.L2_norm();
  CFC_t_norm = CFC_t.L2_norm();
  CFC_p_norm = CFC_p.L2_norm();
  return sqrt(CFC_r_norm*CFC_r_norm + CFC_t_norm*CFC_t_norm +CFC_p_norm*CFC_p_norm);
}
//================================================
//
// Search for apparent horizons
//
//================================================
bool diagnostics::FindHorizon(int timestep, double t, double tau_c, 
			      state * s, curvature *c, 
			      ADM_Source_Terms * adm, Fluxes *fluxes,
			      auxiliary * aux, double adm_mass,
			      double mom_guess, bool force) {
  //
  // check whether it's time to search for horizon
  //
  bool search_for_horizon = horizonfinder->TimeToFindHorizon(timestep,t);
  //
  // if so...
  // 
  if ( (search_for_horizon && adm_mass > 0.0) || force) {
    //
    // first compute curvature invariant
    //
    CurvatureInvariant(s,c,adm,t,tau_c,timestep);
    //
    // then compute derivative of phi
    //
    aux->Compute_phi_derivs(s);
    //
    // now call horizon finder
    //
    bool found_horizon = false; 
    //      adm_mass = 0.5;
    double max_mass = 10.0*adm_mass;
    double delta_mass = adm_mass/10.0;
    double mass_guess = 1.e-4 * adm_mass;
    cout << " DIAGNOSTICS: looking for horizon at time t = " << t << endl;
    int it = 0;
    int it_max = 100;
    while ( (mass_guess <= max_mass) && (!found_horizon) && (it < it_max)) {
      it++;
      found_horizon = horizonfinder->FindHorizon(timestep, t,
						 s, fluxes,
						 mass_guess, mom_guess,
						 (*cosmology).a(t));
      //      mass_guess += delta_mass;
      mass_guess *= 1.3;
    }
    return found_horizon;  
  } else {
    return true;
  }
};
//================================================
//
// Compute ADM mass from surface integral (at radial grid point i)
//    (see eq. (3.139) in Baumgarte & Shapiro, 2010)
//
//================================================
double diagnostics::ADM_Mass_Surface(state * s, curvature * c, auxiliary * aux,
				     int i, double time) {
  //  double mass = 0.0;
  // make sure we're in interior of grid
  if (i >= N_r-N_g) i = N_r-N_g-1;
  if (i < N_g) i = N_g;  // should never evaluate integral at r=0 anyway, but...
  // const double rl = grid->r(i);
  const double a_friedmann = (*cosmology).a(time);
  for (int j = N_g; j < N_t-N_g; j++) {
    // const double stl = grid->sintheta(j);
    //    const double delta_theta = grid->delta_theta(j);
    for (int k = N_g; k < N_p-N_g; k++) {
      // const double delta_phi = grid->delta_phi(k);
      tensor gup(c->gup_rr[i][j][k],c->gup_rt[i][j][k],c->gup_rp[i][j][k],
		 c->gup_tt[i][j][k],c->gup_tp[i][j][k],c->gup_pp[i][j][k]);
      rank3tens Delta_Gam(c->DG_r_rr(i,j,k),c->DG_r_rt(i,j,k),c->DG_r_rp(i,j,k),
			  c->DG_r_tt(i,j,k),c->DG_r_tp(i,j,k),c->DG_r_pp(i,j,k),
			  c->DG_t_rr(i,j,k),c->DG_t_rt(i,j,k),c->DG_t_rp(i,j,k),
			  c->DG_t_tt(i,j,k),c->DG_t_tp(i,j,k),c->DG_t_pp(i,j,k),
			  c->DG_p_rr(i,j,k),c->DG_p_rt(i,j,k),c->DG_p_rp(i,j,k),
			  c->DG_p_tt(i,j,k),c->DG_p_tp(i,j,k),c->DG_p_pp(i,j,k));
      vect dphi(s->phi.dr(i,j,k), s->phi.dtheta(i,j,k), s->phi.dphi(i,j,k));
      double Gamma = c->DG_r(i,j,k);
      double Dpsi_up = 0.0;
      for (int a = 0; a < 3; a++) {
	  // CHECK cosmological term
	Dpsi_up += gup[0][a] * exp(s->phi(i,j,k)) * sqrt(a_friedmann) * dphi[a];
	for (int b = 0; b < 3; b++)
	  Gamma -= gup[0][a] * Delta_Gam[b][a][b];
      }
      double term1 = Gamma/(16.0*PI);
      double term2 = Dpsi_up/(2.0*PI);
      // mass += (term1 - term2)*rl*rl*stl*delta_theta*delta_phi;
      aux->integrand[i][j][k] = (term1 - term2);
    }
  }
  return aux->integrand.surface_integral(i);
};
//================================================
//
// Compute Komar Angular Momentum from surface integral (at radial grid point i)
//    (see eq. (3.187) in Baumgarte & Shapiro, 2010)
//
//================================================
double diagnostics::Angular_Momentum(state * s, curvature * c, auxiliary * aux, 
				     int i, double time) {
  // make sure we're in interior of grid
  if (i >= N_r-N_g) i = N_r-N_g-1;
  if (i < N_g) i = N_g;  // should never evaluate integral at r=0 anyway, but...  
  const double a_friedmann = (*cosmology).a(time);
  //
  // compute integrand
  //
  const double rl = grid->r(i);
  const double r2 = rl*rl;
  for (int j = N_g; j < N_t-N_g; j++) {
    const double sintheta = grid->sintheta(j);
    const double sin2theta = sintheta*sintheta;
    for (int k = N_g; k < N_p-N_g; k++) {
      const double psil = exp(s->phi(i,j,k)) * sqrt(a_friedmann);
      const double psi6 = psil*psil*psil*psil*psil*psil;
      //
      // \bar \gamma^{rk} \tilde A_{kphi}
      //
      const double temp = 
	c->gup_rr(i,j,k) * ( rl * sintheta  ) * s->a_rp(i,j,k) +
	c->gup_rt(i,j,k) * ( r2 * sintheta  ) * s->a_tp(i,j,k) +
	c->gup_rp(i,j,k) * ( r2 * sin2theta ) * s->a_pp(i,j,k); 
      //
      // factor for integrand
      // r^2 sin theta part taken care off in gridfunction.surface_integral...
      //
      const double det_surf = (1.0 + s->h_tt(i,j,k)) * (1.0 + s->h_pp(i,j,k)) - s->h_tp(i,j,k)*s->h_tp(i,j,k);
      aux->integrand[i][j][k] = psi6 * det_surf * temp / sqrt(c->gup_rr(i,j,k));
    }
  }
  return aux->integrand.surface_integral(i) / (8.0*PI);  
};
//================================================
//
// Compute Komar Linear Momentum from surface integral (at radial grid point i)
//    (see eq. (3.195) in Baumgarte & Shapiro, 2010)
// NOTE: this implementation assumes momentum in z-direction, i.e. along symmetry axis
//================================================
double diagnostics::Linear_Momentum(state * s, curvature * c, auxiliary * aux, 
				    int i, double time) {
#ifndef AXISYMMETRY
  cout << " DIAGNOSTICS: CAREFUL: Linear momentum evaluated in z-direction only!! " << endl; 
#endif
  //
  // assume that, in spherical symmetry, p = 0...
  //
  if (grid->N_theta_int() == 2)
    return 0.0;
  // make sure we're in interior of grid
  if (i >= N_r-N_g) i = N_r-N_g-1;
  if (i < N_g) i = N_g;  // should never evaluate integral at r=0 anyway, but...  
  const double a_friedmann = (*cosmology).a(time);
  //
  // compute integrand
  //
  const double rl = grid->r(i);
  const double r2 = rl*rl;
  for (int j = N_g; j < N_t-N_g; j++) {
    const double sintheta = grid->sintheta(j);
    const double costheta = grid->costheta(j);
    const double sin2theta = sintheta*sintheta;
    for (int k = N_g; k < N_p-N_g; k++) {
      const double psil = exp(s->phi(i,j,k)) * sqrt(a_friedmann);
      const double psi6 = psil*psil*psil*psil*psil*psil;
      //
      // K^r_{~r} = \bar \gamma^{rk} \tilde A_{kr} + K / 3
      //
      const double Krr = 
	c->gup_rr(i,j,k) *                     s->a_rr(i,j,k) +
	c->gup_rt(i,j,k) * ( rl            ) * s->a_rt(i,j,k) +
	c->gup_rp(i,j,k) * ( rl * sintheta ) * s->a_rp(i,j,k) +
	s->K(i,j,k) / 3.0; 
      //
      //
      // K^r_{~theta} = \bar \gamma^{rk} \tilde A_{ktheta}
      //
      const double Krt = 
	c->gup_rr(i,j,k) * ( rl            ) * s->a_rt(i,j,k) +
	c->gup_rt(i,j,k) * ( r2            ) * s->a_tt(i,j,k) +
	c->gup_rp(i,j,k) * ( r2 * sintheta ) * s->a_tp(i,j,k); 
      //
      //
      // factor for integrand
      // r^2 sin theta part taken care off in gridfunction.surface_integral...
      //
      const double temp = costheta * ( Krr - s->K(i,j,k) ) - sintheta * Krt / rl;
      const double det_surf = (1.0 + s->h_tt(i,j,k)) * (1.0 + s->h_pp(i,j,k)) - s->h_tp(i,j,k)*s->h_tp(i,j,k);
      aux->integrand[i][j][k] = psi6 * det_surf * temp / sqrt(c->gup_rr(i,j,k));      
    }
  }
  return aux->integrand.surface_integral(i) / (8.0*PI);  
};
  
//================================================
//
// Compute proper radius, as well as s_axial
//
//================================================
void diagnostics::Compute_Proper_Radius(state *s, curvature *c, auxiliary *aux, int sigma){
  //
  // find proper distance from center
  //
  for (int j = N_g; j < N_t - N_g; j++)
    for (int k = N_g; k < N_p - N_g; k++) {
      // for each ray integrate outwards: start with proper distance
      // of first interior grid point from origin
      double R = 0.5*exp(2.0*s->phi(0.0,j,k))*
	sqrt(1.0 + s->h_rr(0.0,j,k))*grid->delta_r(N_g);
      R_prop[N_g][j][k] = R;
      for (int i = N_g+1; i < N_r; i++) {
	const double r_mid = 0.5 * ( grid->r(i) + grid->r(i-1) );
	R += exp(2.0*s->phi(r_mid,j,k))*sqrt(1.0 + s->h_rr(r_mid,j,k))
	  *grid->delta_r(i);
	R_prop[i][j][k] = R;
      }
    }
  R_prop.fill_ghosts();
  //
  // compute axial radius (invariant in axisymmetry)
  // NOTICE: s_axial = \xi_phi / sin(theta)
  // Take out sin(theta) so that value is accurate in spherical symmetry
  // (i.e. don't introduce artificial theta dependence, which then 
  // introduces numerical error when interpolating to equator)
  // NO - don't take out sin(theta) for now...
  //
  for (int i = N_g; i < N_r; i++) {
    const double rl = grid->r(i);
    for (int j = N_g; j < N_t - N_g; j++) {
      const double stl = grid->sintheta(j);
      // const double stl = 1.0;
      for (int k = N_g; k < N_p - N_g; k++) {
	s_axial[i][j][k] = exp(2.0 * s->phi(i,j,k)) * rl * stl * sqrt(1.0 + s->h_pp(i,j,k));
      }
    }
  }
  s_axial.fill_ghosts();
  //
  // Finally, compute curvature invariant zeta as defined in Ledvinka and Khirnov, 2021
  //
  zeta_max = 0.0;
  aux->DivShift(s, c);
  aux->Lie_Metric(s->shift_r, s->shift_t, s->shift_p,
                  s->h_rr, s->h_rt, s->h_rp, s->h_tt, s->h_tp, s->h_pp);
  const double twothirds = 2./3.;
  for (int i = N_g; i < N_r - N_g; i++) {
    const double rl = grid->r(i);
    for (int j = N_g; j < N_t - N_g; j++) {
      const double stl = grid->sintheta(j);
      const double ctl = grid->costheta(j);
      const double cotl = ctl/stl;
      for (int k = N_g; k < N_p - N_g; k++) {
	const double lapse_l = s->lapse(i,j,k);
	const double psim4 = exp(-4.0*s->phi(i,j,k)) ;
	const double sl = s_axial(i,j,k);
	const double s_r = s_axial.dr(i,j,k);
	const double s_theta = (2.0*s->phi.dtheta(i,j,k) + cotl +
				0.5*s->h_pp.dtheta(i,j,k)/(1.0 + s->h_pp(i,j,k)))*sl;
	// const double h_pp_t = - 2.0 * lapse_l * s->a_pp(i,j,k) +
        //   aux->Lt_pp(i,j,k) - sigma * twothirds * aux->div_shift(i,j,k) * (1.0 + s->h_pp(i,j,k));
	// const double phi_t = - lapse_l * (s->K(i,j,k) + 2.0 * s->Theta(i,j,k))/6.0 + 
        //   aux->Lie(s->shift_r,s->shift_t,s->shift_p,s->phi,i,j,k) +
	//   sigma * aux->div_shift(i,j,k)/6.0;
	// const double s_t = (2.*phi_t + 0.5*h_pp_t / (1 + s->h_pp(i,j,k))) * sl;
	vect grad_s(s_r, s_theta, 0.0);   // assume axisymmetry!!
	tensor gup(c->gup_rr(i,j,k), c->gup_rt(i,j,k), c->gup_rp(i,j,k),
		   c->gup_tt(i,j,k), c->gup_tp(i,j,k), c->gup_pp(i,j,k));
	vect shift(s->shift_r(i,j,k), s->shift_t(i,j,k)/rl, 0.0);   // assume axisymmetry!!
	// double grad_s2 = (- s_t*s_t + 2.0*shift[0]*s_r*s_t + 2.0*shift[1]*s_theta*s_t)
	//  / (lapse_l * lapse_l);
	const double Kpp = exp(4.0*s->phi(i,j,k))*rl*rl*stl*stl*
	  ( s->a_pp(i,j,k) + (1.0 + s->h_pp(i,j,k))*s->K(i,j,k)/3.0 );
	double grad_s2 = - Kpp*Kpp / (sl*sl); 
	for (int a = 0; a < 3; a++) 
	  for (int b = 0; b < 3; b++) 
	    // grad_s2 += ( psim4 * gup[a][b] - shift[a]*shift[b]/(lapse_l*lapse_l) )
	    //   * grad_s[a]*grad_s[b];
	    grad_s2 += psim4 * gup[a][b] * grad_s[a]*grad_s[b];
	zeta[i][j][k] = (1.0 - grad_s2)/(sl*sl);
	if (fabs(zeta(i,j,k)) > fabs(zeta_max))
	  zeta_max = zeta(i,j,k);
      }
    }
  }
  zeta.fill_ghosts();
}
//====================================================================
//
// Compute electric and magnetic parts of Weyl tensor; curvature invariant
//
//====================================================================
bool diagnostics::CurvatureInvariant(state *s, curvature *c, ADM_Source_Terms *adm,
				     double t, double tau_c, int step) {
  I_Re_max = 0.0;
  I_Im_max = 0.0;
  J_Re_max = 0.0;
  J_Im_max = 0.0;
  const double a_friedmann = (*cosmology).a(t);
  //
  // As warm-up exercise compute *physical* extrinsic curvature and store in tilde values a_ij_t
  // 
  // for (int i = N_g; i < N_r-N_g; i++)    
  //   for (int j = N_g; j < N_theta-N_g; j++)
  //     for (int k = N_g; k < N_phi-N_g; k++) {
  // 	//
  // 	// define conformal factor and fourth power
  // 	//
  // 	const double psil = exp(phi_o(i,j,k));
  // 	const double psi4 = psil*psil*psil*psil * a_friedmann*a_friedmann;
  // 	//
  // 	// and here we go:
  // 	// 
  // 	a_rr_t[i][j][k] = psi4 * ( a_rr_o(i,j,k) + K_o(i,j,k) * (1.0 + h_rr_o(i,j,k)) / 3.0 );
  // 	a_rt_t[i][j][k] = psi4 * ( a_rt_o(i,j,k) + K_o(i,j,k) * (      h_rt_o(i,j,k)) / 3.0 ) * r[i];
  // 	a_rp_t[i][j][k] = psi4 * ( a_rp_o(i,j,k) + K_o(i,j,k) * (      h_rp_o(i,j,k)) / 3.0 ) * r[i] * sintheta[j];
  // 	a_tt_t[i][j][k] = psi4 * ( a_tt_o(i,j,k) + K_o(i,j,k) * (1.0 + h_tt_o(i,j,k)) / 3.0 ) * r2[i];
  // 	a_tp_t[i][j][k] = psi4 * ( a_tp_o(i,j,k) + K_o(i,j,k) * (      h_tp_o(i,j,k)) / 3.0 ) * r2[i] * sintheta[j];
  // 	a_pp_t[i][j][k] = psi4 * ( a_pp_o(i,j,k) + K_o(i,j,k) * (1.0 + h_pp_o(i,j,k)) / 3.0 ) * r2[i] * sin2theta[j];
  //     } 
  // a_rr_t.fill_ghosts();
  // a_rt_t.fill_ghosts();
  // a_rp_t.fill_ghosts();
  // a_tt_t.fill_ghosts();
  // a_tp_t.fill_ghosts();
  // a_pp_t.fill_ghosts();
  //
  //
  //
#pragma omp parallel for collapse(3)  
  for (int i = N_g; i < N_r-N_g; i++) {     
    for (int j = N_g; j < N_t-N_g; j++) {
      for (int k = N_g; k < N_p-N_g; k++) {
	const double rl = grid->r(i);
	const double r2 = rl * rl;
	const double stl = grid->sintheta(j);
	const double st2 = stl * stl;
	const double ctl = grid->costheta(j);
	//
	// will need lots of quantities, let's start with partial
	// derivatives of conformal factor...
	//
	vect D_phi(s->phi.dr(i,j,k),s->phi.dtheta(i,j,k),s->phi.dphi(i,j,k));
	tensor DD_phi(s->phi.ddr(i,j,k), s->phi.drdtheta(i,j,k), s->phi.drdphi(i,j,k),
		      s->phi.ddtheta(i,j,k), s->phi.dthetadphi(i,j,k), s->phi.ddphi(i,j,k));
	//
	// now assemble *physical* spatial Christoffel symbols.  Start with \Delta \Gamma's...
	//
	rank3tens Gam(c->DG_r_rr(i,j,k),c->DG_r_rt(i,j,k),c->DG_r_rp(i,j,k),
		      c->DG_r_tt(i,j,k),c->DG_r_tp(i,j,k),c->DG_r_pp(i,j,k),
		      c->DG_t_rr(i,j,k),c->DG_t_rt(i,j,k),c->DG_t_rp(i,j,k),
		      c->DG_t_tt(i,j,k),c->DG_t_tp(i,j,k),c->DG_t_pp(i,j,k),
		      c->DG_p_rr(i,j,k),c->DG_p_rt(i,j,k),c->DG_p_rp(i,j,k),
		      c->DG_p_tt(i,j,k),c->DG_p_tp(i,j,k),c->DG_p_pp(i,j,k));
	//
	// ... now add \tilde \Gamma's to get \bar \Gamma's
	//
	Gam[0][1][1] -= rl;
	Gam[0][2][2] -= rl * st2;
	Gam[1][2][2] -= stl*ctl;
	Gam[1][0][1] += 1.0/rl;
	Gam[1][1][0] += 1.0/rl;
	Gam[2][0][2] += 1.0/rl;
	Gam[2][2][0] += 1.0/rl;
	Gam[2][1][2] += ctl/stl;
	Gam[2][2][1] += ctl/stl;
	//
	// Now compute covariant second derivatives of conformal factor
	//
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++)
	    for (int c = 0; c < 3; c++) 
	      DD_phi[a][b] -= D_phi[c] * Gam[c][a][b];
	//
	// Store conformally related metric in g and inverse metric in gup
	//
	tensor g(1.0 + s->h_rr(i,j,k), rl*s->h_rt(i,j,k), rl*stl*s->h_rp(i,j,k),
		 r2*(1.0 + s->h_tt(i,j,k)), r2*stl*s->h_tp(i,j,k), r2*st2*(1.0 + s->h_pp(i,j,k)));
	tensor gup(c->gup_rr[i][j][k],c->gup_rt[i][j][k],c->gup_rp[i][j][k],
		   c->gup_tt[i][j][k],c->gup_tp[i][j][k],c->gup_pp[i][j][k]);
	//
	// store \bar R_{ij} into R...
	//
	tensor R(c->R_rr(i,j,k),c->R_rt(i,j,k),c->R_rp(i,j,k),
		 c->R_tt(i,j,k),c->R_tp(i,j,k),c->R_pp(i,j,k));
	//
	// and add conformal factor terms to \bar R_{ij} to get physical (spatial) R_{ij} 
	// (e.g. eq. (10a) in Brown, 2009)
	//
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++) {
	    R[a][b] += - 2.0 * DD_phi[a][b] + 4.0 * D_phi[a] * D_phi[b];
	    for (int c = 0; c < 3; c++)
	      for (int d = 0; d < 3; d++) 
		R[a][b] += -2.0 * g[a][b] * gup[c][d] * ( DD_phi[c][d] + 2.0 * D_phi[c] * D_phi[d] );
	  }
	//
	// define conformal factor and fourth power
	//
	const double psil = exp(s->phi(i,j,k));
	const double psi4 = psil*psil*psil*psil * a_friedmann*a_friedmann;
	//
	// assign physical extrinsic curvature
	//
	tensor K(psi4 * ( s->a_rr(i,j,k) + s->K(i,j,k) * (1.0 + s->h_rr(i,j,k)) / 3.0 ),
		 psi4 * ( s->a_rt(i,j,k) + s->K(i,j,k) * (      s->h_rt(i,j,k)) / 3.0 ) * rl,
		 psi4 * ( s->a_rp(i,j,k) + s->K(i,j,k) * (      s->h_rp(i,j,k)) / 3.0 ) * rl * stl,
		 psi4 * ( s->a_tt(i,j,k) + s->K(i,j,k) * (1.0 + s->h_tt(i,j,k)) / 3.0 ) * r2,
		 psi4 * ( s->a_tp(i,j,k) + s->K(i,j,k) * (      s->h_tp(i,j,k)) / 3.0 ) * r2 * stl,
		 psi4 * ( s->a_pp(i,j,k) + s->K(i,j,k) * (1.0 + s->h_pp(i,j,k)) / 3.0 ) * r2 * st2);
	//	  tensor K(a_rr_t(i,j,k),a_rt_t(i,j,k),a_rp_t(i,j,k),a_tt_t(i,j,k),a_tp_t(i,j,k),a_pp_t(i,j,k));
	//
	// compute physical metric
	//
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++) {
	    g[a][b] *= psi4;
	    gup[a][b] /= psi4;
	  }
	//
	// now compute electric part of Weyl tensor, including matter 
	//   (define stress-tensor - not rescaled, so can use as is...)
	//
	tensor E;
	tensor Stress(adm->S_rr(i,j,k), adm->S_rt(i,j,k), adm->S_rp(i,j,k), 
		      adm->S_tt(i,j,k), adm->S_tp(i,j,k), adm->S_pp(i,j,k));
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++) {
	    E[a][b] = R[a][b] + s->K(i,j,k) * K[a][b]
	      + 4.0*PI*( - Stress[a][b] + g[a][b] * ( adm->trace_S(i,j,k) - 4.0*adm->rho_ADM(i,j,k) ) / 3.0 );
	    for (int c = 0; c < 3; c++)
	      for (int d = 0; d < 3; d++) 
		E[a][b] -= gup[c][d] * K[a][c] * K[d][b];
	  }
	//
	// Now on to the magnetic part.  First compute partial derivative of (physical) K_{ij}.
	// To do so, first compute partial derivatives of \bar \gamma_{ij} and \bar A_{ij} 
	//
	rank3tens D_A(s->a_rr.dr(i,j,k),
		      rl * s->a_rt.dr(i,j,k) + s->a_rt(i,j,k),
		      rl * stl * s->a_rp.dr(i,j,k) + stl * s->a_rp(i,j,k),
		      r2 * s->a_tt.dr(i,j,k) + 2.0 * rl * s->a_tt(i,j,k),
		      r2 * stl * s->a_tp.dr(i,j,k) + 2.0 * rl * stl * s->a_tp(i,j,k),
		      r2 * st2 * s->a_pp.dr(i,j,k) + 2.0 * rl * st2 * s->a_pp(i,j,k),
		      s->a_rr.dtheta(i,j,k),
		      rl * s->a_rt.dtheta(i,j,k),
		      rl * stl * s->a_rp.dtheta(i,j,k) + rl * ctl * s->a_rp(i,j,k),
		      r2 * s->a_tt.dtheta(i,j,k),
		      r2 * stl * s->a_tp.dtheta(i,j,k) + r2 * ctl * s->a_tp(i,j,k),
		      r2 * st2 * s->a_pp.dtheta(i,j,k) + 2.0 * r2 * ctl * stl * s->a_pp(i,j,k),
		      s->a_rr.dphi(i,j,k),
		      rl * s->a_rt.dphi(i,j,k),
		      rl * stl * s->a_rp.dphi(i,j,k),
		      r2 * s->a_tt.dphi(i,j,k),
		      r2 * stl * s->a_tp.dphi(i,j,k),
		      r2 * st2 * s->a_pp.dphi(i,j,k));
	rank3tens D_gamma(s->h_rr.dr(i,j,k),
			  rl * s->h_rt.dr(i,j,k) + s->h_rt(i,j,k),
			  rl * stl * s->h_rp.dr(i,j,k) + stl * s->h_rp(i,j,k),
			  r2 * s->h_tt.dr(i,j,k) + 2.0 * rl * (1.0 + s->h_tt(i,j,k)),
			  r2 * stl * s->h_tp.dr(i,j,k) + 2.0 * rl * stl * s->h_tp(i,j,k),
			  r2 * st2 * s->h_pp.dr(i,j,k) + 2.0 * rl * st2 * (1.0 + s->h_pp(i,j,k)),
			  s->h_rr.dtheta(i,j,k),
			  rl * s->h_rt.dtheta(i,j,k),
			  rl * stl * s->h_rp.dtheta(i,j,k) + rl * ctl * s->h_rp(i,j,k),
			  r2 * s->h_tt.dtheta(i,j,k),
			  r2 * stl * s->h_tp.dtheta(i,j,k) + r2 * ctl * s->h_tp(i,j,k),
			  r2 * st2 * s->h_pp.dtheta(i,j,k) + 2.0 * r2 * ctl * stl * (1.0 + s->h_pp(i,j,k)),
			  s->h_rr.dphi(i,j,k),
			  rl * s->h_rt.dphi(i,j,k),
			  rl * stl * s->h_rp.dphi(i,j,k),
			  r2 * s->h_tt.dphi(i,j,k),
			  r2 * stl * s->h_tp.dphi(i,j,k),
			  r2 * st2 * s->h_pp.dphi(i,j,k));
	vect D_K_trace(s->K.dr(i,j,k),s->K.dtheta(i,j,k),s->K.dphi(i,j,k));
	rank3tens D_K;
	//
	//
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++)
	    for (int c = 0; c < 3; c++) {
	      // do we need derivative of metric here? In the end, D_a \gamma_{bc} = 0...
	      D_K[a][b][c] = 4.0 * D_phi[a] * K[b][c] + psi4 * D_A[a][b][c] + 
		psi4 * s->K(i,j,k) * D_gamma[a][b][c] / 3.0 + g[b][c] * D_K_trace[a] / 3.0;
	    }
	//
	// notation: D_r K_ij = D_K[0][i][j]
	// //
	// D_K[0][0][0] = a_rr_t.dr(i,j,k);
	// D_K[1][0][0] = a_rr_t.dtheta(i,j,k);
	// D_K[2][0][0] = a_rr_t.dphi(i,j,k);
	// //
	// D_K[0][1][0] = D_K[0][0][1] = a_rt_t.dr(i,j,k);
	// D_K[1][1][0] = D_K[1][0][1] = a_rt_t.dtheta(i,j,k);
	// D_K[2][1][0] = D_K[2][0][1] = a_rt_t.dphi(i,j,k);
	// //
	// D_K[0][2][0] = D_K[0][0][2] = a_rp_t.dr(i,j,k);
	// D_K[1][2][0] = D_K[1][0][2] = a_rp_t.dtheta(i,j,k);
	// D_K[2][2][0] = D_K[2][0][2] = a_rp_t.dphi(i,j,k);
	// //
	// D_K[0][1][1] = a_tt_t.dr(i,j,k);
	// D_K[1][1][1] = a_tt_t.dtheta(i,j,k);
	// D_K[2][1][1] = a_tt_t.dphi(i,j,k);
	// //
	// D_K[0][2][1] = D_K[0][1][2] = a_tp_t.dr(i,j,k);
	// D_K[1][2][1] = D_K[1][1][2] = a_tp_t.dtheta(i,j,k);
	// D_K[2][2][1] = D_K[2][1][2] = a_tp_t.dphi(i,j,k);
	// //
	// D_K[0][2][2] = a_pp_t.dr(i,j,k);
	// D_K[1][2][2] = a_pp_t.dtheta(i,j,k);
	// D_K[2][2][2] = a_pp_t.dphi(i,j,k);
	  
	//
	// add conformal factor terms to \bar \Gamma's to get physical \Gamma's
	// 
	tensor delta(1.0,0.0,0.0,1.0,0.0,1.0);
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++)
	    for (int c = 0; c < 3; c++) {
	      Gam[a][b][c] += 2.0 * ( delta[a][b] * D_phi[c] + delta[a][c] * D_phi[b] );
	      for (int d = 0; d < 3; d++)
		Gam[a][b][c] -= 2.0 * g[b][c] * gup[a][d] * D_phi[d]; 
	    }
	//
	// now add connection terms to compute covariant derivative
	//
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++)
	    for (int c = 0; c < 3; c++)
	      for (int d = 0; d < 3; d++)
		D_K[a][b][c] -= K[b][d] * Gam[d][a][c] + K[d][c] * Gam[d][a][b];
	//	
	// NOTE: can't use rank3tens to define Levi-Civita because rank3tens assumes 
	// symmetry on last two indices...
	// So: add terms by hand (is faster anyway)
	//
	double sqrtdet = sqrt(c->det(i,j,k))*r2*stl*psi4*psil*psil*a_friedmann;
	// NOTE: need extra factor of friedmann, since only psi4 is rescaled above...
	//
	// NOTE: No need to add matter term, if in the end we compute symmetric part
	//
	tensor B(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for (int c = 0; c < 3; c++)
	  for (int d = 0; d < 3; d++) {
	    B[0][0] += sqrtdet * ( gup[1][c] * gup[2][d] - gup[2][c] * gup[1][d] ) * D_K[d][c][0];
	    B[1][0] += sqrtdet * ( gup[1][c] * gup[2][d] - gup[2][c] * gup[1][d] ) * D_K[d][c][1];
	    B[0][1] += sqrtdet * ( gup[2][c] * gup[0][d] - gup[0][c] * gup[2][d] ) * D_K[d][c][0];
	    B[2][0] += sqrtdet * ( gup[1][c] * gup[2][d] - gup[2][c] * gup[1][d] ) * D_K[d][c][2];
	    B[0][2] += sqrtdet * ( gup[0][c] * gup[1][d] - gup[1][c] * gup[0][d] ) * D_K[d][c][0];
	    B[1][1] += sqrtdet * ( gup[2][c] * gup[0][d] - gup[0][c] * gup[2][d] ) * D_K[d][c][1];
	    B[2][1] += sqrtdet * ( gup[2][c] * gup[0][d] - gup[0][c] * gup[2][d] ) * D_K[d][c][2];
	    B[1][2] += sqrtdet * ( gup[0][c] * gup[1][d] - gup[1][c] * gup[0][d] ) * D_K[d][c][1];
	    B[2][2] += sqrtdet * ( gup[0][c] * gup[1][d] - gup[1][c] * gup[0][d] ) * D_K[d][c][2];
	  }
	//
	// fill in symmetries
	//
	B[1][0] = 0.5 * ( B[0][1] + B[1][0] );  B[0][1] = B[1][0];
	B[2][0] = 0.5 * ( B[0][2] + B[2][0] );  B[0][2] = B[2][0];
	B[2][1] = 0.5 * ( B[1][2] + B[2][1] );  B[1][2] = B[2][1];
	// 
	// fill gridfunctions
	//
	E_rr[i][j][k] = E[0][0];
	E_rt[i][j][k] = E[0][1];
	E_rp[i][j][k] = E[0][2];
	E_tt[i][j][k] = E[1][1];
	E_tp[i][j][k] = E[1][2];
	E_pp[i][j][k] = E[2][2];
	B_rr[i][j][k] = B[0][0];
	B_rt[i][j][k] = B[0][1];
	B_rp[i][j][k] = B[0][2];
	B_tt[i][j][k] = B[1][1];
	B_tp[i][j][k] = B[1][2];
	B_pp[i][j][k] = B[2][2];
	//
	// finally compute curvature scalars electric and magnetic parts
	//
	tensor E_mixed(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);	
	tensor B_mixed(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++)
	    for (int c = 0; c < 3; c++) {
	      E_mixed[a][b] += gup[a][c] * E[c][b];
	      B_mixed[a][b] += gup[a][c] * B[c][b];
	    }
	
	// to get started, raise first index of 
	double i_re = 0.0;
	double i_im = 0.0;
	double j_re = 0.0;
	double j_im = 0.0;
	for (int a = 0; a < 3; a++)
	  for (int b = 0; b < 3; b++) {
	    i_re += E_mixed[a][b] * E_mixed[b][a] - B_mixed[a][b] * B_mixed[b][a];
	    i_im += E_mixed[a][b] * B_mixed[b][a] + B_mixed[a][b] * E_mixed[b][a]; 
	    for (int c = 0; c < 3; c++) {
	      j_re -= E_mixed[a][b] * E_mixed[b][c] * E_mixed[c][a] 
		- E_mixed[a][b] * B_mixed[b][c] * B_mixed[c][a]  
		- B_mixed[a][b] * E_mixed[b][c] * B_mixed[c][a] 
		- B_mixed[a][b] * B_mixed[b][c] * E_mixed[c][a]; 
	      j_im += B_mixed[a][b] * B_mixed[b][c] * B_mixed[c][a] 
		- E_mixed[a][b] * E_mixed[b][c] * B_mixed[c][a]  
		- E_mixed[a][b] * B_mixed[b][c] * E_mixed[c][a] 
		- B_mixed[a][b] * E_mixed[b][c] * E_mixed[c][a]; 
	    }
	    // for (int d = 0; d < 3; d++) {
	    // 	i_re += gup[a][b] * gup[c][d] * ( E[a][c] * E[b][d] - B[a][c] * B[b][d] );
	    // 	i_im += gup[a][b] * gup[c][d] * ( E[a][c] * B[b][d] + B[a][c] * E[b][d] );
	  }
	i_re /= 2.0;
	i_im /= 2.0;
	j_re /= 6.0;
	j_im /= 6.0;
	I_Re[i][j][k] = i_re;
	I_Im[i][j][k] = i_im;
	J_Re[i][j][k] = j_re;
	J_Im[i][j][k] = j_im;
	//
	// CHECK this hack: disregard innermost gridpoint in computation of maximum...
	//
	if (fabs(i_re) > fabs(I_Re_max) && i > N_g) I_Re_max = i_re;
	if (fabs(i_im) > fabs(I_Im_max) && i > N_g) I_Im_max = i_im;
	if (fabs(j_re) > fabs(J_Re_max) && i > N_g) J_Re_max = j_re;
	if (fabs(j_im) > fabs(J_Im_max) && i > N_g) J_Im_max = j_im;
	psi4_Re[i][j][k] =   0.5*(E_tt(i,j,k) - 2.0*B_tp(i,j,k)/stl - E_pp(i,j,k)/st2) / r2;
	psi4_Im[i][j][k] = - 0.5*(B_tt(i,j,k) + 2.0*E_tp(i,j,k)/stl - B_pp(i,j,k)/st2) / r2;
      }
    }
  } 
  E_rr.fill_ghosts();
  E_rt.fill_ghosts();
  E_rp.fill_ghosts();
  E_tt.fill_ghosts();
  E_tp.fill_ghosts();
  E_pp.fill_ghosts();
  B_rr.fill_ghosts();
  B_rt.fill_ghosts();
  B_rp.fill_ghosts();
  B_tt.fill_ghosts();
  B_tp.fill_ghosts();
  B_pp.fill_ghosts();
  I_Re.fill_ghosts();
  I_Im.fill_ghosts();
  J_Re.fill_ghosts();
  J_Im.fill_ghosts();
  psi4_Re.fill_ghosts();
  psi4_Im.fill_ghosts();
  //
  // finally, update global maximum, if needed...
  //
  bool print = false;
  if (fabs(I_Re_max) >= fabs(I_Re_max_max)) {
    I_Re_max_max = I_Re_max;
    new_max_last_timestep = true;
  } else {
    if (new_max_last_timestep) print = true;
    new_max_last_timestep = false;
  }
  if (fabs(I_Im_max) > fabs(I_Im_max_max)) I_Im_max_max = I_Im_max;
  if (fabs(J_Re_max) > fabs(J_Re_max_max)) J_Re_max_max = J_Re_max;
  if (fabs(J_Im_max) > fabs(J_Im_max_max)) J_Im_max_max = J_Im_max;
  //  cout << " I_Re_max = " << I_Re_max << " I_Re_max_max = " << I_Re_max_max << " new_max_last_timestep = " <<
  //    new_max_last_timestep << " print = " << print << endl;
  if (print) {
    dump_fcts(t, tau_c, step, "_max");
    s->dump_fcts(t, tau_c, step, "_max");
  }
  //
  // compute integral
  //
  I_Re_int = 0.0;
  I_Im_int = 0.0;
  J_Re_int = 0.0;
  J_Im_int = 0.0;
  for (int i = N_g; i < N_r-N_g; i++) {     
    const double rl = grid->r(i);
    const double r2 = rl * rl;
    const double r_up = 0.5*( grid->r(i+1) + rl );
    const double r_low = 0.5*( grid->r(i-1) + rl );
    //      const double dr = delta_x / dxdr(i);
    const double dr = r_up - r_low;
    for (int j = N_g; j < N_t-N_g; j++) {
      const double stl = grid->sintheta(j);
      for (int k = N_g; k < N_p-N_g; k++) {
	I_Re_int += fabs(I_Re(i,j,k)) * r2 * stl * sqrt(c->det(i,j,k))*exp(6.0*s->phi(i,j,k))
	  *dr*grid->delta_theta(j)*grid->delta_phi(k);
	I_Im_int += fabs(I_Im(i,j,k)) * r2 * stl * sqrt(c->det(i,j,k))*exp(6.0*s->phi(i,j,k))
	  *dr*grid->delta_theta(j)*grid->delta_phi(k);
	J_Re_int += fabs(J_Re(i,j,k)) * r2 * stl * sqrt(c->det(i,j,k))*exp(6.0*s->phi(i,j,k))
	  *dr*grid->delta_theta(j)*grid->delta_phi(k);
	J_Im_int += fabs(J_Im(i,j,k)) * r2 * stl * sqrt(c->det(i,j,k))*exp(6.0*s->phi(i,j,k))
	  *dr*grid->delta_theta(j)*grid->delta_phi(k);
      }
    }
  }



  return print;
};
