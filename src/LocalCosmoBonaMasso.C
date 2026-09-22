#include "nr3.h"
#include "tensors.h"
#include "Cosmology.h"
#include "Grid.h"
#include "State.h"
#include "Bona_Masso_function.h"
#include "Matter.h"
#include "Slicing.h"

inline void LocalCosmoBonaMasso::dot_lapse(state * c, state * derivs, double t) {
    const double K_0 = (*cosmology).K0(t);
    for (int i = N_g; i < N_r-N_g; i++) {
      const double rl = grid->r(i);
      for (int j = N_g; j < N_t-N_g; j++) {
	const double sintheta = grid->sintheta(j);
	for (int k = N_g; k < N_p-N_g; k++) {
	  const double lapse_l = c->lapse(i,j,k);
	  const double f = (*bonamasso_f)(lapse_l);
	  const double shift_r = c->shift_r(i,j,k);
	  const double shift_t = c->shift_t(i,j,k) / rl;
	  const double shift_p = c->shift_p(i,j,k) / (rl * sintheta);
	  // Compute A^2, need to add the non-conformally flat term!!!!!!!
	  const double a_rr = c->a_rr(i,j,k);
	  const double a_rt = c->a_rt(i,j,k);
	  const double a_rp = c->a_rp(i,j,k);
	  const double a_tt = c->a_tt(i,j,k);
	  const double a_tp = c->a_tp(i,j,k);
	  const double a_pp = c->a_pp(i,j,k);
	  const double h_rr = c->h_rr(i,j,k);
	  const double h_rt = c->h_rt(i,j,k);
	  const double h_rp = c->h_rp(i,j,k);
	  const double h_tt = c->h_tt(i,j,k);
	  const double h_tp = c->h_tp(i,j,k);
	  const double h_pp = c->h_pp(i,j,k);
	  tensor a_l(a_rr, a_rt, a_rp, a_tt, a_tp, a_pp);
	  tensor g_down(1. + h_rr, h_rt, h_rp, 1. + h_tt, h_tp, 1. + h_pp);
	  tensor g_up = g_down.inverse();
	  double a2 = 0.0;
	  for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			for (int c = 0; c < 3; c++) {
				for (int d = 0; d<3; d++){
					a2 += g_up[a][c] * g_up[b][d] * a_l[c][d] * a_l[a][b];
				}
			}
		}
	  }
	  derivs->lapse[i][j][k] = - lapse_l * lapse_l * f *
	    ( c->K(i,j,k) + sqrt(24*PI * matter->adm_sources->rho_ADM(i,j,k)  + (3./2.) * a2))
		 + eta_KO * c->lapse.KO(i,j,k) +
	    + shift_r * c->lapse.dr(i,j,k,shift_r) 
	    + shift_t * c->lapse.dtheta(i,j,k,shift_t) 
	    + shift_p * c->lapse.dphi(i,j,k,shift_p);
		//cout << " LOCALCOSMOBONAMASSO: dot_lapse = " << derivs->lapse[i][j][k] << endl;
 	}
      }
    }
  };