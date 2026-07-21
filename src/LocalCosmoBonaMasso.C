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
	  derivs->lapse[i][j][k] = - lapse_l * lapse_l * f *
	    ( c->K(i,j,k) + sqrt(24*PI * matter->adm_sources->rho_ADM(i,j,k) ) ) + eta_KO * c->lapse.KO(i,j,k) +
	    + shift_r * c->lapse.dr(i,j,k,shift_r) 
	    + shift_t * c->lapse.dtheta(i,j,k,shift_t) 
	    + shift_p * c->lapse.dphi(i,j,k,shift_p);
 	}
      }
    }
  };