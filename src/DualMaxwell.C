//================================================
//
// File that contains routines for class DualMaxwell
//
//================================================
#include "Grid.h"
#include "Matter.h"
#include <ctime>
#include "tensors.h"
//
//================================================
// Initialize electromagnetic fields
//================================================
// #define FLAT

void DualMaxwell::Initialize(state* s, curvature* c, diagnostics* d) {
    //
    // set up initial data 
    //
    // NOTE: we'll assume that time derivative and shift vanish initially
    // 
    indata->Initialize_DualMaxwell(last->a_r, last->a_t, last->a_p,
        last->as_r, last->as_t, last->as_p);
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
};
//===============================================
// Compute RHS sides for matter equations
//===============================================
void DualMaxwell::Compute_RHS(state* s, curvature* c, double time) {
    //
    // compute time derivative of vector potentials A and A*
    // 
    dot_a_as(inter, s, time);
    //
    // Outer boundaries:
    //
    for (int i = 0; i < derivs->N_fcts; i++) {
        derivs->fct_list[i]->derivs_outerboundary(inter->fct_list[i]);
    }

    //   dump->slice(0.0,0.0,12473, derivs->a_r.Address());
    //   dump->slice(0.0,0.0,12473, derivs->a_t.Address());
    //   dump->slice(0.0,0.0,12473, derivs->a_p.Address());
    //   dump->slice(0.0,0.0,12473, derivs->as_r.Address());
    //   dump->slice(0.0,0.0,12473, derivs->as_t.Address());
    //   dump->slice(0.0,0.0,12473, derivs->as_p.Address());

    //   cout << "Just Dumped. Goodbye !!" << endl;
    //   exit(0);

};
//
//===============================================
// dot A and dot A*  (Note: will assume that gauge variables phi = phi* = 0)
//===============================================
//
void DualMaxwell::dot_a_as(dualmaxwell_state* m, state* s, double time) {
    //
    // start with Lie derivatives
    //
    Lie_covariant(s->shift_r, s->shift_t, s->shift_p,
        m->a_r, m->a_t, m->a_p,
        derivs->a_r, derivs->a_t, derivs->a_p);
    Lie_covariant(s->shift_r, s->shift_t, s->shift_p,
        m->as_r, m->as_t, m->as_p,
        derivs->as_r, derivs->as_t, derivs->as_p);
    //
    // now go to each gridpoint
    //
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = grid->r(i);
        for (int j = N_g; j < N_t - N_g; j++) {
            const double sinthetal = grid->sintheta(j);
            const double costhetal = grid->costheta(j);
            const double cotthetal = costhetal / sinthetal;
            // const double rst = rl*sinthetal;
            for (int k = N_g; k < N_p - N_g; k++) {
                const double psil = exp(s->phi(i, j, k));
                const double psi2 = psil * psil;
                const double psi4 = psil * psil * psil * psil;
                const double lapsel = s->lapse(i, j, k);
                //
                // compute curl of a and as
                //
                double curl_a_r, curl_a_t, curl_a_p;   // rescaled, indices upstairs
                curl(inter->a_r, inter->a_t, inter->a_p,
                    curl_a_r, curl_a_t, curl_a_p, i, j, k, s);
                vect curl_a(curl_a_r, curl_a_t, curl_a_p);
                double curl_as_r, curl_as_t, curl_as_p;   // rescaled, indices upstairs
                curl(inter->as_r, inter->as_t, inter->as_p,
                    curl_as_r, curl_as_t, curl_as_p, i, j, k, s);
                vect curl_as(curl_as_r, curl_as_t, curl_as_p);
                //
                // assign physical (but rescaled) metric 
                //
                tensor g_down(psi4 * (1.0 + s->h_rr(i, j, k)), psi4 * s->h_rt(i, j, k), psi4 * s->h_rp(i, j, k),
                    psi4 * (1.0 + s->h_tt(i, j, k)), psi4 * s->h_tp(i, j, k), psi4 * (1.0 + s->h_pp(i, j, k)));
                //
                // now lower index on curl
                // 
                vect curl_a_down(0.0, 0.0, 0.0);
                vect curl_as_down(0.0, 0.0, 0.0);
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        curl_a_down[a] += g_down[a][b] * curl_a[b];
                        curl_as_down[a] += g_down[a][b] * curl_as[b];
                    }
                //
                // now add to right-hand sides:
                // 
                derivs->a_r[i][j][k] -= lapsel * curl_as_down[0];
                derivs->a_t[i][j][k] -= lapsel * curl_as_down[1];
                derivs->a_p[i][j][k] -= lapsel * curl_as_down[2];
                derivs->as_r[i][j][k] += lapsel * curl_a_down[0];
                derivs->as_t[i][j][k] += lapsel * curl_a_down[1];
                derivs->as_p[i][j][k] += lapsel * curl_a_down[2];
                //
                // finally add Kreiss-Oliger terms
                // 
                derivs->a_r[i][j][k] += eta_KO * m->a_r.KO(i, j, k);
                derivs->a_t[i][j][k] += eta_KO * m->a_t.KO(i, j, k);
                derivs->a_p[i][j][k] += eta_KO * m->a_p.KO(i, j, k);
                derivs->as_r[i][j][k] += eta_KO * m->as_r.KO(i, j, k);
                derivs->as_t[i][j][k] += eta_KO * m->as_t.KO(i, j, k);
                derivs->as_p[i][j][k] += eta_KO * m->as_p.KO(i, j, k);
            }
        }
    }
    derivs->a_r.fill_ghosts();
    derivs->a_t.fill_ghosts();
    derivs->a_p.fill_ghosts();
    derivs->as_r.fill_ghosts();
    derivs->as_t.fill_ghosts();
    derivs->as_p.fill_ghosts();

};


//================================================
//
// Lie derivative of covariant rank-1 tensor 
// NOTE: uses upwind differencing for advective derivatives of a (but not of shift)
//
//================================================  
void DualMaxwell::Lie_covariant(gf3d& b_r, gf3d& b_t, gf3d& b_p,
    gf3d& a_r, gf3d& a_t, gf3d& a_p,
    gf3d& Ll_r, gf3d& Ll_t, gf3d& Ll_p) {
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = grid->r(i);
        const double r2 = rl * rl;
        for (int j = N_g; j < N_t - N_g; j++) {
            const double stl = grid->sintheta(j);
            const double st2 = stl * stl;
            const double ctl = grid->costheta(j);
            for (int k = N_g; k < N_p - N_g; k++) {
                // 
                // construct physical variables from vector variables
                // 
                const double beta_r = b_r(i, j, k);
                const double beta_t = b_t(i, j, k) / rl;
                const double beta_p = b_p(i, j, k) / (rl * stl);
                const double A_r = a_r(i, j, k);
                const double A_t = a_t(i, j, k) * rl;
                const double A_p = a_p(i, j, k) * (rl * stl);
                //
                // compute partial derivatives of beta and A in terms of code
                // variables b and a.  Notation:  b^i_{,j} = b_i_j
                //
                const double beta_r_r = b_r.dr(i, j, k);
                const double beta_r_t = b_r.dtheta(i, j, k);
                const double beta_r_p = b_r.dphi(i, j, k);
                const double beta_t_r = b_t.dr(i, j, k) / rl - b_t(i, j, k) / r2;
                const double beta_t_t = b_t.dtheta(i, j, k) / rl;
                const double beta_t_p = b_t.dphi(i, j, k) / rl;
                const double beta_p_r = b_p.dr(i, j, k) / (rl * stl) - b_p(i, j, k) / (r2 * stl);
                const double beta_p_t = b_p.dtheta(i, j, k) / (rl * stl) - b_p(i, j, k) * ctl / (rl * st2);
                const double beta_p_p = b_p.dphi(i, j, k) / (rl * stl);

                const double A_r_r = a_r.dr(i, j, k, beta_r);
                const double A_r_t = a_r.dtheta(i, j, k, beta_t);
                const double A_r_p = a_r.dphi(i, j, k, beta_p);
                const double A_t_r = a_t.dr(i, j, k, beta_r) * rl + a_t(i, j, k);
                const double A_t_t = a_t.dtheta(i, j, k, beta_t) * rl;
                const double A_t_p = a_t.dphi(i, j, k, beta_p) * rl;
                const double A_p_r = a_p.dr(i, j, k, beta_r) * (rl * stl) + a_p(i, j, k) * stl;
                const double A_p_t = a_p.dtheta(i, j, k, beta_t) * (rl * stl) + a_p(i, j, k) * ctl * rl;
                const double A_p_p = a_p.dphi(i, j, k, beta_p) * (rl * stl);
                //
                // now compute Lie derivative
                //
                Ll_r[i][j][k] =
                    beta_r * A_r_r + beta_t * A_r_t + beta_p * A_r_p +
                    A_r * beta_r_r + A_t * beta_t_r + A_p * beta_p_r;

                Ll_t[i][j][k] =
                    beta_r * A_t_r + beta_t * A_t_t + beta_p * A_t_p +
                    A_r * beta_t_t + A_t * beta_t_t + A_p * beta_p_t;

                Ll_p[i][j][k] =
                    beta_r * A_p_r + beta_t * A_p_t + beta_p * A_p_p +
                    A_r * beta_r_p + A_t * beta_t_p + A_p * beta_p_p;
                //
                // finally rescale coefficients
                //
                Ll_t[i][j][k] /= rl;
                Ll_p[i][j][k] /= rl * stl;
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
void DualMaxwell::ADM_Sources(state* s, curvature* c) {
    //
    // define 1 / ( 4 \pi )...
    // 
    const double oo4p = 1.0 / (4.0 * PI);
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = grid->r(i);
        const double r2 = rl * rl;
        for (int j = N_g; j < N_t - N_g; j++) {
            const double sinthetal = grid->sintheta(j);
            const double costhetal = grid->costheta(j);
            const double rst = rl * sinthetal;
            for (int k = N_g; k < N_p - N_g; k++) {
                const double psil = exp(s->phi(i, j, k));
                const double psi2 = psil * psil;
                const double psi4 = psil * psil * psil * psil;
                const double psi6 = psi2 * psi4;
                const double psim4 = 1. / psi4;
                //
                // get B^i and E^i from curl of A and *A
                // 
                double b_r, b_t, b_p;   // rescaled, indices upstairs
                curl(inter->a_r, inter->a_t, inter->a_p, b_r, b_t, b_p, i, j, k, s);
                aux->B_r[i][j][k] = b_r;
                aux->B_t[i][j][k] = b_t;
                aux->B_p[i][j][k] = b_p;
                double e_r, e_t, e_p;   // rescaled, indices upstairs
                curl(inter->as_r, inter->as_t, inter->as_p, e_r, e_t, e_p, i, j, k, s);
                aux->E_r[i][j][k] = e_r;
                aux->E_t[i][j][k] = e_t;
                aux->E_p[i][j][k] = e_p;
                //
                // assign physical (but rescaled) metric 
                //
                tensor g_down(psi4 * (1.0 + s->h_rr(i, j, k)), psi4 * s->h_rt(i, j, k), psi4 * s->h_rp(i, j, k),
                    psi4 * (1.0 + s->h_tt(i, j, k)), psi4 * s->h_tp(i, j, k), psi4 * (1.0 + s->h_pp(i, j, k)));
                //
                // now compute dot products and lower components
                //
                vect b_up(b_r, b_t, b_p);
                vect e_up(e_r, e_t, e_p);
                vect b_down(0.0, 0.0, 0.0);
                vect e_down(0.0, 0.0, 0.0);
                double b2 = 0.0;
                double e2 = 0.0;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        b2 += g_down[a][b] * b_up[a] * b_up[b];
                        e2 += g_down[a][b] * e_up[a] * e_up[b];
                        b_down[a] += g_down[a][b] * b_up[b];
                        e_down[a] += g_down[a][b] * e_up[b];
                    }
                //
                // compute energy density
                //
                adm_sources->rho_ADM[i][j][k] = oo4p * (e2 + b2) / 2.0;
                //
                // compute *unrescaled* downstairs momentum densities
                //
                adm_sources->S_r[i][j][k] = psi6 * oo4p * (e_t * b_p - e_p * b_t);
                adm_sources->S_t[i][j][k] = psi6 * oo4p * (e_p * b_r - e_r * b_p) * rl;
                adm_sources->S_p[i][j][k] = psi6 * oo4p * (e_r * b_t - e_t * b_r) * rst;
                //
                // compute *unrescaled* downstairs components of the stress tensor
                //
                adm_sources->S_rr[i][j][k] = oo4p * (0.5 * g_down[0][0] * (e2 + b2)
                    - e_down[0] * e_down[0] - b_down[0] * b_down[0]);
                adm_sources->S_rt[i][j][k] = oo4p * (0.5 * g_down[0][1] * (e2 + b2)
                    - e_down[0] * e_down[1] - b_down[0] * b_down[1]) * rl;
                adm_sources->S_rp[i][j][k] = oo4p * (0.5 * g_down[0][2] * (e2 + b2)
                    - e_down[0] * e_down[2] - b_down[0] * b_down[2]) * rst;
                adm_sources->S_tt[i][j][k] = oo4p * (0.5 * g_down[1][1] * (e2 + b2)
                    - e_down[1] * e_down[1] - b_down[1] * b_down[1]) * r2;
                adm_sources->S_tp[i][j][k] = oo4p * (0.5 * g_down[1][2] * (e2 + b2)
                    - e_down[1] * e_down[2] - b_down[1] * b_down[2]) * rl * rst;
                adm_sources->S_pp[i][j][k] = oo4p * (0.5 * g_down[2][2] * (e2 + b2)
                    - e_down[2] * e_down[2] - b_down[2] * b_down[2]) * rst * rst;
                //
                // ... and trace ...
                //
                adm_sources->trace_S[i][j][k] = adm_sources->rho_ADM[i][j][k];
                //
                // sanity check:
                //
                // if (i == N_g && j == N_g && k == N_g) {
                //   tensor gup(psim4 * c->gup_rr(i,j,k), psim4 * c->gup_rt(i,j,k), psim4 * c->gup_rp(i,j,k),
                // 	     psim4 * c->gup_tt(i,j,k), psim4 * c->gup_tp(i,j,k), psim4 * c->gup_pp(i,j,k));
                //   tensor Sdown(adm_sources->S_rr(i,j,k), adm_sources->S_rt(i,j,k), adm_sources->S_rp(i,j,k),
                // 	       adm_sources->S_tt(i,j,k), adm_sources->S_tp(i,j,k), adm_sources->S_pp(i,j,k)); 
                //   double alt_trace = 0.0;
                //   for (int a = 0; a < 3; a++)
                //     for (int b = 0; b < 3; b++) {
                //       alt_trace += gup[a][b] * Sdown[a][b];
                //     }
                //   cout << " DUALMAXWELL: CHECK Trace S : " << adm_sources->trace_S(i,j,k)
                //        << " or " << alt_trace << endl;
                // }
            }
        }
    }
    //
    adm_sources->S_r.fill_ghosts();
    adm_sources->S_t.fill_ghosts();
    adm_sources->S_p.fill_ghosts();
    adm_sources->rho_ADM.fill_ghosts();
    rho_center = adm_sources->rho_ADM(0.0, N_g, N_g);
    if (rho_center > rho_c_max) rho_c_max = rho_center;
    rho_max = adm_sources->rho_ADM.max(rho_i, rho_j, rho_k);
    if (rho_max > rho_max_MAX) rho_max_MAX = rho_max;
    drhoddr = adm_sources->rho_ADM.ddr(N_g, N_g, N_g);
    s->lapse.min(lapse_i, lapse_j, lapse_k);
    //
    // compute diagnostic function Omega
    //
    // for (int i = N_g; i < N_r; i++) {   
    //   double r_l = rho_ADM.r(i);
    //   for (int j = N_g; j < N_theta - N_g; j++) 
    //     for (int k = N_g; k < N_phi - N_g; k++) {
    // 	double psi_l = exp((*phi_o)(i,j,k));
    // 	double h_pp_l = (*h_pp_o)(i,j,k);
    // 	double r_areal = r_l*psi_l*psi_l*sqrt(1.0 + h_pp_l);
    // 	Omega[i][j][k] = 4.0*PI*r_areal*r_areal*rho_ADM(i,j,k);
    //     }
    // }
    // Omega.fill_ghosts();
    // 
};


//===============================================================
// Compute physical curl: returns scaled eps^ijk ( D_j A_k - D_k A_j )
// where D_i is covariant derivative with respect to reference metric,
// but epsilon is physical (rescaled) 3D epsilon: psi^{-6}[ijk]
//===============================================================
void DualMaxwell::curl(gf3d& a_r, gf3d& a_t, gf3d& a_p,
    double& curl_a_r, double& curl_a_t, double& curl_a_p,
    int i, int j, int k, state* s) {
    const double rl = a_r.r(i);
    const double stl = a_r.sintheta(j);
    const double ctl = a_r.costheta(j);
    const double cot = ctl / stl;
    const double psil = exp(s->phi(i, j, k));
    const double psim6 = 1. / (psil * psil * psil * psil * psil * psil);
    // const double D_r_A_t = rl*a_t.dr(i,j,k);
    // const double D_r_A_p = rl*stl*a_p.dr(i,j,k);
    // const double D_t_A_r = a_r.dtheta(i,j,k) - a_t(i,j,k);
    // const double D_t_A_p = rl*stl*a_p.dtheta(i,j,k);
    // const double D_p_A_r = - stl*a_p(i,j,k);
    // const double D_p_A_t = - rl*ctl*a_p(i,j,k);
    curl_a_r = psim6 * (a_p.dtheta(i, j, k) + cot * a_p(i, j, k)) / rl;
    curl_a_t = -psim6 * (a_p.dr(i, j, k) + a_p(i, j, k) / rl);
    curl_a_p = psim6 * (a_t.dr(i, j, k) - (a_r.dtheta(i, j, k) - a_t(i, j, k)) / rl);
}
