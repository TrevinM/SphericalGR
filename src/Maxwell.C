//================================================
//
// File that contains routines for class Maxwell
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

void Maxwell::Initialize(state* s, curvature* c, diagnostics* d) {
    //
    // set up initial data 
    //
    // NOTE: we'll assume that time derivative and shift vanish initially
    // 
    indata->Initialize_Maxwell(last->e_p, last->a_p);
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
void Maxwell::Compute_RHS(state* s, curvature* c, double time) {
    //
    // compute time derivative of vector potential A
    // 
    dot_a_p(inter, s, time);
    //
    // compute time derivative of electric field E
    // 
    dot_e_p(inter, s, c, time);
    //
    // Outer boundaries:
    //
    for (int i = 0; i < derivs->N_fcts; i++) {
        derivs->fct_list[i]->derivs_outerboundary(inter->fct_list[i]);
    }
};
//
//===============================================
// dot A  (Note: will assume that gauge variable phi = 0)
//===============================================
//
void Maxwell::dot_a_p(maxwell_state* m, state* s, double time) {
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = grid->r(i);
        for (int j = N_g; j < N_t - N_g; j++) {
            const double sinthetal = grid->sintheta(j);
            const double costhetal = grid->costheta(j);
            const double cotthetal = costhetal / sinthetal;
            // const double rst = rl*sinthetal;
            for (int k = N_g; k < N_p - N_g; k++) {
                //
                // first Lie derivative of  along shift
                //
                const double br = s->shift_r(i, j, k); // use rescaled vectors
                const double bt = s->shift_t(i, j, k);
                //	const double bp = s->shift_p(i,j,k);
#ifdef FLAT
                derivs->a_p[i][j][k] = 0.0;
#else
                derivs->a_p[i][j][k] = br * (m->a_p.dr(i, j, k, br) + m->a_p(i, j, k) / rl) +
                    bt * (m->a_p.dtheta(i, j, k, bt) + cotthetal * m->a_p(i, j, k)) / rl;
#endif /* FLAT */
                //
                // and compute RHS...
                //
                double gamma_pp = 1.0;
#ifndef FLAT
                const double psi4 = exp(4.0 * s->phi(i, j, k));
                gamma_pp = psi4 * (1.0 + s->h_pp(i, j, k));
#endif
                const double e_p_low = gamma_pp * m->e_p(i, j, k);
#ifdef FLAT
                derivs->a_p[i][j][k] -= e_p_low;
#else
                derivs->a_p[i][j][k] -= s->lapse(i, j, k) * e_p_low;
#endif  /* FLAT */
            }
        }
    }
};
//===============================================
// dot E
//===============================================
void Maxwell::dot_e_p(maxwell_state* m, state* s, curvature* c, double time) {
    //
    // just to make sure...
    //    
    c->det.fill_ghosts();
    //
    // and now...
    //
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = grid->r(i);
        for (int j = N_g; j < N_t - N_g; j++) {
            const double sinthetal = grid->sintheta(j);
            const double costhetal = grid->costheta(j);
            //      const double rst = rl*sinthetal;
            const double sin2theta = sinthetal * sinthetal;
            const double cotthetal = costhetal / sinthetal;
            for (int k = N_g; k < N_p - N_g; k++) {
#ifdef FLAT
                const double lapsel = 1.0;
#else
                const double lapsel = s->lapse(i, j, k);
#endif
                //
                // Compute Lie-derivative along shift...
                // 
                const double br = s->shift_r(i, j, k); // use rescaled vectors
                const double bt = s->shift_t(i, j, k);
                //	const double bp = s->shift_p(i,j,k);
                //
#ifdef FLAT
                derivs->e_p[i][j][k] = 0.0;
#else
                derivs->e_p[i][j][k] = br * (m->e_p.dr(i, j, k, br) - m->e_p(i, j, k) / rl) +
                    bt * (m->e_p.dtheta(i, j, k, bt) - cotthetal * m->e_p(i, j, k)) / rl
                    //
                    // ... and add extrinsic curvature term
                    // 
                    +lapsel * s->K(i, j, k) * m->e_p(i, j, k);
#endif
                // define physical (not conformal!) inverse metric, *not* rescaled
#ifdef FLAT
                tensor gup(1.0, 0.0, 0.0,
                    1.0 / (rl * rl), 0.0,
                    1.0 / (rl * rl * sin2theta));
#else
                const Doub psim4 = exp(-4.0 * s->phi(i, j, k));
                tensor gup(psim4 * c->gup_rr(i, j, k), psim4 * c->gup_rt(i, j, k), psim4 * c->gup_rp(i, j, k),
                    psim4 * c->gup_tt(i, j, k), psim4 * c->gup_tp(i, j, k), psim4 * c->gup_pp(i, j, k));
#endif
                // 
                // now define tensor \hat D_j \hat D_m A_l; also not rescaled
                //
                // IMPORTANT: object rank3tens assumes symmetry on last two indices,
                // so read this as 
                //            A_l;jm
                // with semicolon representing hatted cov. deriv.
                //
                // also: will assume axisymmetry, and only phi component of A 
                // is non-zero...
                //
                const double a_pl = m->a_p(i, j, k);
                const double a_p_r = m->a_p.dr(i, j, k);
                const double a_p_t = m->a_p.dtheta(i, j, k);
                const double a_p_rr = m->a_p.ddr(i, j, k);
                const double a_p_rt = m->a_p.drdtheta(i, j, k);
                const double a_p_tt = m->a_p.ddtheta(i, j, k);
                rank3tens DDA(0.0,                              // A_r;rr
                    0.0,
                    -sinthetal * (a_p_r - a_pl / rl),   // A_r;rp
                    0.0,
                    -sinthetal * (a_p_t - cotthetal * a_pl), // A_r,tp
                    0.0,
                    0.0,
                    0.0,
                    -rl * costhetal * (a_p_r - a_pl / rl),    // A_t,rp
                    0.0,
                    -rl * costhetal * (a_p_t - cotthetal * a_pl), // A_t,tp
                    0.0,                               // A_t,pp
                    rl * sinthetal * a_p_rr,               // A_p,rr
                    rl * sinthetal * (a_p_rt - a_p_t / rl), // A_p,rt
                    0.0,
                    rl * sinthetal * (a_p_tt + rl * a_p_r), // A_p,tt
                    0.0,
                    rl * sinthetal * (rl * sin2theta * a_p_r +   // A_p,pp
                        sinthetal * costhetal * a_p_t -
                        a_pl));
                //
                // now add second-derivative terms; rescale result by multiplying with r sin(theta)
                //
                Doub term1 = 0.0;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        for (int c = 0; c < 3; c++) {
                            term1 +=
                                (gup[2][a] * gup[c][b] - gup[c][a] * gup[2][b]) * DDA[b][c][a];
                        }
                derivs->e_p[i][j][k] += lapsel * term1 * rl * sinthetal;   // add rescaled version to e_p_RHS
                //
                // Now add all other terms..
                //
#ifndef FLAT
    //
    // define first derivatives of A_p; NOTE: ordered differently from above: DA[0][2] = D_r A_phi
    // (not rescaled)
    // 
                tensor DA(0.0, 0.0, rl * sinthetal * a_p_r,
                    0.0, 0.0, rl * sinthetal * a_p_t,
                    -sinthetal * a_pl, -rl * costhetal * a_pl, 0.0);
                //
                // Assign derivatives of physical, spatial metric to rank-3 tensor D_gamma
                //    example: D_gamma[0][1][2] = D_r \gamma_{tp} 
                //                              = e^{4 \phi} ( D_r \bar \gamma_{tp} + 4 \bar \gamma_{tp} D_r \phi ) 
                //                              = e^{4 \phi} D_r \bar \gamma_{tp} + 4 \gamma_{tp} D_r \phi 
                //
                tensor g_conf = gup.inverse();  // physical spatial metric, not conformally related
                const double e4p = exp(4.0 * s->phi(i, j, k));
                rank3tens D_gamma(e4p * c->Dr_e_rr(i, j, k) + 4.0 * g_conf[0][0] * s->phi.dr(i, j, k),
                    e4p * c->Dr_e_rt(i, j, k) + 4.0 * g_conf[0][1] * s->phi.dr(i, j, k),
                    e4p * c->Dr_e_rp(i, j, k) + 4.0 * g_conf[0][2] * s->phi.dr(i, j, k),
                    e4p * c->Dr_e_tt(i, j, k) + 4.0 * g_conf[1][1] * s->phi.dr(i, j, k),
                    e4p * c->Dr_e_tp(i, j, k) + 4.0 * g_conf[1][2] * s->phi.dr(i, j, k),
                    e4p * c->Dr_e_pp(i, j, k) + 4.0 * g_conf[2][2] * s->phi.dr(i, j, k),
                    e4p * c->Dt_e_rr(i, j, k) + 4.0 * g_conf[0][0] * s->phi.dtheta(i, j, k),
                    e4p * c->Dt_e_rt(i, j, k) + 4.0 * g_conf[0][1] * s->phi.dtheta(i, j, k),
                    e4p * c->Dt_e_rp(i, j, k) + 4.0 * g_conf[0][2] * s->phi.dtheta(i, j, k),
                    e4p * c->Dt_e_tt(i, j, k) + 4.0 * g_conf[1][1] * s->phi.dtheta(i, j, k),
                    e4p * c->Dt_e_tp(i, j, k) + 4.0 * g_conf[1][2] * s->phi.dtheta(i, j, k),
                    e4p * c->Dt_e_pp(i, j, k) + 4.0 * g_conf[2][2] * s->phi.dtheta(i, j, k),
                    e4p * c->Dp_e_rr(i, j, k) + 4.0 * g_conf[0][0] * s->phi.dphi(i, j, k),
                    e4p * c->Dp_e_rt(i, j, k) + 4.0 * g_conf[0][1] * s->phi.dphi(i, j, k),
                    e4p * c->Dp_e_rp(i, j, k) + 4.0 * g_conf[0][2] * s->phi.dphi(i, j, k),
                    e4p * c->Dp_e_tt(i, j, k) + 4.0 * g_conf[1][1] * s->phi.dphi(i, j, k),
                    e4p * c->Dp_e_tp(i, j, k) + 4.0 * g_conf[1][2] * s->phi.dphi(i, j, k),
                    e4p * c->Dp_e_pp(i, j, k) + 4.0 * g_conf[2][2] * s->phi.dphi(i, j, k));
                //
                // now compute derivatives of upper metric, ie.
                //     D_a \gamma^{bc} = - \gamma^{bd} \gamma^{ce} D_a \gamma_{de} 
                // 	
                rank3tens D_g_up;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        for (int c = 0; c < 3; c++) {
                            D_g_up[a][b][c] = 0.0;
                            for (int d = 0; d < 3; d++)
                                for (int e = 0; e < 3; e++)
                                    D_g_up[a][b][c] -= gup[b][d] * gup[c][e] * D_gamma[a][d][e];
                        }
                //
                // also set up derivatives of lapse, determinant, and conformal factor 
                // RECALL: code variable det is actually \bar \gamma / \hat \gamma 
                //
                vect dlapse(s->lapse.dr(i, j, k), s->lapse.dtheta(i, j, k), s->lapse.dphi(i, j, k));
                vect ddet(c->det.dr(i, j, k), c->det.dtheta(i, j, k), c->det.dphi(i, j, k));
                vect dphi(s->phi.dr(i, j, k), s->phi.dtheta(i, j, k), s->phi.dphi(i, j, k));
                //
                // combine derivatives to find 
                //  (\gamma / \hat \gamma)^{-1/2} D_j (\alpha \sqrt(\gamma / \hat \gamma)) = 
                //
                //	const double det_3D = exp(6.0 * (*phi_c)(i,j,k)) * sqrt((*det)(i,j,k)); 
                vect D_lapsedet;
                for (int a = 0; a < 3; a++)
                    D_lapsedet[a] = dlapse[a] + 6.0 * lapsel * dphi[a]
                    + 0.5 * lapsel * ddet[a] / c->det(i, j, k);
                Doub term2 = 0.0;
                Doub term3 = 0.0;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        for (int c = 0; c < 3; c++) {
                            term2 += (gup[a][b] * D_g_up[a][2][c] +
                                gup[2][c] * D_g_up[a][a][b] -
                                gup[2][b] * D_g_up[a][a][c] -
                                gup[a][c] * D_g_up[a][2][b]) * DA[c][b];
                            term3 += (gup[2][c] * gup[a][b] - gup[a][c] * gup[2][b]) * DA[c][b] * D_lapsedet[a];

                        }
                derivs->e_p[i][j][k] += (lapsel * term2 + term3) * rl * sinthetal;
#endif  /* FLAT */
                // finally add Kreiss-Oliger...
                derivs->e_p[i][j][k] += eta_KO * m->e_p.KO(i, j, k);
                //       }
                //     }
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
void Maxwell::ADM_Sources(state* s, curvature* c) {
    //
    // define 1 / ( 4 \pi )...
    // 
    const double oo4p = 1.0 / (4.0 * PI);
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = grid->r(i);
        for (int j = N_g; j < N_t - N_g; j++) {
            const double sinthetal = grid->sintheta(j);
            const double costhetal = grid->costheta(j);
            const double rst = rl * sinthetal;
            for (int k = N_g; k < N_p - N_g; k++) {
                //
                // define *physical* inverse metric...
                //
                const double em4p = exp(-4.0 * (s->phi)(i, j, k));
                tensor gup(em4p * (c->gup_rr)(i, j, k), em4p * (c->gup_rt)(i, j, k), em4p * (c->gup_rp)(i, j, k),
                    em4p * (c->gup_tt)(i, j, k), em4p * (c->gup_tp)(i, j, k), em4p * (c->gup_pp)(i, j, k));
                //
                // ... and physical metric:
                // 
                tensor gdown = gup.inverse();
                //
                // unrescaled version of eletric field:
                // 
                const double E_phi_up = inter->e_p(i, j, k) / rst;
                const double E_phi_down = gdown[2][2] * E_phi_up;
                const double E2 = E_phi_down * E_phi_up;
                //
                // define first derivatives of A_p; NOTE: ordered as in DA[0][2] = D_r A_phi
                // (not rescaled)
                // 
                const double a_p_r = inter->a_p.dr(i, j, k);
                const double a_p_t = inter->a_p.dtheta(i, j, k);
                const double a_pl = inter->a_p(i, j, k);
                tensor DA(0.0, 0.0, rl * sinthetal * a_p_r,
                    0.0, 0.0, rl * sinthetal * a_p_t,
                    -sinthetal * a_pl, -rl * costhetal * a_pl, 0.0);
                //
                // compute F_{ab} F^{ab}
                //
                double F2 = 0.0;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        for (int c = 0; c < 3; c++)
                            for (int d = 0; d < 3; d++)
                                F2 += gup[a][c] * gup[b][d] *
                                (DA[a][b] - DA[b][a]) * (DA[c][d] - DA[d][c]);
                F2 -= 2.0 * E2;
                //
                // Now compute ADM density...
                // 
                adm_sources->rho_ADM[i][j][k] = oo4p * (E2 + 0.25 * F2);
                //
                // ... fluxes...
                //
                adm_sources->S_r[i][j][k] = oo4p * E_phi_up * (DA[0][2] - DA[2][0]);
                adm_sources->S_t[i][j][k] = oo4p * E_phi_up * (DA[1][2] - DA[2][1]);
                adm_sources->S_p[i][j][k] = 0.0;     // Axisymmetry...
                //
                // ... stresses...
                // 
                tensor BB;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        BB[a][b] = 0.0;
                        for (int c = 0; c < 3; c++)
                            for (int d = 0; d < 3; d++)
                                BB[a][b] += gup[c][d] * (DA[a][d] - DA[d][a]) * (DA[b][c] - DA[c][b]);
                    }
                adm_sources->S_rr[i][j][k] = oo4p * (BB[0][0] - 0.25 * gdown[0][0] * F2);
                adm_sources->S_rt[i][j][k] = oo4p * (BB[0][1] - 0.25 * gdown[0][1] * F2);
                adm_sources->S_rp[i][j][k] = oo4p * (BB[0][2] - 0.25 * gdown[0][2] * F2);
                adm_sources->S_tt[i][j][k] = oo4p * (BB[1][1] - 0.25 * gdown[1][1] * F2);
                adm_sources->S_tp[i][j][k] = oo4p * (BB[1][2] - 0.25 * gdown[1][2] * F2);
                adm_sources->S_pp[i][j][k] = oo4p * (BB[2][2] - 0.25 * gdown[2][2] * F2 - E_phi_down * E_phi_down);
                //
                // ... and trace ...
                //
                adm_sources->trace_S[i][j][k] = adm_sources->rho_ADM[i][j][k];
                //
            }
        }
    }
    //
    adm_sources->rho_ADM.fill_ghosts();
    rho_center = adm_sources->rho_ADM(0.0, N_g, N_g);
    if (rho_center > rho_c_max) rho_c_max = rho_center;
    drhoddr = adm_sources->rho_ADM.ddr(N_g, N_g, N_g);
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
