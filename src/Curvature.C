//
//==================================================================================
// Routines that belong to class Curvature and -- guess what -- compute curvature
//==================================================================================
//
#include "tensors.h"
#include "Manager.h"

/*  Notation...

We split the metric into the form

\bar \gamma_{ij} = \eta_{ij} + \epsilon_{ij}

where \eta_{ij} is the flat metric (in spherical polar coordinates), and where we write
\epsilon_{ij} in the form

                /  h_rr                  r h_rt                  r sin(theta) h_rp     \
\epsilon_{ij} = |  r h_rt                r^2 h_tt                r^2 sin(theta) h_tp   |
                \  r sin(theta) h_rp     r^2 sin(theta) h_tp     r^2 sin^2(theta) h_pp /

In the first routine we compute the covariant derivative of the metric

D_i \bar \gamma_{jk} = D_i \epsilon_{jk},

where the covariant derivative is that associated with the flat metric \eta_{ij}.

We write the (traceless, conformally rescaled) extrinsic curvature in a similar way:

                /  a_rr                  r a_rt                  r sin(theta) a_rp     \
\bar A_{ij}   = |  r a_rt                r^2 a_tt                r^2 sin(theta) a_tp   |
                \  r sin(theta) a_rp     r^2 sin(theta) a_tp     r^2 sin^2(theta) a_pp /

Finally, we write the conformal connection functions \Lambda^i in the form

              /  lam_r                  \
\Lambda^i  =  |  lam_t / r              |
              \  lam_p / (r sin(theta)) /

*/

//=========================================================================================
// Compute (covariant) derivatives of metric from current state
//=========================================================================================

void curvature::Compute_Metric_Derivatives(state* c) {
#ifdef _VERBOSE_
    cout << " Entering Compute_Metric_Derivatives..." << endl;
#endif
    // Covariant derivatives with respect to the flat metric in spherical polar coordinates...
    //
#pragma omp parallel for collapse(3)  
    for (int i = N_g; i < N_r - N_g; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                const double rl = grid->r(i);
                const double r2 = rl * rl;
                const double stl = grid->sintheta(j);
                const double st2 = stl * stl;
                const double ctl = grid->costheta(j);
                Dr_e_rr[i][j][k] = c->h_rr.dr(i, j, k);
                Dr_e_rt[i][j][k] = rl * c->h_rt.dr(i, j, k);
                Dr_e_rp[i][j][k] = rl * stl * c->h_rp.dr(i, j, k);
                Dr_e_tt[i][j][k] = r2 * c->h_tt.dr(i, j, k);
                Dr_e_tp[i][j][k] = r2 * stl * c->h_tp.dr(i, j, k);
                Dr_e_pp[i][j][k] = r2 * st2 * c->h_pp.dr(i, j, k);

                Dt_e_rr[i][j][k] = c->h_rr.dtheta(i, j, k) - 2.0 * c->h_rt(i, j, k);
                Dt_e_rt[i][j][k] = rl * (c->h_rt.dtheta(i, j, k) + c->h_rr(i, j, k) - c->h_tt(i, j, k));
                Dt_e_rp[i][j][k] = rl * stl * (c->h_rp.dtheta(i, j, k) - c->h_tp(i, j, k));
                Dt_e_tt[i][j][k] = r2 * (c->h_tt.dtheta(i, j, k) + 2.0 * c->h_rt(i, j, k));
                Dt_e_tp[i][j][k] = r2 * stl * (c->h_tp.dtheta(i, j, k) + c->h_rp(i, j, k));
                Dt_e_pp[i][j][k] = r2 * st2 * c->h_pp.dtheta(i, j, k);

                Dp_e_rr[i][j][k] = c->h_rr.dphi(i, j, k) - 2.0 * stl * c->h_rp(i, j, k);
                Dp_e_rt[i][j][k] = rl * (c->h_rt.dphi(i, j, k) - ctl * c->h_rp(i, j, k) - stl * c->h_tp(i, j, k));
                Dp_e_rp[i][j][k] = rl * stl * (c->h_rp.dphi(i, j, k) + stl * c->h_rr(i, j, k) +
                    ctl * c->h_rt(i, j, k) - stl * c->h_pp(i, j, k));
                Dp_e_tt[i][j][k] = r2 * (c->h_tt.dphi(i, j, k) - 2.0 * ctl * c->h_tp(i, j, k));
                Dp_e_tp[i][j][k] = r2 * stl * (c->h_tp.dphi(i, j, k) + stl * c->h_rt(i, j, k) +
                    ctl * c->h_tt(i, j, k) - ctl * c->h_pp(i, j, k));
                Dp_e_pp[i][j][k] = r2 * st2 * (c->h_pp.dphi(i, j, k) + 2.0 * stl * c->h_rp(i, j, k) +
                    2.0 * ctl * c->h_tp(i, j, k));
            }
        }
    }
};

//=========================================================================================
// Now compute inverse of metric from current state
//=========================================================================================
void curvature::Compute_Inverse_Metric(state* c) {
#ifdef _VERBOSE_
    cout << " Entering Compute_Inverse_Metric..." << endl;
#endif
    //
    // first compute determinant (without factors of r^4 sin^2 theta - will take care of that for individual
    //                            metric components below)
    //
    Compute_Determinant(c);
    //
#pragma omp parallel for collapse(3)
    for (int i = 0; i < N_r; i++) {
        // CHECK: now include ghost zones
        for (int j = 0; j < N_t; j++) {
            for (int k = 0; k < N_p; k++) {
                const double rl = grid->r(i);
                const double r2 = rl * rl;
                const double stl = grid->sintheta(j);
                const double st2 = stl * stl;
                const double det_l = det(i, j, k);
                gup_rr[i][j][k] = ((1.0 + c->h_tt(i, j, k)) * (1.0 + c->h_pp(i, j, k)) - c->h_tp(i, j, k) * c->h_tp(i, j, k)) / det_l;
                gup_rt[i][j][k] = -(c->h_rt(i, j, k) * (1.0 + c->h_pp(i, j, k)) - c->h_tp(i, j, k) * c->h_rp(i, j, k)) / rl / det_l;
                gup_rp[i][j][k] = (c->h_rt(i, j, k) * c->h_tp(i, j, k) - (1.0 + c->h_tt(i, j, k)) * c->h_rp(i, j, k)) / (rl * stl) / det_l;
                gup_tt[i][j][k] = ((1.0 + c->h_rr(i, j, k)) * (1.0 + c->h_pp(i, j, k)) - c->h_rp(i, j, k) * c->h_rp(i, j, k)) / r2 / det_l;
                gup_tp[i][j][k] = -((1.0 + c->h_rr(i, j, k)) * c->h_tp(i, j, k) - c->h_rt(i, j, k) * c->h_rp(i, j, k)) / (r2 * stl) / det_l;
                gup_pp[i][j][k] = ((1.0 + c->h_rr(i, j, k)) * (1.0 + c->h_tt(i, j, k)) - c->h_rt(i, j, k) * c->h_rt(i, j, k)) / (r2 * st2) / det_l;
            }
        }
    }
};

//============================================================================
// Compute determinant of metric (divided by r^4 sin^2 theta)
//============================================================================
int curvature::Compute_Determinant(state* c) {
    int error = 0;
#ifdef _VERBOSE_
    cout << " Entering Compute_Determinant()..." << endl;
#endif
#pragma omp parallel for collapse(3)
    for (int i = 0; i < N_r; i++)
        for (int j = 0; j < N_t; j++)
            for (int k = 0; k < N_p; k++) {
                //
                // Note: this is determinant divided by r^4 sin^2 theta. 
                //
                det[i][j][k] = ((1.0 + c->h_rr(i, j, k)) * (1.0 + c->h_tt(i, j, k)) * (1.0 + c->h_pp(i, j, k)) +
                    c->h_rt(i, j, k) * c->h_tp(i, j, k) * c->h_rp(i, j, k) +
                    c->h_rp(i, j, k) * c->h_rt(i, j, k) * c->h_tp(i, j, k) -
                    (1.0 + c->h_rr(i, j, k)) * c->h_tp(i, j, k) * c->h_tp(i, j, k) -
                    c->h_rt(i, j, k) * c->h_rt(i, j, k) * (1.0 + c->h_pp(i, j, k)) -
                    c->h_rp(i, j, k) * (1.0 + c->h_tt(i, j, k)) * c->h_rp(i, j, k));
                if (det(i, j, k) <= 0.0) error = 1;
                if (!isfinite(det(i, j, k))) {
                    // cout << " Compute_Determinant: determinant = " 
                    // << det(i,j,k) << " at (" 
                    // << i << "," << j << "," << k << ")" << endl;
                    error = 2;
                }
            }
    return error;
};
//============================================================================
// Compute trace of (traceless part of) extrinsic curvature
//============================================================================
void curvature::Compute_Trace(state* c) {
#ifdef _VERBOSE_
    cout << " Entering Compute_Trace()..." << endl;
#endif
    for (int i = 0; i < N_r; i++) {
        const double rl = grid->r(i);
        const double r2 = rl * rl;
        for (int j = 0; j < N_t; j++) {
            const double stl = grid->sintheta(j);
            const double st2 = stl * stl;
            for (int k = 0; k < N_p; k++) {
                trace_A[i][j][k] =
                    gup_rr(i, j, k) * c->a_rr(i, j, k) +
                    gup_rt(i, j, k) * c->a_rt(i, j, k) * rl * 2.0 +
                    gup_rp(i, j, k) * c->a_rp(i, j, k) * rl * stl * 2.0 +
                    gup_tt(i, j, k) * c->a_tt(i, j, k) * r2 +
                    gup_tp(i, j, k) * c->a_tp(i, j, k) * r2 * stl * 2.0 +
                    gup_pp(i, j, k) * c->a_pp(i, j, k) * r2 * st2;
            }
        }
    }
};
/*

Next we compute the connection

\Delta \Gamma^i_{jk} \equiv \bar \Gamma^i_{jk} - \Gamma0^i_{jk}
                     = (1/2) * \bar \gamma^{il} ( D_j \bar \gamma_{kl} + D_k \bar \gamma_{jl} - D_l \bar \gamma_{jk} )
                     = (1/2) * \bar \gamma^{il} ( D_j \epsilon_{kl} + D_k \epsilon_{jl} - D_l \epsilon_{jk} )

where D_i is the covariant derivative associated with the flat metric \eta_{ij}.

*/

void curvature::Compute_Connection() {
#ifdef _VERBOSE_
    cout << " Entering Compute_Connection..." << endl;
#endif
#pragma omp parallel for collapse(3)
    for (int ig = N_g; ig < N_r - N_g; ig++)
        for (int jg = N_g; jg < N_t - N_g; jg++)
            for (int kg = N_g; kg < N_p - N_g; kg++) {
                //
                // set up tensors...
                //
                // Assign upper metric to tensor
                tensor gup(gup_rr[ig][jg][kg], gup_rt[ig][jg][kg], gup_rp[ig][jg][kg], gup_tt[ig][jg][kg], gup_tp[ig][jg][kg], gup_pp[ig][jg][kg]);
                // Assign derivatives of metric to rank-3 tensor
                //    example: D_eps[0][1][2] = D_r \epsilon_{tp}
                rank3tens D_eps(Dr_e_rr[ig][jg][kg], Dr_e_rt[ig][jg][kg], Dr_e_rp[ig][jg][kg],
                    Dr_e_tt[ig][jg][kg], Dr_e_tp[ig][jg][kg], Dr_e_pp[ig][jg][kg],
                    Dt_e_rr[ig][jg][kg], Dt_e_rt[ig][jg][kg], Dt_e_rp[ig][jg][kg],
                    Dt_e_tt[ig][jg][kg], Dt_e_tp[ig][jg][kg], Dt_e_pp[ig][jg][kg],
                    Dp_e_rr[ig][jg][kg], Dp_e_rt[ig][jg][kg], Dp_e_rp[ig][jg][kg],
                    Dp_e_tt[ig][jg][kg], Dp_e_tp[ig][jg][kg], Dp_e_pp[ig][jg][kg]);
                // Now loop over indices
                rank3tens Delta_Gam;
                int nn = 3;
                for (int i = 0; i < nn; i++)
                    for (int j = 0; j < nn; j++)
                        for (int k = 0; k < nn; k++) {
                            Delta_Gam[i][j][k] = 0.0;
                            for (int l = 0; l < nn; l++)
                                // save factor of 1/2 for later...
                                Delta_Gam[i][j][k] += gup[i][l] * (D_eps[j][k][l] + D_eps[k][j][l] - D_eps[l][j][k]);
                        }
                //
                DG_r_rr[ig][jg][kg] = 0.5 * Delta_Gam[0][0][0];
                DG_r_rt[ig][jg][kg] = 0.5 * Delta_Gam[0][0][1];
                DG_r_rp[ig][jg][kg] = 0.5 * Delta_Gam[0][0][2];
                DG_r_tt[ig][jg][kg] = 0.5 * Delta_Gam[0][1][1];
                DG_r_tp[ig][jg][kg] = 0.5 * Delta_Gam[0][1][2];
                DG_r_pp[ig][jg][kg] = 0.5 * Delta_Gam[0][2][2];
                //
                DG_t_rr[ig][jg][kg] = 0.5 * Delta_Gam[1][0][0];
                DG_t_rt[ig][jg][kg] = 0.5 * Delta_Gam[1][0][1];
                DG_t_rp[ig][jg][kg] = 0.5 * Delta_Gam[1][0][2];
                DG_t_tt[ig][jg][kg] = 0.5 * Delta_Gam[1][1][1];
                DG_t_tp[ig][jg][kg] = 0.5 * Delta_Gam[1][1][2];
                DG_t_pp[ig][jg][kg] = 0.5 * Delta_Gam[1][2][2];
                //
                DG_p_rr[ig][jg][kg] = 0.5 * Delta_Gam[2][0][0];
                DG_p_rt[ig][jg][kg] = 0.5 * Delta_Gam[2][0][1];
                DG_p_rp[ig][jg][kg] = 0.5 * Delta_Gam[2][0][2];
                DG_p_tt[ig][jg][kg] = 0.5 * Delta_Gam[2][1][1];
                DG_p_tp[ig][jg][kg] = 0.5 * Delta_Gam[2][1][2];
                DG_p_pp[ig][jg][kg] = 0.5 * Delta_Gam[2][2][2];
                //
                // now compute contractions 
                //
                vect DG;
                for (int i = 0; i < nn; i++) {
                    DG[i] = 0.0;
                    for (int j = 0; j < nn; j++)
                        for (int k = 0; k < nn; k++)
                            DG[i] += gup[j][k] * Delta_Gam[i][j][k];
                }
                DG_r[ig][jg][kg] = 0.5 * DG[0];
                DG_t[ig][jg][kg] = 0.5 * DG[1] * grid->r(ig);
                DG_p[ig][jg][kg] = 0.5 * DG[2] * grid->r(ig) * grid->sintheta(jg);
            };
};

//====================================================================================
//
// Now compute Ricci tensor, following eq. (19) in Brown (2009) (arXiv:0902.3652v2)
//
//====================================================================================
void curvature::Compute_Ricci(state* c) {
#ifdef _VERBOSE_
    cout << " Entering Compute_Ricci..." << endl;
#endif
    //
    // first compute auxiliary quantity X = exp(-2 phi)
    //
    //
#pragma omp parallel for collapse(3)
    for (int ig = N_g; ig < N_r - N_g; ig++) {
        for (int jg = N_g; jg < N_t - N_g; jg++) {
            for (int kg = N_g; kg < N_p - N_g; kg++) {
                double rl = grid->r(ig);
                double r2 = rl * rl;
                double stl = grid->sintheta(jg);
                double st2 = stl * stl;
                double ctl = grid->costheta(jg);
                //
                // Start with term g^{ij} \partial_i D_j g_{lm} "by hand"
                //
                // NOTE: we're assuming that the second derivatives of g_{lm} commute, i.e. D_i D_j g_{lm} = D_j D_i g_{lm}.  This is justified since
                // the Riemann tensor associated with the flat background metric vanishes, and there is no torsion (see, e.g., eq. (3.114) in Carroll).
                //
                R_rr[ig][jg][kg] =
                    // g^{rr} \partial_r D_r g_{rr}
                    gup_rr(ig, jg, kg) * (c->h_rr.ddr(ig, jg, kg)) +
                    // g^{rt} \partial_r D_t g_{rr}
                    gup_rt(ig, jg, kg) * (c->h_rr.drdtheta(ig, jg, kg) - 2.0 * c->h_rt.dr(ig, jg, kg)) * 2.0 +
                    // g^{rp} \partial_r D_p g_{rr}
                    gup_rp(ig, jg, kg) * (c->h_rr.drdphi(ig, jg, kg) - 2.0 * stl * c->h_rp.dr(ig, jg, kg)) * 2.0 +
                    // g^{tt} \partial_t D_t g_{rr}
                    gup_tt(ig, jg, kg) * (c->h_rr.ddtheta(ig, jg, kg) - 2.0 * c->h_rt.dtheta(ig, jg, kg)) +
                    // g^{tp} \partial_t D_p g_{rr}
                    gup_tp(ig, jg, kg) * (c->h_rr.dthetadphi(ig, jg, kg) - 2.0 * ctl * c->h_rp(ig, jg, kg)
                        - 2.0 * stl * c->h_rp.dtheta(ig, jg, kg)) * 2.0 +
                    // g^{pp} \partial_p D_p g_{rr}
                    gup_pp(ig, jg, kg) * (c->h_rr.ddphi(ig, jg, kg) - 2.0 * stl * c->h_rp.dphi(ig, jg, kg));
                //
                //
                R_rt[ig][jg][kg] =
                    // g^{rr} \partial_r D_r g_{rt}
                    gup_rr(ig, jg, kg) * (rl * c->h_rt.ddr(ig, jg, kg) + c->h_rt.dr(ig, jg, kg)) +
                    // g^{rt} \partial_r D_t g_{rt}
                    gup_rt(ig, jg, kg) * (rl * (c->h_rt.drdtheta(ig, jg, kg) + c->h_rr.dr(ig, jg, kg) - c->h_tt.dr(ig, jg, kg)) +
                        (c->h_rt.dtheta(ig, jg, kg) + c->h_rr(ig, jg, kg) - c->h_tt(ig, jg, kg))) * 2.0 +
                    // g^{rp} \partial_r D_p g_{rt}
                    gup_rp(ig, jg, kg) * (rl * (c->h_rt.drdphi(ig, jg, kg) - ctl * c->h_rp.dr(ig, jg, kg) - stl * c->h_tp.dr(ig, jg, kg)) +
                        (c->h_rt.dphi(ig, jg, kg) - ctl * c->h_rp(ig, jg, kg) - stl * c->h_tp(ig, jg, kg))) * 2.0 +
                    // g^{tt} \partial_t D_t g_{rt}
                    gup_tt(ig, jg, kg) * (rl * (c->h_rt.ddtheta(ig, jg, kg) + c->h_rr.dtheta(ig, jg, kg) - c->h_tt.dtheta(ig, jg, kg))) +
                    // g^{tp} \partial_t D_p g_{rt}
                    gup_tp(ig, jg, kg) * (rl * (c->h_rt.dthetadphi(ig, jg, kg) + stl * c->h_rp(ig, jg, kg) - ctl * c->h_rp.dtheta(ig, jg, kg) -
                        ctl * c->h_tp(ig, jg, kg) - stl * c->h_tp.dtheta(ig, jg, kg))) * 2.0 +
                    // g^{pp} \partial_p D_p g_{rt}
                    gup_pp(ig, jg, kg) * (rl * (c->h_rt.ddphi(ig, jg, kg) - ctl * c->h_rp.dphi(ig, jg, kg) - stl * c->h_tp.dphi(ig, jg, kg)));
                //
                //
                R_rp[ig][jg][kg] =
                    // g^{rr} \partial_r D_r g_{rp}
                    gup_rr(ig, jg, kg) * (rl * stl * c->h_rp.ddr(ig, jg, kg) + stl * c->h_rp.dr(ig, jg, kg)) +
                    // g^{rt} \partial_r D_t g_{rp}
                    gup_rt(ig, jg, kg) * (rl * stl * (c->h_rp.drdtheta(ig, jg, kg) - c->h_tp.dr(ig, jg, kg)) +
                        stl * (c->h_rp.dtheta(ig, jg, kg) - c->h_tp(ig, jg, kg))) * 2.0 +
                    // g^{rp} \partial_r D_p g_{rp}
                    gup_rp(ig, jg, kg) * (rl * stl * (c->h_rp.drdphi(ig, jg, kg) + stl * c->h_rr.dr(ig, jg, kg) +
                        ctl * c->h_rt.dr(ig, jg, kg) - stl * c->h_pp.dr(ig, jg, kg)) +
                        stl * (c->h_rp.dphi(ig, jg, kg) + stl * c->h_rr(ig, jg, kg) +
                            ctl * c->h_rt(ig, jg, kg) - stl * c->h_pp(ig, jg, kg))) * 2.0 +
                    // g^{tt} \partial_t D_t g_{rp}
                    gup_tt(ig, jg, kg) * (rl * stl * (c->h_rp.ddtheta(ig, jg, kg) - c->h_tp.dtheta(ig, jg, kg)) +
                        rl * ctl * (c->h_rp.dtheta(ig, jg, kg) - c->h_tp(ig, jg, kg))) +
                    // g^{tp} \partial_t D_p g_{rp}
                    gup_tp(ig, jg, kg) * (rl * stl * (c->h_rp.dthetadphi(ig, jg, kg) + ctl * c->h_rr(ig, jg, kg) +
                        stl * c->h_rr.dtheta(ig, jg, kg) - stl * c->h_rt(ig, jg, kg) +
                        ctl * c->h_rt.dtheta(ig, jg, kg) - ctl * c->h_pp(ig, jg, kg) -
                        stl * c->h_pp.dtheta(ig, jg, kg)) +
                        rl * ctl * (c->h_rp.dphi(ig, jg, kg) + stl * c->h_rr(ig, jg, kg) +
                            ctl * c->h_rt(ig, jg, kg) - stl * c->h_pp(ig, jg, kg))) * 2.0 +
                    // g^{pp} \partial_p D_p g_{rp}
                    gup_pp(ig, jg, kg) * (rl * stl * (c->h_rp.ddphi(ig, jg, kg) + stl * c->h_rr.dphi(ig, jg, kg) +
                        ctl * c->h_rt.dphi(ig, jg, kg) - stl * c->h_pp.dphi(ig, jg, kg)));
                //
                //
                R_tt[ig][jg][kg] =
                    // g^{rr} \partial_r D_r g_{tt}
                    gup_rr(ig, jg, kg) * (r2 * c->h_tt.ddr(ig, jg, kg) + 2.0 * rl * c->h_tt.dr(ig, jg, kg)) +
                    // g^{rt} \partial_r D_t g_{tt}
                    gup_rt(ig, jg, kg) * (r2 * (c->h_tt.drdtheta(ig, jg, kg) + 2.0 * c->h_rt.dr(ig, jg, kg)) +
                        2.0 * rl * (c->h_tt.dtheta(ig, jg, kg) + 2.0 * c->h_rt(ig, jg, kg))) * 2.0 +
                    // g^{rp} \partial_r D_p g_{tt}
                    gup_rp(ig, jg, kg) * (r2 * (c->h_tt.drdphi(ig, jg, kg) - 2.0 * ctl * c->h_tp.dr(ig, jg, kg)) +
                        2.0 * rl * (c->h_tt.dphi(ig, jg, kg) - 2.0 * ctl * c->h_tp(ig, jg, kg))) * 2.0 +
                    // g^{tt} \partial_t D_t g_{tt}
                    gup_tt(ig, jg, kg) * (r2 * (c->h_tt.ddtheta(ig, jg, kg) + 2.0 * c->h_rt.dtheta(ig, jg, kg))) +
                    // g^{tp} \partial_t D_p g_{tt}
                    gup_tp(ig, jg, kg) * (r2 * (c->h_tt.dthetadphi(ig, jg, kg) + 2.0 * stl * c->h_tp(ig, jg, kg)
                        - 2.0 * ctl * c->h_tp.dtheta(ig, jg, kg))) * 2.0 +
                    // g^{pp} \partial_p D_p g_{tt}
                    gup_pp(ig, jg, kg) * (r2 * (c->h_tt.ddphi(ig, jg, kg) - 2.0 * ctl * c->h_tp.dphi(ig, jg, kg)));
                //
                //
                R_tp[ig][jg][kg] =
                    // g^{rr} \partial_r D_r g_{tp}
                    gup_rr(ig, jg, kg) * (r2 * stl * c->h_tp.ddr(ig, jg, kg) + 2.0 * rl * stl * c->h_tp.dr(ig, jg, kg)) +
                    // g^{rt} \partial_r D_t g_{tp}
                    gup_rt(ig, jg, kg) * (r2 * stl * (c->h_tp.drdtheta(ig, jg, kg) + c->h_rp.dr(ig, jg, kg)) +
                        2.0 * rl * stl * (c->h_tp.dtheta(ig, jg, kg) + c->h_rp(ig, jg, kg))) * 2.0 +
                    // g^{rp} \partial_r D_p g_{tp}
                    gup_rp(ig, jg, kg) * (r2 * stl * (c->h_tp.drdphi(ig, jg, kg) + stl * c->h_rt.dr(ig, jg, kg) +
                        ctl * c->h_tt.dr(ig, jg, kg) - ctl * c->h_pp.dr(ig, jg, kg)) +
                        2.0 * rl * stl * (c->h_tp.dphi(ig, jg, kg) + stl * c->h_rt(ig, jg, kg) +
                            ctl * c->h_tt(ig, jg, kg) - ctl * c->h_pp(ig, jg, kg))) * 2.0 +
                    // g^{tt} \partial_t D_t g_{tp}
                    gup_tt(ig, jg, kg) * (r2 * stl * (c->h_tp.ddtheta(ig, jg, kg) + c->h_rp.dtheta(ig, jg, kg)) +
                        r2 * ctl * (c->h_tp.dtheta(ig, jg, kg) + c->h_rp(ig, jg, kg))) +
                    // g^{tp} \partial_t D_p g_{tp}
                    gup_tp(ig, jg, kg) * (r2 * stl * (c->h_tp.dthetadphi(ig, jg, kg)
                        + ctl * c->h_rt(ig, jg, kg) + stl * c->h_rt.dtheta(ig, jg, kg)
                        - stl * c->h_tt(ig, jg, kg) + ctl * c->h_tt.dtheta(ig, jg, kg)
                        + stl * c->h_pp(ig, jg, kg) - ctl * c->h_pp.dtheta(ig, jg, kg)) +
                        r2 * ctl * (c->h_tp.dphi(ig, jg, kg) + stl * c->h_rt(ig, jg, kg) +
                            ctl * c->h_tt(ig, jg, kg) - ctl * c->h_pp(ig, jg, kg))) * 2.0 +
                    // g^{pp} \partial_p D_p g_{tp}
                    gup_pp(ig, jg, kg) * (r2 * stl * (c->h_tp.ddphi(ig, jg, kg) + stl * c->h_rt.dphi(ig, jg, kg) +
                        ctl * c->h_tt.dphi(ig, jg, kg) - ctl * c->h_pp.dphi(ig, jg, kg)));
                //
                //
                R_pp[ig][jg][kg] =
                    // g^{rr} \partial_r D_r g_{pp}
                    gup_rr(ig, jg, kg) * (r2 * st2 * c->h_pp.ddr(ig, jg, kg) + 2.0 * rl * st2 * c->h_pp.dr(ig, jg, kg)) +
                    // g^{rt} \partial_r D_t g_{pp}
                    gup_rt(ig, jg, kg) * (r2 * st2 * c->h_pp.drdtheta(ig, jg, kg) + 2.0 * rl * st2 * c->h_pp.dtheta(ig, jg, kg)) * 2.0 +
                    // g^{rp} \partial_r D_p g_{pp}
                    gup_rp(ig, jg, kg) * (r2 * st2 * (c->h_pp.drdphi(ig, jg, kg)
                        + 2.0 * stl * c->h_rp.dr(ig, jg, kg) + 2.0 * ctl * c->h_tp.dr(ig, jg, kg)) +
                        2.0 * rl * st2 * (c->h_pp.dphi(ig, jg, kg)
                            + 2.0 * stl * c->h_rp(ig, jg, kg) + 2.0 * ctl * c->h_tp(ig, jg, kg))) * 2.0 +
                    // g^{tt} \partial_t D_t g_{pp}
                    gup_tt(ig, jg, kg) * (r2 * st2 * c->h_pp.ddtheta(ig, jg, kg) + 2.0 * r2 * stl * ctl * c->h_pp.dtheta(ig, jg, kg)) +
                    // g^{tp} \partial_t D_p g_{pp}
                    gup_tp(ig, jg, kg) * (r2 * st2 * (c->h_pp.dthetadphi(ig, jg, kg)
                        + 2.0 * ctl * c->h_rp(ig, jg, kg) + 2.0 * stl * c->h_rp.dtheta(ig, jg, kg)
                        - 2.0 * stl * c->h_tp(ig, jg, kg) + 2.0 * ctl * c->h_tp.dtheta(ig, jg, kg)) +
                        2.0 * r2 * stl * ctl * (c->h_pp.dphi(ig, jg, kg)
                            + 2.0 * stl * c->h_rp(ig, jg, kg)
                            + 2.0 * ctl * c->h_tp(ig, jg, kg))) * 2.0 +
                    // g^{pp} \partial_p D_p g_{pp}
                    gup_pp(ig, jg, kg) * (r2 * st2 * (c->h_pp.ddphi(ig, jg, kg)
                        + 2.0 * stl * c->h_rp.dphi(ig, jg, kg) + 2.0 * ctl * c->h_tp.dphi(ig, jg, kg)));
                /*
                Now add connection terms

                - g^{ij} D_k g_{lm} Gam^k_{ij} - g^{ij} D_j g_{km} Gam^k_{il} - g^{ij} D_j g_{lk} Gam^k_{im}
                   = - D_k g_{lm} Gam^k - g^{ij} D_j g_{km} Gam^k_{il} - g^{ij} D_j g_{lk} Gam^k_{im}

                where

                   Gam^k = g^{ij} Gam^k_{ij}

                    Recall that Gam is connection of flat metric in spherical polar coordinates here.
                */
                // first compute Conn^k = gup^{ij} Gam^k_{ij
                //
                const double Conn_r = -gup_tt(ig, jg, kg) * rl - gup_pp(ig, jg, kg) * rl * st2;
                const double Conn_t = 2.0 * gup_rt(ig, jg, kg) / rl - gup_pp(ig, jg, kg) * stl * ctl;
                const double Conn_p = 2.0 * gup_rp(ig, jg, kg) / rl + 2.0 * gup_tp(ig, jg, kg) * ctl / stl;
                //
                // ... then add connection terms
                /*
                R_lm[ig][jg][kg] -= Dr_e_lm(ig,jg,kg) * Conn_r + Dt_e_lm(ig,jg,kg) * Conn_t + Dp_e_lm(ig,jg,kg) * Conn_p +
                  gup_rr(ig,jg,kg) * ( Dr_e_rm(ig,jg,kg) * Conn^r_{rl} + Dr_e_tm(ig,jg,kg) * Conn^t_{rl} + Dr_e_pm(ig,jg,kg) * Conn^p_{rl} +
                               Dr_e_rl(ig,jg,kg) * Conn^r_{rm} + Dr_e_tl(ig,jg,kg) * Conn^t_{rm} + Dr_e_pl(ig,jg,kg) * Conn^p_{rm} ) +
                  gup_rt(ig,jg,kg) * ( Dt_e_rm(ig,jg,kg) * Conn^r_{rl} + Dt_e_tm(ig,jg,kg) * Conn^t_{rl} + Dt_e_pm(ig,jg,kg) * Conn^p_{rl} +
                               Dt_e_rl(ig,jg,kg) * Conn^r_{rm} + Dt_e_tl(ig,jg,kg) * Conn^t_{rm} + Dt_e_pl(ig,jg,kg) * Conn^p_{rm} ) * 2.0 +
                  gup_rp(ig,jg,kg) * ( Dp_e_rm(ig,jg,kg) * Conn^r_{rl} + Dp_e_tm(ig,jg,kg) * Conn^t_{rl} + Dp_e_pm(ig,jg,kg) * Conn^p_{rl} +
                               Dp_e_rl(ig,jg,kg) * Conn^r_{rm} + Dp_e_tl(ig,jg,kg) * Conn^t_{rm} + Dp_e_pl(ig,jg,kg) * Conn^p_{rm} ) * 2.0+
                  gup_tt(ig,jg,kg) * ( Dt_e_rm(ig,jg,kg) * Conn^r_{tl} + Dt_e_tm(ig,jg,kg) * Conn^t_{tl} + Dt_e_pm(ig,jg,kg) * Conn^p_{tl} +
                               Dt_e_rl(ig,jg,kg) * Conn^r_{tm} + Dt_e_tl(ig,jg,kg) * Conn^t_{tm} + Dt_e_pl(ig,jg,kg) * Conn^p_{tm} ) +
                  gup_tp(ig,jg,kg) * ( Dp_e_rm(ig,jg,kg) * Conn^r_{tl} + Dp_e_tm(ig,jg,kg) * Conn^t_{tl} + Dp_e_pm(ig,jg,kg) * Conn^p_{tl} +
                               Dp_e_rl(ig,jg,kg) * Conn^r_{tm} + Dp_e_tl(ig,jg,kg) * Conn^t_{tm} + Dp_e_pl(ig,jg,kg) * Conn^p_{tm} ) * 2.0+
                  gup_pp(ig,jg,kg) * ( Dp_e_rm(ig,jg,kg) * Conn^r_{pl} + Dp_e_tm(ig,jg,kg) * Conn^t_{pl} + Dp_e_pm(ig,jg,kg) * Conn^p_{pl} +
                               Dp_e_rl(ig,jg,kg) * Conn^r_{pm} + Dp_e_tl(ig,jg,kg) * Conn^t_{pm} + Dp_e_pl(ig,jg,kg) * Conn^p_{pm} );
                */

                R_rr[ig][jg][kg] -= Dr_e_rr(ig, jg, kg) * Conn_r + Dt_e_rr(ig, jg, kg) * Conn_t + Dp_e_rr(ig, jg, kg) * Conn_p +
                    2.0 * (gup_tt(ig, jg, kg) * Dt_e_rt(ig, jg, kg) +
                        2.0 * gup_tp(ig, jg, kg) * Dp_e_rt(ig, jg, kg) +
                        gup_pp(ig, jg, kg) * Dp_e_rp(ig, jg, kg)) / rl;
                //
                //
                R_rt[ig][jg][kg] -= Dr_e_rt(ig, jg, kg) * Conn_r + Dt_e_rt(ig, jg, kg) * Conn_t + Dp_e_rt(ig, jg, kg) * Conn_p +
                    gup_rr(ig, jg, kg) * (Dr_e_rt(ig, jg, kg) / rl) +
                    gup_rt(ig, jg, kg) * (Dt_e_rt(ig, jg, kg) / rl) * 2.0 +
                    gup_rp(ig, jg, kg) * (Dp_e_rt(ig, jg, kg) / rl) * 2.0 +
                    gup_tt(ig, jg, kg) * (Dt_e_tt(ig, jg, kg) / rl -
                        Dt_e_rr(ig, jg, kg) * rl) +
                    gup_tp(ig, jg, kg) * (Dp_e_tt(ig, jg, kg) / rl -
                        Dp_e_rr(ig, jg, kg) * rl) * 2.0 +
                    gup_pp(ig, jg, kg) * (Dp_e_tp(ig, jg, kg) / rl +
                        Dp_e_rp(ig, jg, kg) * ctl / stl);
                //
                //
                R_rp[ig][jg][kg] -= Dr_e_rp(ig, jg, kg) * Conn_r + Dt_e_rp(ig, jg, kg) * Conn_t + Dp_e_rp(ig, jg, kg) * Conn_p +
                    gup_rr(ig, jg, kg) * (Dr_e_rp(ig, jg, kg) / rl) +
                    gup_rt(ig, jg, kg) * (Dt_e_rp(ig, jg, kg) / rl) * 2.0 +
                    gup_rp(ig, jg, kg) * (Dp_e_rp(ig, jg, kg) / rl) * 2.0 +
                    gup_tt(ig, jg, kg) * (Dt_e_tp(ig, jg, kg) / rl +
                        Dt_e_rp(ig, jg, kg) * ctl / stl) +
                    gup_tp(ig, jg, kg) * (Dp_e_tp(ig, jg, kg) / rl +
                        Dp_e_rp(ig, jg, kg) * ctl / stl) * 2.0 +
                    gup_pp(ig, jg, kg) * (Dp_e_pp(ig, jg, kg) / rl -
                        Dp_e_rr(ig, jg, kg) * rl * st2 - Dp_e_rt(ig, jg, kg) * stl * ctl);
                //
                //
                R_tt[ig][jg][kg] -= Dr_e_tt(ig, jg, kg) * Conn_r + Dt_e_tt(ig, jg, kg) * Conn_t + Dp_e_tt(ig, jg, kg) * Conn_p +
                    2.0 * (gup_rr(ig, jg, kg) * Dr_e_tt(ig, jg, kg) / rl +
                        2.0 * gup_rt(ig, jg, kg) * Dt_e_tt(ig, jg, kg) / rl +
                        2.0 * gup_rp(ig, jg, kg) * Dp_e_tt(ig, jg, kg) / rl -
                        gup_tt(ig, jg, kg) * Dt_e_rt(ig, jg, kg) * rl -
                        2.0 * gup_tp(ig, jg, kg) * Dp_e_rt(ig, jg, kg) * rl +
                        gup_pp(ig, jg, kg) * Dp_e_tp(ig, jg, kg) * ctl / stl);
                //
                //
                R_tp[ig][jg][kg] -= Dr_e_tp(ig, jg, kg) * Conn_r + Dt_e_tp(ig, jg, kg) * Conn_t + Dp_e_tp(ig, jg, kg) * Conn_p +
                    gup_rr(ig, jg, kg) * (Dr_e_tp(ig, jg, kg) / rl +
                        Dr_e_tp(ig, jg, kg) / rl) +
                    gup_rt(ig, jg, kg) * (Dt_e_tp(ig, jg, kg) / rl +
                        Dt_e_tp(ig, jg, kg) / rl) * 2.0 +
                    gup_rp(ig, jg, kg) * (Dp_e_tp(ig, jg, kg) / rl +
                        Dp_e_tp(ig, jg, kg) / rl) * 2.0 +
                    gup_tt(ig, jg, kg) * (-Dt_e_rp(ig, jg, kg) * rl +
                        Dt_e_tp(ig, jg, kg) * ctl / stl) +
                    gup_tp(ig, jg, kg) * (-Dp_e_rp(ig, jg, kg) * rl +
                        Dp_e_tp(ig, jg, kg) * ctl / stl) * 2.0 +
                    gup_pp(ig, jg, kg) * (Dp_e_pp(ig, jg, kg) * ctl / stl -
                        Dp_e_rt(ig, jg, kg) * rl * st2 - Dp_e_tt(ig, jg, kg) * stl * ctl);
                //
                R_pp[ig][jg][kg] -= Dr_e_pp(ig, jg, kg) * Conn_r + Dt_e_pp(ig, jg, kg) * Conn_t + Dp_e_pp(ig, jg, kg) * Conn_p +
                    2.0 * (gup_rr(ig, jg, kg) * Dr_e_pp(ig, jg, kg) / rl +
                        2.0 * gup_rt(ig, jg, kg) * Dt_e_pp(ig, jg, kg) / rl +
                        2.0 * gup_rp(ig, jg, kg) * Dp_e_pp(ig, jg, kg) / rl +
                        gup_tt(ig, jg, kg) * Dt_e_pp(ig, jg, kg) * ctl / stl +
                        2.0 * gup_tp(ig, jg, kg) * Dp_e_pp(ig, jg, kg) * ctl / stl +
                        gup_pp(ig, jg, kg) * (-Dp_e_rp(ig, jg, kg) * rl * st2 -
                            Dp_e_tp(ig, jg, kg) * stl * ctl));
                //
                // now multiply R_{lm} with - (1/2)...
                //
                R_rr[ig][jg][kg] *= -0.5;
                R_rt[ig][jg][kg] *= -0.5;
                R_rp[ig][jg][kg] *= -0.5;
                R_tt[ig][jg][kg] *= -0.5;
                R_tp[ig][jg][kg] *= -0.5;
                R_pp[ig][jg][kg] *= -0.5;
                //
                // compute D_a \lambda^b and store in D_lam
                //
                const double Dr_lam_r = c->lam_r.dr(ig, jg, kg);
                const double Dt_lam_r = c->lam_r.dtheta(ig, jg, kg) - c->lam_t(ig, jg, kg);
                const double Dp_lam_r = c->lam_r.dphi(ig, jg, kg) - stl * c->lam_p(ig, jg, kg);
                //
                const double Dr_lam_t = c->lam_t.dr(ig, jg, kg) / rl;
                const double Dt_lam_t = (c->lam_t.dtheta(ig, jg, kg) + c->lam_r(ig, jg, kg)) / rl;
                const double Dp_lam_t = (c->lam_t.dphi(ig, jg, kg) - ctl * c->lam_p(ig, jg, kg)) / rl;
                //
                const double Dr_lam_p = c->lam_p.dr(ig, jg, kg) / (rl * stl);
                const double Dt_lam_p = c->lam_p.dtheta(ig, jg, kg) / (rl * stl);
                const double Dp_lam_p = (c->lam_p.dphi(ig, jg, kg) + stl * c->lam_r(ig, jg, kg) + ctl * c->lam_t(ig, jg, kg)) / (rl * stl);
                //
                tensor D_lam(Dr_lam_r, Dr_lam_t, Dr_lam_p,
                    Dt_lam_r, Dt_lam_t, Dt_lam_p,
                    Dp_lam_r, Dp_lam_t, Dp_lam_p);
                //
                // store conformal metric in g
                //
                tensor g(1.0 + c->h_rr(ig, jg, kg), rl * c->h_rt(ig, jg, kg), rl * stl * c->h_rp(ig, jg, kg),
                    r2 * (1.0 + c->h_tt(ig, jg, kg)), r2 * stl * c->h_tp(ig, jg, kg), r2 * st2 * (1.0 + c->h_pp(ig, jg, kg)));
                //
                // store Ricci tensor in R
                //
                tensor R(R_rr(ig, jg, kg), R_rt(ig, jg, kg), R_rp(ig, jg, kg),
                    R_tt(ig, jg, kg), R_tp(ig, jg, kg), R_pp(ig, jg, kg));
                //
                // store Delta Gamma^i_{jk} in Delta_Gam
                //
                rank3tens Delta_Gam(DG_r_rr(ig, jg, kg), DG_r_rt(ig, jg, kg), DG_r_rp(ig, jg, kg), DG_r_tt(ig, jg, kg), DG_r_tp(ig, jg, kg), DG_r_pp(ig, jg, kg),
                    DG_t_rr(ig, jg, kg), DG_t_rt(ig, jg, kg), DG_t_rp(ig, jg, kg), DG_t_tt(ig, jg, kg), DG_t_tp(ig, jg, kg), DG_t_pp(ig, jg, kg),
                    DG_p_rr(ig, jg, kg), DG_p_rt(ig, jg, kg), DG_p_rp(ig, jg, kg), DG_p_tt(ig, jg, kg), DG_p_tp(ig, jg, kg), DG_p_pp(ig, jg, kg));
                //
                // assign storage for Gamma's with first index lowered
                //
                rank3tens D_Gam_low;
                //
                // IMPORTANT NOTE: I'm using independent variable Lam here -- in Brown arXiv:0902.3652 this term is instead computed from 
                //     \bar g^{de} \Delta \Gamma^c_{de}
                // ANOTHER IMPORTANT NOTE: Now I am using \Delta Gamma^a here...
                // 
                vect Lam(c->lam_r(ig, jg, kg), c->lam_t(ig, jg, kg) / rl, c->lam_p(ig, jg, kg) / (rl * stl));
                vect DG(DG_r(ig, jg, kg), DG_t(ig, jg, kg) / rl, DG_p(ig, jg, kg) / (rl * stl));
                //
                // Assign upper metric to tensor
                //
                tensor gup(gup_rr[ig][jg][kg], gup_rt[ig][jg][kg], gup_rp[ig][jg][kg], gup_tt[ig][jg][kg], gup_tp[ig][jg][kg], gup_pp[ig][jg][kg]);
                //
                // first lower first index on Delta_Gam and store in D_Gam_low
                // 
                // loop over three free indices
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        for (int c = 0; c < 3; c++) {
                            D_Gam_low[a][b][c] = 0.0;
                            for (int d = 0; d < 3; d++) {
                                D_Gam_low[a][b][c] += g[a][d] * Delta_Gam[d][b][c];
                            }
                        }
                //
                // now loop over free indices...
                // 
                for (int a = 0; a < 3; a++)
                    for (int b = a; b < 3; b++) {
                        //
                        // loop over dummies
                        //
                        for (int c = 0; c < 3; c++) {
                            R[a][b] += 0.5 * (g[c][a] * D_lam[b][c] + g[c][b] * D_lam[a][c] + DG[c] * D_Gam_low[a][b][c] + DG[c] * D_Gam_low[b][a][c]);
                            // R[a][b] += 0.5 * ( g[c][a]*D_lam[b][c] + g[c][b]*D_lam[a][c] + Lam[c]*D_Gam_low[a][b][c] + Lam[c]*D_Gam_low[b][a][c] ); 
                            for (int d = 0; d < 3; d++) {
                                for (int e = 0; e < 3; e++) {
                                    R[a][b] += gup[c][d] * (Delta_Gam[e][c][a] * D_Gam_low[b][e][d] +
                                        Delta_Gam[e][c][b] * D_Gam_low[a][e][d] +
                                        Delta_Gam[e][a][c] * D_Gam_low[e][b][d]);
                                }
                            }
                        }
                    }
                R_rr[ig][jg][kg] = R[0][0];
                R_rt[ig][jg][kg] = R[0][1];
                R_rp[ig][jg][kg] = R[0][2];
                R_tt[ig][jg][kg] = R[1][1];
                R_tp[ig][jg][kg] = R[1][2];
                R_pp[ig][jg][kg] = R[2][2];
                //
                // one more thing: compute trace...
                //
                trace_R[ig][jg][kg] = 0.0;
                R[1][0] = R[0][1];
                R[2][0] = R[0][2];
                R[2][1] = R[1][2];
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        trace_R[ig][jg][kg] += gup[a][b] * R[a][b];
            }
        }
    }
};


