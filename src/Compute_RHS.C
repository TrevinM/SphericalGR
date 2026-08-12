#include "Manager.h"
#include "Container.h"

void Manager::Compute_RHS(double time) {
    if (cowling == 1)
        derivs->equals(0.0);
    else {
        //
        // Compute curvature from intermediate state inter
        //
        //    curve->Update(inter);
        //
        // Compute divergence of shift
        // 
        aux->DivShift(inter, curve);
        //
        // Compute time derivative of metric
        // 
        dot_metric(inter, time);
        //
        // Compute time derivative of extrinsic curvature
        // 
        dot_ext_curv(inter, time);
        //
        // Compute time derivative of connection coefficients (Lambda's)
        // 
        dot_connection(inter, time);
        //
        // Compute time derivatives of lapse and shift
        //
        // NOTE: call dot_shift *after* dot_Gamma, since time
        // derivatives of Gamma are needed for some gauge conditions.
        // 
        Container::slicing->dot_lapse(inter, derivs);
        Container::gauge->dot_shift(inter, derivs);
        //
        // Outer boundaries:
        //
        if (char_OB == 0)
            for (int i = 0; i < derivs->N_fcts; i++) {
                derivs->fct_list[i]->derivs_outerboundary(inter->fct_list[i]);
            }
    }
    derivs->fill_ghosts();
};
//
//==============================================================
// time derivative of metric 
//==============================================================
//
void Manager::dot_metric(state* s, double time) {
    //
    // first compute Lie derivative of metric
    // 
    aux->Lie_Metric(s->shift_r, s->shift_t, s->shift_p,
        s->h_rr, s->h_rt, s->h_rp, s->h_tt, s->h_tp, s->h_pp);
    //
    // compute trace of A_ij (CHECK: is this necessary?)
    // 
    curve->Compute_Trace(s);
    //
    // now compute time derivatives
    // 
    const double twothirds = 2.0 / 3.0;
#pragma omp parallel for collapse(3)  
    for (int i = N_g; i < N_r - N_g; i++)
        for (int j = N_g; j < N_t - N_g; j++)
            for (int k = N_g; k < N_p - N_g; k++) {
                const double lapse_l = s->lapse(i, j, k);
                const double A = curve->trace_A(i, j, k);
                derivs->h_rr[i][j][k] = lapse_l * (-2.0 * s->a_rr(i, j, k) + twothirds * (1.0 + s->h_rr(i, j, k)) * A) +
                    aux->Lt_rr(i, j, k) + eta_KO * s->h_rr.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * (1.0 + s->h_rr(i, j, k));
                derivs->h_rt[i][j][k] = lapse_l * (-2.0 * s->a_rt(i, j, k) + twothirds * (s->h_rt(i, j, k)) * A) +
                    aux->Lt_rt(i, j, k) + eta_KO * s->h_rt.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * (s->h_rt(i, j, k));
                derivs->h_rp[i][j][k] = lapse_l * (-2.0 * s->a_rp(i, j, k) + twothirds * (s->h_rp(i, j, k)) * A) +
                    aux->Lt_rp(i, j, k) + eta_KO * s->h_rp.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * (s->h_rp(i, j, k));
                derivs->h_tt[i][j][k] = lapse_l * (-2.0 * s->a_tt(i, j, k) + twothirds * (1.0 + s->h_tt(i, j, k)) * A) +
                    aux->Lt_tt(i, j, k) + eta_KO * s->h_tt.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * (1.0 + s->h_tt(i, j, k));
                derivs->h_tp[i][j][k] = lapse_l * (-2.0 * s->a_tp(i, j, k) + twothirds * (s->h_tp(i, j, k)) * A) +
                    aux->Lt_tp(i, j, k) + eta_KO * s->h_tp.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * (s->h_tp(i, j, k));
                derivs->h_pp[i][j][k] = lapse_l * (-2.0 * s->a_pp(i, j, k) + twothirds * (1.0 + s->h_pp(i, j, k)) * A) +
                    aux->Lt_pp(i, j, k) + eta_KO * s->h_pp.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * (1.0 + s->h_pp(i, j, k));
                //
                derivs->phi[i][j][k] = -lapse_l * (s->K(i, j, k) + 2.0 * s->Theta(i, j, k)) / 6.0 +
                    aux->Lie(s->shift_r, s->shift_t, s->shift_p, s->phi, i, j, k) +
                    eta_KO * s->phi.KO(i, j, k) + sigma * aux->div_shift(i, j, k) / 6.0
                    // add cosmological term
                    - Container::cosmology->Hubble(time) / 2.0;
            }
};
//
//==============================================================
// time derivative of extrinsic curvature
//==============================================================
//
void Manager::dot_ext_curv(state* s, double time) {
    const double lambda = Container::cosmology->Lambda();
    const double a_friedmann = Container::cosmology->a(time);
    //
    // first compute Lie derivative of extrinsic curvature...
    // 
    aux->Lie_ExCurv(s->shift_r, s->shift_t, s->shift_p,
        s->a_rr, s->a_rt, s->a_rp, s->a_tt, s->a_tp, s->a_pp);
    //
    // ... and auxiliary quantity X
    // 
    aux->Compute_X(s, a_friedmann);
    //
    // now compute time derivatives
    //
    const double twothirds = 2.0 / 3.0;
#pragma omp parallel for collapse(3)
    for (int i = N_g; i < N_r - N_g; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                const double rl = Container::grid->r(i);
                const double r2 = rl * rl;
                const double stl = Container::grid->sintheta(j);
                const double st2 = stl * stl;
                const double ctl = Container::grid->costheta(j);
                //
                // First the terms that, in PIRK scheme, were "well-behaved" L3 terms...
                    //
                // Assign upper metric to tensor
                //
                tensor gup(curve->gup_rr[i][j][k], curve->gup_rt[i][j][k], curve->gup_rp[i][j][k],
                    curve->gup_tt[i][j][k], curve->gup_tp[i][j][k], curve->gup_pp[i][j][k]);
                //
                // store extrinsic curvature (lower indices) in A (not rescaled)
                //
                tensor A(s->a_rr(i, j, k), rl * s->a_rt(i, j, k), rl * stl * s->a_rp(i, j, k),
                    r2 * s->a_tt(i, j, k), r2 * stl * s->a_tp(i, j, k), r2 * st2 * s->a_pp(i, j, k));
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
                const double lapse_l = s->lapse(i, j, k);
                const double K_l = s->K(i, j, k);
                const double T_l = s->Theta(i, j, k);
                const double K_plus_T = K_l + kappa_ric * T_l;
                derivs->a_rr[i][j][k] = aux->Lt_rr(i, j, k) + lapse_l * (s->a_rr(i, j, k) * K_plus_T - 2.0 * AdotA[0][0]);
                derivs->a_rt[i][j][k] = aux->Lt_rt(i, j, k) + lapse_l * (s->a_rt(i, j, k) * K_plus_T - 2.0 * AdotA[0][1] / rl);
                derivs->a_rp[i][j][k] = aux->Lt_rp(i, j, k) + lapse_l * (s->a_rp(i, j, k) * K_plus_T - 2.0 * AdotA[0][2] / (rl * stl));
                derivs->a_tt[i][j][k] = aux->Lt_tt(i, j, k) + lapse_l * (s->a_tt(i, j, k) * K_plus_T - 2.0 * AdotA[1][1] / r2);
                derivs->a_tp[i][j][k] = aux->Lt_tp(i, j, k) + lapse_l * (s->a_tp(i, j, k) * K_plus_T - 2.0 * AdotA[1][2] / (r2 * stl));
                derivs->a_pp[i][j][k] = aux->Lt_pp(i, j, k) + lapse_l * (s->a_pp(i, j, k) * K_plus_T - 2.0 * AdotA[2][2] / (r2 * st2));
                const double kappa_1 = kappa_11 + kappa_12 * T_l * T_l;
                derivs->K[i][j][k] = aux->Lie(s->shift_r, s->shift_t, s->shift_p, s->K, i, j, k)
                    // note cosmological term
                    + lapse_l * (K_plus_T * K_plus_T / 3.0 + A2_trace
                        + 4.0 * PI * (Container::matter->adm_sources->rho_ADM(i, j, k) + Container::matter->adm_sources->trace_S(i, j, k)) - lambda
                        + kappa_1 * (1.0 - kappa_2) * s->Theta(i, j, k));
                // NOTE: absence of "Riccatti" term: replaced K_plus_2T with K_l;  also need to fix lambda term if present.
                if (lambda != 0.0) cerr << " FIX lambda term in Theta equation! " << endl;
                if (z4) derivs->Theta[i][j][k] = aux->Lie(s->shift_r, s->shift_t, s->shift_p, s->Theta, i, j, k)
                    + lapse_l * (K_plus_T * K_plus_T / 3.0 - 0.5 * A2_trace - 8.0 * PI * Container::matter->adm_sources->rho_ADM(i, j, k)
                        - kappa_1 * (2.0 + kappa_2) * s->Theta(i, j, k));

                //
                // add Kreiss-Oliger terms as well as divergence of shift for Lagrangian formulation
                //
                derivs->a_rr[i][j][k] += eta_KO * s->a_rr.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * s->a_rr(i, j, k);
                derivs->a_rt[i][j][k] += eta_KO * s->a_rt.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * s->a_rt(i, j, k);
                derivs->a_rp[i][j][k] += eta_KO * s->a_rp.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * s->a_rp(i, j, k);
                derivs->a_tt[i][j][k] += eta_KO * s->a_tt.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * s->a_tt(i, j, k);
                derivs->a_tp[i][j][k] += eta_KO * s->a_tp.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * s->a_tp(i, j, k);
                derivs->a_pp[i][j][k] += eta_KO * s->a_pp.KO(i, j, k) - sigma * twothirds * aux->div_shift(i, j, k) * s->a_pp(i, j, k);
                derivs->K[i][j][k] += eta_KO * s->K.KO(i, j, k);
                if (z4) derivs->Theta[i][j][k] += eta_KO * s->Theta.KO(i, j, k);
                //
                // Now the terms that, in PIRK scheme, were "stiff" L2 terms...
                // 
                // store Delta Gamma^i_{jk} in Delta_Gam
                rank3tens Delta_Gam(curve->DG_r_rr(i, j, k), curve->DG_r_rt(i, j, k), curve->DG_r_rp(i, j, k),
                    curve->DG_r_tt(i, j, k), curve->DG_r_tp(i, j, k), curve->DG_r_pp(i, j, k),
                    curve->DG_t_rr(i, j, k), curve->DG_t_rt(i, j, k), curve->DG_t_rp(i, j, k),
                    curve->DG_t_tt(i, j, k), curve->DG_t_tp(i, j, k), curve->DG_t_pp(i, j, k),
                    curve->DG_p_rr(i, j, k), curve->DG_p_rt(i, j, k), curve->DG_p_rp(i, j, k),
                    curve->DG_p_tt(i, j, k), curve->DG_p_tp(i, j, k), curve->DG_p_pp(i, j, k));
                // store derivatives in tensors
                vect D_lapse(s->lapse.dr(i, j, k), s->lapse.dtheta(i, j, k), s->lapse.dphi(i, j, k));
                vect D_phi(s->phi.dr(i, j, k), s->phi.dtheta(i, j, k), s->phi.dphi(i, j, k));
                tensor DD_lapse(s->lapse.ddr(i, j, k), s->lapse.drdtheta(i, j, k), s->lapse.drdphi(i, j, k),
                    s->lapse.ddtheta(i, j, k), s->lapse.dthetadphi(i, j, k), s->lapse.ddphi(i, j, k));
                tensor DD_phi(s->phi.ddr(i, j, k), s->phi.drdtheta(i, j, k), s->phi.drdphi(i, j, k),
                    s->phi.ddtheta(i, j, k), s->phi.dthetadphi(i, j, k), s->phi.ddphi(i, j, k));
                //
                // alternative: compute derivative of phi from X
                //
                vect D_X(aux->X.dr(i, j, k), aux->X.dtheta(i, j, k), aux->X.dphi(i, j, k));
                tensor DD_X(aux->X.ddr(i, j, k), aux->X.drdtheta(i, j, k), aux->X.drdphi(i, j, k),
                    aux->X.ddtheta(i, j, k), aux->X.dthetadphi(i, j, k), aux->X.ddphi(i, j, k));
                const double factor1 = a_friedmann * exp(2.0 * s->phi(i, j, k)) / 2.0;
                const double factor2 = a_friedmann * a_friedmann * exp(4.0 * s->phi(i, j, k)) / 2.0;
                //
                D_phi[0] = -factor1 * D_X[0];
                D_phi[1] = -factor1 * D_X[1];
                D_phi[2] = -factor1 * D_X[2];
                //
                DD_phi[0][0] = -factor1 * DD_X[0][0] + factor2 * D_X[0] * D_X[0];
                DD_phi[0][1] = -factor1 * DD_X[0][1] + factor2 * D_X[0] * D_X[1];
                DD_phi[1][0] = -factor1 * DD_X[1][0] + factor2 * D_X[1] * D_X[0];
                DD_phi[0][2] = -factor1 * DD_X[0][2] + factor2 * D_X[0] * D_X[2];
                DD_phi[2][0] = -factor1 * DD_X[2][0] + factor2 * D_X[2] * D_X[0];
                DD_phi[1][1] = -factor1 * DD_X[1][1] + factor2 * D_X[1] * D_X[1];
                DD_phi[1][2] = -factor1 * DD_X[1][2] + factor2 * D_X[1] * D_X[2];
                DD_phi[2][1] = -factor1 * DD_X[2][1] + factor2 * D_X[2] * D_X[1];
                DD_phi[2][2] = -factor1 * DD_X[2][2] + factor2 * D_X[2] * D_X[2];
                //
                // add Delta connection terms to second derivatives
                //
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        for (int c = 0; c < 3; c++) {
                            DD_lapse[a][b] -= D_lapse[c] * Delta_Gam[c][a][b];
                            DD_phi[a][b] -= D_phi[c] * Delta_Gam[c][a][b];
                        }
                // add background connection terms to second derivatives
                DD_lapse[0][1] -= D_lapse[1] / rl;
                DD_lapse[1][0] -= D_lapse[1] / rl;
                DD_lapse[0][2] -= D_lapse[2] / rl;
                DD_lapse[2][0] -= D_lapse[2] / rl;
                DD_lapse[1][1] += D_lapse[0] * rl;
                DD_lapse[1][2] -= D_lapse[2] * ctl / stl;
                DD_lapse[2][1] -= D_lapse[2] * ctl / stl;
                DD_lapse[2][2] += D_lapse[0] * rl * st2 + D_lapse[1] * stl * ctl;

                DD_phi[0][1] -= D_phi[1] / rl;
                DD_phi[1][0] -= D_phi[1] / rl;
                DD_phi[0][2] -= D_phi[2] / rl;
                DD_phi[2][0] -= D_phi[2] / rl;
                DD_phi[1][1] += D_phi[0] * rl;
                DD_phi[1][2] -= D_phi[2] * ctl / stl;
                DD_phi[2][1] -= D_phi[2] * ctl / stl;
                DD_phi[2][2] += D_phi[0] * rl * st2 + D_phi[1] * stl * ctl;
                // define two more tensors as well as metric for later purposes...
                tensor R(curve->R_rr(i, j, k), curve->R_rt(i, j, k), curve->R_rp(i, j, k),
                    curve->R_tt(i, j, k), curve->R_tp(i, j, k), curve->R_pp(i, j, k));
                // NOTE: S_ij is assumed *not* to be rescaled (like R_ij)
                tensor Stress(Container::matter->adm_sources->S_rr(i, j, k), Container::matter->adm_sources->S_rt(i, j, k), Container::matter->adm_sources->S_rp(i, j, k),
                    Container::matter->adm_sources->S_tt(i, j, k), Container::matter->adm_sources->S_tp(i, j, k), Container::matter->adm_sources->S_pp(i, j, k));
                tensor L2;
                tensor g(1.0 + s->h_rr(i, j, k), rl * s->h_rt(i, j, k), rl * stl * s->h_rp(i, j, k),
                    r2 * (1.0 + s->h_tt(i, j, k)), r2 * stl * s->h_tp(i, j, k), r2 * st2 * (1.0 + s->h_pp(i, j, k)));
                //
                // see eq. (22b) in B09
                // 
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        L2[a][b] = -2.0 * lapse_l * DD_phi[a][b] + 4.0 * lapse_l * D_phi[a] * D_phi[b] +
                        2.0 * (D_lapse[a] * D_phi[b] + D_phi[a] * D_lapse[b]) - DD_lapse[a][b] +
                        lapse_l * (R[a][b] - 8.0 * PI * Stress[a][b]);
                //
                // Now remove trace of L2 (NOTE: could speed this up by doing it by hand)
                // 
                L2.remove_trace(gup, g);
                //
                // also need Laplace operator of lapse (could use .trace tensor function...)
                // simultaneously compute Laplace operator phi (needed for Hamiltonian constraint)
                //
                double L2_K = 0.0;
                double Lap_phi = 0.0;
                double DphiDphi = 0.0;
                //
                // see eq. (22d) in B09
                //
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        L2_K += gup[a][b] * (DD_lapse[a][b] + 2.0 * D_lapse[a] * D_phi[b]);
                        Lap_phi += gup[a][b] * DD_phi[a][b];
                        DphiDphi += gup[a][b] * D_phi[a] * D_phi[b];
                    }
                // finalize...  Note cosmological expansion term
                double e4p = exp(-4.0 * s->phi(i, j, k)) / (a_friedmann * a_friedmann);
                derivs->a_rr[i][j][k] += e4p * L2[0][0];
                derivs->a_rt[i][j][k] += e4p * L2[0][1] / rl;
                derivs->a_rp[i][j][k] += e4p * L2[0][2] / (rl * stl);
                derivs->a_tt[i][j][k] += e4p * L2[1][1] / r2;
                derivs->a_tp[i][j][k] += e4p * L2[1][2] / (r2 * stl);
                derivs->a_pp[i][j][k] += e4p * L2[2][2] / (r2 * st2);
                derivs->K[i][j][k] -= e4p * L2_K;
                if (z4) derivs->Theta[i][j][k] += 0.5 * lapse_l * e4p * (curve->trace_R(i, j, k) - 8.0 * DphiDphi -
                    8.0 * Lap_phi);
            }
        }
    }
    //  if (z4 != 1) derivs->Theta.equals(0.0); 
};
//
//==============================================================
// time derivative of connection coefficients
//==============================================================
//
void Manager::dot_connection(state* s, double time) {
    const double a_friedmann = Container::cosmology->a(time);
    //
    // Comment added 5/2/22: careful in this routine with rescaling...
    //
    //
    // first compute Lie derivative... - returns *rescaled* quantities
    //
    aux->Lie(s->shift_r, s->shift_t, s->shift_p,
        s->lam_r, s->lam_t, s->lam_p,
        derivs->lam_r, derivs->lam_t, derivs->lam_p);
    //
    // and make sure det has its ghosts filled
    // 
    curve->det.fill_ghosts();
    //
    // then compute time derivatives
    //
    const double twothirds = 2.0 / 3.0;
#pragma omp parallel for collapse(3)  
    for (int i = N_g; i < N_r - N_g; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                const double rl = Container::grid->r(i);
                const double r2 = rl * rl;
                const double stl = Container::grid->sintheta(j);
                const double st2 = stl * stl;
                const double rst = rl * stl;
                const double ctl = Container::grid->costheta(j);
                //
                // "L3" terms in old PIRK scheme
                //
                // All terms *rescaled*...
                //
                derivs->lam_r[i][j][k] += eta_KO * s->lam_r.KO(i, j, k) + sigma * twothirds * curve->DG_r(i, j, k) * aux->div_shift(i, j, k);
                derivs->lam_t[i][j][k] += eta_KO * s->lam_t.KO(i, j, k) + sigma * twothirds * curve->DG_t(i, j, k) * aux->div_shift(i, j, k);
                derivs->lam_p[i][j][k] += eta_KO * s->lam_p.KO(i, j, k) + sigma * twothirds * curve->DG_p(i, j, k) * aux->div_shift(i, j, k);
                //
                // add constraint damping Z4c terms
                // 
                const double lapse_l = s->lapse(i, j, k);
                const double CFC_r = s->lam_r(i, j, k) - curve->DG_r(i, j, k);
                const double CFC_t = s->lam_t(i, j, k) - curve->DG_t(i, j, k);
                const double CFC_p = s->lam_p(i, j, k) - curve->DG_p(i, j, k);
                derivs->lam_r[i][j][k] -= 2.0 * lapse_l * (kappa_11 + kappa_12 * CFC_r * CFC_r) * CFC_r;
                derivs->lam_t[i][j][k] -= 2.0 * lapse_l * (kappa_11 + kappa_12 * CFC_t * CFC_t) * CFC_t;
                derivs->lam_p[i][j][k] -= 2.0 * lapse_l * (kappa_11 + kappa_12 * CFC_p * CFC_p) * CFC_p;
                //	derivs->lam_r[i][j][k] -= 2.0*lapse_l * kappa_1 * ( s->lam_r(i,j,k) - curve->DG_r(i,j,k) ); 
                //	derivs->lam_t[i][j][k] -= 2.0*lapse_l * kappa_1 * ( s->lam_t(i,j,k) - curve->DG_t(i,j,k) ); 
                //	derivs->lam_p[i][j][k] -= 2.0*lapse_l * kappa_1 * ( s->lam_p(i,j,k) - curve->DG_p(i,j,k) ); 
                //
                // now "L2" terms...  Will initially compute them *not* rescaled, and will therefore
                // store them in vector L2 first.
                //
                vect L2;
                // compute first derivatives of shift in terms of code variables  
                // Notation: D_i \beta^j = beta_j_i  (where D_i is *flat* covariant derivative)
                //
                const double beta_r_r = s->shift_r.dr(i, j, k);
                const double beta_r_t = s->shift_r.dtheta(i, j, k) - s->shift_t(i, j, k);
                const double beta_r_p = s->shift_r.dphi(i, j, k) - stl * s->shift_p(i, j, k);
                const double beta_t_r = (s->shift_t.dr(i, j, k)) / rl;
                const double beta_t_t = (s->shift_t.dtheta(i, j, k) + s->shift_r(i, j, k)) / rl;
                const double beta_t_p = (s->shift_t.dphi(i, j, k) - ctl * s->shift_p(i, j, k)) / rl;
                const double beta_p_r = (s->shift_p.dr(i, j, k)) / (rst);
                const double beta_p_t = (s->shift_p.dtheta(i, j, k)) / (rst);
                const double beta_p_p = (s->shift_p.dphi(i, j, k) + stl * s->shift_r(i, j, k) + ctl * s->shift_t(i, j, k)) /
                    (rst);
                //
                // define partial_k D_i beta^j = p_k_beta_j_i
                //
                const double p_r_beta_r_r = s->shift_r.ddr(i, j, k);
                const double p_r_beta_r_t = s->shift_r.drdtheta(i, j, k) - s->shift_t.dr(i, j, k);
                const double p_r_beta_r_p = s->shift_r.drdphi(i, j, k) - stl * s->shift_p.dr(i, j, k);
                const double p_r_beta_t_r = s->shift_t.ddr(i, j, k) / rl - beta_t_r / rl;
                const double p_r_beta_t_t = (s->shift_t.drdtheta(i, j, k) + s->shift_r.dr(i, j, k)) / rl - beta_t_t / rl;
                const double p_r_beta_t_p = (s->shift_t.drdphi(i, j, k) - ctl * s->shift_p.dr(i, j, k)) / rl - beta_t_p / rl;
                const double p_r_beta_p_r = s->shift_p.ddr(i, j, k) / rst - beta_p_r / rl;
                const double p_r_beta_p_t = s->shift_p.drdtheta(i, j, k) / rst - beta_p_t / rl;
                const double p_r_beta_p_p = (s->shift_p.drdphi(i, j, k) + stl * s->shift_r.dr(i, j, k) + ctl * s->shift_t.dr(i, j, k)) / rst
                    - beta_p_p / rl;
                //
                const double p_t_beta_r_r = s->shift_r.drdtheta(i, j, k);
                const double p_t_beta_r_t = s->shift_r.ddtheta(i, j, k) - s->shift_t.dtheta(i, j, k);
                const double p_t_beta_r_p = s->shift_r.dthetadphi(i, j, k) - stl * s->shift_p.dtheta(i, j, k) - ctl * s->shift_p(i, j, k);
                const double p_t_beta_t_r = s->shift_t.drdtheta(i, j, k) / rl;
                const double p_t_beta_t_t = (s->shift_t.ddtheta(i, j, k) + s->shift_r.dtheta(i, j, k)) / rl;
                const double p_t_beta_t_p = (s->shift_t.dthetadphi(i, j, k) - ctl * s->shift_p.dtheta(i, j, k)
                    + stl * s->shift_p(i, j, k)) / rl;
                const double p_t_beta_p_r = s->shift_p.drdtheta(i, j, k) / rst - beta_p_r * ctl / stl;
                const double p_t_beta_p_t = s->shift_p.ddtheta(i, j, k) / rst - beta_p_t * ctl / stl;
                const double p_t_beta_p_p = (s->shift_p.dthetadphi(i, j, k) + stl * s->shift_r.dtheta(i, j, k) + ctl * s->shift_r(i, j, k)
                    + ctl * s->shift_t.dtheta(i, j, k) - stl * s->shift_t(i, j, k)) / rst
                    - beta_p_p * ctl / stl;
                //
                const double p_p_beta_r_r = s->shift_r.drdphi(i, j, k);
                const double p_p_beta_r_t = s->shift_r.dthetadphi(i, j, k) - s->shift_t.dphi(i, j, k);
                const double p_p_beta_r_p = s->shift_r.ddphi(i, j, k) - stl * s->shift_p.dphi(i, j, k);
                const double p_p_beta_t_r = s->shift_t.drdphi(i, j, k) / rl;
                const double p_p_beta_t_t = (s->shift_t.dthetadphi(i, j, k) + s->shift_r.dphi(i, j, k)) / rl;
                const double p_p_beta_t_p = (s->shift_t.ddphi(i, j, k) - ctl * s->shift_p.dphi(i, j, k)) / rl;
                const double p_p_beta_p_r = s->shift_p.drdphi(i, j, k) / rst;
                const double p_p_beta_p_t = s->shift_p.dthetadphi(i, j, k) / rst;
                const double p_p_beta_p_p = (s->shift_p.ddphi(i, j, k) + stl * s->shift_r.dphi(i, j, k) + ctl * s->shift_t.dphi(i, j, k)) / rst;
                //
                // Now compute gup^{lm} \partial_l D_m \beta^i
                //
                //derivs->lam_r[i][j][k] +=
                L2[0] =
                    //   gup^rr \partial_r D_r \beta^r + 
                    curve->gup_rr(i, j, k) * p_r_beta_r_r +
                    // 2 gup^rt \partial_r D_t \beta^r +   // rather: symmetrize term!!!
                    curve->gup_rt(i, j, k) * (p_r_beta_r_t + p_t_beta_r_r) +
                    // 2 gup^rp \partial_r D_p \beta^r + 
                    curve->gup_rp(i, j, k) * (p_r_beta_r_p + p_p_beta_r_r) +
                    //   gup^tt \partial_t D_t \beta^r + 
                    curve->gup_tt(i, j, k) * p_t_beta_r_t +
                    // 2 gup^tp \partial_t D_p \beta^r + 
                    curve->gup_tp(i, j, k) * (p_t_beta_r_p + p_p_beta_r_t) +
                    //   gup^pp \partial_p D_p \beta^r 
                    curve->gup_pp(i, j, k) * p_p_beta_r_p;
                //
                //	derivs->lam_t[i][j][k] += 
                L2[1] =
                    curve->gup_rr(i, j, k) * p_r_beta_t_r +
                    curve->gup_rt(i, j, k) * (p_r_beta_t_t + p_t_beta_t_r) +
                    curve->gup_rp(i, j, k) * (p_r_beta_t_p + p_p_beta_t_r) +
                    curve->gup_tt(i, j, k) * p_t_beta_t_t +
                    curve->gup_tp(i, j, k) * (p_t_beta_t_p + p_p_beta_t_t) +
                    curve->gup_pp(i, j, k) * p_p_beta_t_p;
                //
                // 	derivs->lam_p[i][j][k] +=
                L2[2] =
                    curve->gup_rr(i, j, k) * p_r_beta_p_r +
                    curve->gup_rt(i, j, k) * (p_r_beta_p_t + p_t_beta_p_r) +
                    curve->gup_rp(i, j, k) * (p_r_beta_p_p + p_p_beta_p_r) +
                    curve->gup_tt(i, j, k) * p_t_beta_p_t +
                    curve->gup_tp(i, j, k) * (p_t_beta_p_p + p_p_beta_p_t) +
                    curve->gup_pp(i, j, k) * p_p_beta_p_p;
                // 
                // Now add connection terms
                // 	
                // first compute Conn^k = gup^{ij} Gam^k_{ij
                //
                // (where Gam^k_{ij} are Christoffel symbols associated with *flat* metric)
                const double Conn_r = -curve->gup_tt(i, j, k) * rl - curve->gup_pp(i, j, k) * rl * st2;
                const double Conn_t = 2.0 * curve->gup_rt(i, j, k) / rl - curve->gup_pp(i, j, k) * stl * ctl;
                const double Conn_p = 2.0 * curve->gup_rp(i, j, k) / rl + 2.0 * curve->gup_tp(i, j, k) * ctl / stl;
                //
                // Fixed version - old version commented out
                //
                const double Gh_r_tt = -rl;
                const double Gh_r_pp = -rl * st2;
                const double Gh_t_pp = -stl * ctl;
                const double Gh_t_rt = 1.0 / rl;
                const double Gh_p_rp = 1.0 / rl;
                const double Gh_p_tp = ctl / stl;
                // derivs->lam_r[i][j][k]
                L2[0] += -beta_r_r * Conn_r - beta_r_t * Conn_t - beta_r_p * Conn_p +
                    curve->gup_rt(i, j, k) * (Gh_r_tt * beta_t_r) +
                    curve->gup_rp(i, j, k) * (Gh_r_pp * beta_p_r) +
                    curve->gup_tt(i, j, k) * (Gh_r_tt * beta_t_t) +
                    curve->gup_tp(i, j, k) * (Gh_r_tt * beta_t_p + Gh_r_pp * beta_p_t) +
                    curve->gup_pp(i, j, k) * (Gh_r_pp * beta_p_p);
                // curve->gup_tt(i,j,k) * ( beta_t_t * (- rl) ) +  
                // curve->gup_tp(i,j,k) * ( beta_t_p * (- rl) ) * 2.0 + 
                // curve->gup_pp(i,j,k) * ( beta_p_p * (- rl * st2 ) );  
              //
              //	derivs->lam_t[i][j][k]
                L2[1] += -beta_t_r * Conn_r - beta_t_t * Conn_t - beta_t_p * Conn_p +
                    curve->gup_rr(i, j, k) * (Gh_t_rt * beta_t_r) +
                    curve->gup_rt(i, j, k) * (Gh_t_rt * beta_t_t + Gh_t_rt * beta_r_r) +
                    curve->gup_rp(i, j, k) * (Gh_t_rt * beta_t_p + Gh_t_pp * beta_p_r) +
                    curve->gup_tt(i, j, k) * (Gh_t_rt * beta_r_t) +
                    curve->gup_tp(i, j, k) * (Gh_t_rt * beta_r_p + Gh_t_pp * beta_p_t) +
                    curve->gup_pp(i, j, k) * (Gh_t_pp * beta_p_p);
                //
                // derivs->lam_p[i][j][k]
                L2[2] += -beta_p_r * Conn_r - beta_p_t * Conn_t - beta_p_p * Conn_p +
                    curve->gup_rr(i, j, k) * (Gh_p_rp * beta_p_r) +
                    curve->gup_rt(i, j, k) * (Gh_p_rp * beta_p_t + Gh_p_tp * beta_p_r) +
                    curve->gup_rp(i, j, k) * (Gh_p_rp * beta_p_p + Gh_p_rp * beta_r_r + Gh_p_tp * beta_t_r) +
                    curve->gup_tt(i, j, k) * (Gh_p_tp * beta_p_t) +
                    curve->gup_tp(i, j, k) * (Gh_p_rp * beta_r_t + Gh_p_tp * beta_t_t + Gh_p_tp * beta_p_p) +
                    curve->gup_pp(i, j, k) * (Gh_p_rp * beta_r_p + Gh_p_tp * beta_t_p);
                //
                // rescale L2 terms and add to derivs->lam_i - ACTUALLY: will do later!
                //
                //	derivs->lam_r[i][j][k] += L2[0]
                //	derivs->lam_t[i][j][k] += L2[1] * rl;
                //	derivs->lam_p[i][j][k] += L2[2] * rl * stl;	
                // 
                // L2[i] now contains term gup^{lm} D_l D_m beta^i  (where D_l is *flat* covariant derivative)
                //
                //====================================================================================================
                //
                // Now work on term 1/3 \bar D^i \bar D_j beta^j.
                //
                // Will start with partial derivatives of divergence of shift, \partial_l \bar D_j beta^j
                //   (take partial derivatives of expression for divergence of shift in routine DivShift()...
                //
                vect partialDivShift(0.0, 0.0, 0.0);
                if (sigma == 1) {  // don't bother doing this for Eulerian formulation...
                    partialDivShift[0] =
                        s->shift_r.ddr(i, j, k) +
                        s->shift_t.drdtheta(i, j, k) / rl - s->shift_t.dtheta(i, j, k) / r2 +
                        s->shift_p.drdphi(i, j, k) / rst - s->shift_p.dphi(i, j, k) / (r2 * stl) +
                        2.0 * s->shift_r.dr(i, j, k) / rl - 2.0 * s->shift_r(i, j, k) / r2 +
                        s->shift_t.dr(i, j, k) * ctl / rst - s->shift_t(i, j, k) * ctl / (r2 * stl) +
                        (s->shift_r.dr(i, j, k) * curve->det.dr(i, j, k) + s->shift_r(i, j, k) * curve->det.ddr(i, j, k) +
                            s->shift_t.dr(i, j, k) * curve->det.dtheta(i, j, k) / rl + s->shift_t(i, j, k) * curve->det.drdtheta(i, j, k) / rl -
                            s->shift_t(i, j, k) * curve->det.dtheta(i, j, k) / r2 +
                            s->shift_p.dr(i, j, k) * curve->det.dphi(i, j, k) / rst + s->shift_p(i, j, k) * curve->det.drdphi(i, j, k) / rst
                            - s->shift_p(i, j, k) * curve->det.dphi(i, j, k) / (r2 * stl)) /
                        (2.0 * curve->det(i, j, k)) -
                        (s->shift_r(i, j, k) * curve->det.dr(i, j, k) + s->shift_t(i, j, k) * curve->det.dtheta(i, j, k) / rl + s->shift_p(i, j, k) * curve->det.dphi(i, j, k) / rst) /
                        (2.0 * curve->det(i, j, k) * curve->det(i, j, k)) * curve->det.dr(i, j, k);
                    //
                    partialDivShift[1] =
                        s->shift_r.drdtheta(i, j, k) +
                        s->shift_t.ddtheta(i, j, k) / rl +
                        s->shift_p.dthetadphi(i, j, k) / rst - s->shift_p.dphi(i, j, k) * ctl / (rl * st2) +
                        2.0 * s->shift_r.dtheta(i, j, k) / rl +
                        s->shift_t.dtheta(i, j, k) * ctl / rst - s->shift_t(i, j, k) / (rl * st2) +
                        (s->shift_r.dtheta(i, j, k) * curve->det.dr(i, j, k) + s->shift_r(i, j, k) * curve->det.drdtheta(i, j, k) +
                            s->shift_t.dtheta(i, j, k) * curve->det.dtheta(i, j, k) / rl + s->shift_t(i, j, k) * curve->det.ddtheta(i, j, k) / rl +
                            s->shift_p.dtheta(i, j, k) * curve->det.dphi(i, j, k) / rst + s->shift_p(i, j, k) * curve->det.dthetadphi(i, j, k) / rst -
                            s->shift_p(i, j, k) * curve->det.dphi(i, j, k) * ctl / (rl * st2)) / (2.0 * curve->det(i, j, k)) -
                        (s->shift_r(i, j, k) * curve->det.dr(i, j, k) + s->shift_t(i, j, k) * curve->det.dtheta(i, j, k) / rl + s->shift_p(i, j, k) * curve->det.dphi(i, j, k) / rst) /
                        (2.0 * curve->det(i, j, k) * curve->det(i, j, k)) * curve->det.dtheta(i, j, k);
                    //
                    partialDivShift[2] =
                        s->shift_r.drdphi(i, j, k) +
                        s->shift_t.dthetadphi(i, j, k) / rl +
                        s->shift_p.ddphi(i, j, k) / rst +
                        2.0 * s->shift_r.dphi(i, j, k) / rl + s->shift_t.dphi(i, j, k) * ctl / rst +
                        (s->shift_r.dphi(i, j, k) * curve->det.dr(i, j, k) + s->shift_r(i, j, k) * curve->det.drdphi(i, j, k) +
                            s->shift_t.dphi(i, j, k) * curve->det.dtheta(i, j, k) / rl + s->shift_t(i, j, k) * curve->det.dthetadphi(i, j, k) / rl +
                            s->shift_p.dphi(i, j, k) * curve->det.dphi(i, j, k) / rst + s->shift_p(i, j, k) * curve->det.ddphi(i, j, k) / rst) / (2.0 * curve->det(i, j, k)) -
                        (s->shift_r(i, j, k) * curve->det.dr(i, j, k) + s->shift_t(i, j, k) * curve->det.dtheta(i, j, k) / rl + s->shift_p(i, j, k) * curve->det.dphi(i, j, k) / rst) /
                        (2.0 * curve->det(i, j, k) * curve->det(i, j, k)) * curve->det.dphi(i, j, k);
                }
                //
                // Assign upper metric to tensor
                //
                tensor gup(curve->gup_rr[i][j][k], curve->gup_rt[i][j][k], curve->gup_rp[i][j][k],
                    curve->gup_tt[i][j][k], curve->gup_tp[i][j][k], curve->gup_pp[i][j][k]);
                //
                // raise indices of A_{ij}...
                //
                // store extrinsic curvature (lower indices) in A
                //
                tensor A(s->a_rr(i, j, k), rl * s->a_rt(i, j, k), rst * s->a_rp(i, j, k),
                    r2 * s->a_tt(i, j, k), r2 * stl * s->a_tp(i, j, k), r2 * st2 * s->a_pp(i, j, k));
                tensor Aup;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        Aup[a][b] = 0.0;
                        for (int c = 0; c < 3; c++)
                            for (int d = 0; d < 3; d++) {
                                Aup[a][b] += gup[a][c] * gup[b][d] * A[c][d];
                            }
                    }
                //
                // Define vectors with derivatives of lapse, phi and K
                //
                vect D_lapse(s->lapse.dr(i, j, k), s->lapse.dtheta(i, j, k), s->lapse.dphi(i, j, k));
                vect D_phi(s->phi.dr(i, j, k), s->phi.dtheta(i, j, k), s->phi.dphi(i, j, k));
                vect D_K(s->K.dr(i, j, k), s->K.dtheta(i, j, k), s->K.dphi(i, j, k));
                //
                // Define vector with Container::matter source term S_i (defined *downstairs*!), assumed *not* to be rescaled
                //
                vect S_matter(Container::matter->adm_sources->S_r(i, j, k), Container::matter->adm_sources->S_t(i, j, k), Container::matter->adm_sources->S_p(i, j, k));
                //
                // alternative: compute derivative of phi from X - note cosmological factor
                //
                vect D_X(aux->X.dr(i, j, k), aux->X.dtheta(i, j, k), aux->X.dphi(i, j, k));
                const double factor = -exp(2.0 * s->phi(i, j, k)) / 2.0 * a_friedmann;
                D_phi[0] = factor * D_X[0];
                D_phi[1] = factor * D_X[1];
                D_phi[2] = factor * D_X[2];
                //
                // store Delta Gamma^i_{jk} in Delta_Gam
                //
                rank3tens Delta_Gam(curve->DG_r_rr(i, j, k), curve->DG_r_rt(i, j, k), curve->DG_r_rp(i, j, k),
                    curve->DG_r_tt(i, j, k), curve->DG_r_tp(i, j, k), curve->DG_r_pp(i, j, k),
                    curve->DG_t_rr(i, j, k), curve->DG_t_rt(i, j, k), curve->DG_t_rp(i, j, k),
                    curve->DG_t_tt(i, j, k), curve->DG_t_tp(i, j, k), curve->DG_t_pp(i, j, k),
                    curve->DG_p_rr(i, j, k), curve->DG_p_rt(i, j, k), curve->DG_p_rp(i, j, k),
                    curve->DG_p_tt(i, j, k), curve->DG_p_tp(i, j, k), curve->DG_p_pp(i, j, k));
                //
                // Define and fill vector L2
                // 
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        L2[a] -= 2.0 * Aup[a][b] * (D_lapse[b] - 6.0 * lapse_l * D_phi[b])
                            + lapse_l * gup[a][b] * (4.0 / 3.0 * D_K[b] + 16.0 * PI * S_matter[b]);
                        L2[a] += sigma * gup[a][b] * partialDivShift[b] / 3.0;
                        for (int c = 0; c < 3; c++)
                            L2[a] += 2.0 * lapse_l * Aup[b][c] * Delta_Gam[a][b][c];
                    }
                }
                //	
                // store result in gridfunction, taking into account proper scaling
                //
                derivs->lam_r[i][j][k] += L2[0];
                derivs->lam_t[i][j][k] += L2[1] * rl;
                derivs->lam_p[i][j][k] += L2[2] * rl * stl;
                //
                // done...
                //
            }
        }
    }
};

