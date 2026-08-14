//================================================
//
// NOTE: In this version all hydro variables are properly rescaled.
//
// REFERENCE: Thierfelder, Bernuzzi & Bruegmann, arXiv:1104.4751 (TBB)
//      Duez, Liu, Shapiro & Stephens, arXiv:astro-ph/0503420 (DLSS)
//================================================
//
// NOTE: use state vectors for conservative variables and fluxes 
// as defined in (21) and (24) of TBB - i.e. *with* factors of determinant of 
// metric and lapse
//
// i.e.: gridfunction D, e.g., is D-component of q_D = gamma * D
//
//================================================
//
// File that includes routines for class Hydro
//
//================================================
#include "Grid.h"
#include "Hydro.h"
#include "tensors.h"
//================================================
//
// Read input from file "Hydro_Input"
//
//================================================
int Hydro::ReadInput(int& N_particles, double& R_max) {
    int error = 0;
    ifstream infile;
    infile.open("Hydro_Input");
    if (!infile) {
        cerr << " HYDRO: Can't open input file Hydro_Input " << endl;
        cerr << " HYDRO: Will use default values " << endl;
        error = 1;
        return error;
    }
    char buf[500], c;
    infile.get(buf, 500, '='); infile.get(c); infile >> slope_limiter_type;
    infile.get(buf, 500, '='); infile.get(c); infile >> riemann_type;
    infile.get(buf, 500, '='); infile.get(c); infile >> partial;
    infile.get(buf, 500, '='); infile.get(c); infile >> f_atm;
    infile.get(buf, 500, '='); infile.get(c); infile >> f_thr;
    infile.get(buf, 500, '='); infile.get(c); infile >> fix_hydro;
    infile.get(buf, 500, '='); infile.get(c); infile >> n_fixed;
    infile.get(buf, 500, '='); infile.get(c); infile >> N_particles;
    infile.get(buf, 500, '='); infile.get(c); infile >> R_max;
    if (infile.eof()) {
        cerr << " HYDRO: Error reading input file HYDRO_Input " << endl;
        cerr << " HYDRO: Will use default values " << endl;
        error = 2;
    }
    return error;
}
//================================================
//
// Read grid function from checkpoint file
//
//================================================
bool Hydro::read_gf(gf3d& gf, int timestep) {
    bool found = false;
    ifstream infile;
    stringstream filename;
    filename << gf.Name() << "_" << setfill('0') << setw(8)
        << timestep << ".cpt" << ends;
    infile.open(filename.str().c_str());
    // cout << " INDATA: looking for checkpoint file "
    //	 << filename.str().c_str() << endl;
    if (infile) {
        cout << " HYDRO: reading initial data from checkpoint file "
            << filename.str().c_str() << endl;
        //
        // NOTE: this logic has to match that in gridfunction.checkpoint!
        //
        double value = 0.0;
        for (int i = N_g; i < N_r; i++)
            for (int j = N_g; j < N_t - N_g; j++)
                for (int k = N_g; k < N_p - N_g; k++) {
                    infile >> value;
                    gf[i][j][k] = value;
                }
        found = true;
        infile.close();
    }
    return found;
};

//================================================
//
// Initialize fluid
//
// NOTE: for hydro, we handle checkpoints differently from other functions.  The point is that for
// hydro we usually initialize using "cold_eps" - i.e. we set T = 0 or s = const or something like
// that.  That will no longer hold at later times.  We therefore initialize from the conserved
// quantities for checkpoints.
//
//================================================
void Hydro::Initialize(state* s, curvature* c, diagnostics* d) {
    aux->hydro_errors.equals(0.0);
    if (checkpoint->ReadFromChkpt()) {
        //
        // initialization for reading from checkpoint: read in conserved quantities
        // 
        int timestep = checkpoint->TimeStepStart();
        bool found_chkpt_file;
        found_chkpt_file = read_gf(last->D, timestep);
        found_chkpt_file = read_gf(last->tau, timestep);
        found_chkpt_file = read_gf(last->S_r, timestep);
        found_chkpt_file = read_gf(last->S_t, timestep);
        found_chkpt_file = read_gf(last->S_p, timestep);
        last->D.fill_ghosts();
        last->tau.fill_ghosts();
        last->S_r.fill_ghosts();
        last->S_t.fill_ghosts();
        last->S_p.fill_ghosts();
        Recovery(last, s, c);
        // rho_atm = f_atm * aux->rho_0.max();
        rho_atm = f_atm;
        rho_thr = f_thr * rho_atm;
    } else {
        //
        // or otherwise provide analytical values for primitive variables
        // 
        indata->Initialize_Matter(aux->rho_0, aux->v_r, aux->v_t, aux->v_p);
        // rho_atm = f_atm * aux->rho_0.max();
        rho_atm = f_atm;
        rho_thr = f_thr * rho_atm;
        for (int i = N_g; i < N_r; i++)
            for (int j = N_g; j < N_t - N_g; j++)
                for (int k = N_g; k < N_p - N_g; k++) {
                    if (aux->rho_0(i, j, k) < rho_thr)
                        aux->rho_0[i][j][k] = rho_atm;
                    if (i == N_g && j == N_g && k == N_g)
                        cout << " P = " << aux->p[i][j][k] << endl;
                    aux->eps[i][j][k] = eos->cold_eps(aux->rho_0(i, j, k));
                    aux->p[i][j][k] = eos->P(aux->rho_0(i, j, k), aux->eps(i, j, k));
                    if (i == N_g && j == N_g && k == N_g)
                        cout << " P = " << aux->p[i][j][k] << endl;
                }
        //
        // Now compute conserved quantities on initial slice
        //
        Compute_Conserved_Vars(last, s, c);
    }
    //
    // "Recover" primitive variables - deals with atmosphere for t=0
    //
    Recovery(last, s, c);
    //
    //  Test_Recovery(s, c);
    //
    // Finally: add field variables to particle dump functions:
    //
    if (particles != NULL) {
        particles->add_to_dump_list("Particle_Fct_List", s);
        particles->add_to_dump_list("Particle_Fct_List", c);
        particles->add_to_dump_list("Particle_Fct_List", d);
    }
};

//================================================
//
// Compute right-hand sides...
//
//================================================
void Hydro::Compute_RHS(state* s, curvature* c, double time) {
    //
    // inverse Cowling approximation: hydro-without-hydro...
    //
    if (cowling == 2)
        derivs->equals(0.0);
    else {
        //
        // Otherwise: compute primitives from matter state inter
        //
        Recovery(inter, s, c);
        //
        // Construct fluxes; involves
        // - reconstruction at cell interfaces (in Reconstruct_i)
        // - computation of fluxes at interfaces (in Compute_Fluxes)
        // - solution of local Riemann problem (in riemann-solver->flux)
        //
        Construct_Fluxes(s, c);
        //
        // then compute time derivatives of conserved quantities
        //
        dot_conserved_vars(inter, s, c);
        //
        // Outer boundaries (CHECK: may have to fix...)
        //
        // for (int i = 0; i < derivs->N_fcts; i++)
        //   derivs->fct_list[i]->derivs_outerboundary(inter->fct_list[i]);
        //
        // finally: compute right-hand sides for particles
        // 
        if (particles != NULL)
            particles->Compute_RHS(&s->lapse, &aux->W,
                &s->shift_r, &s->shift_t, &s->shift_p,
                &aux->v_r, &aux->v_t, &aux->v_p);
    }
};

//================================================
// Compute conserved variables (multiplied with determinant of metric, 
//    i.e. members of q vector!) from primitive ones
//
// Either on entire grid, or locally.  
// First global version which calls local version
//
//================================================
int Hydro::Compute_Conserved_Vars(hydro_state* m, state* s, curvature* c) {
    int error = 0;
#pragma omp parallel for collapse(3)
    for (int i = N_g; i < N_r; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                double r_l = grid->r(i);
                double sintheta_l = grid->sintheta(j);
                double h_rr_l = s->h_rr(i, j, k);
                double h_rt_l = s->h_rt(i, j, k);
                double h_rp_l = s->h_rp(i, j, k);
                double h_tt_l = s->h_tt(i, j, k);
                double h_tp_l = s->h_tp(i, j, k);
                double h_pp_l = s->h_pp(i, j, k);
                double phi_l = s->phi(i, j, k);
                double lapse_l = s->lapse(i, j, k);
                double det_l = c->det(i, j, k);
                double D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l;
                const double rho_l = aux->rho_0(i, j, k);
                const double eps_l = aux->eps(i, j, k);
                double p_l = eos->P(rho_l, eps_l);
                error = Compute_Conserved_Vars(rho_l, eps_l, p_l,
                    aux->v_r(i, j, k), aux->v_t(i, j, k), aux->v_p(i, j, k),
                    h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l,
                    phi_l, det_l, r_l, sintheta_l,
                    D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l);
                aux->hydro_errors[i][j][k] += error;
                m->D[i][j][k] = D_l;
                m->S_r[i][j][k] = S_r_l;
                m->S_t[i][j][k] = S_t_l;
                m->S_p[i][j][k] = S_p_l;
                m->tau[i][j][k] = tau_l;
            }
        }
    }
    m->D.fill_ghosts();
    m->S_r.fill_ghosts();
    m->S_t.fill_ghosts();
    m->S_p.fill_ghosts();
    m->tau.fill_ghosts();
    return error;
};
//
//================================================
// Now local version
//================================================
//
int Hydro::Compute_Conserved_Vars(double rho_0, double eps, double p,
    double v_r, double v_t, double v_p,
    double h_rr, double h_rt, double h_rp,
    double h_tt, double h_tp, double h_pp,
    double phi, double det,
    double r, double sintheta,
    double& D,
    double& S_r, double& S_t, double& S_p,
    double& tau, double& W, double& v2) {
    int error = 0;
    // CHECK: assuming that v^i *is* rescaled!
#ifdef SR
    const double det_3D = 1.0;
    const double e4p = 1.0;
    tensor g_conf(1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
#else
    const double det_3D = exp(6.0 * phi) * sqrt(det);
    const double e4p = exp(4.0 * phi);
    tensor g_conf(1.0 + h_rr, h_rt, h_rp, 1.0 + h_tt, h_tp, 1.0 + h_pp);
#endif
    vect v_up(v_r, v_t, v_p);
    vect v_low;
    double vdotv = 0.0;
    // lower index on v and compute dot product
    for (int a = 0; a < 3; a++) {
        v_low[a] = 0.0;
        for (int b = 0; b < 3; b++) {
            v_low[a] += e4p * g_conf[a][b] * v_up[b];
            vdotv += e4p * g_conf[a][b] * v_up[a] * v_up[b];
        }
    }
    v2 = vdotv;
    if (v2 >= 1.0) {
        // cout << " Trouble in Compute_Conserved_Vars: v2 = " << v2 << " r = " << r << endl;
        // exit(0);
        error = 1;
        v2 = vdotv = 0.95;
    }
    //
    // compute Lorentz factor W
    //
    W = 1.0 / sqrt(1.0 - vdotv);
    //
    // compute matter density D
    //
    if (partial)
        D = det_3D * r * r * sintheta * W * rho_0;
    else
        D = det_3D * W * rho_0;
    //
    // compute momentum density S_i (indices down-stairs!) (again: *rescaled!*)
    //
    const double W2 = W * W;
    const double hl = eos->h(rho_0, eps);
    S_r = det_3D * W2 * rho_0 * hl * v_low[0];
    S_t = det_3D * W2 * rho_0 * hl * v_low[1];
    S_p = det_3D * W2 * rho_0 * hl * v_low[2];
    //
    // compute energy density tau
    //
    //  tau = W2 * rho_0 * hl - eos->P(rho_0,eps) - D;
    if (partial)
        tau = det_3D * r * r * sintheta * (W2 * rho_0 * hl - p) - D;
    else
        tau = det_3D * (W2 * rho_0 * hl - p) - D;
    return error;
};
//================================================
// Compute primitive variables from conserved ones:
//
// pass in current time-level of conserved variables, computes
// primitive variables
//================================================
int Hydro::Recovery(hydro_state* m, state* s, curvature* c) {
    int error = 0;
    //
    // set up atmosphere treatment...
    // 
    //  const double rho_atm = f_atm * aux->rho_0.max();
    //  const double rho_atm = f_atm * 0.4;
    //  const double rho_thr = f_thr * rho_atm;
    const double eps_atm = eos->cold_eps(rho_atm);
    const double p_atm = eos->P(rho_atm, eps_atm);
    //
    // now go to each grid-point
    //
    int i_test = -65;
    int j_test = 3;
    int k_test = 3;
    // NOTE TOLERANCE! - consider adjusting...
    const double tol = 1.e-8;
#pragma omp parallel for collapse(3)
    for (int i = N_g; i < N_r; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                double r_l = grid->r(i);
                double sintheta_l = grid->sintheta(j);
#ifdef SR
                const double det_3D = 1.0;
#else
                const double det_3D = exp(6.0 * s->phi(i, j, k)) * sqrt(c->det(i, j, k));
#endif
                //
                // define local conserved quantities
                // RECALL: gridfunctions are components for "q" vector; see eq. (50) in
                //         Montero, Baumgarte & Mueller
                //
                double Dl = m->D[i][j][k] / det_3D;
                double taul = m->tau[i][j][k] / det_3D;
                if (partial) {
                    Dl /= r_l * r_l * sintheta_l;
                    taul /= r_l * r_l * sintheta_l;
                }
                //
                // store S_i and metric in tensors to compute S2
                //
                vect S(m->S_r[i][j][k] / det_3D, m->S_t[i][j][k] / det_3D, m->S_p[i][j][k] / det_3D);
                vect S_up;
#ifdef SR
                const double em4p = 1.0;
                tensor g_conf(1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
#else
                const double em4p = exp(-4.0 * s->phi(i, j, k));
                tensor g_conf(1.0 + s->h_rr(i, j, k), s->h_rt(i, j, k), s->h_rp(i, j, k),
                    1.0 + s->h_tt(i, j, k), s->h_tp(i, j, k), 1.0 + s->h_pp(i, j, k));
#endif
                tensor gup = g_conf.inverse();

                double S2 = 0.0;
                for (int a = 0; a < 3; a++) {
                    S_up[a] = 0.0;
                    for (int b = 0; b < 3; b++) {
                        S_up[a] += em4p * gup[a][b] * S[b];
                        S2 += em4p * gup[a][b] * S[a] * S[b];
                    }
                }
                //
                // Now check condition (34) in Etienne et.al., PRD 77, 084002 (2008) - 
                // CHECK: put in a "safety" factor of 0.99 or so??
                //
                const double S2_max = taul * (taul + 2.0 * Dl);
                // CHECK!!!
                // if (S2 > 0.98 * S2_max) {
                if (S2 > 1.0 * S2_max) {
                    aux->hydro_errors[i][j][k] += 10;
                    // cout << " Oooops at r = " << r_l << " ratio = " << S2 / S2_max << endl;
                    // check: change back to 0.98??
                    const double factor = sqrt(fabs(0.999 * S2_max / S2));
                    m->S_r[i][j][k] *= factor;
                    m->S_t[i][j][k] *= factor;
                    m->S_p[i][j][k] *= factor;
                    for (int a = 0; a < 3; a++) {
                        S[a] *= factor;
                        S_up[a] *= factor;
                    }
                    S2 *= factor * factor;
                    //  rho_0_guess = 0.0;  // will send this into atmosphere...
                }
                double p_guess = aux->p[i][j][k];  // take old pressure value as starting value for iteration
                //
                // compute "guess" value for rho_0 to decide whether it's atmosphere or not
                //
                double rho_0_guess = 0.0;
                double W_guess = 1.0;
                const double tauDp_guess = taul + p_guess + Dl;
                if (tauDp_guess > 0.0) {
                    W_guess = tauDp_guess / sqrt(fabs(tauDp_guess * tauDp_guess - S2));
                    rho_0_guess = Dl / W_guess;
                }
                //
                // Now, if rho_0_guess too small...
                //
                if (rho_0_guess <= rho_thr || taul < 0.0 || Dl < 0.0) {
                    //===========================================================
                    // ATMOSPHERE
                    //===========================================================
                    aux->rho_0[i][j][k] = rho_atm;
                    aux->v_r[i][j][k] = 0.0;
                    aux->v_t[i][j][k] = 0.0;
                    aux->v_p[i][j][k] = 0.0;
                    aux->eps[i][j][k] = eps_atm;
                    aux->p[i][j][k] = p_atm;
                    //
                    // update conserved variables accordingly
                    //
                    double h_rr_l = s->h_rr(i, j, k);
                    double h_rt_l = s->h_rt(i, j, k);
                    double h_rp_l = s->h_rp(i, j, k);
                    double h_tt_l = s->h_tt(i, j, k);
                    double h_tp_l = s->h_tp(i, j, k);
                    double h_pp_l = s->h_pp(i, j, k);
                    double phi_l = s->phi(i, j, k);
                    double lapse_l = s->lapse(i, j, k);
                    double det_l = c->det(i, j, k);
                    //
                    double D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l;
                    Compute_Conserved_Vars(rho_atm, eps_atm, p_atm, 0.0, 0.0, 0.0,
                        h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l,
                        phi_l, det_l,
                        r_l, sintheta_l,
                        D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l);
                    m->D[i][j][k] = D_l;
                    m->S_r[i][j][k] = S_r_l;
                    m->S_t[i][j][k] = S_t_l;
                    m->S_p[i][j][k] = S_p_l;
                    m->tau[i][j][k] = tau_l;
                    aux->W[i][j][k] = W_l;
                } else {
                    //===========================================================
                    // Not in atmosphere - now use Newton-Raphson to iterate...
                    //===========================================================
                    double dfdp = 0.0;
                    double f_guess = f(p_guess, taul, Dl, S2, dfdp);
                    int it = 0;
                    if (i == i_test && j == j_test && k == k_test)
                        cout << " Recovery: p_guess = " << setprecision(16) << p_guess << " after "
                        << it << " iterations, f_guess = " << f_guess << endl;
                    int it_max = 50;
                    // NOTE: tolerance for relative error in p...
                    while (it < it_max && abs(f_guess) > tol * p_guess) {
                        it++;
                        // 
                        // improved value of p:
                        //
                        p_guess -= f_guess / dfdp;
                        f_guess = f(p_guess, taul, Dl, S2, dfdp);
                        if (i == i_test && j == j_test && k == k_test)
                            cout << " Recovery: p_guess = " << p_guess << " after "
                            << it << " iterations, f_guess = " << f_guess
                            << " tau*(tau + 2 D) - S2 = " << taul * (taul + 2.0 * Dl) - S2 << endl;
                    }
                    // NOTE: again...
                    if (abs(f_guess) < tol * p_guess) {
                        //
                        // iteration converged... compute primitive variables
                        //
                        const double taupD = taul + p_guess + Dl;
                        const double taupD2 = taupD * taupD;
                        const double Wl = taupD / sqrt(taupD2 - S2);
                        // if (Wl != 1.0) cout << " at r = " << r_l << " Wl = " << setprecision(16)<< Wl 
                        //  			<< " taul = " << taul << " p_guess = " << p_guess 
                        //  			<< " Dl = " << Dl << " i = " << i << endl;
                        aux->W[i][j][k] = Wl;
                        aux->rho_0[i][j][k] = Dl / Wl;
                        aux->eps[i][j][k] = (sqrt(taupD2 - S2) - Wl * p_guess - Dl) / Dl;
                        // if (aux->eps(i,j,k) < 0.0) 
                        if (i == i_test && j == j_test && k == k_test)
                            cout << " In Recovery : taupD = " << taupD << " p_guess = "
                            << p_guess << " P " << eos->P(aux->rho_0(i, j, k), aux->eps(i, j, k)) << " Dl = " << Dl
                            << " S2 = " << S2 << " i " << i << " rho_0 "
                            << aux->rho_0(i, j, k) << " eps " << aux->eps(i, j, k) << endl;
                        aux->v_r[i][j][k] = S_up[0] / taupD;
                        aux->v_t[i][j][k] = S_up[1] / taupD;
                        aux->v_p[i][j][k] = S_up[2] / taupD;
                        aux->p[i][j][k] = eos->P(aux->rho_0(i, j, k), aux->eps(i, j, k));
                        //
                        // in atmosphere after all...
                        //
                        if (aux->rho_0[i][j][k] <= rho_thr || aux->eps[i][j][k] <= eps_atm) {
                            //===========================================================
                            // ATMOSPHERE
                            //===========================================================
                            aux->rho_0[i][j][k] = rho_atm;
                            aux->v_r[i][j][k] = 0.0;
                            aux->v_t[i][j][k] = 0.0;
                            aux->v_p[i][j][k] = 0.0;
                            aux->eps[i][j][k] = eps_atm;
                            aux->p[i][j][k] = p_atm;
                            //
                            // update conserved variables accordingly
                            //
                            double h_rr_l = s->h_rr(i, j, k);
                            double h_rt_l = s->h_rt(i, j, k);
                            double h_rp_l = s->h_rp(i, j, k);
                            double h_tt_l = s->h_tt(i, j, k);
                            double h_tp_l = s->h_tp(i, j, k);
                            double h_pp_l = s->h_pp(i, j, k);
                            double phi_l = s->phi(i, j, k);
                            double lapse_l = s->lapse(i, j, k);
                            double det_l = c->det(i, j, k);
                            //
                            double D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l;
                            Compute_Conserved_Vars(rho_atm, eps_atm, p_atm, 0.0, 0.0, 0.0,
                                h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l,
                                phi_l, det_l,
                                r_l, sintheta_l,
                                D_l, S_r_l, S_t_l, S_p_l, tau_l, W_l, v2_l);
                            m->D[i][j][k] = D_l;
                            m->S_r[i][j][k] = S_r_l;
                            m->S_t[i][j][k] = S_t_l;
                            m->S_p[i][j][k] = S_p_l;
                            m->tau[i][j][k] = tau_l;
                            aux->W[i][j][k] = W_l;
                        }
                    } else {
                        aux->hydro_errors[i][j][k] += 100;
                        //
                        // iteration did not converge
                        //
                        const double tauDp = taul + p_guess + Dl;
                        const double W_guess = tauDp / sqrt(tauDp * tauDp - S2);
                        // cout << " In Recovery: iteration did not converge at r = " << r_l << " i,j,k = "  
                        //	 << i << ", " << j << ", " << k << " S2 = " << S2 << " tau(tau+2D) " 
                        //	 << taul*(taul+ 2.0*Dl) << endl;
                        //
                        // take some emergency measures...
                        //
                        aux->rho_0[i][j][k] = rho_0_guess;
                        aux->eps[i][j][k] = eos->cold_eps(rho_0_guess);
                        aux->p[i][j][k] = eos->P(rho_0_guess, aux->eps(i, j, k));
                        aux->W[i][j][k] = 1.0;
                        aux->v_r[i][j][k] = 0.0;
                        aux->v_t[i][j][k] = 0.0;
                        aux->v_p[i][j][k] = 0.0;
                        error++;
                    }
                }
            }
        }
    };
    aux->rho_0.fill_ghosts();
    aux->eps.fill_ghosts();
    aux->v_r.fill_ghosts();
    aux->v_t.fill_ghosts();
    aux->v_p.fill_ghosts();
    aux->p.fill_ghosts();
    aux->W.fill_ghosts();
    m->D.fill_ghosts();
    m->S_r.fill_ghosts();
    m->S_t.fill_ghosts();
    m->S_p.fill_ghosts();
    m->tau.fill_ghosts();
    recovery_error_counter += error;
    return error;
};
//============================================================
// compute function f(p) (see eq. (A.5) in TBB)
//============================================================
double Hydro::f(double p, double tau, double D, double S2, double& dfdp) {
    const double taupD = tau + p + D;
    const double taupD2 = taupD * taupD;
    const double number = taupD2 - S2;
    if (number < 0.0) {
        cout << " trouble in HYDRO::f : number = " << taupD2 - S2 << " not supposed to be negative! " << endl;
    }
    // CHECK hack...
    const double SQR = sqrt(fabs(number));
    //
    // compute W(p) (eq. (A.2))
    //
    const double W = taupD / SQR;
    //
    // compute rho_0(p) and D(p) (eqs. (A.3) and (A.4))
    // 
    const double rho_0 = D / W;
    const double eps = (SQR - W * p - D) / D;
    //
    // compute partial rho / partial p (eq. (A.8))
    //
    const double drhodp = D * S2 / (taupD2 * SQR);
    //
    // compute partial eps / partial p (eq. (A.9))
    //
    const double depsdp = p * S2 / (D * (taupD2 - S2) * SQR);
    //
    // compute P and its derivatives from rho_0 and eps using eos 
    //
    const double P = eos->P(rho_0, eps);
    const double chi = eos->dPdrho_0(rho_0, eps);
    const double kappa = eos->dPdepsilon(rho_0, eps);
    //
    // compute dfdp  (eq. (A.7))
    //
    dfdp = 1.0 - chi * drhodp - kappa * depsdp;
    //
    //return f = p - P
    //
    return p - P;
};

//
//================================================
// Reconstruct primitive variables at grid interfaces
//
// Fluxes are stored at cell interfaces,
//
//      p_L_i = p_i-1/2-eps
//      p_R_i = p_i-1/2+eps
//
//  Will loop over interfaces, and compute
// left and right fluxes at each interface
//
//     |   x   |   x   |   x   |   x   |   x   |
//                         i
//                   i-1/2
//           df_m    df_c    df_p
//              sigma_m sigma_p     
//                    L R
//
// For logarithmic grid: assume that cell interfaces are half-way between 
// grid-points
//
//================================================
//
int Hydro::Reconstruct_r(gf3d& fct, gf3d& fct_L, gf3d& fct_R) {
    int error = 0;
    //
    // loop over interior grid interfaces
    //
    for (int i = N_g; i < N_r - N_g + 1; i++) {
        const double rmt = grid->r(i - 2);
        const double rmo = grid->r(i - 1);
        const double rl = grid->r(i);
        const double rpo = grid->r(i + 1);
        for (int j = N_g; j < N_t - N_g; j++)
            for (int k = N_g; k < N_p - N_g; k++) {
                // recall: index i refers to interface i-1/2...
                const double df_m = (fct[i - 1][j][k] - fct[i - 2][j][k]) / (rmo - rmt);
                const double df_c = (fct[i][j][k] - fct[i - 1][j][k]) / (rl - rmo);
                const double df_p = (fct[i + 1][j][k] - fct[i][j][k]) / (rpo - rl);
                bool minmod;
                // use minmod no matter what in innermost 2 gridcells? 
                if (i < N_g + 2)
                    minmod = true;
                else
                    minmod = false;
                minmod = false;   // no...
                const double sigma_m = slope->limiter(df_m, df_c, minmod);
                const double sigma_p = slope->limiter(df_c, df_p, minmod);
                // reconstruct at i-1/2
                fct_L[i][j][k] = fct[i - 1][j][k] + 0.5 * sigma_m * (rl - rmo);
                fct_R[i][j][k] = fct[i][j][k] - 0.5 * sigma_p * (rl - rmo);
            }
    }
    return error;
};
int Hydro::Reconstruct_t(gf3d& fct, gf3d& fct_L, gf3d& fct_R) {
    int error = 0;
    //
    // loop over interior grid interfaces
    //
    for (int i = N_g; i < N_r - N_g; i++)
        for (int j = N_g; j < N_t - N_g + 1; j++) {
            const double thetamt = grid->theta(j - 2);
            const double thetamo = grid->theta(j - 1);
            const double thetal = grid->theta(j);
            const double thetapo = grid->theta(j + 1);
            for (int k = N_g; k < N_p - N_g; k++) {
                const double df_m = (fct[i][j - 1][k] - fct[i][j - 2][k]) / (thetamo - thetamt);
                const double df_c = (fct[i][j][k] - fct[i][j - 1][k]) / (thetal - thetamo);
                const double df_p = (fct[i][j + 1][k] - fct[i][j][k]) / (thetapo - thetal);
                bool minmod = false;
                const double sigma_m = slope->limiter(df_m, df_c, minmod);
                const double sigma_p = slope->limiter(df_c, df_p, minmod);
                fct_L[i][j][k] = fct[i][j - 1][k] + 0.5 * sigma_m * (thetal - thetamo);
                fct_R[i][j][k] = fct[i][j][k] - 0.5 * sigma_p * (thetal - thetamo);
            }
        }
    return error;
};
int Hydro::Reconstruct_p(gf3d& fct, gf3d& fct_L, gf3d& fct_R) {
    int error = 0;
    //
    // loop over interior grid interfaces
    //
    const double dphi = grid->delta_phi();    // so far grid in phi is uniform...
    for (int i = N_g; i < N_r - N_g; i++)
        for (int j = N_g; j < N_t - N_g; j++)
            for (int k = N_g; k < N_p - N_g + 1; k++) {
                const double df_m = (fct[i][j][k - 1] - fct[i][j][k - 2]) / dphi;
                const double df_c = (fct[i][j][k] - fct[i][j][k - 1]) / dphi;
                const double df_p = (fct[i][j][k + 1] - fct[i][j][k]) / dphi;
                bool minmod = false;
                const double sigma_m = slope->limiter(df_m, df_c, minmod);
                const double sigma_p = slope->limiter(df_c, df_p, minmod);
                fct_L[i][j][k] = fct[i][j][k - 1] + 0.5 * sigma_m * dphi;
                fct_R[i][j][k] = fct[i][j][k] - 0.5 * sigma_p * dphi;
            }
    return error;
};
//
//================================================
// Compute eigenvalues (see eqs. (35) and (36) in TBB)
//================================================
//
int Hydro::lambda(double lapse, double v, double v2,
    double shift, double gamma, double cs,
    double& lambda_0, double& lambda_p, double& lambda_m) {
    int error = 0;
    if (v2 >= 1.0) {
        // cout << " Trouble in max_lambda: v2 = " << v2 << endl;
        error++;
    }
    if (cs >= 1.0) cs_error += 1;
    // if (cs >= 1.0 || cs < 0.0) {
    //   cout << " Trouble in max_lambda: cs = " << cs << endl;
    //   error++;
    // }
    const double cs2 = cs * cs;
    const double omv2cs2 = 1.0 - v2 * cs2;
    const double omcs2 = 1.0 - cs2;
    double number = (1.0 - v2) * (gamma * omv2cs2 - v * v * omcs2);
    //  if (number < 0.0) cout << " More trouble in max_lambda... " << endl;
    const double help = cs * sqrt(fabs(number));
    lambda_0 = lapse * v - shift;
    lambda_p = lapse / omv2cs2 * (v * omcs2 + help) - shift;
    lambda_m = lapse / omv2cs2 * (v * omcs2 - help) - shift;
    // if (number != fabs(number) ) {
    //   cout << " number = " << number << "   " << fabs(number) << endl;
    //   cout << "          " << v2 << "  "   << omv2cs2 << "   " <<  endl; 
    //   cout << "          " << lambda_p << "  "   << lambda_m << "   " <<  endl; 
    //   lambda_p = 0.0;
    //   lambda_m = 0.0;
    // }
    return error;
};

//
//================================================
// Construct fluxes at interfaces
//================================================
//
int Hydro::Construct_Fluxes(state* s, curvature* c) {
    int error = 0;
    //
    //==============================================
    // deal with r-direction
    //==============================================
    //
    // first reconstruct primitive variables at cell interfaces
    //
    Reconstruct_r(aux->rho_0, aux->rho_L, aux->rho_R);
    Reconstruct_r(aux->v_r, aux->v_r_L, aux->v_r_R);
    Reconstruct_r(aux->v_t, aux->v_t_L, aux->v_t_R);
    Reconstruct_r(aux->v_p, aux->v_p_L, aux->v_p_R);
    Reconstruct_r(aux->eps, aux->eps_L, aux->eps_R);
    //
    // now go to each interior cell interface
    //
    int i_test = -66;
    int j_test = 3;
    int k_test = 3;
#pragma omp parallel for collapse(3)
    for (int i = N_g; i < N_r - N_g + 1; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                // NOTE: interface half-way between grid points:
                double r_l = 0.5 * (grid->r(i - 1) + grid->r(i));
                double sintheta_l = grid->sintheta(j);
                //
                // interpolate grid functions to left interfaces
                //
                double h_rr_l = s->h_rr(r_l, j, k);
                double h_rt_l = s->h_rt(r_l, j, k);
                double h_rp_l = s->h_rp(r_l, j, k);
                double h_tt_l = s->h_tt(r_l, j, k);
                double h_tp_l = s->h_tp(r_l, j, k);
                double h_pp_l = s->h_pp(r_l, j, k);
                double phi_l = s->phi(r_l, j, k);
                double lapse_l = s->lapse(r_l, j, k);
                double shift_r_l = s->shift_r(r_l, j, k);
                double det_l = c->det(r_l, j, k);
#ifdef SR
                double g_up_rr_l = 1.0;
#else
                double g_up_rr_l = ((1.0 + h_tt_l) * (1.0 + h_pp_l) - h_tp_l * h_tp_l) / det_l;
                g_up_rr_l *= exp(-4.0 * phi_l);
#endif
                //
                // compute conserved variables, fluxes and max eigenvalues from left states
                //
                double D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L;
                const double rho_L_l = aux->rho_L(i, j, k);
                const double eps_L_l = aux->eps_L(i, j, k);
                double p_L = eos->P(rho_L_l, eps_L_l);
                Compute_Conserved_Vars(rho_L_l, eps_L_l, p_L,
                    aux->v_r_L(i, j, k), aux->v_t_L(i, j, k), aux->v_p_L(i, j, k),
                    h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, phi_l, det_l,
                    r_l, sintheta_l,
                    D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L);
                double f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L;
                Compute_Fluxes(0, r_l, sintheta_l, lapse_l, phi_l, det_l,
                    D_L, aux->v_r_L(i, j, k), shift_r_l, S_r_L, S_t_L, S_p_L, p_L, tau_L,
                    f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L);
                double cs = eos->sound_speed(rho_L_l, eps_L_l);
                double lam_0_L, lam_p_L, lam_m_L;
                error += lambda(lapse_l, aux->v_r_L(i, j, k), v2_L, shift_r_l, g_up_rr_l, cs,
                    lam_0_L, lam_p_L, lam_m_L);
                //
                // compute conserved variables, fluxes and max eigenvalues from right states
                //
                double D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R;
                const double rho_R_l = aux->rho_R(i, j, k);
                const double eps_R_l = aux->eps_R(i, j, k);
                double p_R = eos->P(rho_R_l, eps_R_l);
                Compute_Conserved_Vars(rho_R_l, eps_R_l, p_R,
                    aux->v_r_R(i, j, k), aux->v_t_R(i, j, k), aux->v_p_R(i, j, k),
                    h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, phi_l, det_l,
                    r_l, sintheta_l,
                    D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R);
                double f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R;
                Compute_Fluxes(0, r_l, sintheta_l, lapse_l, phi_l, det_l,
                    D_R, aux->v_r_R(i, j, k), shift_r_l, S_r_R, S_t_R, S_p_R, p_R, tau_R,
                    f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R);
                cs = eos->sound_speed(rho_R_l, eps_R_l);
                double lam_0_R, lam_p_R, lam_m_R;
                error += lambda(lapse_l, aux->v_r_R(i, j, k), v2_R, shift_r_l, g_up_rr_l, cs,
                    lam_0_R, lam_p_R, lam_m_R);
                //
                // now compute fluxes using approximate Riemann solver
                //
                aux->f_D_r[i][j][k] = riemann_solver->flux(f_D_L, f_D_R, D_L, D_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_r_r[i][j][k] = riemann_solver->flux(f_S_r_L, f_S_r_R, S_r_L, S_r_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_t_r[i][j][k] = riemann_solver->flux(f_S_t_L, f_S_t_R, S_t_L, S_t_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_p_r[i][j][k] = riemann_solver->flux(f_S_p_L, f_S_p_R, S_p_L, S_p_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_tau_r[i][j][k] = riemann_solver->flux(f_tau_L, f_tau_R, tau_L, tau_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                if (i == i_test && j == j_test && k == k_test)
                    cout << " in Construct_Flux : " << setprecision(16) << aux->f_tau_r(i, j, k)
                    << " tau = " << inter->tau(i, j, k)
                    << " tau_R = " << tau_R << " tau_L = " << tau_L
                    << endl;
            }
        }
    }
    //
    //==============================================
    // deal with theta-direction
    //==============================================
    //
    // first reconstruct primite variables at cell interfaces
    //
    Reconstruct_t(aux->rho_0, aux->rho_L, aux->rho_R);
    Reconstruct_t(aux->v_r, aux->v_r_L, aux->v_r_R);
    Reconstruct_t(aux->v_t, aux->v_t_L, aux->v_t_R);
    Reconstruct_t(aux->v_p, aux->v_p_L, aux->v_p_R);
    Reconstruct_t(aux->eps, aux->eps_L, aux->eps_R);
    //
    // now go to each interior cell interface
    //
#pragma omp parallel for collapse(3)
    for (int i = N_g; i < N_r - N_g; i++) {
        for (int j = N_g; j < N_t - N_g + 1; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                double r_l = grid->r(i);
                double theta_l = 0.5 * (grid->theta(j - 1) + grid->theta(j));
                double sintheta_l = sin(theta_l);
                //
                // interpolate grid functions to left interfaces
                //
                double h_rr_l = s->h_rr(i, theta_l, k);
                double h_rt_l = s->h_rt(i, theta_l, k);
                double h_rp_l = s->h_rp(i, theta_l, k);
                double h_tt_l = s->h_tt(i, theta_l, k);
                double h_tp_l = s->h_tp(i, theta_l, k);
                double h_pp_l = s->h_pp(i, theta_l, k);
                double phi_l = s->phi(i, theta_l, k);
                double lapse_l = s->lapse(i, theta_l, k);
                double shift_t_l = s->shift_t(i, theta_l, k);
                double det_l = c->det(i, theta_l, k);
#ifdef SR
                double g_up_tt_l = 1.0;
#else
                double g_up_tt_l = ((1.0 + h_rr_l) * (1.0 + h_pp_l) - h_rp_l * h_rp_l) / det_l;
                g_up_tt_l *= exp(-4.0 * phi_l);
#endif
                //
                // compute conserved variables, fluxes and max eigenvalues from left states
                //
                double D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L;
                const double rho_L_l = aux->rho_L(i, j, k);
                const double eps_L_l = aux->eps_L(i, j, k);
                double p_L = eos->P(rho_L_l, eps_L_l);
                Compute_Conserved_Vars(rho_L_l, eps_L_l, p_L,
                    aux->v_r_L(i, j, k), aux->v_t_L(i, j, k), aux->v_p_L(i, j, k),
                    h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, phi_l, det_l,
                    r_l, sintheta_l,
                    D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L);
                double f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L;
                Compute_Fluxes(1, r_l, sintheta_l, lapse_l, phi_l, det_l,
                    D_L, aux->v_t_L(i, j, k), shift_t_l, S_r_L, S_t_L, S_p_L, p_L, tau_L,
                    f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L);
                double cs = eos->sound_speed(rho_L_l, eps_L_l);
                double lam_0_L, lam_p_L, lam_m_L;
                error += lambda(lapse_l, aux->v_t_L(i, j, k), v2_L, shift_t_l, g_up_tt_l, cs, lam_0_L, lam_p_L, lam_m_L);
                //
                // compute conserved variables, fluxes and max eigenvalues from right states
                //
                double D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R;
                const double rho_R_l = aux->rho_R(i, j, k);
                const double eps_R_l = aux->eps_R(i, j, k);
                double p_R = eos->P(rho_R_l, eps_R_l);
                Compute_Conserved_Vars(rho_R_l, eps_R_l, p_R,
                    aux->v_r_R(i, j, k), aux->v_t_R(i, j, k), aux->v_p_R(i, j, k),
                    h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, phi_l, det_l,
                    r_l, sintheta_l,
                    D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R);
                double f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R;
                Compute_Fluxes(1, r_l, sintheta_l, lapse_l, phi_l, det_l,
                    D_R, aux->v_t_R(i, j, k), shift_t_l, S_r_R, S_t_R, S_p_R, p_R, tau_R,
                    f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R);
                cs = eos->sound_speed(rho_R_l, eps_R_l);
                double lam_0_R, lam_p_R, lam_m_R;
                error += lambda(lapse_l, aux->v_t_R(i, j, k), v2_R, shift_t_l, g_up_tt_l, cs, lam_0_R, lam_p_R, lam_m_R);
                //
                // now compute fluxes using approximate Riemann solver
                //
                aux->f_D_t[i][j][k] = riemann_solver->flux(f_D_L, f_D_R, D_L, D_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_r_t[i][j][k] = riemann_solver->flux(f_S_r_L, f_S_r_R, S_r_L, S_r_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_t_t[i][j][k] = riemann_solver->flux(f_S_t_L, f_S_t_R, S_t_L, S_t_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_p_t[i][j][k] = riemann_solver->flux(f_S_p_L, f_S_p_R, S_p_L, S_p_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_tau_t[i][j][k] = riemann_solver->flux(f_tau_L, f_tau_R, tau_L, tau_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
            }
        }
    }
    //
    //==============================================
    // deal with phi-direction
    //==============================================
    //
    // first reconstruct primite variables at cell interfaces
    //
    Reconstruct_p(aux->rho_0, aux->rho_L, aux->rho_R);
    Reconstruct_p(aux->v_r, aux->v_r_L, aux->v_r_R);
    Reconstruct_p(aux->v_t, aux->v_t_L, aux->v_t_R);
    Reconstruct_p(aux->v_p, aux->v_p_L, aux->v_p_R);
    Reconstruct_p(aux->eps, aux->eps_L, aux->eps_R);
    //
    // now go to each interior cell interface
    //
#pragma omp parallel for collapse(3)
    for (int i = N_g; i < N_r - N_g; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g + 1; k++) {
                double r_l = grid->r(i);
                double sintheta_l = grid->sintheta(j);
                double p_l = 0.5 * (grid->phi(k - 1) + grid->phi(k));
                //
                // interpolate grid functions to left interfaces
                //
                double h_rr_l = s->h_rr(i, j, p_l);
                double h_rt_l = s->h_rt(i, j, p_l);
                double h_rp_l = s->h_rp(i, j, p_l);
                double h_tt_l = s->h_tt(i, j, p_l);
                double h_tp_l = s->h_tp(i, j, p_l);
                double h_pp_l = s->h_pp(i, j, p_l);
                double phi_l = s->phi(i, j, p_l);
                double lapse_l = s->lapse(i, j, p_l);
                double shift_p_l = s->shift_p(i, j, p_l);
                double det_l = c->det(i, j, p_l);
#ifdef SR
                double g_up_pp_l = 1.0;
#else
                double g_up_pp_l = ((1.0 + h_rr_l) * (1.0 + h_tt_l) - h_rt_l * h_rt_l) / det_l;
                g_up_pp_l *= exp(-4.0 * phi_l);
#endif
                //
                // compute conserved variables, fluxes and max eigenvalues from left states
                //
                double D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L;
                const double rho_L_l = aux->rho_L(i, j, k);
                const double eps_L_l = aux->eps_L(i, j, k);
                double p_L = eos->P(rho_L_l, eps_L_l);
                Compute_Conserved_Vars(rho_L_l, eps_L_l, p_L,
                    aux->v_r_L(i, j, k), aux->v_t_L(i, j, k), aux->v_p_L(i, j, k),
                    h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, phi_l, det_l,
                    r_l, sintheta_l,
                    D_L, S_r_L, S_t_L, S_p_L, tau_L, W_L, v2_L);
                double f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L;
                Compute_Fluxes(2, r_l, sintheta_l, lapse_l, phi_l, det_l,
                    D_L, aux->v_p_L(i, j, k), shift_p_l, S_r_L, S_t_L, S_p_L, p_L, tau_L,
                    f_D_L, f_S_r_L, f_S_t_L, f_S_p_L, f_tau_L);
                double cs = eos->sound_speed(rho_L_l, eps_L_l);
                double lam_0_L, lam_p_L, lam_m_L;
                error += lambda(lapse_l, aux->v_p_L(i, j, k), v2_L, shift_p_l, g_up_pp_l, cs,
                    lam_0_L, lam_p_L, lam_m_L);
                //
                // compute conserved variables, fluxes and max eigenvalues from right states
                //
                double D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R;
                const double rho_R_l = aux->rho_R(i, j, k);
                const double eps_R_l = aux->eps_R(i, j, k);
                double p_R = eos->P(rho_R_l, eps_R_l);
                Compute_Conserved_Vars(rho_R_l, eps_R_l, p_R,
                    aux->v_r_R(i, j, k), aux->v_t_R(i, j, k), aux->v_p_R(i, j, k),
                    h_rr_l, h_rt_l, h_rp_l, h_tt_l, h_tp_l, h_pp_l, phi_l, det_l,
                    r_l, sintheta_l,
                    D_R, S_r_R, S_t_R, S_p_R, tau_R, W_R, v2_R);
                double f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R;
                Compute_Fluxes(2, r_l, sintheta_l, lapse_l, phi_l, det_l,
                    D_R, aux->v_p_R(i, j, k), shift_p_l, S_r_R, S_t_R, S_p_R, p_R, tau_R,
                    f_D_R, f_S_r_R, f_S_t_R, f_S_p_R, f_tau_R);
                cs = eos->sound_speed(rho_R_l, eps_R_l);
                double lam_0_R, lam_p_R, lam_m_R;
                error += lambda(lapse_l, aux->v_p_R(i, j, k), v2_R, shift_p_l, g_up_pp_l, cs,
                    lam_0_R, lam_p_R, lam_m_R);
                //
                // now compute fluxes using approximate Riemann solver
                //
                aux->f_D_p[i][j][k] = riemann_solver->flux(f_D_L, f_D_R, D_L, D_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_r_p[i][j][k] = riemann_solver->flux(f_S_r_L, f_S_r_R, S_r_L, S_r_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_t_p[i][j][k] = riemann_solver->flux(f_S_t_L, f_S_t_R, S_t_L, S_t_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_S_p_p[i][j][k] = riemann_solver->flux(f_S_p_L, f_S_p_R, S_p_L, S_p_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
                aux->f_tau_p[i][j][k] = riemann_solver->flux(f_tau_L, f_tau_R, tau_L, tau_R,
                    lam_0_L, lam_p_L, lam_m_L,
                    lam_0_R, lam_p_R, lam_m_R);
            }
        }
    }
    //
    return error;
};
//
//================================================
// Compute Fluxes from conserved hydro variables 
//
// Get correct component by 
// (1) passing in corresponding component of v and shift  -- i.e., pass 
//     in v^r and \beta^r, get r components of fluxes)
// (2) pass in correct index (0 for r, 1 for theta, 2 for phi)
//
//================================================
//
int Hydro::Compute_Fluxes(int index, double r, double sintheta,
    double lapse, double phi, double det, double D,
    double v, double shift, double S_r, double S_t, double S_p,
    double p, double tau,
    double& f_D,
    double& f_S_r, double& f_S_t, double& f_S_p,
    double& f_tau) {
    int error = 0;
#ifdef SR
    double det_3D = 1.0;
    const double vel = v;
#else
    double det_3D = exp(6.0 * phi) * sqrt(det);
    const double vel = v - shift / lapse;
#endif
    f_D = lapse * D * vel;
    f_S_r = lapse * (S_r * vel + det_3D * (index == 0) * p);
    f_S_t = lapse * (S_t * vel + det_3D * (index == 1) * p);
    f_S_p = lapse * (S_p * vel + det_3D * (index == 2) * p);
    if (partial)
        f_tau = lapse * (tau * vel + det_3D * r * r * sintheta * p * v);
    else
        f_tau = lapse * (tau * vel + det_3D * p * v);
    return error;
};
//
//================================================
// 
// Compute right-hand sides...
//
// Three parts:
//   - divergence of fluxes
//   - flat connection terms
//   - source terms
//
//================================================
//
int Hydro::dot_conserved_vars(hydro_state* m, state* s, curvature* c) {
    int error = 0;
    int i_test = -66;
    int j_test = 3;
    int k_test = 3;
    //
    // loop over interior grid
    //
#pragma omp parallel for collapse(3)
    for (int i = N_g; i < N_r - N_g; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                const double rl = grid->r(i);
                const double dr = 0.5 * (grid->r(i + 1) - grid->r(i - 1));
                const double sintheta = grid->sintheta(j);
                const double rst = rl * sintheta;
                const double sin2theta = sintheta * sintheta;
                const double costheta = grid->costheta(j);
                const double cottheta = costheta / sintheta;
                const double dtheta = 0.5 * (grid->theta(j + 1) - grid->theta(j - 1));
                //
                // compute dr : distance between cell-interfaces (different
                // from distances between grid points!
                //
                const double dphi = grid->delta_phi();
                //================================================
                // first compute divergence of fluxes - note extra factors 
                // and terms due to rescaling
                //================================================
                derivs->D[i][j][k] = -
                    (aux->f_D_r[i + 1][j][k] - aux->f_D_r[i][j][k]) / dr -
                    (aux->f_D_t[i][j + 1][k] - aux->f_D_t[i][j][k]) / dtheta / rl -
                    (aux->f_D_p[i][j][k + 1] - aux->f_D_p[i][j][k]) / dphi / rst;
                derivs->S_r[i][j][k] = -
                    (aux->f_S_r_r[i + 1][j][k] - aux->f_S_r_r[i][j][k]) / dr -
                    (aux->f_S_r_t[i][j + 1][k] - aux->f_S_r_t[i][j][k]) / dtheta / rl -
                    (aux->f_S_r_p[i][j][k + 1] - aux->f_S_r_p[i][j][k]) / dphi / rst;
                derivs->S_t[i][j][k] = -
                    (aux->f_S_t_r[i + 1][j][k] - aux->f_S_t_r[i][j][k]) / dr -
                    (aux->f_S_t_t[i][j + 1][k] - aux->f_S_t_t[i][j][k]) / dtheta / rl -
                    (aux->f_S_t_p[i][j][k + 1] - aux->f_S_t_p[i][j][k]) / dphi / rst;
                derivs->S_p[i][j][k] = -
                    (aux->f_S_p_r[i + 1][j][k] - aux->f_S_p_r[i][j][k]) / dr -
                    (aux->f_S_p_t[i][j + 1][k] - aux->f_S_p_t[i][j][k]) / dtheta / rl -
                    (aux->f_S_p_p[i][j][k + 1] - aux->f_S_p_p[i][j][k]) / dphi / rst;
                derivs->tau[i][j][k] = -
                    (aux->f_tau_r[i + 1][j][k] - aux->f_tau_r[i][j][k]) / dr -
                    (aux->f_tau_t[i][j + 1][k] - aux->f_tau_t[i][j][k]) / dtheta / rl -
                    (aux->f_tau_p[i][j][k + 1] - aux->f_tau_p[i][j][k]) / dphi / rst;
                if (i == i_test && j == j_test && k == k_test) cout << " 1 " << setprecision(16) << derivs->tau(i, j, k) << " f_tau_r = " << aux->f_tau_r(i, j, k) << "   " << aux->f_tau_r(i + 1, j, k) << endl;
                //================================================
                // compute fluxes at cell centers for flat connection terms
                // as well as extra terms due to rescaling
                //================================================
                double f_D_r_c, f_D_t_c, f_D_p_c;
                double f_S_r_r_c, f_S_r_t_c, f_S_r_p_c;
                double f_S_t_r_c, f_S_t_t_c, f_S_t_p_c;
                double f_S_p_r_c, f_S_p_t_c, f_S_p_p_c;
                double f_tau_r_c, f_tau_t_c, f_tau_p_c;
                // radial fluxes
                Compute_Fluxes(0, rl, sintheta, s->lapse[i][j][k], s->phi[i][j][k], c->det[i][j][k],
                    m->D[i][j][k], aux->v_r[i][j][k], s->shift_r[i][j][k],
                    m->S_r[i][j][k], m->S_t[i][j][k], m->S_p[i][j][k],
                    aux->p[i][j][k], m->tau[i][j][k],
                    f_D_r_c, f_S_r_r_c, f_S_t_r_c, f_S_p_r_c, f_tau_r_c);
                // theta fluxes
                Compute_Fluxes(1, rl, sintheta, s->lapse[i][j][k], s->phi[i][j][k], c->det[i][j][k],
                    m->D[i][j][k], aux->v_t[i][j][k], s->shift_t[i][j][k],
                    m->S_r[i][j][k], m->S_t[i][j][k], m->S_p[i][j][k],
                    aux->p[i][j][k], m->tau[i][j][k],
                    f_D_t_c, f_S_r_t_c, f_S_t_t_c, f_S_p_t_c, f_tau_t_c);
                // // phi fluxes
                Compute_Fluxes(2, rl, sintheta, s->lapse[i][j][k], s->phi[i][j][k], c->det[i][j][k],
                    m->D[i][j][k], aux->v_p[i][j][k], s->shift_p[i][j][k],
                    m->S_r[i][j][k], m->S_t[i][j][k], m->S_p[i][j][k],
                    aux->p[i][j][k], m->tau[i][j][k],
                    f_D_p_c, f_S_r_p_c, f_S_t_p_c, f_S_p_p_c, f_tau_p_c);
                if (!partial)
                    derivs->D[i][j][k] += -f_D_r_c * 2.0 / rl - f_D_t_c * cottheta / rl;
                derivs->S_r[i][j][k] += (f_S_t_t_c + f_S_p_p_c) / rl
                    - f_S_r_r_c * 2.0 / rl - f_S_r_t_c * cottheta / rl;  // term due to rescaling
                derivs->S_t[i][j][k] += -f_S_r_t_c / rl + f_S_t_r_c / rl + f_S_p_p_c * cottheta / rl
                    - f_S_t_r_c * 2.0 / rl - f_S_t_t_c * cottheta / rl  // term due to rescaling 
                    - f_S_t_r_c / rl;   // term due to rescaling
                derivs->S_p[i][j][k] += -f_S_r_p_c / rl - f_S_t_p_c * cottheta / rl
                    + f_S_p_r_c / rl + f_S_p_t_c * cottheta / rl
                    - f_S_p_r_c * 2.0 / rl - f_S_p_t_c * cottheta / rl
                    - f_S_p_r_c / rl - cottheta / rl * f_S_p_t_c;  // terms due to rescaling
                if (!partial)
                    derivs->tau[i][j][k] += -f_tau_r_c * 2.0 / rl - f_tau_t_c * cottheta / rl;
                if (i == i_test && j == j_test && k == k_test) cout << " 2 " << derivs->tau(i, j, k) << " D = " << derivs->D(i, j, k) << endl;
                //
                //================================================
                // Finally compute source terms on right-hand side
                //================================================
                //
                // first compute components of stress-energy tensor (see eqs. (27) - (30) in TBB)
                //
#ifndef SR
                const double enth = eos->h(aux->rho_0(i, j, k), aux->eps(i, j, k));
                const double p = eos->P(aux->rho_0(i, j, k), aux->eps(i, j, k));   // CHECK: do we need this call?
                const double rhohW2 = aux->rho_0(i, j, k) * enth * aux->W(i, j, k) * aux->W(i, j, k);
                const double lapsel = s->lapse(i, j, k);
                //
                // NOTE: We will do the following calculation in terms of physical, i.e. *not* rescaled
                // variables.  Why?  D_i \gamma_jk in rest of code is not rescaled,
                // and it is easier to use that variable here, too.
                //
                //
                // physical velocity
                vect v_up(aux->v_r(i, j, k), aux->v_t(i, j, k) / rl, aux->v_p(i, j, k) / rst);
                vect v_down;
                tensor g_conf(1.0 + s->h_rr(i, j, k),
                    rl * s->h_rt(i, j, k),
                    rst * s->h_rp(i, j, k),
                    rl * rl * (1.0 + s->h_tt(i, j, k)),
                    rl * rst * s->h_tp(i, j, k),
                    rst * rst * (1.0 + s->h_pp(i, j, k)));
                tensor g_up = g_conf.inverse();
                const double e4p = exp(4.0 * s->phi(i, j, k));
                //
                // lower index of v^i (result is physical velocity, not rescaled)
                //
                for (int a = 0; a < 3; a++) {
                    v_down[a] = 0.0;
                    for (int b = 0; b < 3; b++) {
                        v_down[a] += e4p * g_conf[a][b] * v_up[b];
                    }
                }
                //
                // compute components of stress-energy tensor (use *unrescaled* shift)
                //
                vect shift(s->shift_r(i, j, k), s->shift_t(i, j, k) / rl, s->shift_p(i, j, k) / rst);
                const double T00 = (rhohW2 - p) / (lapsel * lapsel);
                vect T0i_up;  // none of these are rescaled!
                vect T0i_down;
                tensor Tij_up;
                for (int a = 0; a < 3; a++) {
                    T0i_up[a] = (rhohW2 * (v_up[a] - shift[a] / lapsel) + p * shift[a] / lapsel) / lapsel;
                    T0i_down[a] = rhohW2 * v_down[a] / lapsel;
                    for (int b = 0; b < 3; b++) {
                        Tij_up[a][b] = rhohW2 * (v_up[a] - shift[a] / lapsel) * (v_up[b] - shift[b] / lapsel)
                            + p * (g_up[a][b] / e4p - shift[a] * shift[b] / (lapsel * lapsel));
                    }
                }
                //
                // Assign derivatives of physical, spatial metric to rank-3 tensor D_gamma
                //    example: D_gamma[0][1][2] = D_r \gamma_{tp} 
                //                              = e^{4 \phi} ( D_r \bar \gamma_{tp} + 4 \bar \gamma_{tp} D_r \phi ) 
                //
                rank3tens D_gamma(e4p * (c->Dr_e_rr(i, j, k) + 4.0 * g_conf[0][0] * s->phi.dr(i, j, k)),
                    e4p * (c->Dr_e_rt(i, j, k) + 4.0 * g_conf[0][1] * s->phi.dr(i, j, k)),
                    e4p * (c->Dr_e_rp(i, j, k) + 4.0 * g_conf[0][2] * s->phi.dr(i, j, k)),
                    e4p * (c->Dr_e_tt(i, j, k) + 4.0 * g_conf[1][1] * s->phi.dr(i, j, k)),
                    e4p * (c->Dr_e_tp(i, j, k) + 4.0 * g_conf[1][2] * s->phi.dr(i, j, k)),
                    e4p * (c->Dr_e_pp(i, j, k) + 4.0 * g_conf[2][2] * s->phi.dr(i, j, k)),
                    e4p * (c->Dt_e_rr(i, j, k) + 4.0 * g_conf[0][0] * s->phi.dtheta(i, j, k)),
                    e4p * (c->Dt_e_rt(i, j, k) + 4.0 * g_conf[0][1] * s->phi.dtheta(i, j, k)),
                    e4p * (c->Dt_e_rp(i, j, k) + 4.0 * g_conf[0][2] * s->phi.dtheta(i, j, k)),
                    e4p * (c->Dt_e_tt(i, j, k) + 4.0 * g_conf[1][1] * s->phi.dtheta(i, j, k)),
                    e4p * (c->Dt_e_tp(i, j, k) + 4.0 * g_conf[1][2] * s->phi.dtheta(i, j, k)),
                    e4p * (c->Dt_e_pp(i, j, k) + 4.0 * g_conf[2][2] * s->phi.dtheta(i, j, k)),
                    e4p * (c->Dp_e_rr(i, j, k) + 4.0 * g_conf[0][0] * s->phi.dphi(i, j, k)),
                    e4p * (c->Dp_e_rt(i, j, k) + 4.0 * g_conf[0][1] * s->phi.dphi(i, j, k)),
                    e4p * (c->Dp_e_rp(i, j, k) + 4.0 * g_conf[0][2] * s->phi.dphi(i, j, k)),
                    e4p * (c->Dp_e_tt(i, j, k) + 4.0 * g_conf[1][1] * s->phi.dphi(i, j, k)),
                    e4p * (c->Dp_e_tp(i, j, k) + 4.0 * g_conf[1][2] * s->phi.dphi(i, j, k)),
                    e4p * (c->Dp_e_pp(i, j, k) + 4.0 * g_conf[2][2] * s->phi.dphi(i, j, k)));
                //
                // assign components of physical extrinsic cuvature to tensor K 
                //
                tensor Aij(s->a_rr(i, j, k),
                    rl * s->a_rt(i, j, k),
                    rl * sintheta * s->a_rp(i, j, k),
                    rl * rl * s->a_tt(i, j, k),
                    rl * rl * sintheta * s->a_tp(i, j, k),
                    rl * rl * sin2theta * s->a_pp(i, j, k));
                tensor Kij;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        Kij[a][b] = e4p * (Aij[a][b] + g_conf[a][b] * s->K(i, j, k) / 3.0);
                //
                // compute covariant derivative of shift with respect to background metric: \hat D_i \beta^j
                // result is physical scaling (i.e. not rescaled)
                //
                const double shift_r_r = s->shift_r.dr(i, j, k);
                const double shift_r_t = s->shift_r.dtheta(i, j, k) - s->shift_t(i, j, k);
                const double shift_r_p = s->shift_r.dphi(i, j, k) - sintheta * s->shift_p(i, j, k);
                const double shift_t_r = (s->shift_t.dr(i, j, k)) / rl;
                const double shift_t_t = (s->shift_t.dtheta(i, j, k) + s->shift_r(i, j, k)) / rl;
                const double shift_t_p = (s->shift_t.dphi(i, j, k) - costheta * s->shift_p(i, j, k)) / rl;
                const double shift_p_r = (s->shift_p.dr(i, j, k)) / (rst);
                const double shift_p_t = (s->shift_p.dtheta(i, j, k)) / (rst);
                const double shift_p_p = (s->shift_p.dphi(i, j, k) + sintheta * s->shift_r(i, j, k)
                    + costheta * s->shift_t(i, j, k)) / (rst);
                // CHECK that indices are right...
                tensor D_shift(shift_r_r, shift_r_t, shift_r_p,
                    shift_t_r, shift_t_t, shift_t_p,
                    shift_p_r, shift_p_t, shift_p_p);
                vect D_lapse(s->lapse.dr(i, j, k), s->lapse.dtheta(i, j, k), s->lapse.dphi(i, j, k));
                //
                // square root of spacetime determinant
                //
                const double det_4D = s->lapse(i, j, k) * exp(6.0 * s->phi(i, j, k)) * sqrt(c->det(i, j, k));
                double det_4D_tau = det_4D;
                if (partial) det_4D_tau *= rl * rl * sintheta;
                //
                // now compute source terms, and take care of rescaling
                //
                derivs->S_r[i][j][k] += -det_4D * T00 * lapsel * D_lapse[0];
                derivs->S_t[i][j][k] += -det_4D * T00 * lapsel * D_lapse[1] / rl;
                derivs->S_p[i][j][k] += -det_4D * T00 * lapsel * D_lapse[2] / rst;
                for (int a = 0; a < 3; a++) {
                    derivs->S_r[i][j][k] += det_4D * T0i_down[a] * D_shift[a][0];  // CHECK ordering of indices
                    derivs->S_t[i][j][k] += det_4D * T0i_down[a] * D_shift[a][1] / rl;
                    derivs->S_p[i][j][k] += det_4D * T0i_down[a] * D_shift[a][2] / rst;
                    derivs->tau[i][j][k] += det_4D_tau * (-T00 * shift[a] - T0i_up[a]) * D_lapse[a];
                    for (int b = 0; b < 3; b++) {
                        const double term_ab = T00 * shift[a] * shift[b] / 2.0 + T0i_up[a] * shift[b] + Tij_up[a][b] / 2.0;
                        derivs->S_r[i][j][k] += det_4D * term_ab * D_gamma[0][a][b];
                        derivs->S_t[i][j][k] += det_4D * term_ab * D_gamma[1][a][b] / rl;
                        derivs->S_p[i][j][k] += det_4D * term_ab * D_gamma[2][a][b] / rst;
                        derivs->tau[i][j][k] += det_4D_tau * 2.0 * term_ab * Kij[a][b];
                    }
                }
                if (i == i_test && j == j_test && k == k_test) cout << " 3 " << setprecision(16) << derivs->tau(i, j, k) << " D = " << derivs->D(i, j, k) << endl;
#endif  /* SR */
            }
        }
    }
    if (cowling == 3) {
        derivs->S_r.equals(0.0);
        derivs->S_t.equals(0.0);
        derivs->S_p.equals(0.0);
    }
    // 
    // don't evolve density variables for n_fixed innermost radial gridpoints
    // if flag fix_hydro is set:
    //
    if (fix_hydro)
        for (int i = N_g; i < N_g + n_fixed; i++)
            for (int j = N_g; j < N_t - N_g; j++)
                for (int k = N_g; k < N_p - N_g; k++) {
                    derivs->D[i][j][k] = 0.0;
                    derivs->tau[i][j][k] = 0.0;
                }
    //
    return error;
};
//
//================================================
//
// Compute ADM sources
//
// NOTE: will compute unrescaled sources, as expected in field equations
//================================================
//
void Hydro::ADM_Sources(state* s, curvature* c) {
#pragma omp parallel for collapse(3)  
    for (int i = N_g; i < N_r - N_g; i++) {
        for (int j = N_g; j < N_t - N_g; j++) {
            for (int k = N_g; k < N_p - N_g; k++) {
                const double rl = grid->r(i);
                const double sinthetal = grid->sintheta(j);
                const double rst = rl * sinthetal;
                const double e4p = exp(4.0 * s->phi(i, j, k));
                tensor gup(c->gup_rr(i, j, k), c->gup_rt(i, j, k), c->gup_rp(i, j, k),
                    c->gup_tt(i, j, k), c->gup_tp(i, j, k), c->gup_pp(i, j, k));
                tensor g_conf = gup.inverse();
                //
                // adopt "Valencia" definition of v^i, use unrescaled velocity now
                //
                vect v_up(aux->v_r(i, j, k), aux->v_t(i, j, k) / rl, aux->v_p(i, j, k) / rst);
                double vdotv = 0.0;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        vdotv += e4p * g_conf[a][b] * v_up[a] * v_up[b];
                //
                // compute Lorentz factor W
                //
                double Wl = aux->W(i, j, k);
                //
                // compute u_i (lower indices!)
                //
                vect u;
                for (int a = 0; a < 3; a++) {
                    u[a] = 0.0;
                    for (int b = 0; b < 3; b++)
                        u[a] += Wl * e4p * g_conf[a][b] * v_up[b];
                }
                //
                // find some thermodynamic quantities
                //	
                const double rho_rest = aux->rho_0(i, j, k);
                const double epsilon = aux->eps(i, j, k);
                const double enth = eos->h(rho_rest, epsilon);
                const double P = eos->P(rho_rest, epsilon);
                //
                // Now compute ADM density...
                // 
                adm_sources->rho_ADM[i][j][k] = rho_rest * enth * Wl * Wl - P;
                //
                // ... fluxes ...
                //
                adm_sources->S_r[i][j][k] = rho_rest * enth * Wl * u[0];
                adm_sources->S_t[i][j][k] = rho_rest * enth * Wl * u[1];
                adm_sources->S_p[i][j][k] = rho_rest * enth * Wl * u[2];
                //
                // ... stresses ...
                //
                adm_sources->S_rr[i][j][k] = P * e4p * g_conf[0][0] + rho_rest * enth * u[0] * u[0];
                adm_sources->S_rt[i][j][k] = P * e4p * g_conf[0][1] + rho_rest * enth * u[0] * u[1];
                adm_sources->S_rp[i][j][k] = P * e4p * g_conf[0][2] + rho_rest * enth * u[0] * u[2];
                adm_sources->S_tt[i][j][k] = P * e4p * g_conf[1][1] + rho_rest * enth * u[1] * u[1];
                adm_sources->S_tp[i][j][k] = P * e4p * g_conf[1][2] + rho_rest * enth * u[1] * u[2];
                adm_sources->S_pp[i][j][k] = P * e4p * g_conf[2][2] + rho_rest * enth * u[2] * u[2];
                //
                // ... and trace of stress:
                //
                adm_sources->trace_S[i][j][k] = 3.0 * P + rho_rest * enth * (Wl * Wl - 1.0);
                //
                //=====================================
                // now compute (unrescaled) fluxes \rho_0 u^i (index up!)
                //=====================================
                //
                double lapsel = s->lapse(i, j, k);
                vect shift(s->shift_r(i, j, k),
                    s->shift_t(i, j, k) / rl,
                    s->shift_p(i, j, k) / rst);
                fluxes->j_r[i][j][k] = rho_rest * Wl * (v_up[0] - shift[0] / lapsel);
                fluxes->j_t[i][j][k] = rho_rest * Wl * (v_up[1] - shift[1] / lapsel);
                fluxes->j_p[i][j][k] = rho_rest * Wl * (v_up[2] - shift[2] / lapsel);
            }
        }
    }
};
//
//================================================
// 
// Rest Mass
//
//================================================
//
double Hydro::Rest_Mass(hydro_state* m) {
    double rest_mass = 0.0;
    if (N_t == 6) {
        // assuming spherical symmetry...
        const double PI = acos(-1.0);
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double dr = grid->delta_r(i);
            rest_mass += m->D(i, 2, 2) * rl * rl * 4.0 * PI * dr;
        }
    } else {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double dr = grid->delta_r(i);
            for (int j = N_g; j < N_t - N_g; j++) {
                const double dtheta = grid->delta_theta(j);
                const double sintheta = grid->sintheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    if (!horizonfinder->InsideHorizon(rl, j, k) && aux->rho_0(i, j, k) > f_thr * rho_atm) {
                        const double dphi = grid->delta_phi(k);
                        if (partial)
                            rest_mass += m->D(i, j, k) * dr * dtheta * dphi;
                        else
                            rest_mass += m->D(i, j, k) * rl * rl * sintheta * dr * dtheta * dphi;
                    }
                }
            }
        }
    }
#ifdef EQSYMMETRY
    return 2.0 * rest_mass;
#else
    return rest_mass;
#endif
}
//
//================================================
// 
// Unbound (rest) Mass
//
//================================================
//
double Hydro::Unbound_Mass(hydro_state* m) {
    double rest_mass = 0.0;
    if (N_t == 6) {
        // assuming spherical symmetry...
        const double PI = acos(-1.0);
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double dr = grid->delta_r(i);
            if (aux->ut_low(i, 2, 2) < -1.0)
                rest_mass += m->D(i, 2, 2) * rl * rl * 4.0 * PI * dr;
        }
    } else {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double dr = grid->delta_r(i);
            for (int j = N_g; j < N_t - N_g; j++) {
                const double dtheta = grid->delta_theta(j);
                const double sintheta = grid->sintheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    if (!horizonfinder->InsideHorizon(rl, j, k) &&
                        aux->ut_low(i, j, k) < -1.0 && aux->v_r(i, j, k) > 0.0) {
                        const double dphi = grid->delta_phi(k);
                        if (partial)
                            rest_mass += m->D(i, j, k) * dr * dtheta * dphi;
                        else
                            rest_mass += m->D(i, j, k) * rl * rl * sintheta * dr * dtheta * dphi;
                    }
                }
            }
        }
    }
#ifdef EQSYMMETRY
    return 2.0 * rest_mass;
#else
    return rest_mass;
#endif
}
//
//================================================
// 
// Test recovery
//
//================================================
//
int Hydro::Test_Recovery(state* s, curvature* c) {
    int i_test = 2;
    int j_test = 2;
    int k_test = 2;
    cout << " HYDRO: rho_0 initially = " << aux->rho_0(i_test, j_test, k_test) << endl;
    cout << " HYDRO: Computing conserveds from primitives..." << endl;
    Compute_Conserved_Vars(last, s, c);
    // const double Dl =  4.04157e-11;
    // const double factor = 0.01;
    // const double taul = 4.69192e-12;
    // const double Srl = - factor * 3.3e-9;
    // const double pl_max = 1.e-11;
    // double dfdp = 0.0;
    // for (double pl = - pl_max; pl <= pl_max; pl += pl_max/20) {
    //   cout << " pl = " << setw(16) << pl << " f = " << f(pl, taul, Dl, 1.e-21, dfdp) << " dfdp = " << dfdp << endl;
    // }
    for (int i = 0; i < N_r; i++)
        for (int j = 0; j < N_t; j++)
            for (int k = 0; k < N_p; k++) {
                aux->rho_0[i][j][k] = 0.0;
            }
    // 	const double det_3D = exp(6.0*s->phi(i,j,k))*sqrt( c->det(i,j,k) );
    // 	double r_l = grid->r(i);
    // 	double sintheta_l = grid->sintheta(j);
    // 	last->D[i][j][k] = Dl * det_3D;
    // 	last->tau[i][j][k] = taul * det_3D;
    // 	if (partial) {
    // 	  last->D[i][j][k] *= r_l * r_l * sintheta_l;
    // 	  last->tau[i][j][k] *= r_l * r_l * sintheta_l;
    // 	}	
    // 	last->S_r[i][j][k] = - Srl * det_3D;
    // 	last->S_t[i][j][k] = 0.0;
    //     }
    int error = Recovery(last, s, c);
    cout << " HYDRO: recovered primitives with " << error << " errors " << endl;
    cout << " HYDRO: rho_0 after recovery = " << aux->rho_0(i_test, j_test, k_test) << endl;
    cout << " HYDRO: v_r after recovery = " << aux->v_r(i_test, j_test, k_test) << endl;
    return error;
};
