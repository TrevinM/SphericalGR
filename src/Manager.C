// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing all the integration stuff...
//================================================
#ifndef MANAGER_H
#define MANAGER_H

#include "Manager.h"

//===========================================
// Constructor
//===========================================
Manager::Manager(int matter_type, int sigma_i, int cowling_i, double eta_i,
    int z4_i, double kappa_11_i, double kappa_12_i, double kappa_2_i,
    double kappa_ric_i, int RK_order_i,
    int char_OB_i, int solve_constraints_i) :
    sigma(sigma_i), cowling(cowling_i), eta(eta_i),
    z4(z4_i), kappa_11(kappa_11_i), kappa_12(kappa_12_i),
    kappa_2(kappa_2_i), kappa_ric(kappa_ric_i),
    RK_order(RK_order_i), char_OB(char_OB_i),
    solve_constraints(solve_constraints_i),
    r(grid->N_r_tot()), r2(grid->N_r_tot()), theta(grid->N_theta_tot()),
    sintheta(grid->N_theta_tot()), sin2theta(grid->N_theta_tot()),
    costheta(grid->N_theta_tot()), phi(grid->N_phi_tot()) {
    //
    cout << " MANAGER: Constructing manager... " << endl;
    //
    // initial times and step counter (provided by CheckPoint.h, which
    // returns zero unless data are read in from checkpoint files)
    //
    t = checkpoint->TStart();
    cout << " MANAGER: t = " << t << endl;
    tau_c = checkpoint->TauStart();
    step = checkpoint->TimeStepStart();
    grid->Setup_Grid(r, r2, theta, sintheta, sin2theta, costheta, phi);
    //
    // Sanity check
    //
#ifdef EIGHTHORDERTHETA
    if (grid->N_ghosts() < 4) {
        cerr << " MANAGER: Need N_ghosts >= 4 when using eighth-order differencing! " << endl;
        exit(0);
    }
#endif
    //
    // create state for dynamical variables
    //
    cout << " MANAGER: setting up states... " << endl;
    last = new state(grid, dump, "last");
    derivs = new state(grid, dump, "derivs");
    inter = new state(grid, dump, "inter");
    updates = new state(grid, dump, "updates");
    // create dump list for state last:
    last->assemble_dump_list("Dump_List");
    // let checkpointer know about last:
    checkpoint->CollectDynVariables(last);

    // 
    // create curvature class
    // 
    curve = new curvature(grid, dump, "curve");
    curve->assemble_dump_list("Dump_List");
    //
    // create class with auxiliary functions
    //
    aux = new auxiliary(grid, dump, "auxiliaries");
    aux->assemble_dump_list("Dump_List");
    //
    // create horizon finder
    //
    horizonfinder = new HorizonFinder(grid, last, curve, aux,
        monitor->Filestem());
    //
    // create diagnostics class
    //
    cout << " MANAGER: setting up diagnostics... " << endl;
    constraints = new diagnostics(grid, dump, monitor, cosmology,
        horizonfinder, "constraints");
    constraints->assemble_dump_list("Dump_List");
    //
    // set time step
    //
    SetTimeStep();
    eta_KO = eta;
    cout << " MANAGER: using Kreiss-Oliger coefficient " << eta_KO << endl;
    gauge->Set_eta(eta_KO);
    //
    // create matter
    //
    if (matter_type == 1) {
        matter = new Vacuum(grid, dump, indata, cowling, cosmology);
    } else if (matter_type == 3) {
        if (!eos) {
            cerr << " MANAGER: Can't set up Hydro without eos! " << endl;
            exit(1);
        }
        matter = new Hydro(grid, dump, indata, eos, cowling, cosmology,
            monitor, horizonfinder, checkpoint);
    } else if (matter_type == 5) {
        matter = new ScalarField(grid, dump, indata, cowling, cosmology,
            monitor, eta_KO, checkpoint);
    } else if (matter_type == 6) {
        matter = new RadHydro(grid, dump, indata, eos, cowling, cosmology,
            monitor, checkpoint);
    } else if (matter_type == 7) {
        matter = new Maxwell(grid, dump, indata, cowling, char_OB, cosmology,
            monitor, eta_KO, checkpoint);
    } else if (matter_type == 8) {
        matter = new DualMaxwell(grid, dump, indata, cowling, char_OB, cosmology,
            monitor, eta_KO, checkpoint);
    } else {
        cerr << " MANAGER: Unknown matter type!!! " << endl;
    }
    //===========================================
    // set up profiles
    //===========================================
    if (!(profiles == NULL)) {
        profiles->add_to_profile_list("Profile_List", last);
        profiles->add_to_profile_list("Profile_List", curve);
        profiles->add_to_profile_list("Profile_List", aux);
        profiles->add_to_profile_list("Profile_List", constraints);
        profiles->initialize_file(constraints->R_prop.Address());
    }
    //
    //===========================================
    // initialize wave extraction
    //===========================================
    if (!(waves == NULL)) {
        waves->Initialize(constraints->psi4_Re.Address(),
            constraints->psi4_Im.Address());
    }
    N_g = inter->grid->N_ghosts();
    N_r = inter->grid->N_r_tot();
    N_t = inter->grid->N_theta_tot();
    N_p = inter->grid->N_phi_tot();
    //
    steps_between_regrids = 100;
    timestep_last_regrid = -steps_between_regrids;  // to allow regrid at t=0
    //
    PI = acos(-1.0);
}

//===========================================
// Find time step
//===========================================
Manager::SetTimeStep() {
    int N_g = grid->N_ghosts();
    dt = grid->courant_factor() * 0.5 * grid->delta_r(N_g)
        * grid->delta_theta(N_g);
    cout << " MANAGER: using dt = " << dt << endl;
};
Manager::SetTimeStep(state* s) {
    int N_g = grid->N_ghosts();
    const double gamma_tt = (1.0 + s->h_tt(N_g, N_g, N_g)) * exp(4.0 * s->phi(N_g, N_g, N_g));
    // const double factor = sqrt(gamma_tt) / s->lapse(N_g, N_g, N_g); // + abs(s->shift_r(N_g, N_g, N_g));
    const double factor = sqrt(gamma_tt); // + abs(s->shift_r(N_g, N_g, N_g));
    dt = grid->courant_factor() * 0.5 * grid->delta_r(N_g) * grid->delta_theta(N_g) * factor;
};
//===========================================
// Initialize
//===========================================
Manager::Set_t_max(double set_t) { t_max = set_t; }

//================================================
//
// Rescale metric
// 
//================================================
Manager::Rescale_Metric(state* s) {
    curve->Compute_Determinant(s);
#pragma omp parallel for collapse(3)
    for (int i = 0; i < N_r; i++)
        for (int j = 0; j < N_t; j++)
            for (int k = 0; k < N_p; k++) {
                const double factor = pow(curve->det_init(i, j, k) / curve->det(i, j, k), 1.0 / 3.0);
                s->h_rr[i][j][k] = factor * (1.0 + s->h_rr(i, j, k)) - 1.0;
                s->h_rt[i][j][k] *= factor;
                s->h_rp[i][j][k] *= factor;
                s->h_tt[i][j][k] = factor * (1.0 + s->h_tt(i, j, k)) - 1.0;
                s->h_tp[i][j][k] *= factor;
                s->h_pp[i][j][k] = factor * (1.0 + s->h_pp(i, j, k)) - 1.0;
                s->phi[i][j][k] -= log(factor);
            }
    curve->Compute_Determinant(s);
    curve->Compute_Inverse_Metric(s);
};
//================================================
//
// Remove trace
// 
//================================================
Manager::Remove_Trace(state* s) {
    curve->Compute_Trace(s);
    const double onethird = 1.0 / 3.0;
#pragma omp parallel for collapse(3)
    for (int i = 0; i < N_r; i++)
        for (int j = 0; j < N_t; j++)
            for (int k = 0; k < N_p; k++) {
                s->a_rr[i][j][k] -= onethird * (1.0 + s->h_rr(i, j, k)) * curve->trace_A(i, j, k);
                s->a_rt[i][j][k] -= onethird * (s->h_rt(i, j, k)) * curve->trace_A(i, j, k);
                s->a_rp[i][j][k] -= onethird * (s->h_rp(i, j, k)) * curve->trace_A(i, j, k);
                s->a_tt[i][j][k] -= onethird * (1.0 + s->h_tt(i, j, k)) * curve->trace_A(i, j, k);
                s->a_tp[i][j][k] -= onethird * (s->h_tp(i, j, k)) * curve->trace_A(i, j, k);
                s->a_pp[i][j][k] -= onethird * (1.0 + s->h_pp(i, j, k)) * curve->trace_A(i, j, k);
                // K[i][j][k] += curve->trace_A(i,j,k);
            }
    curve->Compute_Trace(s);
};


#endif  /* MANAGER_H */
