// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing all the integration stuff...
//================================================
#ifndef MANAGER_H
#define MANAGER_H

#include "Manager.h"
#include "GR.h"

//===========================================
// Constructor
//===========================================
Manager::Manager(int matter_type, int sigma_i, int cowling_i, double eta_i,
    int z4_i, double kappa_11_i, double kappa_12_i, double kappa_2_i,
    double kappa_ric_i, int RK_order_i,
    int char_OB_i, int solve_constraints_i) :
    sigma(sigma_i), cowling(cowling_i), eta(eta_i), z4(z4_i), kappa_11(kappa_11_i),
    kappa_12(kappa_12_i), kappa_2(kappa_2_i), kappa_ric(kappa_ric_i), RK_order(RK_order_i),
    char_OB(char_OB_i), solve_constraints(solve_constraints_i), r(GR::grid->N_r_tot()),
    r2(GR::grid->N_r_tot()), theta(GR::grid->N_theta_tot()), sintheta(GR::grid->N_theta_tot()),
    sin2theta(GR::grid->N_theta_tot()), costheta(GR::grid->N_theta_tot()), phi(GR::grid->N_phi_tot()) {
    //
    cout << " MANAGER: Constructing manager... " << endl;
    //
    // initial times and step counter (provided by CheckPoint.h, which
    // returns zero unless data are read in from checkpoint files)
    //
    t = GR::checkpoint->TStart();
    cout << " MANAGER: t = " << t << endl;
    tau_c = GR::checkpoint->TauStart();
    step = GR::checkpoint->TimeStepStart();
    GR::grid->Setup_Grid(r, r2, theta, sintheta, sin2theta, costheta, phi);
    //
    // Sanity check
    //
#ifdef EIGHTHORDERTHETA
    if (GR::grid->N_ghosts() < 4) {
        cerr << " MANAGER: Need N_ghosts >= 4 when using eighth-order differencing! " << endl;
        exit(0);
    }
#endif
    //
    // create state for dynamical variables
    //
    cout << " MANAGER: setting up states... " << endl;
    last = new state(GR::grid, GR::dump, "last");
    derivs = new state(GR::grid, GR::dump, "derivs");
    inter = new state(GR::grid, GR::dump, "inter");
    updates = new state(GR::grid, GR::dump, "updates");
    // create dump list for state last:
    last->assemble_dump_list("Dump_List");
    // let checkpointer know about last:
    GR::checkpoint->CollectDynVariables(last);

    // 
    // create curvature class
    // 
    curve = new curvature(GR::grid, GR::dump, "curve");
    curve->assemble_dump_list("Dump_List");
    //
    // create class with auxiliary functions
    //
    aux = new auxiliary(GR::grid, GR::dump, "auxiliaries");
    aux->assemble_dump_list("Dump_List");
    //
    // create horizon finder
    //
    horizonfinder = new HorizonFinder(GR::grid, last, curve, aux,
        GR::monitor->Filestem());
    //
    // create diagnostics class
    //
    cout << " MANAGER: setting up diagnostics... " << endl;
    constraints = new diagnostics(GR::grid, GR::dump, GR::monitor, GR::cosmology,
        horizonfinder, "constraints");
    constraints->assemble_dump_list("Dump_List");
    //
    // set time step
    //
    SetTimeStep();
    eta_KO = eta;
    cout << " MANAGER: using Kreiss-Oliger coefficient " << eta_KO << endl;
    GR::gauge->Set_eta(eta_KO);
    //
    // create matter
    //
    if (matter_type == 1) {
        matter = new Vacuum(GR::grid, GR::dump, GR::indata, cowling, GR::cosmology);
    } else if (matter_type == 3) {
        if (!GR::eos) {
            cerr << " MANAGER: Can't set up Hydro without eos! " << endl;
            exit(1);
        }
        matter = new Hydro(GR::grid, GR::dump, GR::indata, GR::eos, cowling, GR::cosmology,
            GR::monitor, horizonfinder, GR::checkpoint);
    } else if (matter_type == 5) {
        matter = new ScalarField(GR::grid, GR::dump, GR::indata, cowling, GR::cosmology,
            GR::monitor, eta_KO, GR::checkpoint);
    } else if (matter_type == 6) {
        matter = new RadHydro(GR::grid, GR::dump, GR::indata, GR::eos, cowling, GR::cosmology,
            GR::monitor, GR::checkpoint);
    } else if (matter_type == 7) {
        matter = new Maxwell(GR::grid, GR::dump, GR::indata, cowling, char_OB, GR::cosmology,
            GR::monitor, eta_KO, GR::checkpoint);
    } else if (matter_type == 8) {
        matter = new DualMaxwell(GR::grid, GR::dump, GR::indata, cowling, char_OB, GR::cosmology,
            GR::monitor, eta_KO, GR::checkpoint);
    } else {
        cerr << " MANAGER: Unknown matter type!!! " << endl;
    }
    //===========================================
    // set up profiles
    //===========================================
    if (!(GR::profiles == NULL)) {
        GR::profiles->add_to_profile_list("Profile_List", last);
        GR::profiles->add_to_profile_list("Profile_List", curve);
        GR::profiles->add_to_profile_list("Profile_List", aux);
        GR::profiles->add_to_profile_list("Profile_List", constraints);
        GR::profiles->initialize_file(constraints->R_prop.Address());
    }
    //
    //===========================================
    // initialize wave extraction
    //===========================================
    if (!(GR::waves == NULL)) {
        GR::waves->Initialize(constraints->psi4_Re.Address(),
            constraints->psi4_Im.Address());
    }
    N_g = inter->GR::grid->N_ghosts();
    N_r = inter->GR::grid->N_r_tot();
    N_t = inter->GR::grid->N_theta_tot();
    N_p = inter->GR::grid->N_phi_tot();
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
    int N_g = GR::grid->N_ghosts();
    dt = GR::grid->courant_factor() * 0.5 * GR : grid->delta_r(N_g)
        * GR : grid->delta_theta(N_g);
    cout << " MANAGER: using dt = " << dt << endl;
};
Manager::SetTimeStep(state* s) {
    int N_g = GR::grid->N_ghosts();
    const double gamma_tt = (1.0 + s->h_tt(N_g, N_g, N_g)) * exp(4.0 * s->phi(N_g, N_g, N_g));
    // const double factor = sqrt(gamma_tt) / s->lapse(N_g, N_g, N_g); // + abs(s->shift_r(N_g, N_g, N_g));
    const double factor = sqrt(gamma_tt); // + abs(s->shift_r(N_g, N_g, N_g));
    dt = GR::grid->courant_factor() * 0.5 * GR::grid->delta_r(N_g) * GR::grid->delta_theta(N_g) * factor;
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
