// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing all the integration stuff...
//================================================

#include "Manager.h"
#include "Container.h"

//===========================================
// Constructor
//===========================================
state* Manager::last = nullptr;
state* Manager::derivs = nullptr;
state* Manager::inter = nullptr;
state* Manager::updates = nullptr;
curvature* Manager::curve = nullptr;
auxiliary* Manager::aux = nullptr;
HorizonFinder* Manager::horizonfinder = nullptr;
diagnostics* Manager::constraints = nullptr;
Matter* Manager::matter = nullptr;


Manager::Manager(int matter_type, int sigma_i, int cowling_i, double eta_i,
    int z4_i, double kappa_11_i, double kappa_12_i, double kappa_2_i,
    double kappa_ric_i, int RK_order_i,
    int char_OB_i, int solve_constraints_i) :
    sigma(sigma_i), cowling(cowling_i), eta(eta_i), z4(z4_i), kappa_11(kappa_11_i),
    kappa_12(kappa_12_i), kappa_2(kappa_2_i), kappa_ric(kappa_ric_i), RK_order(RK_order_i),
    char_OB(char_OB_i), solve_constraints(solve_constraints_i), r(Container::grid->N_r_tot()),
    r2(Container::grid->N_r_tot()), theta(Container::grid->N_theta_tot()), sintheta(Container::grid->N_theta_tot()),
    sin2theta(Container::grid->N_theta_tot()), costheta(Container::grid->N_theta_tot()), phi(Container::grid->N_phi_tot()) {
    //
    cout << " MANAGER: Constructing manager... " << endl;
    //
    // initial times and step counter (provided by CheckPoint.h, which
    // returns zero unless data are read in from checkpoint files)
    //
    t = Container::checkpoint->TStart();
    cout << " MANAGER: t = " << t << endl;
    tau_c = Container::checkpoint->TauStart();
    step = Container::checkpoint->TimeStepStart();
    Container::grid->Setup_Grid(r, r2, theta, sintheta, sin2theta, costheta, phi);
    //
    // Sanity check
    //
#ifdef EIGHTHORDERTHETA
    if (Container::grid->N_ghosts() < 4) {
        cerr << " MANAGER: Need N_ghosts >= 4 when using eighth-order differencing! " << endl;
        exit(0);
    }
#endif
    //
    // create state for dynamical variables
    //
    cout << " MANAGER: setting up states... " << endl;
    last = new state(Container::grid, Container::dump, "last");
    derivs = new state(Container::grid, Container::dump, "derivs");
    inter = new state(Container::grid, Container::dump, "inter");
    updates = new state(Container::grid, Container::dump, "updates");
    // create dump list for state last:
    last->assemble_dump_list("Dump_List");
    // let checkpointer know about last:
    Container::checkpoint->CollectDynVariables(last);

    // 
    // create curvature class
    // 
    curve = new curvature(Container::grid, Container::dump, "curve");
    curve->assemble_dump_list("Dump_List");
    //
    // create class with auxiliary functions
    //
    aux = new auxiliary(Container::grid, Container::dump, "auxiliaries");
    aux->assemble_dump_list("Dump_List");
    //
    // create horizon finder
    //
    horizonfinder = new HorizonFinder(Container::grid, last, curve, aux,
        Container::monitor->Filestem());
    //
    // create diagnostics class
    //
    cout << " MANAGER: setting up diagnostics... " << endl;
    constraints = new diagnostics(Container::grid, Container::dump, Container::monitor, Container::cosmology,
        horizonfinder, "constraints");
    constraints->assemble_dump_list("Dump_List");
    //
    // set time step
    //
    SetTimeStep();
    eta_KO = eta;
    cout << " MANAGER: using Kreiss-Oliger coefficient " << eta_KO << endl;
    Container::gauge->Set_eta(eta_KO);
    //
    // create matter
    //
    if (matter_type == 1) {
        matter = new Vacuum(Container::grid, Container::dump, Container::indata, cowling, Container::cosmology);
    } else if (matter_type == 3) {
        if (!Container::eos) {
            cerr << " MANAGER: Can't set up Hydro without eos! " << endl;
            exit(1);
        }
        matter = new Hydro(Container::grid, Container::dump, Container::indata, Container::eos, cowling, Container::cosmology,
            Container::monitor, horizonfinder, Container::checkpoint);
    } else if (matter_type == 5) {
        matter = new ScalarField(Container::grid, Container::dump, Container::indata, cowling, Container::cosmology,
            Container::monitor, eta_KO, Container::checkpoint);
    } else if (matter_type == 6) {
        matter = new RadHydro(Container::grid, Container::dump, Container::indata, Container::eos, cowling, Container::cosmology,
            Container::monitor, Container::checkpoint);
    } else if (matter_type == 7) {
        matter = new Maxwell(Container::grid, Container::dump, Container::indata, cowling, char_OB, Container::cosmology,
            Container::monitor, eta_KO, Container::checkpoint);
    } else if (matter_type == 8) {
        matter = new DualMaxwell(Container::grid, Container::dump, Container::indata, cowling, char_OB, Container::cosmology,
            Container::monitor, eta_KO, Container::checkpoint);
    } else {
        cerr << " MANAGER: Unknown matter type!!! " << endl;
    }
    //===========================================
    // set up profiles
    //===========================================
    if (!(Container::profiles == NULL)) {
        Container::profiles->add_to_profile_list("Profile_List", last);
        Container::profiles->add_to_profile_list("Profile_List", curve);
        Container::profiles->add_to_profile_list("Profile_List", aux);
        Container::profiles->add_to_profile_list("Profile_List", constraints);
        Container::profiles->initialize_file(constraints->R_prop.Address());
    }
    //
    //===========================================
    // initialize wave extraction
    //===========================================
    if (!(Container::waves == NULL)) {
        Container::waves->Initialize(constraints->psi4_Re.Address(),
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
void Manager::SetTimeStep() {
    int N_g = Container::grid->N_ghosts();
    dt = Container::grid->courant_factor() * 0.5 * Container::grid->delta_r(N_g)
        * Container::grid->delta_theta(N_g);
    cout << " MANAGER: using dt = " << dt << endl;
};
void Manager::SetTimeStep(state* s) {
    int N_g = Container::grid->N_ghosts();
    const double gamma_tt = (1.0 + s->h_tt(N_g, N_g, N_g)) * exp(4.0 * s->phi(N_g, N_g, N_g));
    // const double factor = sqrt(gamma_tt) / s->lapse(N_g, N_g, N_g); // + abs(s->shift_r(N_g, N_g, N_g));
    const double factor = sqrt(gamma_tt); // + abs(s->shift_r(N_g, N_g, N_g));
    dt = Container::grid->courant_factor() * 0.5 * Container::grid->delta_r(N_g) * Container::grid->delta_theta(N_g) * factor;
};

//===========================================
// Initialize
//===========================================
void Manager::Set_t_max(double set_t) { t_max = set_t; }

//================================================
//
// Rescale metric
// 
//================================================
void Manager::Rescale_Metric(state* s) {
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
void Manager::Remove_Trace(state* s) {
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

//===========================================
// Destructor
//===========================================
Manager::~Manager() {
    delete last;
    delete derivs;
    delete inter;
    delete updates;
    delete curve;
    delete aux;
    delete constraints;
    delete matter;
    delete horizonfinder;
    cout << " MANAGER: Destructing manager - bye! " << endl;
}