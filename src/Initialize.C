//================================================
//
// Initialize fields, routine of manager.h
//
// Note: will initialize State last
//================================================
#include "Manager.h"
#include "Container.h"

bool Manager::Initialize() {
    //
    // initialize gravitational fields
    // 
    cout << " INITIALIZE: t = " << t << endl;
    bool success = Container::indata->Initialize_Metric(last);
    if (!success) return false;
    Container::indata->Initialize_Lapse(last);
    Container::indata->Initialize_Shift(last);
    Container::indata->Initialize_ExCurvature(last);

    curve->Update(last);
    Container::indata->Initialize_ConfConn(last, curve);
    // if (indata->Type() == brill || indata->Type() == kerr || indata->Type() == linwave) {
    // //  if (indata->Type() == brill || indata->Type() == kerr ) {
    //   indata->Initialize_ConfConn(last, curve);
    // } else {
    //   indata->Initialize_ConfConn(last); 
    // }
    //  last->lapse.slicemax(3);
    //
    // evaluate curvature and remember initial determinant 
    // CHECK: can we do that above instead, before calling Initialize_ConfConn?
    curve->Update(last);
    curve->InitDet();
    //
    // initialize matter
    //
    Container::matter->Initialize(last, curve, constraints);
    Container::matter->ADM_Sources(last, curve);
    //
    // if desired, solve constraints
    // 
    if (solve_constraints) {
        cout << " INITIALIZE: Hamiltonian constraint before solving: "
            << constraints->Hamiltonian(last, curve, aux, Container::matter->adm_sources)
            << endl;
        ConstraintSolver* constraint_solver;
        constraint_solver = new ConstraintSolver(Container::grid, last, curve, Container::matter);
        constraint_solver->SolveHamiltonian();
        delete constraint_solver;
        cout << " INITIALIZE: Hamiltonian constraint after solving: "
            << constraints->Hamiltonian(last, curve, aux, Container::matter->adm_sources)
            << endl;
    }
    //
    // find mass and angular momentum
    //
    int i = 7.5 / 8. * (N_r - N_g);
    double mass = constraints->ADM_Mass_Surface(last, curve, aux, i);
    double ang_mom = constraints->Angular_Momentum(last, curve, aux, i);
    double lin_mom = constraints->Linear_Momentum(last, curve, aux, i);
    cout << " MANAGER: Found M = " << mass << ", J = " << ang_mom
        << " and P_z = " << lin_mom
        << " at radius " << Container::grid->r(i) << endl;
    Container::monitor->note(step, t, tau_c, mass, ang_mom,
        lin_mom, last->phi(0.0, N_g, N_g),
        last->lapse(0.0, N_g, N_g), last->lapse.min(),
        last->K(0.0, N_g, N_g), Container::grid->RegridCriterion());
    //
    // look for horizons
    //
    bool force = false;
    constraints->FindHorizon(step, t, tau_c, last, curve,
        Container::matter->adm_sources, Container::matter->fluxes,
        aux, mass, lin_mom, force);
    Container::matter->Note(step, t, tau_c);
    //
    // evaluate constraints
    //
    double Ham_norm, Ham_norm_ex, Mom_r_norm, Mom_t_norm, Mom_p_norm;
    double CFC_r_norm, CFC_t_norm, CFC_p_norm;
    Ham_norm_ex =
        constraints->Hamiltonian(last, curve, aux, Container::matter->adm_sources, true);
    Ham_norm = constraints->Hamiltonian(last, curve, aux, Container::matter->adm_sources);
    constraints->MomentumConstraint(last, curve, inter,
        Container::matter->adm_sources,
        Mom_r_norm, Mom_t_norm,
        Mom_p_norm);
    constraints->ConnectionFunctionConstraint(last, curve,
        CFC_r_norm,
        CFC_t_norm,
        CFC_p_norm);
    constraints->Compute_Proper_Radius(last, curve, aux, sigma);
    Container::monitor->note_constraints(step, t, tau_c, Ham_norm, Ham_norm_ex,
        Mom_r_norm, Mom_t_norm, Mom_p_norm,
        CFC_r_norm, CFC_t_norm, CFC_p_norm);
    //
    // compute curvature invariants
    //
    double I_max = constraints->CurvatureInvariant(last, curve,
        Container::matter->adm_sources, t, tau_c, step);
    constraints->Note_Invariants(step, t, tau_c, Container::monitor);
    if (Container::waves != NULL) {
        Container::waves->Update(0.0);   // needed in order to project psi4 data...
        Container::waves->write_wave_data(t, tau_c);
    }
    Container::matter->Compute_Diagnostics(last, curve);
    //
    // check whether shift needs to be computed (for self-similar shift)
    //
    // if (gauge->GaugeType() == self_sim) {
    //    Compute_RHS(0.0);
    //    Container::gauge->overwrite_shift(last, derivs, 0.0);
    //  };
    //
    // now dump initial data
    //
    last->dump_fcts(t, tau_c, step);
    curve->dump_fcts(t, tau_c, step);
    aux->dump_fcts(t, tau_c, step);
    Container::matter->dump_fcts(t, tau_c, step);
    constraints->dump_fcts(t, tau_c, step);
    if (!(Container::profiles == NULL)) {
        Container::profiles->write_profile(t, tau_c);
    }
    //
    return success;
}
