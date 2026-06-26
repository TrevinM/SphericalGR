//================================================
//
// Initialize fields, routine of manager.h
//
// Note: will initialize State last
//================================================
#include "Manager.h"

bool Manager::Initialize() {
  //
  // initialize gravitational fields
  // 
  cout << " INITIALIZE: t = " << t << endl;
  bool success = indata->Initialize_Metric(last);
  if (!success) return false;
  indata->Initialize_Lapse(last);
  indata->Initialize_Shift(last);
  indata->Initialize_ExCurvature(last);
  
  curve->Update(last);
  indata->Initialize_ConfConn(last, curve);
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
  matter->Initialize(last, curve, constraints);
  matter->ADM_Sources(last, curve);
  //
  // if desired, solve constraints
  // 
  if (solve_constraints) {
    cout << " INITIALIZE: Hamiltonian constraint before solving: " 
	 << constraints->Hamiltonian(last,curve,aux,matter->adm_sources)
	 << endl;
    ConstraintSolver * constraint_solver;
    constraint_solver = new ConstraintSolver(grid,last,curve,matter);
    constraint_solver->SolveHamiltonian();
    delete constraint_solver;
    cout << " INITIALIZE: Hamiltonian constraint after solving: " 
	 << constraints->Hamiltonian(last,curve,aux,matter->adm_sources) 
	 << endl;
  }
  //
  // find mass and angular momentum
  //
  int i = 7.5 / 8. * (N_r - N_g);
  double mass = constraints->ADM_Mass_Surface(last,curve,aux,i);
  double ang_mom = constraints->Angular_Momentum(last,curve,aux,i);
  double lin_mom = constraints->Linear_Momentum(last,curve,aux,i);
  cout << " MANAGER: Found M = " << mass << ", J = " << ang_mom
       << " and P_z = " << lin_mom 
       << " at radius " << grid->r(i) << endl;
  monitor->note(step, t, tau_c, mass, ang_mom, 
		lin_mom, last->phi(0.0, N_g, N_g), 
		last->lapse(0.0, N_g, N_g), last->lapse.min(), 
		last->K(0.0, N_g, N_g), RegridCriterion());
  //
  // look for horizons
  //
  bool force = false;
  constraints->FindHorizon(step, t, tau_c, last, curve,
			   matter->adm_sources, matter->fluxes,
			   aux, mass, lin_mom, force);
  matter->Note(step, t, tau_c);
  //
  // evaluate constraints
  //
  double Ham_norm, Ham_norm_ex, Mom_r_norm, Mom_t_norm, Mom_p_norm;
  double CFC_r_norm, CFC_t_norm, CFC_p_norm;
  Ham_norm_ex = 
    constraints->Hamiltonian(last,curve,aux,matter->adm_sources, true);
  Ham_norm = constraints->Hamiltonian(last,curve,aux,matter->adm_sources);
  constraints->MomentumConstraint(last,curve,inter,
				  matter->adm_sources,
				  Mom_r_norm, Mom_t_norm, 
				  Mom_p_norm);
  constraints->ConnectionFunctionConstraint(last,curve,
					    CFC_r_norm,
					    CFC_t_norm,
					    CFC_p_norm);
  constraints->Compute_Proper_Radius(last);
  monitor->note_constraints(step, t, tau_c, Ham_norm, Ham_norm_ex,
			    Mom_r_norm, Mom_t_norm, Mom_p_norm, 
			    CFC_r_norm, CFC_t_norm, CFC_p_norm);
  //
  // compute curvature invariants
  //
  double I_max = constraints->CurvatureInvariant(last,curve,
						 matter->adm_sources, t, tau_c, step);
  constraints->Note_Invariants(step, t, tau_c, monitor);
  if (waves != NULL) {
    waves->Update(0.0);   // needed in order to project psi4 data...
    waves->write_wave_data(t, tau_c);
  }
  matter->Compute_Diagnostics(last,curve);
  //
  // check whether shift needs to be computed (for self-similar shift)
  //
  // if (gauge->GaugeType() == self_sim) {
  //    Compute_RHS(0.0);
  //    gauge->overwrite_shift(last, derivs, 0.0);
  //  };
  //
  // now dump initial data
  //
  last->dump_fcts(t, tau_c, step);
  curve->dump_fcts(t, tau_c, step);
  aux->dump_fcts(t, tau_c, step);
  matter->dump_fcts(t, tau_c, step);
  constraints->dump_fcts(t, tau_c, step);
  if (!(profiles == NULL)) {
    profiles->write_profile(t, tau_c);
  }
  //
  return success;
}
