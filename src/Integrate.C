//================================================
//
// Integrate fields, routine of manager.h
//
//================================================
#include "Manager.h"

bool Manager::Integrate(double t_max) {
  while (t < t_max) {
    // SetTimeStep(last);
    // for self-similar shift update shift here:
    if (gauge->GaugeType() == self_sim) {
      Compute_RHS(t);
      gauge->overwrite_shift(last, derivs, t);
    };

    //================================================
    // carry out one time step
    //================================================
    if (RK_order == 4) Step_RK4(t, tau_c);
    else if (RK_order == 3) Step_RK3(t, tau_c);
    else if (RK_order == 2) Step_ICN(t, tau_c);
    else {
      cerr << " Runge Kutta integrator of order " << RK_order
	   << " is not supported! " << endl;
      return false;
    }
    step++; 
    //================================================
    // check for regridding
    //================================================
    Regrid(t, tau_c, step, t_max);
    //================================================
    // update photons
    //================================================
    //      inter->equals(last);
    curve->Update(last);
    if (!(photons == NULL)) {
      aux->Compute_r_derivs(last, curve);
      photons->Update_Photons(t, tau_c, dt, last, curve, aux, 
			      constraints->I_Re.Address(),
			      constraints->J_Re.Address(),
			      grid->r_max());
    }    
    //================================================
    // extract waves
    //================================================
    if (waves != NULL) {
      double I_Re_max = constraints->CurvatureInvariant(last,curve,
							matter->adm_sources,
							t, tau_c, step);
      waves->Update(dt);
    }
    //================================================
    // check whether it's time to write check point
    //================================================ 
    if (step % checkpoint->CheckPointStep() == 0)
      checkpoint->WriteCheckPoint(step, t, tau_c);
    //================================================
    // check whether it's time to dump
    //================================================ 
    // (mass needed for horizons below)
    int i = 7.5 / 8. * (N_r - N_g);
    double mass = constraints->ADM_Mass_Surface(last,curve,aux,i);
    if (dump->time_to_dump(step)) {
      last->dump_fcts(t, tau_c, step);
      curve->dump_fcts(t, tau_c, step);
      aux->dump_fcts(t, tau_c, step);
      //================================================
      // Compute matter diagnostics
      //================================================ 
      constraints->Compute_Proper_Radius(last, curve, aux, sigma);
      matter->Compute_Diagnostics(last,curve);
      matter->dump_fcts(t, tau_c, step);
      //================================================
      // Also force a search for horizon, so that we have
      // horizon data at the same time
      //================================================ 
      bool force = true;
      double lin_mom_guess = 0.0;
      constraints->FindHorizon(step, t, tau_c, last, curve,
			       matter->adm_sources, matter->fluxes,
			       aux, mass,  lin_mom_guess, force);
      //================================================
      // evaluate constraints
      //================================================ 
      double Ham_norm, Ham_norm_ex;
      double Mom_r_norm, Mom_t_norm, Mom_p_norm;
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
      constraints->dump_fcts(t, tau_c, step);
      monitor->note_constraints(step, t, tau_c, Ham_norm, Ham_norm_ex,
				Mom_r_norm, Mom_t_norm, Mom_p_norm, 
				CFC_r_norm, CFC_t_norm, CFC_p_norm);
      //================================================
      // evaluate curvature invariants
      //================================================ 
      if (waves == NULL) // otherwise they are computed above already...
	double I_Re_max = constraints->CurvatureInvariant(last,curve,
							  matter->adm_sources,
							  t, tau_c, step);
      constraints->Note_Invariants(step, t, tau_c, monitor);      
    }
    //================================================
    // check whether it's time to note...
    //================================================
    if (monitor->time_to_note(step)) {
      double ang_mom = constraints->Angular_Momentum(last,curve,aux,i);
      double lin_mom = constraints->Linear_Momentum(last,curve,aux,i);
      monitor->note(step, t, tau_c, mass, ang_mom, 
		    lin_mom, last->phi(0.0, N_g, N_g), 
		    last->lapse(0.0, N_g, N_g), last->lapse.min(),
		    last->K(0.0, N_g, N_g), RegridCriterion());
      //================================================
      // evaluate constraints
      //================================================ 
      double Ham_norm, Ham_norm_ex;
      double Mom_r_norm, Mom_t_norm, Mom_p_norm;
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
      monitor->note_constraints(step, t, tau_c, Ham_norm, Ham_norm_ex,
				Mom_r_norm, Mom_t_norm, Mom_p_norm, 
				CFC_r_norm, CFC_t_norm, CFC_p_norm);
      //================================================
      // evaluate curvature invariants
      //================================================ 
      double I_Re_max = constraints->CurvatureInvariant(last,curve,
							matter->adm_sources, t, tau_c, step);
      constraints->Note_Invariants(step, t, tau_c, monitor);
      constraints->Compute_Proper_Radius(last, curve, aux, sigma);
      if (profiles != NULL) profiles->write_profile(t, tau_c);
      if (photons != NULL) photons->Monitor(t, tau_c);
      if (waves != NULL) waves->write_wave_data(t, tau_c);
      matter->Compute_Diagnostics(last,curve);
      matter->Note(step, t, tau_c);
    }
    //================================================
    // look for horizons
    //================================================
    constraints->FindHorizon(step, t, tau_c, last, curve,
			     matter->adm_sources, matter->fluxes,
			     aux, mass);
    //================================================
    // Finally: check for NaN's...
    //================================================
    if (!last->FINITE() || !aux->FINITE() || !curve->FINITE() ) {
      last->dump_fcts(t, tau_c, step, "_crash");
      curve->dump_fcts(t, tau_c, step, "_crash");
      aux->dump_fcts(t, tau_c, step, "_crash");
      constraints->dump_fcts(t, tau_c, step, "_crash");
      matter->dump_fcts(t, tau_c, step, "_crash");
      return false;
    }
  }
  return true;
}

//===========================================
// Carry out one RK4 time step
//===========================================
void Manager::Step_RK4(double & t, double & tau_c) {
  //
  // before getting started...
  //
  inter->equals(last);
  matter->Start_RK();  // sets matter_inter = matter_last
  //
  // evaluate k_1 at current time
  // 
  Compute_RHS(t);              // computes derivs from inter
  matter->Compute_RHS(inter, curve, t);
  updates->add(last, dt / 6.0, derivs);
  matter->Update(dt / 6.0);
  tau_c += inter->lapse(0.0,N_g,N_g) * dt / 6.0;
  //
  // evaluate k_2 at middle point
  // 
  inter->add(last, 0.5 * dt, derivs);
  inter->fill_ghosts();
  if (char_OB == 1) inter->char_OB(last, 0.5 *dt);
  curve->Update(inter);
  matter->Compute_inter(0.5 * dt);
  matter->ADM_Sources(inter, curve);
  Compute_RHS(t + 0.5 * dt);    // always compute derivs from inter
  matter->Compute_RHS(inter, curve, t + 0.5 * dt);
  updates->add(dt / 3.0, derivs);
  matter->Update(dt / 3.0);
  tau_c += inter->lapse(0.0,N_g,N_g) * dt / 3.0;
  //
  // evaluate k_3 at middle point
  // 
  inter->add(last, 0.5 * dt, derivs);
  inter->fill_ghosts();
  if (char_OB == 1) inter->char_OB(last, 0.5 *dt);
  curve->Update(inter);
  matter->Compute_inter(0.5 * dt);
  matter->ADM_Sources(inter, curve);
  Compute_RHS(t + 0.5 * dt);    // always compute derivs from inter
  matter->Compute_RHS(inter, curve, t + 0.5 * dt);
  updates->add(dt / 3.0, derivs);
  matter->Update(dt / 3.0);
  tau_c += inter->lapse(0.0,N_g,N_g) * dt / 3.0;
  //
  // evaluate k_4 at last point
  //
  inter->add(last, dt, derivs);
  inter->fill_ghosts();
  if (char_OB == 1) inter->char_OB(last, dt);
  curve->Update(inter);
  matter->Compute_inter(dt);
  matter->ADM_Sources(inter, curve);
  Compute_RHS(t + dt);    // always compute derivs from inter
  matter->Compute_RHS(inter, curve, t + dt);    
  updates->add(dt / 6.0, derivs);
  matter->Update(dt / 6.0);
  tau_c += inter->lapse(0.0,N_g,N_g) * dt / 6.0;
  //
  //
  // if using Lagrangian formulation, rescale metric
  //  
  if (sigma == 1) {
    Rescale_Metric(updates);
    Remove_Trace(updates);
  }
  if (char_OB == 1) updates->char_OB(last, dt);
  last->equals(updates);
  last->fill_ghosts();
  curve->Update(last);
  matter->Finish_RK(last, curve, dt);
  //
  t += dt;
}


//===========================================
// Carry out one RK3 time step (Ralston's third-order method)
// (see, e.g., Appendix A in Calabrese, Hinder & Husa, gr-qc/0503056)
//===========================================
void Manager::Step_RK3(double & t, double & tau_c) {
  //
  // before getting started...
  //
  inter->equals(last);
  matter->Start_RK();  // sets matter_inter = matter_last
  //
  // evaluate k_1 at current time
  // 
  Compute_RHS(t);              // computes derivs from inter
  matter->Compute_RHS(inter, curve, t);
  updates->add(last, 2.0 * dt / 9.0, derivs);
  matter->Update(2.0 * dt / 9.0);
  tau_c += inter->lapse(0.0,N_g,N_g) * 2.0 * dt / 9.0;
  //
  // evaluate k_2 at middle point
  // 
  inter->add(last, 0.5 * dt, derivs);
  inter->fill_ghosts();
  if (char_OB == 1) inter->char_OB(last, 0.5 *dt);
  curve->Update(inter);
  matter->Compute_inter(0.5 * dt);
  matter->ADM_Sources(inter, curve);
  Compute_RHS(t + 0.5 * dt);    // always compute derivs from inter
  matter->Compute_RHS(inter, curve, t + 0.5 * dt);
  updates->add(dt / 3.0, derivs);
  matter->Update(dt / 3.0);
  tau_c += inter->lapse(0.0,N_g,N_g) * dt / 3.0;
  //
  // evaluate k_3 at point 3/4 down...
  // 
  inter->add(last, 0.75 * dt, derivs);
  inter->fill_ghosts();
  if (char_OB == 1) inter->char_OB(last, 0.75 *dt);
  curve->Update(inter);
  matter->Compute_inter(0.75 * dt);
  matter->ADM_Sources(inter, curve);
  Compute_RHS(t + 0.75 * dt);    // always compute derivs from inter
  matter->Compute_RHS(inter, curve, t + 0.75 * dt);
  updates->add(4.0 * dt / 9.0, derivs);
  matter->Update(4.0 * dt / 9.0);
  tau_c += inter->lapse(0.0,N_g,N_g) * 4.0 * dt / 9.0;
  //
  //
  // if using Lagrangian formulation, rescale metric
  //  
  if (sigma == 1) {
    Rescale_Metric(last);
    Remove_Trace(last);
  }
  if (char_OB == 1) updates->char_OB(last, dt);
  last->equals(updates);
  last->fill_ghosts();
  curve->Update(last);
  matter->Finish_RK(last, curve, dt);
  //
  t += dt;
}


//===========================================
// Carry out one iterative Crank-Nicholson time step 
// (see, e.g., Teukolsky, PRD 61, 087501 (2000)
//===========================================
void Manager::Step_ICN(double & t, double & tau_c) {
  //
  // before getting started...
  //
  inter->equals(last);
  matter->Start_RK();  // sets matter_inter = matter_last
  //
  // evaluate k_1 at current time
  // 
  Compute_RHS(t);              // computes derivs from inter
  matter->Compute_RHS(inter, curve, t);
  updates->add(last, 0.5 * dt, derivs);
  matter->Update(0.5 * dt);
  tau_c += inter->lapse(0.0,N_g,N_g) * 0.5 * dt;
  //
  // evaluate k_2 at end point, but don't update
  // 
  inter->add(last, dt, derivs);
  inter->fill_ghosts();
  if (char_OB == 1) inter->char_OB(last, dt);
  curve->Update(inter);
  matter->Compute_inter(dt);
  matter->ADM_Sources(inter, curve);
  Compute_RHS(t + dt);    // always compute derivs from inter
  matter->Compute_RHS(inter, curve, t + dt);
  //  updates->add(dt / 3.0, derivs);
  //  matter->Update(dt / 3.0);
  //  tau_c += inter->lapse(0.0,N_g,N_g) * dt / 3.0;
  //
  // evaluate k_3 at end point again
  // 
  inter->add(last, dt, derivs);
  inter->fill_ghosts();
  if (char_OB == 1) inter->char_OB(last, dt);
  curve->Update(inter);
  matter->Compute_inter(dt);
  matter->ADM_Sources(inter, curve);
  Compute_RHS(t + dt);    // always compute derivs from inter
  matter->Compute_RHS(inter, curve, t + dt);
  updates->add(0.5 * dt, derivs);
  matter->Update(0.5 * dt);
  tau_c += inter->lapse(0.0,N_g,N_g) * 0.5 * dt;
  //
  //
  // if using Lagrangian formulation, rescale metric
  //  
  if (sigma == 1) {
    Rescale_Metric(last);
    Remove_Trace(last);
  }
  if (char_OB == 1) updates->char_OB(last, dt);
  last->equals(updates);
  last->fill_ghosts();
  curve->Update(last);
  matter->Finish_RK(last, curve, dt);
  //
  t += dt;
}
