//================================================
//
// Member function of Manager.h that handles regridding
//
//================================================
#include "Manager.h"

double Manager::Regrid(double t, double tau_c, int timestep, double & t_max) {
  //  cout << " Manager: in Regrid()... " << endl;
  //
  // check whether it's time to regrid
  // 
  double criterion = RegridCriterion();
  bool regrid = grid->TimeToRegrid(criterion);
  if (timestep - timestep_last_regrid < steps_between_regrids) 
    regrid = false;
  //
  // if so, go ahead...
  //
  if (regrid) {
    //
    // take note before...
    //
    int i = 7.5 / 8. * (N_r - N_g);
    double mass = constraints->ADM_Mass_Surface(last,curve,aux,i);
    double ang_mom = constraints->Angular_Momentum(last,curve,aux,i);
    double lin_mom = constraints->Linear_Momentum(last,curve,aux,i);
    monitor->note(step, t, tau_c, mass, ang_mom, 
		  lin_mom, last->phi(0.0, N_g, N_g), 
		  last->lapse(0.0, N_g, N_g), last->lapse.min(),
		  last->K(0.0, N_g, N_g), RegridCriterion(), true);
    //
    timestep_last_regrid = timestep;
    //
    //    cout << " MANAGER: regridding..." << endl;
    last->dump_fcts(t,tau_c,timestep,"_before_regrid");
    matter->dump_fcts(t,tau_c,timestep, "_before_regrid");
    constraints->dump_fcts(t,tau_c,timestep, "_before_regrid");
    //
    // Compute_Inverse_Metric(h_rr_o, h_rt_o, h_rp_o,h_tt_o, h_tp_o, h_pp_o);
    // Compute_Metric_Derivatives(h_rr_o, h_rt_o, h_rp_o,h_tt_o, h_tp_o, h_pp_o);
    // Compute_Connection();
    // Compute_Ricci(h_rr_o, h_rt_o, h_rp_o,h_tt_o, h_tp_o, h_pp_o, lam_r_o,
    // 		  lam_t_o, lam_p_o);
    // Compute_Trace(a_rr_o,a_rt_o,a_rp_o,a_tt_o,a_tp_o,a_pp_o);
    // cout << " REGRID: Hamiltonian before regrid = " << Hamiltonian() << endl;
    // dump->cond_dump(t,tau_c,timestep,&Ham,"_before_regrid");
    //
    // NOTE: don't overwrite old r's yet: need these for 
    // interpolations of gridfunctions (done in gridfunction.Regrid() ).
    // Overwrite old r's with new r's only once all gridfunctions have
    // been interpolated to new grid
    //
    VecDoub r_new(N_r);
    //
    // returns new grid points only -- need to finish up regrid with call
    // to Setup_Radial_Grid() below.
    //
    bool regrid = grid->Regrid(r_new, criterion);
    //
    // now make sure to regrid every dynamical variable
    //
    last->Regrid(r_new);
    curve->det_init.Regrid(r_new);
    // 
    // regrid matter variables
    //
    matter->Regrid(r_new);
    //
    // update regrid ADM matter sources
    //
    // rho_ADM.Regrid(r_new);
    // S_r.Regrid(r_new);
    // S_t.Regrid(r_new);
    // S_p.Regrid(r_new);
    // S_rr.Regrid(r_new);
    // S_rt.Regrid(r_new);
    // S_rp.Regrid(r_new);
    // S_tt.Regrid(r_new);
    // S_tp.Regrid(r_new);
    // S_pp.Regrid(r_new);
    // S.Regrid(r_new);
    //
    // finally update radius
    // 
    grid->Setup_Radial_Grid(r, r2);
    //
    last->dump_fcts(t,tau_c,timestep,"_after_regrid");
    matter->dump_fcts(t,tau_c,timestep,"_after_regrid");
    constraints->dump_fcts(t,tau_c,timestep,"_after_regrid");
    monitor->regrid(grid->r_max());
    horizonfinder->SetRMaxMin();
    SetTimeStep();
    //
    // compute curvature quantities
    //
    // Compute_Inverse_Metric(h_rr_o, h_rt_o, h_rp_o,h_tt_o, h_tp_o, h_pp_o);
    // Compute_Metric_Derivatives(h_rr_o, h_rt_o, h_rp_o,h_tt_o, h_tp_o, h_pp_o);
    // Compute_Connection();
    // Compute_Ricci(h_rr_o, h_rt_o, h_rp_o,h_tt_o, h_tp_o, h_pp_o, lam_r_o,
    // 		  lam_t_o, lam_p_o);
    //
    // Also, just to compute functions A2 and Dphi that are used in Hamiltonian,
    // compute some right-hand sides of evolution equations
    // NOTE: we'll assume that a=1
    // CHECK: generalize for a \ne 1!
    // 
    //
    //
    criterion = RegridCriterion();
    //
    // finally adjust t_max if necessary
    //
    if (t + grid->r_max() < t_max) {
      t_max = t + grid->r_max();
      cout << " REGRID: set t_max to " << t_max << endl;
      monitor->set_t_max(t_max);
    }
    //
    // ... and take note right after regrid...
    //
    mass = constraints->ADM_Mass_Surface(last,curve,aux,i);
    ang_mom = constraints->Angular_Momentum(last,curve,aux,i);
    lin_mom = constraints->Linear_Momentum(last,curve,aux,i);
    monitor->note(step, t, tau_c, mass, ang_mom, 
		  lin_mom, last->phi(0.0, N_g, N_g), 
		  last->lapse(0.0, N_g, N_g), last->lapse.min(),
		  last->K(0.0, N_g, N_g), criterion, true);
  }
  return criterion;
}

//================================================
// evaluates some criterion so that if this number becomes larger than 
// the cutoff specified in Grid_Input, we regrid
//================================================
double Manager::RegridCriterion() {
  if (grid->Regrid_Type() == 0) {
    if (!strcmp(matter->Name(),"vacuum")) {  // vacuum...
      double diff = 0.0; 
      double max_diff = 0.0;
      // first compute maximum of I_Re
      double I_Re_max = 0.0;
      // for (int i = N_g; i < N_r - N_g; i++) 
      //   for (int j = N_g; j < N_t - N_g; j++) 
      // 	for (int k = N_g; k < N_p - N_g; k++) {
      // 	  if (fabs(constraints->I_Re(i,j,k)) > I_Re_max)
      // 	    I_Re_max = fabs(constraints->I_Re(i,j,k));
      // 	}
      //
      for (int i = N_g; i < N_r - N_g; i++) 
        for (int j = N_g; j < N_t - N_g; j++) 
          for (int k = N_g; k < N_p - N_g; k++) {
            diff = last->lapse(i+1,j,k) - last->lapse(i,j,k);
            // diff = (constraints->I_Re(i+1,j,k) - constraints->I_Re(i,j,k)) / I_Re_max;
            if (abs(diff) > max_diff) max_diff = abs(diff);
	        }
      return max_diff;
    } else {
      return matter->RegridCriterion();
    }
  }
  else if (grid->Regrid_Type() == 1) {
    return grid->RegridCriterion(tau_c);
  }
  else if (grid->Regrid_Type() == 2) {
    return (grid->r_max() - grid->r_max_final()) - (t_max - t);
  }
  else {
    return 0;
  }
}





