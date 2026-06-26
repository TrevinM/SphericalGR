// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing all the integration stuff...
//================================================
#ifndef MANAGER_H
#define MANAGER_H

#include "Grid.h"
#include "dumper.h"
#include "Monitor.h"
#include "CheckPoint.h"
#include "Slicing.h"
#include "Gauge.h"
#include "Cosmology.h"
#include "State.h"
#include "Curvature.h"
#include "Auxiliary.h"
#include "Diagnostics.h"
#include "InData.h"
#include "Matter.h"
#include "EOS.h"
#include "HorizonFinder.h"
#include "Photons.h"
#include "ConstraintSolver.h"
#include "Profiles.h"
#include "WaveExtraction.h"
#include "CheckPoint.h"

class Manager {
 private:
  double t, tau_c; // coordinate time and proper time at center
  int step;        // step counter
  double dt;
  Grid *grid;
  dumper *dump;
  Monitor *monitor;
  CheckPoint *checkpoint;
  Profiles *profiles;
  WaveExtraction *waves;
  InData *indata;
  EOS *eos;
  Slicing *slicing;
  Gauge *gauge;
  Cosmology *cosmology;
  diagnostics * constraints;
  state *last;
  state *derivs;
  state *inter;
  state *updates;
  curvature *curve;
  auxiliary *aux;
  Matter * matter; 
  HorizonFinder * horizonfinder;
  Photons * photons;
  int sigma;      // decides between Eulerian and Lagrangian formalism
  int cowling;    // Cowling approximation (0: no, 1: fix gravity, 2: fix matter)
  double eta, eta_KO;   // Kreiss-Oliger coefficients
  int z4;         // switch between BSSN (0) and Z4c (1)
  double kappa_11, kappa_12, kappa_2, kappa_ric;   // parameters for Z4c
  int RK_order;    // order of Runge-Kutta integrator
  int char_OB;     // switch for implementation of Sommerfeld BCs
  int solve_constraints; // whether or not to solve constraints 
  VecDoub r, r2, theta, sintheta, sin2theta, costheta, phi;
  int N_g, N_r, N_t, N_p;   // number of ghost and grid points
  //
  int steps_between_regrids, timestep_last_regrid;
  double PI;
 public:
  //===========================================
  // Constructor
  //===========================================
  Manager(Grid *grid_i, dumper *dump_i, Monitor *monitor_i,
	  CheckPoint *checkpoint_i,
	  Profiles *profiles_i, WaveExtraction *waves_i, InData *indata_i,
	  EOS *eos_i, Slicing *slicing_i, Gauge *gauge_i, 
	  Cosmology * cosmology_i, Photons * photons_i,
	  int matter_type, int sigma_i, int cowling_i, double eta_i,
	  int z4_i, double kappa_11_i, double kappa_12_i, double kappa_2_i,
	  double kappa_ric_i, int RK_order_i,
	  int char_OB_i, int solve_constraints_i) :
    grid(grid_i),
    dump(dump_i), monitor(monitor_i), checkpoint(checkpoint_i),
    profiles(profiles_i),
    waves(waves_i), indata(indata_i), eos(eos_i),
    slicing(slicing_i), gauge(gauge_i), cosmology(cosmology_i), 
    photons(photons_i),
    sigma(sigma_i), cowling(cowling_i), eta(eta_i),
    z4(z4_i), kappa_11(kappa_11_i), kappa_12(kappa_12_i), 
    kappa_2(kappa_2_i), kappa_ric(kappa_ric_i), 
    RK_order(RK_order_i), char_OB(char_OB_i), 
    solve_constraints(solve_constraints_i),
    r(grid_i->N_r_tot()),r2(grid_i->N_r_tot()), theta(grid_i->N_theta_tot()),
    sintheta(grid_i->N_theta_tot()),sin2theta(grid_i->N_theta_tot()),
    costheta(grid_i->N_theta_tot()),phi(grid_i->N_phi_tot())
  {
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
				      monitor->Filestem() );
    //
    // create diagnostics class
    //
    cout << " MANAGER: setting up diagnostics... " << endl;
    constraints = new diagnostics(grid,dump,monitor,cosmology,
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
  // Integrate (in Integrate.C)
  //===========================================
  bool Integrate(double t_max);
  //===========================================
  // Carry out one Runge-Kutta time step (in Integrate.C)
  //===========================================
  void Step_RK4(double & t, double & tau_c);
  void Step_RK3(double & t, double & tau_c);
  void Step_ICN(double & t, double & tau_c);
  //===========================================
  // Find time step
  //===========================================
  void SetTimeStep() {
    int N_g = grid->N_ghosts(); 
    dt = grid->courant_factor() * 0.5 * grid->delta_r(N_g)
      * grid->delta_theta(N_g); 
    cout << " MANAGER: using dt = " << dt << endl;
  };
  void SetTimeStep(state *s) {
    int N_g = grid->N_ghosts();
    const double gamma_tt = (1.0 + s->h_tt(N_g, N_g, N_g))*exp(4.0 * s->phi(N_g, N_g, N_g));
    // const double factor = sqrt(gamma_tt) / s->lapse(N_g, N_g, N_g); // + abs(s->shift_r(N_g, N_g, N_g));
    const double factor = sqrt(gamma_tt); // + abs(s->shift_r(N_g, N_g, N_g));
    dt = grid->courant_factor() * 0.5 * grid->delta_r(N_g) * grid->delta_theta(N_g) * factor;
  }; 
  //===========================================
  // Initialize
  //===========================================
  bool Initialize();
  //===========================================
  // Compute right-hand sides of equations...
  //===========================================
  void Compute_RHS(double time = 0.0);
  void dot_metric(state * c, double time = 0.0);
  void dot_ext_curv(state * c, double time = 0.0);
  void dot_connection(state * c, double time = 0.0);
  //===========================================
  // Tests (all in Tests.C)
  //===========================================
  void Test_Indices();
  void Test_Ricci_for_Schwarzschild();
  void Test_Flat_Metric();
  void Test_EOS();
  void Test_Derivatives();
  //===========================================
  // Regrid
  //===========================================
  double Regrid(double t, double tau_c, int timestep, double & t_max);
  double RegridCriterion();
  //================================================
  //
  // Rescale metric
  // 
  //================================================
  void Rescale_Metric(state *s) {
    curve->Compute_Determinant(s);
#pragma omp parallel for collapse(3)
    for (int i = 0; i < N_r; i++)    
      for (int j = 0; j < N_t; j++)
	for (int k = 0; k < N_p; k++) {
	  const double factor = pow(curve->det_init(i,j,k)/curve->det(i,j,k),1.0/3.0);
	  s->h_rr[i][j][k] = factor*(1.0 + s->h_rr(i,j,k)) - 1.0;
	  s->h_rt[i][j][k] *= factor;
	  s->h_rp[i][j][k] *= factor;
	  s->h_tt[i][j][k] = factor*(1.0 + s->h_tt(i,j,k)) - 1.0;
	  s->h_tp[i][j][k] *= factor;
	  s->h_pp[i][j][k] = factor*(1.0 + s->h_pp(i,j,k)) - 1.0;
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
  void Remove_Trace(state *s) {
    curve->Compute_Trace(s);
    const double onethird = 1.0/3.0;
#pragma omp parallel for collapse(3)
    for (int i = 0; i < N_r; i++)    
      for (int j = 0; j < N_t; j++)
	for (int k = 0; k < N_p; k++) {
	  s->a_rr[i][j][k] -= onethird * (1.0 + s->h_rr(i,j,k)) * curve->trace_A(i,j,k);
	  s->a_rt[i][j][k] -= onethird * (      s->h_rt(i,j,k)) * curve->trace_A(i,j,k);
	  s->a_rp[i][j][k] -= onethird * (      s->h_rp(i,j,k)) * curve->trace_A(i,j,k);
	  s->a_tt[i][j][k] -= onethird * (1.0 + s->h_tt(i,j,k)) * curve->trace_A(i,j,k);
	  s->a_tp[i][j][k] -= onethird * (      s->h_tp(i,j,k)) * curve->trace_A(i,j,k);
	  s->a_pp[i][j][k] -= onethird * (1.0 + s->h_pp(i,j,k)) * curve->trace_A(i,j,k);
	  // K[i][j][k] += curve->trace_A(i,j,k);
	}
    curve->Compute_Trace(s);
  };  
  //===========================================
  // Destructor
  //===========================================
  ~Manager() {
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
};
  

#endif  /* MANAGER_H */
