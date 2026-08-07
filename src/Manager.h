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
  double t, tau_c, t_max; // coordinate time and proper time at center
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
	  int char_OB_i, int solve_constraints_i);

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
  void SetTimeStep();
  void SetTimeStep(state *s);

  //===========================================
  // Initialize
  //===========================================
  void Set_t_max(double set_t);
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
  void Rescale_Metric(state *s);

    //================================================
    //
    // Remove trace
    // 
    //================================================
  void Remove_Trace(state *s);

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
