// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing monitor stuff (writes out
// diagnostics as function of time)
//================================================
//
#ifndef MONITOR_H
#define MONITOR_H

// #define _DEBUG_

#include <ctime>
#include "nr3.h"
#include "InData.h"
#include "Slicing.h"
#include "Gauge.h"
#include "Cosmology.h"
#include "gridfunction.h"

class Monitor {
private:
  int note_step;
  InData *indata;
  Slicing *slicing;
  Gauge *gauge;
  const char * filestem;
  //   ostringstream name;
  ofstream outfile, constoutfile, invoutfile;
public:
  //================================================
  // Constructor (opens file and prints header)
  //================================================
  Monitor(int step, const char * filename_i, Grid *grid, InData *indata_i,
	  Slicing * slicing_i, Gauge * gauge_i, Cosmology * cosmology,
	  double eta_KO, int sigma, 
	  int z4, double kappa_11, double kappa_12, double kappa_2, 
	  double kappa_ric, int RK_order, int char_OB ) :
    note_step(step), indata(indata_i), 
    slicing(slicing_i), gauge(gauge_i), filestem(filename_i)
  {
    // name << filename_i << ends;
    // cout << name.str().c_str() << endl; 
    //=======================================
    // extract information from grid
    //=======================================
    int N_r = grid->N_r_int();
    int N_theta = grid->N_theta_int();
    int N_phi = grid->N_phi_int();
    double r_param = grid->r_parameter();
    double theta_param = grid->theta_parameter();
    double r_max = grid->r_max();
    double Courant = grid->courant_factor();
    //=======================================
    // get current time
    //=======================================
    time_t clocktime;
    struct tm * currenttime;
    time(&clocktime);
    currenttime=localtime(&clocktime);
    //=======================================
    // create monitor file
    //=======================================    
    ostringstream monitorfilename;
    monitorfilename << filestem << "_" << N_r << "_" << N_theta 
		    << ".mon" << ends;
    outfile.open(monitorfilename.str().c_str());
    if (outfile)
      cout << " Opened monitor file " << monitorfilename.str().c_str() <<
	" for output " << endl;
    else
      cerr << " Could not open file " 
		       << monitorfilename.str().c_str() << endl;
    outfile.setf(ios::left);
    outfile << "# Monitor file created at " << asctime(currenttime);  
#ifdef AXISYMMETRY
    outfile << "# ======================================== " << endl;
    outfile << "#   Running code assuming AXISYMMETRY     " << endl;
    outfile << "# ======================================== " << endl;
#endif /* AXISYMMETRY */
    outfile << "# Monitoring evolution of " << indata->Name() << endl;
    outfile << "# in asymptotic " << cosmology->Name() << " spacetime " << endl; 
    if (z4) outfile << "# Evolving with Z4c formalism with kappa_11 = " 
		    << kappa_11 << ", kappa_12 = " << kappa_12 
		    << ", kappa_2 = " << kappa_2
		    << " and kappa_ric = " << kappa_ric  
		    << endl;
    else outfile << "# Evolving with BSSN formalism " << endl;
    if (RK_order == 3)
      outfile << "# Evolving with third-order Runge-Kutta " << endl;
    else if (RK_order == 4)
      outfile << "# Evolving with fourth-order Runge-Kutta " << endl;
    else
      outfile << "# Evolving with iterative Crank-Nicholson " << endl;
#ifdef EIGHTHORDER
    outfile << "# Evolving with eighth-order finite-differencing stencils " << endl;
#elif SIXTHORDER
    outfile << "# Evolving with sixth-order finite-differencing stencils " << endl;
#else
    outfile << "# Evolving with fourth-order finite-differencing stencils " << endl;
#endif
    if (char_OB == 0) outfile << "# Implementing Sommerfeld BCs with derivatives " << endl;
    else outfile << "# Implementing Sommerfeld BCs with characteristic interpolation " << endl;
    outfile << "# Evolved with (" << N_r << "," << N_theta << "," << N_phi
	    << ") gridpoints (NOT including ghosts)," << endl;
#ifdef EIGHTHORDER
    outfile << "# using eighth-order spatial differencing with " << grid->N_ghosts() << " ghost zones " << endl;
#elif SIXTHORDER
    outfile << "# using sixth-order spatial differencing with " << grid->N_ghosts() << " ghost zones " << endl;
#else
    outfile << "# using fourth-order spatial differencing with " << grid->N_ghosts() << " ghost zones " << endl;
#endif

    outfile << "# Grid extends to r_max = " << r_max
	    << ", Courant factor = " << Courant << endl;
    outfile << "# Parameters for radial and angular grid : " << r_param 
	    << " and " << theta_param << endl;
    outfile << "# Kreiss-Oliger coefficient = " << eta_KO << endl;
    if (sigma == 0) 
      outfile << "# Using EULERIAN form of evolution equations " << endl; 
    else if (sigma == 1)
      outfile << "# Using LAGRANGIAN form of evolution equations " << endl; 

    outfile << "# Evolution with " << slicing->Name() << " and " << endl;
    outfile << "# " << gauge->Name() << endl;
    outfile << "# Noting values every " << note_step << " time steps." << endl;
    outfile << "# " << setw(14) << "time" 
	    << setw(16) << "tau at center"
	    << setw(16) << "ADM mass"
	    << setw(16) << "Ang. Mom."
	    << setw(16) << "Lin. Mom."
	    << setw(16) << "phi at r=0"
	    << setw(16) << "lapse at r=0"
	    << setw(16) << "minimum lapse"
	    << setw(16) << "K at r=0"
	    << setw(16) << "regrid crit."
	    << endl;
    outfile << "#=============================================================================================================================================" << endl;
    //=======================================
    // create constraint monitor file
    //=======================================    
    ostringstream constraintmonitorfilename;
    constraintmonitorfilename << filestem << "_" << N_r << "_" 
			      << N_theta << ".const_mon" << ends;
    constoutfile.open(constraintmonitorfilename.str().c_str());
    if (constoutfile)
      cout << " Opened monitor file " << constraintmonitorfilename.str().c_str() <<
	" for output " << endl;
    else
      cerr << " Could not open file " 
		       << constraintmonitorfilename.str().c_str() << endl;
    constoutfile.setf(ios::left);
    constoutfile << "# Constraint monitor file created at " << asctime(currenttime);  
#ifdef AXISYMMETRY
    constoutfile << "# ======================================== " << endl;
    constoutfile << "#   Running code assuming AXISYMMETRY     " << endl;
    constoutfile << "# ======================================== " << endl;
#endif /* AXISYMMETRY */
    constoutfile << "# Monitoring evolution of " << indata->Name() << " in asymptotic " << cosmology->Name() 
	    << " spacetime " << endl;
    if (z4) constoutfile << "# Evolving with Z4c formalism with kappa_11 = " 
		    << kappa_11 << ", kappa_12 = " << kappa_12 
		    << ", kappa_2 = " << kappa_2
		    << " and kappa_ric = " << kappa_ric  
		    << endl;
    else constoutfile << "# Evolving with BSSN formalism " << endl;
    if (RK_order == 3)
      constoutfile << "# Evolving with third-order Runge-Kutta " << endl;
    else if (RK_order == 4)
      constoutfile << "# Evolving with fourth-order Runge-Kutta " << endl;
    else
      constoutfile << "# Evolving with iterative Crank-Nicholson " << endl;
    constoutfile << "# Evolved with (" << N_r << "," << N_theta << "," << N_phi
	    << ") gridpoints (NOT including ghosts)," << endl;
    constoutfile << "# Grid extends to r_max = " << r_max
	    << ", Courant factor = " << Courant << endl;
    constoutfile << "# Parameters for radial and angular grid : " << r_param 
	    << " and " << theta_param << endl; 
    constoutfile << "# Kreiss-Oliger coefficient = " << eta_KO << endl;
    if (sigma == 0) 
      constoutfile << "# Using EULERIAN form of evolution equations " << endl; 
    else if (sigma == 1)
      constoutfile << "# Using LAGRANGIAN form of evolution equations " << endl; 
    constoutfile << "# Evolution with " << slicing->Name() << " and " << endl;
    constoutfile << "# " << gauge->Name() << endl;
    constoutfile << "# Noting values every " << note_step << " time steps." << endl;
    constoutfile << "# " << setw(14) << "time" 
		 << setw(16) << "tau at center"
		 << setw(16) << "Hamiltonian"
		 << setw(16) << "Ham (excised)"
		 << setw(16) << "Mom_r"
		 << setw(16) << "Mom_t"
		 << setw(16) << "Mom_p"
		 << setw(16) << "CFC_r"
		 << setw(16) << "CFC_t"
		 << setw(16) << "CFC_p"
		 << endl;
    constoutfile << "#============================================================================================================" << endl;

    //=======================================
    // create curvature invariant monitor file
    //=======================================    
    ostringstream invmonitorfilename;
    invmonitorfilename << filestem << "_" << N_r << "_" 
			     << N_theta << ".inv_mon" << ends;
    invoutfile.open(invmonitorfilename.str().c_str());
    if (invoutfile)
      cout << " Opened monitor file " << invmonitorfilename.str().c_str() <<
	" for output " << endl;
    else
      cerr << " Could not open file " 
		       << invmonitorfilename.str().c_str() << endl;
    invoutfile.setf(ios::left);
    invoutfile << "# Constraint monitor file created at " << asctime(currenttime);  
#ifdef AXISYMMETRY
    invoutfile << "# ======================================== " << endl;
    invoutfile << "#   Running code assuming AXISYMMETRY     " << endl;
    invoutfile << "# ======================================== " << endl;
#endif /* AXISYMMETRY */
    invoutfile << "# Monitoring evolution of " << indata->Name() << " in asymptotic " << cosmology->Name() 
	    << " spacetime " << endl;
    if (z4) invoutfile << "# Evolving with Z4c formalism with kappa_11 = " 
		    << kappa_11 << ", kappa_12 = " << kappa_12 
		    << ", kappa_2 = " << kappa_2
		    << " and kappa_ric = " << kappa_ric  
		    << endl;
    else invoutfile << "# Evolving with BSSN formalism " << endl;
    if (RK_order == 3)
      invoutfile << "# Evolving with third-order Runge-Kutta " << endl;
    else if (RK_order == 4)
      invoutfile << "# Evolving with fourth-order Runge-Kutta " << endl;
    else
      invoutfile << "# Evolving with iterative Crank-Nicholson " << endl;
    invoutfile << "# Evolved with (" << N_r << "," << N_theta << "," << N_phi
	    << ") gridpoints (NOT including ghosts)," << endl;
    invoutfile << "# Grid extends to r_max = " << r_max
	    << ", Courant factor = " << Courant << endl;
    invoutfile << "# Parameters for radial and angular grid : " << r_param 
	    << " and " << theta_param << endl; 
    invoutfile << "# Kreiss-Oliger coefficient = " << eta_KO << endl;
    if (sigma == 0) 
      invoutfile << "# Using EULERIAN form of evolution equations " << endl; 
    else if (sigma == 1)
      invoutfile << "# Using LAGRANGIAN form of evolution equations " << endl; 
    invoutfile << "# Evolution with " << slicing->Name() << " and " << endl;
    invoutfile << "# " << gauge->Name() << endl;
    invoutfile << "# Noting values every " << note_step << " time steps." << endl;
    invoutfile << "# " << setw(14) << "time" 
		 << setw(16) << "tau at center"
		 << setw(16) << "I_re (center)"
		 << setw(16) << "I_re (cur max)"
		 << setw(16) << "I_re_int"
		 << setw(16) << "I_im (center)"
		 << setw(16) << "I_im (cur max)"
		 << setw(16) << "I_im_int"
		 << setw(16) << "J_re (center)"
		 << setw(16) << "J_re (cur max)"
		 << setw(16) << "J_re_int"
		 << setw(16) << "J_im (center)"
		 << setw(16) << "J_im (cur max)"
		 << setw(16) << "J_im_int"
		 << endl;
    invoutfile << "#============================================================================================================" << endl;

  };
  //================================================
  // check whether it's time to take note...
  //================================================
  bool time_to_note(int time_step) { return (time_step % note_step == 0); }
  //================================================
  // Note routines
  //================================================
  void note(int time_step, double time, double tau, double mass, 
	    double ang_mom, double P, double phi, double lapse, 
	    double lapse_min, double K, double regrid,
	    bool force_note = false) {
    if (time_to_note(time_step) || force_note) {
      outfile << setw(16) << time 
	      << setw(16) << setprecision(8) << tau 
	      << setw(16) << setprecision(8) << mass
	      << setw(16) << setprecision(8) << ang_mom
	      << setw(16) << setprecision(8) << P
	      << setw(16) << setprecision(8) << phi
	      << setw(16) << setprecision(8) << lapse
	      << setw(16) << setprecision(8) << lapse_min
	      << setw(16) << setprecision(8) << K
	      << setw(16) << setprecision(8) << regrid
	      << endl;
    }
  };
  void note_constraints(int time_step, double time, double tau, 
			double Ham, double Hamex, 
			double Mom_r, double Mom_t, double Mom_p, 
			double CFC_r, double CFC_t, double CFC_p) {
      constoutfile << setw(16) << time 
		   << setw(16) << setprecision(8) << tau 
		   << setw(16) << setprecision(8) << Ham
		   << setw(16) << setprecision(8) << Hamex
		   << setw(16) << setprecision(8) << Mom_r
		   << setw(16) << setprecision(8) << Mom_t
		   << setw(16) << setprecision(8) << Mom_p
		   << setw(16) << setprecision(8) << CFC_r
		   << setw(16) << setprecision(8) << CFC_t
		   << setw(16) << setprecision(8) << CFC_p
		   << endl;
  }
  void note_invariants(int time_step, double time, double tau, 
		       double I_re_center, double I_re_max, double I_re_int, 
		       double I_im_center, double I_im_max, double I_im_int, 
		       double J_re_center, double J_re_max, double J_re_int, 
		       double J_im_center, double J_im_max, double J_im_int) {
      invoutfile << setw(16) << time 
		   << setw(16) << setprecision(8) << tau 
		   << setw(16) << setprecision(8) << I_re_center
		   << setw(16) << setprecision(8) << I_re_max
		   << setw(16) << setprecision(8) << I_re_int
		   << setw(16) << setprecision(8) << I_im_center
		   << setw(16) << setprecision(8) << I_im_max
		   << setw(16) << setprecision(8) << I_im_int
		   << setw(16) << setprecision(8) << J_re_center
		   << setw(16) << setprecision(8) << J_re_max
		   << setw(16) << setprecision(8) << J_re_int
		   << setw(16) << setprecision(8) << J_im_center
		   << setw(16) << setprecision(8) << J_im_max
		   << setw(16) << setprecision(8) << J_im_int
		   << endl;
  };
  //================================================
  // regrid
  //================================================
  void regrid(Doub r_max) {
    outfile << "# Regridding with r_max = " << r_max << endl;
  };
  void set_t_max(Doub t_max) {
    outfile << "# Reset t_max to " << t_max << endl;
  };
  //================================================
  // Warning
  //================================================
  int warning(int error) {
    if (error == 1) {
      outfile << "# Oh no: found Nan's!" << endl;
    }
    return error;
  };
  //================================================
  // return filestem name and note_step
  //================================================
  const char * Filestem() { return filestem; }
  int Note_Step() { return note_step; }
  //================================================
  // Destructor (closes file)
  //================================================
  ~Monitor() { 
    //
    // get current time
    //
    time_t clocktime;
    struct tm * currenttime;
    time(&clocktime);
    currenttime=localtime(&clocktime);
    outfile << "# Monitor file closed at " << asctime(currenttime) 
	    << "# Bye!" << endl;  
    outfile.close(); 
    constoutfile << "# Monitor file closed at " << asctime(currenttime) 
	    << "# Bye!" << endl;  
    constoutfile.close(); 
    invoutfile << "# Monitor file closed at " << asctime(currenttime) 
	    << "# Bye!" << endl;  
    invoutfile.close(); 
  };
};

#endif  /*  MONITOR */
