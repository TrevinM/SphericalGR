// Tell emacs that this is -*-c++-*- mode
#include "nr3.h"
//================================================================
//
// Source File for input
//
//================================================================
int Read_Input(int argc, char * argv[], 
	       int & get_rid_of_output,
	       double & t_max,
	       double & eta_KO, int & sigma,
	       int & cowling,
	       int & z4, double & kappa_11, double & kappa_12,
	       double  & kappa_2, double & kappa_ric,
	       int & RK_order, int & char_OB, int & dump_step,
	       char * monitor_filename, int & note_step,
	       int & read_from_chkpt, int & chkpt_step,
	       int & space_type,
	       int & indata_type, char * indata_input, int & matter_type,
	       int & eos_type, int & slicing_type, int & gauge_type,
	       bool & solve_constraints, bool & rescale_metric,
	       bool & track_photons, bool & write_profiles,
	       bool & extract_waves)
{
  char sourcefiledir[128];
  int error = 0;
  if (argc != 2) {
    cerr << "Syntax: ./SphericalGR Input_file" << endl;
    error = 1;
  } else {
    ifstream infile;
    infile.open(argv[1]);
    if(!infile){
      cerr << "Can't open " << argv[1] << " for input." << endl;
      error = 2;
      return error;
    } else {
      cout << "Reading input from file " << argv[1] << endl;
    }
    char buf[1000],c;
    //    char run_case[8];
    infile.get(buf,1000,'='); infile.get(c); infile.get(c); 
    infile.get(sourcefiledir,128);
    infile.get(buf,1000,'='); infile.get(c); infile >> get_rid_of_output;
    infile.get(buf,1000,'='); infile.get(c); infile >> t_max;
    infile.get(buf,1000,'='); infile.get(c); infile >> eta_KO;
    infile.get(buf,1000,'='); infile.get(c); infile >> sigma;
    infile.get(buf,1000,'='); infile.get(c); infile >> cowling;
    infile.get(buf,1000,'='); infile.get(c); infile >> z4;
    infile.get(buf,1000,'='); infile.get(c); infile >> kappa_11;
    infile.get(buf,1000,'='); infile.get(c); infile >> kappa_12;
    infile.get(buf,1000,'='); infile.get(c); infile >> kappa_2;
    infile.get(buf,1000,'='); infile.get(c); infile >> kappa_ric;
    infile.get(buf,1000,'='); infile.get(c); infile >> RK_order;
    infile.get(buf,1000,'='); infile.get(c); infile >> char_OB;
    infile.get(buf,1000,'='); infile.get(c); infile >> dump_step;
    infile.get(buf,1000,'='); infile.get(c); infile.get(c); 
    infile.get(monitor_filename,64);
    infile.get(buf,1000,'='); infile.get(c); infile >> note_step;
    infile.get(buf,1000,'='); infile.get(c); infile >> read_from_chkpt;
    infile.get(buf,1000,'='); infile.get(c); infile >> chkpt_step;
    infile.get(buf,1000,'='); infile.get(c); infile >> space_type;
    infile.get(buf,1000,'='); infile.get(c); infile >> indata_type;
    infile.get(buf,1000,'='); infile.get(c); infile.get(c); 
    infile.get(indata_input,64);
    infile.get(buf,1000,'='); infile.get(c); infile >> matter_type;
    infile.get(buf,1000,'='); infile.get(c); infile >> eos_type;
    infile.get(buf,1000,'='); infile.get(c); infile >> slicing_type;
    infile.get(buf,1000,'='); infile.get(c); infile >> gauge_type;
    int solve, rescale, photons, profiles, waves;
    infile.get(buf,1000,'='); infile.get(c); infile >> solve;
    infile.get(buf,1000,'='); infile.get(c); infile >> rescale;
    infile.get(buf,1000,'='); infile.get(c); infile >> photons;
    infile.get(buf,1000,'='); infile.get(c); infile >> profiles;
    infile.get(buf,1000,'='); infile.get(c); infile >> waves;
    if(infile.eof()) {
      cerr << "Error reading input file " << argv[1] << endl;
      error = 3;
    }
    cout << "   Maximum time : " << t_max << endl;
    cout << "   Coefficient for Kreiss-Oliger : " << eta_KO << endl;
    if (sigma == 0) 
      cout << "   Using EULERIAN formulation " << endl;
    else if (sigma == 1)
      cout << "   Using LAGRANGIAN formulation " << endl;
    else {
      cerr << "   No formulation for sigma = " << sigma << endl;
      return 6;
    }
    if (cowling == 0) 
      cout << "   Evolving all fields dynamically " << endl;
    else if (cowling == 1)
      cout << "   Using COWLING approximation (grav. fields fixed) " << endl;
    else if (cowling == 2)
      cout << "   Using hydro-without-hydro (matter fields fixed) " << endl;
    else if (cowling == 3)
      cout << "   Using 'shifted' hydro-without-hydro " << endl;
    else {
      cerr << "   No formulation for cowling = " << cowling << endl;
      return 6;
    }
    if (z4 == 0)
      cout << "   Evolving with BSSN " << endl;
    else if (z4 == 1)
      cout << "   Evolving with Z4 with kappa_11 = " << kappa_11
	   << ", kappa_12 = " << kappa_12 
	   << ", kappa_2 = " << kappa_2 
	   << " and kappa_ric = " << kappa_ric 
	   << endl;
    else {
      cerr << "   No formulation for z4 = " << z4 << endl;
      return 7;
    }
    if (RK_order == 3)
      cout << "   Evolving with third-order Runge-Kutta " << endl;
    else if (RK_order == 4)
      cout << "   Evolving with fourth-order Runge-Kutta " << endl;
    else if (RK_order == 2)
      cout << "   Evolving with iterative Crank-Nicholson " << endl;
    else {
      cerr << "   Runge-Kutta order " << RK_order 
	   << " has not been implemented! " << endl;
      return 9;
    }
    if (char_OB == 0)
      cout << "   Applying Sommerfeld conditions via derivative " << endl;
    else if (char_OB == 1)
      cout << "   Applying Sommerfeld conditions via char. interpolation " 
	   << endl;
    else {
      cerr << "   No such outer boundary implementation " << endl;
      return 8;
    }
    cout << "   Dumping output every " << dump_step << " timesteps " << endl;
    cout << "   Monitor writes into file " << monitor_filename << endl;
    cout << "   Monitor will note every " << note_step << " timesteps " << endl;
    if (read_from_chkpt > 0) 
      cout << "   Will start from checkpoint at timestep "
	   << read_from_chkpt << endl;
    else
      cout << "   Will start with initial data (NOT from check point)" << endl;
    cout << "   Will write check point files every " << chkpt_step
	 << " timesteps" << endl;
    if (space_type == 1) 
      cout << "   Asymptotically Minkowski " << endl;
    else if (space_type == 2)
      cout << "   Asymptotically de Sitter " << endl;
    else if (space_type == 3)
      cout << "   Radiation-dominated Universe " << endl;
    else {
      cout << "   NO SUCH ASYMPTOTIC SPACE " << endl;
      error = 5;
    }
    if (indata_type == 1) 
      cout << "   Running with Teukolsky wave initial data " << endl;
    else if (indata_type == 2)
      cout << "   Running with Schwarzschild initial data " << endl;
    else if (indata_type == 3)
      cout << "   Running with flat initial data " << endl;
    else if (indata_type == 4)
      cout << "   Running with TOV initial data " << endl;
    else if (indata_type == 5)
      cout << "   Running with Trumpet initial data " << endl;
    else if (indata_type == 6)
      cout << "   Running with RNS initial data " << endl;
    else if (indata_type == 7)
      cout << "   Running with Brill Wave initial data " << endl;
    else if (indata_type == 8)
      cout << "   Running with Kerr initial data " << endl;
    else if (indata_type == 9)
      cout << "   Running with Bowen-York initial data " << endl;
    else if (indata_type == 10)
      cout << "   Running with Kerr-Schild initial data " << endl;
    else if (indata_type == 11)
      cout << "   Running with Ken's trumpet initial data " << endl;
   else if (indata_type == 12)
      cout << "   Running with Brill_Lindquist initial data " << endl;
   else if (indata_type == 13)
      cout << "   Running with spherical shock-tube initial data " << endl;
   else if (indata_type == 14)
      cout << "   Running with Evans-Coleman initial data for radiation fluid " << endl;
   else if (indata_type == 15)
      cout << "   Running with TOV initial data for radiation fluid " << endl;
   else if (indata_type == 16)
      cout << "   Running with OS initial data " << endl;
   else if (indata_type == 17)
      cout << "   Running with Rotating Perfect Fluid initial data " << endl;
   else if (indata_type == 18)
      cout << "   Running with Bondi accreation initial data " << endl;
   else if (indata_type == 19)
      cout << "   Running with Choptuik initial data " << endl;
   else if (indata_type == 20)
      cout << "   Running with E&M wave initial data " << endl;
   else if (indata_type == 21)
      cout << "   Running with Radiation-Hydro Shock initial data " << endl;
   else if (indata_type == 22)
      cout << "   Running with TOV_BH initial data " << endl;
   else if (indata_type == 23)
      cout << "   Running with SMS_TOV initial data " << endl;
   else if (indata_type == 24)
      cout << "   Running with Gauge Wave initial data " << endl;
   else if (indata_type == 25)
      cout << "   Running with Wind Tunnel initial data " << endl;
   else if (indata_type == 26)
      cout << "   Running with Shibata wave initial data " << endl;
   else if (indata_type == 27)
      cout << "   Running with disk initial data " << endl;
   else if (indata_type == 28)
      cout << "   Running with dual E&M wave initial data " << endl;
   else {
      cout << "   NO SUCH INITIAL DATA " << endl;
      error = 5;
    }
    cout << "   Will read initial data input from file " << indata_input << endl;    
    if (matter_type == 1) 
      cout << "   Running with vacuum " << endl;
    else if (matter_type == 2)
      cout << "   Running with hydro_without_hydro " << endl;
    else if (matter_type == 3)
      cout << "   Running with hydro " << endl;
    else if (matter_type == 4)
      cout << "   Running with perfect fluid " << endl;
    else if (matter_type == 5)
      cout << "   Running with scalar field " << endl;
    else if (matter_type == 6)
      cout << "   Running with radiation hydrodynamics " << endl;
    else if (matter_type == 7)
      cout << "   Running with Maxwell's equations " << endl;
    else if (matter_type == 8)
      cout << "   Running with Dual Maxwell's equations " << endl;
    else {
      cout << "   Unknown matter source! " << endl;
      error = 4;
    }
    if (slicing_type == 1)  
      cout << "   Running with geodesic slicing " << endl;
    else if (slicing_type == 2)
      cout << "   Running with (non-advective) 1+log slicing " << endl;
    else if (slicing_type == 3)
      cout << "   Running with advective 1+log slicing " << endl;
    else if (slicing_type == 4)
      cout << "   Running with harmonic slicing " << endl;
    else if (slicing_type == 5)
      cout << "   Running with maximal slicing " << endl;
    else if (slicing_type == 6)
      cout << "   Running with second-order maximal slicing " << endl;
    else if (slicing_type == 7)
      cout << "   Running with Ken's log slicing " << endl;
    else if (slicing_type == 8)
      cout << "   Running with generalized advective 1+log slicing " << endl;
    else
      error = 4;
    if (gauge_type == 1)
      cout << "   Running with zero shift " << endl;
    else if (gauge_type == 2)
      cout << "   Running with (non-advective) Gamma-driver shift " << endl;
    else if (gauge_type == 3)
      cout << "   Running with advective Gamma-driver shift " << endl;
    else if (gauge_type == 4)
      cout << "   Running with Jena version of (non-advective) Gamma-driver shift " << endl;
    else if (gauge_type == 5)
      cout << "   Running with Jena version of advective Gamma-driver shift " << endl;
    else if (gauge_type == 6)
      cout << "   Running with Modified Gamma-driver shift " << endl;
    else if (gauge_type == 7)
      cout << "   Running with Covariant Jena version of advective Gamma-driver shift " << endl;
    else if (gauge_type == 8)
      cout << "   Running with self-similar shift " << endl;
    else
      error = 4;    
    if (solve == 1) {
      cout << "   Solving Constraints at t = 0 " << endl;
      solve_constraints = true;
    } else {
      solve_constraints = false;
    }
    if (rescale == 1) {
      cout << "   Rescaling initial metric so that det = 1 " << endl;
      rescale_metric = true;
    } else {
      rescale_metric = false;
    }
    if (photons == 1) {
      cout << "   Will track photons emitted from center " << endl;
      track_photons = true;
    } else {
      track_photons = false;
    }
    if (profiles == 1) {
      cout << "   Will write spacetime profiles " << endl;
      write_profiles = true;
    } else {
      write_profiles = false;
    }
    if (waves == 1) {
      cout << "   Will extract gravitational waves " << endl;
      extract_waves = true;
    } else {
      extract_waves = false;
    }
  }  
  //==============================================================
  // Finally, tar up all source file of current version of code
  //============================================================== 
  cout << "  Creating tar-file of source files in directory "
       << sourcefiledir << endl;
  stringstream command;
  command << "tar -czf code.tgz -C " << sourcefiledir << " . "
	  << ends;
  // cout << command.str().c_str() << endl;
  system(command.str().c_str());
  //
  cout << "===================================================" << endl;
  return error;
}
