#include "gridfunction.h"
#include "dumper.h"
#include "InData.h"
#include "Slicing.h"
#include "Gauge.h"

//================================================
// forced dump
//================================================
void dumper::dump(double phys_time, double prop_time, int timestep, gf3d *fct, 
		  const char * suffix)
{
  //================================================
  // find number of gridpoints
  //================================================
  int N_r = fct->dim1();
  int N_theta = fct->dim2();
  int N_phi = fct->dim3();
  int N_g = fct->N_ghosts();
  //================================================
  // make sure that ghosts are filled...
  //================================================
  fct->fill_ghosts();
  //================================================
  // first write out rays
  //================================================
  ofstream outfile;
  ostringstream rayfilename;
  rayfilename << "output/" << fct->Name() << "_rays_" << N_r - 2*N_g 
	      << "_" << N_theta - 2*N_g 
	      << "_" << setfill('0') << setw(8)
	      << timestep << suffix << ends;
  outfile.open(rayfilename.str().c_str());
  if (!outfile) cerr << " Could not open file " 
		     << rayfilename.str().c_str() << endl;
  outfile.setf(ios::right);
  time_t clocktime;
  struct tm * currenttime;
  time(&clocktime);
  currenttime=localtime(&clocktime);
#ifdef _DEBUG_
  cout << " Current time: " << asctime(currenttime);
#endif
  double tl = 0.0;
  double ta = PI/2;  // CHECK: adjust according to preference...  
  int np = N_phi/2;
  //    double tl = fct->theta(nt);
  //    double ta = fct->theta(na);
  double pl = fct->phi(np);
  outfile << "# File created on " << asctime(currenttime);  
#ifdef AXISYMMETRY
  outfile << "# ======================================== " << endl;
  outfile << "# = Running code assuming AXISYMMETRY    = " << endl;
  outfile << "# ======================================== " << endl;
#endif /* AXISYMMETRY */
  outfile << "# Data for " << indata->Name() << " initial data at coord. time " 
  	  << setprecision(16) << setw(24) << phys_time << " and proper central time " << prop_time << endl;
  outfile << "# Evolved with (" << N_r << "," << N_theta << "," << N_phi
  	  << ") gridpoints (including " << N_g << " ghosts)" << endl; 
  outfile << "# Evolution with " << slicing->Name() << " and " 
  	  << gauge->Name() << endl;
#ifdef EIGHTHORDER
  outfile << "# using eighth-order spatial differencing with " << N_g << " ghost zones " << endl;
#elif SIXTHORDER
  outfile << "# using sixth-order spatial differencing with " << N_g << " ghost zones " << endl;
#else
  outfile << "# using fourth-order spatial differencing with " << N_g << " ghost zones " << endl;
#endif
  outfile << "# data 1 at theta = " << tl << " and phi = " << pl << endl;
  outfile << "# data 2 at theta = " << ta << " and phi = " << pl << endl;
  outfile << "# " << setw(22) << "r" << setw(24) << fct->Name() 
	  << setw(24) << fct->Name() << endl;
  outfile << "# " << setw(31) << " (" << setw(6) << tl << "," << pl << ")"
	  << "        (" << setw(6) << ta << "," << pl << ")" << endl;
  outfile << "#========================================================================" << endl;
  outfile.setf(ios::right);
  for (int i = 0; i < N_r; i++) {
    const double rl = fct->r(i);
    outfile << setprecision(16) << setw(24) << rl
	    << setw(24) << (*fct)(i,tl,np) 
	    << setw(24) << (*fct)(i,ta,np); 
    // check whether we should print analytical solution:
    if (indata->Analytical_Solution_Exists()) { 
      if (!strcmp("h_rr", fct->Name()))
	outfile << setprecision(16) << setw(24)
		<< indata->h_rr_analytical(rl,tl,pl,phys_time)
		<< setprecision(16) << setw(24)
		<< indata->h_rr_analytical(rl,ta,pl,phys_time);
      if (!strcmp("E", fct->Name()))
	outfile << setprecision(16) << setw(24)
		<< indata->E_analytical(rl,tl,pl,phys_time)
		<< setprecision(16) << setw(24)
		<< indata->E_analytical(rl,ta,pl,phys_time);
    }
    outfile << endl;
  }
  outfile.close();
};


//================================================
// forced slice dump
//================================================
void dumper::slice(double phys_time, double prop_time, int timestep, gf3d * fct,
		   const char * suffix) 
{
  //================================================
  // find number of gridpoints
  //================================================
  int N_r = fct->dim1();
  int N_theta = fct->dim2();
  int N_phi = fct->dim3();
  int N_g = fct->N_ghosts();
  //================================================
  // now write out slice
  //================================================
  ofstream outfile;
  ostringstream rayfilename;
  rayfilename << "output/" << fct->Name() << "_slice_" << N_r - 2*N_g 
	      << "_" << N_theta - 2*N_g 
	      << "_" << setfill('0') << setw(8)
	      << timestep << suffix << ends;
  outfile.open(rayfilename.str().c_str());
  if (!outfile) cerr << " Could not open file " 
		     << rayfilename.str().c_str() << endl;
  outfile.setf(ios::right);
  time_t clocktime;
  struct tm * currenttime;
  time(&clocktime);
  currenttime=localtime(&clocktime);
#ifdef _DEBUG_
  cout << " Current time: " << asctime(currenttime);
#endif
  int np = N_phi/2;
  double pl = fct->phi(np);
  outfile << "# File created on " << asctime(currenttime);  
  outfile << "# Data for " << indata->Name() << " initial data at coord time " 
  	  << phys_time << " and central proper time " << prop_time << endl;
  outfile << "# Evolved with (" << N_r << "," << N_theta << "," << N_phi
  	  << ") gridpoints (including " << N_g << " ghosts)" << endl;
  outfile << "# Evolution with " << slicing->Name() << " and " 
  	  << gauge->Name() << " shift condition " << endl;
  outfile << "# data for phi = " << pl << endl;
  outfile << "# " << setw(14) << "r" << setw(16) << "theta"
	  << setw(20) << fct->Name() << endl;
  outfile << "#===============================================================" << endl;
  outfile.setf(ios::right);
  for (int i = N_g; i < N_r - N_g; i++) {
  //LATER: put 0 back to N_g and N_theta - N_g
    for (int j = N_g; j < N_theta - N_g; j++) { 
      double rl = fct->r(i);
      double thetal = fct->theta(j);
      outfile << setprecision(8) << setw(16) << rl 
	      << setw(16) << thetal 
	      << setprecision(16)
	      << setw(24) << (*fct)[i][j][np] << endl;
    }
    outfile << endl;
  }
  outfile.close();
};
