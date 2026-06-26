// Tell emacs that this is -*-c++-*- mode
// 
//================================================
//
// Contains all methods for extracting gravitational waves
//
// Note: assuming axisymmetry here, and that only h_+, i.e. only
// real part of psi_4, is non-zero.  Could be generalized easily...
//
//================================================
//
#ifndef WAVEEXTRACTION_H
#define WAVEEXTRACTION_H

#include "Grid.h"
#include "gridfunction.h"
#include "surfacefunction.h"
#include "Slicing.h"
#include "InData.h"
#include "Gauge.h"
#include "Cosmology.h"

class WaveExtraction {
 private:
  Grid *grid;
  InData *indata;
  Slicing *slicing;
  Gauge *gauge;
  Cosmology *cosmology;
  int N_r, N_theta, N_phi, N_g;
  ofstream outfile;
  gf3d * psi4_Re_3d, * psi4_Im_3d;
  double * psi4_Re_20, * psi4_Im_20;
  double * h_plus_20, * h_cross_20;
  double * h_plus_dot_20, * h_cross_dot_20;
  double * r_ext;
  int N_ext;
  const char * filestem;
  int note_step;
  double PI;
public:
  //==========================================
  // Constructor
  //==========================================
  WaveExtraction(const char * filename_i, Grid *grid_i, int note_step_i,
		 InData *indata_i, Slicing *slicing_i,
		 Gauge *gauge_i, Cosmology *cosmology_i) :
    grid(grid_i), indata(indata_i), slicing(slicing_i), gauge(gauge_i),
    cosmology(cosmology_i), note_step(note_step_i), filestem(filename_i)
  {
    N_r = grid->N_r_tot();
    N_theta = grid->N_theta_tot();
    N_phi = grid->N_phi_tot();
    N_g = grid->N_ghosts();
    PI = acos(-1.0);
    psi4_Re_3d = NULL;
    psi4_Im_3d = NULL;
  };
  //==========================================
  // Destructor
  //==========================================
  ~WaveExtraction() {
    if (outfile) {
      outfile.close();
    }
    delete r_ext, psi4_Re_20, psi4_Im_20;
    delete h_plus_20, h_cross_20;
    delete h_plus_dot_20, h_cross_dot_20;
    cout << " WAVEEXTRACTION: destructing waveextraction... " << endl;
  }
  //==========================================
  // initialize wave extraction
  //==========================================
  int Initialize(gf3d * psi4_Re_i, gf3d * psi4_Im_i) {
    psi4_Re_3d = psi4_Re_i;
    psi4_Im_3d = psi4_Im_i;
    //
    // set up arrays for storage...
    //
    N_ext = 3;    // could be changed...
    r_ext = new double[N_ext];
    psi4_Re_20 = new double[N_ext];
    psi4_Im_20 = new double[N_ext];
    h_plus_20 = new double[N_ext];
    h_cross_20 = new double[N_ext];
    h_plus_dot_20 = new double[N_ext];
    h_cross_dot_20 = new double[N_ext];
    //
    // initialize h's...
    //
    for (int i = 0; i < N_ext; i++) {
      h_plus_20[i] = 0.0;
      h_cross_20[i] = 0.0;
      h_plus_dot_20[i] = 0.0;
      h_cross_dot_20[i] = 0.0;
    }
    //
    // Read in extraction radii from input file
    //
    ifstream infile;
    infile.open("WaveExtraction_Input");
    if (!infile) { 
      cerr << " WAVEEXTRACTION: can't open WaveExtraction_Input..." << endl;
      cerr << " WAVEEXTRACTION: will use default values for radii..." << endl;
      r_ext[0] = 5;
      r_ext[1] = 6;
      r_ext[2] = 7;
    } else {
      char buf[1000], c;
      infile.get(buf,1000,'='); infile.get(c); infile >> r_ext[0];
      infile.get(buf,1000,'='); infile.get(c); infile >> r_ext[1];
      infile.get(buf,1000,'='); infile.get(c); infile >> r_ext[2];
    }
    cout << " WAVEEXTRACTION: extraction radii r1 = " << r_ext[0]
	 << ", r2 = " << r_ext[1] << ", r3 = " << r_ext[2] << endl;
    //
    // set up output file...
    //
    ostringstream outfilename;
    outfilename << filestem << "_" << N_r - 2*N_g << "_"
		<< N_theta - 2*N_g << ".wave_extraction" << ends;
    outfile.open(outfilename.str().c_str());
    if (outfile)
      cout << " WAVEEXTRACTION: Opened profile file "
	   << outfilename.str().c_str()
	   << " for output " << endl;
    else
      cerr << " WAVEEXTRACTION: Could not open file "
	   << outfilename.str().c_str() << endl;
    outfile << "# Wave extraction data for " << indata->Name()
	    << " initial data " << endl;
    outfile << "# evolved with (" << N_r << "," << N_theta
	    << "," << N_phi << ") gridpoints (including "
	    << N_g << " ghosts)" << endl; 
    outfile << "# Evolution with " << slicing->Name() << " and " 
	    << gauge->Name() << endl;
    outfile << "# Plotting data every " << note_step << " steps "
	    << endl;
    //    outfile << "# Extraction radii r1 = " << r_ext[0]
    //	    << ", r2 = " << r_ext[1] << ", r3 = " << r_ext[2] << endl;
    outfile << "#                               ";
    for (int i = 0; i < N_ext; i++)
      outfile << "r = " << std::left << setw(50) << r_ext[i];
    outfile << endl;
    outfile << "# " << setw(14) << "coord time"
	    << setw(16) << "proper time";
    for (int i = 0; i < N_ext; i++)
      outfile << setw(18) << "Re(psi4) " << setw(18) << "h_plus_dot"
	      << setw(18) << "h_plus"; 
    outfile << "\n# ==============================";
    for (int i = 0; i < N_ext; i++)
      outfile << "======================================================";
    outfile << endl;
    return 1;
  }
  //==========================================
  // write profile
  //==========================================
  int write_wave_data(double coord_time, double prop_time) {
    outfile << setw(16) << coord_time << setw(16) << prop_time;
    for (int i = 0; i < N_ext; i++) {
      outfile << setw(18) << psi4_Re_20[i]
	      << setw(18) << h_plus_dot_20[i]
	      << setw(18) << h_plus_20[i];
    }
    outfile << endl;
    return 0;
  }
  //==========================================
  // project psi4 to find l=2, m=0 modes
  //==========================================
  int Project() {
#ifndef AXISYMMETRY
    cout << " WAVEEXTRACTION implemented assuming Axisymmetry!!! " << endl;
#endif   /* AXISYMMETRY */
    //
    // contract with -2 Y_{20} = (3/2) \sqrt{5/\pi} \sin^2 \theta to find
    // 20 mode; assuming axisymmetry...
    //
    for (int i = 0; i < N_ext; i++) {
      double rl = r_ext[i];
      psi4_Re_20[i] = 0.0;
      for (int j = N_g; j < N_theta - N_g; j++) {
	double stl = grid->sintheta(j);
	double st2 = stl*stl;
	double dtheta = grid->delta_theta(j);
	double sY20 = st2;     // get to factors later...
	psi4_Re_20[i] += (*psi4_Re_3d)(rl, j, N_g) * sY20 * stl * dtheta; 
      }
      psi4_Re_20[i] *= 2.0 * PI * 1.5 * sqrt(5. / PI);
#ifdef EQSYMMETRY
      psi4_Re_20[i] *= 2.0;
#endif /* EQSYMMETRY */
    }
    return 1;
  }
  //==========================================
  // update functions (using simple first-order method, to avoid
  // having to evaluate curvature invariants at all intermediate steps...)
  //==========================================
  int Update(double dt) {
    //
    // First update (doesn't harm on first call in Initialize...)
    //
    for (int i = 0; i < N_ext; i++) {
      h_plus_dot_20[i] += psi4_Re_20[i] * dt;
      h_plus_20[i] += h_plus_dot_20[i] * dt;
    }
    //
    // then project new values of psi4...
    // 
    Project();
    //
    return 1;
  }
};

#endif  /* WAVEEXTRACTION_H */

