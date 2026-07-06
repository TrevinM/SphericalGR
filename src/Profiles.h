// Tell emacs that this is -*-c++-*- mode
// 
//================================================
//
// Contains all methods for creating spacetime profiles
//
//================================================
//
#ifndef PROFILES_H
#define PROFILES_H

#include "Grid.h"
#include "gridfunction.h"
#include "Slicing.h"
#include "InData.h"
#include "Gauge.h"
#include "Cosmology.h"

class Profiles {
 private:
  Grid *grid;
  InData *indata;
  Slicing *slicing;
  Gauge *gauge;
  Cosmology *cosmology;
  int N_fcts;     // actual number
  int N_fcts_max; // maximum number - set below...
  int N_r, N_theta, N_phi, N_g;
  ofstream outfile_ax, outfile_eq;
  gf3d ** fct_list;
  gf3d * inv_radius;
  const char * filestem;
  int note_step;
  double PI, PI_half;
public:
  //==========================================
  // Constructor
  //==========================================
  Profiles(const char * filename_i, Grid *grid_i, int note_step_i, InData *indata_i, Slicing *slicing_i,
	   Gauge *gauge_i, Cosmology *cosmology_i) :
    grid(grid_i), indata(indata_i), slicing(slicing_i), gauge(gauge_i),
    cosmology(cosmology_i), note_step(note_step_i), filestem(filename_i)
  {
    N_fcts = 0;
    N_fcts_max = 20;
    fct_list = new gf3d*[N_fcts_max];
    N_r = grid->N_r_tot();
    N_theta = grid->N_theta_tot();
    N_phi = grid->N_phi_tot();
    N_g = grid->N_ghosts();
    PI = acos(-1.0);
    PI_half = PI / 2.0;
    inv_radius = NULL;
  };
  //==========================================
  // Destructor
  //==========================================
  ~Profiles() {
    delete fct_list;
    if (outfile_ax) outfile_ax.close();
    if (outfile_eq) outfile_eq.close();
    cout << " PROFILES: destructing profiles... " << endl;
  }
  //==========================================
  // initialize profile files
  //==========================================
  int initialize_file(gf3d * prop_radius) {
    inv_radius = prop_radius;
    ostringstream profilename_ax, profilename_eq;
    //
    // first for axis...
    //
    profilename_ax << "output/" << filestem << "_" << N_r - 2*N_g << "_"
		<< N_theta - 2*N_g << ".axis_profiles" << ends;
    outfile_ax.open(profilename_ax.str().c_str());
    if (outfile_ax)
      cout << " PROFILES: Opened profile file " << profilename_ax.str().c_str()
	   << " for output " << endl;
    else
      cerr << " PROFILES: Could not open file " << profilename_ax.str().c_str()
	   << endl;
    outfile_ax << "# for " << indata->Name() << " initial data " 
	    << endl;
    outfile_ax << "# evolved with (" << N_r << "," << N_theta
	    << "," << N_phi << ") gridpoints (including "
	    << N_g << " ghosts)" << endl; 
    outfile_ax << "# Evolution with " << slicing->Name() << " and " 
	    << gauge->Name() << endl;
    outfile_ax << "# Plotting data every " << note_step << " steps "
	    << endl;
    outfile_ax << "# " << setw(14) << "coord time"
	    << setw(16) << "proper time" << setw(20)
	    << "coord radius" << setw(20) << inv_radius->Name(); 
    for (int l = 0; l < N_fcts; l++)
      outfile_ax << setw(20) << fct_list[l]->Name();
    outfile_ax << endl;
#ifndef EQSYMMETRY
    //      for (int i = N_r - N_g - 1; i >= N_g; i--)
    //	outfile_ax << setprecision(8) << setw(16) << - grid->r(i);
#endif
    //      for (int i = N_g; i < N_r - N_g; i++)
    //	outfile_ax << setprecision(8) << setw(16) << grid->r(i);
    outfile_ax << "\n# ==================================================";
      outfile_ax << "====================";
    for (int l = 0; l < N_fcts; l++)
      outfile_ax << "====================";
    outfile_ax << endl;
    //
    // ... then again for equation:
    //
    profilename_eq << "output/" << filestem << "_" << N_r - 2*N_g << "_"
		<< N_theta - 2*N_g << ".eq_profiles" << ends;
    outfile_eq.open(profilename_eq.str().c_str());
    if (outfile_eq)
      cout << " PROFILES: Opened profile file " << profilename_eq.str().c_str()
	   << " for output " << endl;
    else
      cerr << " PROFILES: Could not open file " << profilename_eq.str().c_str()
	   << endl;
    outfile_eq << "# for " << indata->Name() << " initial data " 
	    << endl;
    outfile_eq << "# evolved with (" << N_r << "," << N_theta
	    << "," << N_phi << ") gridpoints (including "
	    << N_g << " ghosts)" << endl; 
    outfile_eq << "# Evolution with " << slicing->Name() << " and " 
	    << gauge->Name() << endl;
    outfile_eq << "# Plotting data every " << note_step << " steps "
	    << endl;
    outfile_eq << "# " << setw(14) << "coord time"
	    << setw(16) << "proper time" << setw(20)
	    << "coord radius" << setw(20) << inv_radius->Name(); 
    for (int l = 0; l < N_fcts; l++)
      outfile_eq << setw(20) << fct_list[l]->Name();
    outfile_eq << endl;
#ifndef EQSYMMETRY
    //      for (int i = N_r - N_g - 1; i >= N_g; i--)
    //	outfile_eq << setprecision(8) << setw(16) << - grid->r(i);
#endif
    //      for (int i = N_g; i < N_r - N_g; i++)
    //	outfile_eq << setprecision(8) << setw(16) << grid->r(i);
    outfile_eq << "\n# ==================================================";
      outfile_eq << "====================";
    for (int l = 0; l < N_fcts; l++)
      outfile_eq << "====================";
    outfile_eq << endl;
    return 1;
  }
  //==========================================
  // write profile
  //==========================================
  int write_profile(double coord_time, double prop_time) {
      //      outfile_ax << setw(16) << coord_time << setw(16) << prop_time;
    //
    // first on axis... 
    //
#ifndef EQSYMMETRY
    for (int i = N_r - N_g - 1; i >= N_g; i--) {
      outfile_ax << setprecision(8) << setw(16) << coord_time 
	      << setw(16) << prop_time
	      << setw(20) << - grid->r(i)
	      << setw(20) << - (*inv_radius)(i,PI,N_g);
      for (int l = 0; l < N_fcts; l++)
	outfile_ax << setw(20) << (*fct_list[l])(i,PI,N_g);
      outfile_ax << endl;
    }
#endif
    for (int i = N_g; i < N_r - N_g; i++) {
      outfile_ax << setprecision(8) << setw(16) << coord_time
	      << setw(16) << prop_time
	      << setw(20) << grid->r(i)
	      << setw(20) << (*inv_radius)(i,0.0,N_g);
      for (int l = 0; l < N_fcts; l++)
	outfile_ax << setw(20) << (*fct_list[l])(i,0.0,N_g);
      outfile_ax << endl;
    }
    outfile_ax << endl;
    //
    // ... then again on equator:
    //
    for (int i = N_g; i < N_r - N_g; i++) {
      outfile_eq << setprecision(8) << setw(16) << coord_time
	      << setw(16) << prop_time
	      << setw(20) << grid->r(i)
	      << setw(20) << (*inv_radius)(i,PI_half,N_g);
      for (int l = 0; l < N_fcts; l++)
	outfile_eq << setw(20) << (*fct_list[l])(i,PI_half,N_g);
      outfile_eq << endl;
    }
    outfile_eq << endl;
    //
    return 0;
  }
  //==========================================
  // assemble profile list
  //==========================================
  template <class bundle>
  int add_to_profile_list(const char * profile_list_file, bundle *s) {
    ifstream infile;
    infile.open(profile_list_file);
    if (!infile) {
      cerr << " PROFILES: can't open " << profile_list_file
	   << " for input." << endl;
      return 0;
    }
    cout << " PROFILES: looking for functions listed in "
	 << profile_list_file << " in " << s->Name() << " ..." << endl;
    char fct_name[64];
    infile >> fct_name;
    while (!infile.eof() && N_fcts < N_fcts_max) {
      for (int l = 0; l < s->N_fcts; l++) {
	if (!strcmp(fct_name, s->fct_list[l]->Name())) {
	  cout << " ... found match for function name " << fct_name 
	       << " in " << s->Name() << endl;
	  fct_list[N_fcts] = s->fct_list[l]->Address();
	  N_fcts++;
	}
      }
      infile >> fct_name;
    }
    if (N_fcts == N_fcts_max) {
      cout << " PROFILES: ran out of room for functions to track - " 
	   << " increase N_fcts_max! " << endl; 
    }
    return N_fcts;
  }

};

#endif  /* PROFILES_H */

