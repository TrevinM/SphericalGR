// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing fluxes - needed to compute accretion rate onto
// black holes.  Currently implemented for Hydro only!
//
//================================================
#ifndef FLUXES_H
#define FLUXES_H

#include "gridfunction.h"
#include "Grid.h"

class Fluxes {
public:

  // Note: all uppder indices, NOT rescaled!
  // 
  gf3d j_r, j_t, j_p;
  // 
  int N_g, N_r, N_t, N_p;  // ghost zones, *total* number grid points
  Grid *grid;
  dumper *dump;
  const char * name;
  int N_fcts, N_dump;
  gf3d ** fct_list;   // list of all grid functions in fluxes
  gf3d ** dump_list;  // list of all grid functions to be dumped
public:
  //===============================================
  // Constructor
  //===============================================
  Fluxes(Grid *grid_i, dumper *dump_i, const char * name_i) :
    grid(grid_i), dump(dump_i), name(name_i)
  {
    N_fcts = 3;
    fct_list = new gf3d*[N_fcts];
    dump_list = new gf3d*[N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    //
    //
    int gf_counter = 0;
    //
    fct_list[gf_counter] = j_r.setup(grid, 2, "j_r", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = j_t.setup(grid, 2, "j_t", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = j_p.setup(grid, 2, "j_p", gf_counter,-1,-1,+1);
    gf_counter++;
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN FLUXES!!! " << endl;
    //
    N_g = grid->N_ghosts();
    N_r = grid->N_r_tot();
    N_t = grid->N_theta_tot();
    N_p = grid->N_phi_tot();
    
  };
  //===============================================
  // print list of all functions in state
  //===============================================
  void function_names() {
    cout << " FLUXES: List of all functions in Fluxes " << name << ": " << endl;
    for (int i = 0; i < N_fcts; i++)
      cout << "      " << fct_list[i]->Name() << endl;
  };
  //===============================================
  // return name of class
  //===============================================
  const char * Name() { return name; };
  //===============================================
  // Dump grid functions
  //===============================================
  void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
		 const char * suffix = "") {
    if (dump->time_to_dump(timestep) || strcmp(suffix,"") ) {
      cout << " FLUXES: Dumping source terms " << name
	   << " at time t = " << time  << endl;
      for (int i = 0; i < N_dump; i++)
	dump->dump(time,prop_time,timestep,dump_list[i], suffix);
    }
  };
  //===============================================
  // Find all functions to be dumped
  //===============================================
  int assemble_dump_list(const char * dump_list_file) {
    ifstream infile;
    infile.open(dump_list_file);
    if (!infile) {
      cerr << " FLUXES: can't oppen " << dump_list_file << " for input." << endl;
      return 0;
    }
    cout << " FLUXES: reading dump list from file " << dump_list_file << endl;
    N_dump = 0;
    char fct_name[64];
    infile >> fct_name;
    while (!infile.eof()) {
      for (int i = 0; i < N_fcts; i++) {
	if (!strcmp(fct_name, fct_list[i]->Name())) {
	  cout << " ... found function name " << fct_name << endl;
	  dump_list[N_dump] = fct_list[i]->Address();
	  N_dump++;
	}
      }
      infile >> fct_name;
    }
    return N_dump;
  };
  //===============================================
  // Destructor
  //===============================================
  ~Fluxes() {
    delete fct_list;
    delete dump_list;
    cout << " FLUXES: ... closing fluxes ... " << endl;
  };
  //===============================================
  // Updates
  //===============================================
};

#endif /* FLUXES_H */
