// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of matter, i.e. collection
// of grid functions for dynamical matter variables
//
//================================================
#ifndef HYDRO_STATE_H
#define HYDRO_STATE_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class hydro_state {
public:
  //
  // conserved fluid variables
  //  det_3D = \sqrt( \gamma / \hat \gamma )
  //
  //   D = det_3D W \rho_0
  //   S_r = det_3D W^2 \rho_0 h v_r
  //   S_t = det_3D W^2 \rho_0 h v_t / r 
  //   S_p = det_3D W^2 \rho_0 h v_p / (r sin \theta) 
  //   tau = det_3D (W^2 \rho_0 h - p) - W \rho_0)
  //
  gf3d D, tau, S_r, S_t, S_p;
  int N_fcts, N_dump;
  gf3d ** fct_list;   // list of all grid functions in hydro_state
  gf3d ** dump_list;  // list of all grid functions to be dumped
  Grid *grid;
  dumper *dump;
  const char * name;
public:
  //===============================================
  // Constructor
  //===============================================
  hydro_state(Grid *grid_i, dumper *dump_i, const char * name_i) :
    grid(grid_i), dump(dump_i), name(name_i)
  {
    N_fcts = 5;
    fct_list = new gf3d*[N_fcts];
    dump_list = new gf3d*[N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    //
    int gf_counter = 0;
    //
    fct_list[gf_counter] = D.setup(grid, 1, "D", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = tau.setup(grid, 1, "tau", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = S_r.setup(grid, 1, "S_r", gf_counter, -1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = S_t.setup(grid, 1, "S_t", gf_counter, +1, -1, -1);
    gf_counter++;
    fct_list[gf_counter] = S_p.setup(grid, 1, "S_p", gf_counter, -1, -1, +1);
    gf_counter++;
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN HYDRO_STATE!!! " << endl;
  };
  //===============================================
  // check properties
  //===============================================
  void check_props() {
    for (int i = 0; i < N_fcts; i++) {
      (*fct_list)[i].constants();
    }
  };
  //===============================================
  // fill ghosts
  //===============================================
  void fill_ghosts() {
    for (int i = 0; i < N_fcts; i++)
      (*fct_list)[i].fill_ghosts();
  };
  //===============================================
  // addition
  //===============================================
  void add(double factor, hydro_state * rhs) {
    for (int i = 0; i < N_fcts; i++)
      (*fct_list)[i].add(factor, rhs->fct_list[i]);
  };
  //===============================================
  // addition 
  //===============================================
  void add(hydro_state *hydro_state1, double factor, hydro_state * hydro_state2) {
    for (int i = 0; i < N_fcts; i++)
      (*fct_list)[i].add(hydro_state1->fct_list[i], factor, hydro_state2->fct_list[i]);
  };
  //===============================================
  // equals
  //===============================================
  void equals(hydro_state * rhs) {
    for (int i = 0; i < N_fcts; i++) {
      fct_list[i]->equals(rhs->fct_list[i]);
    }
  };
  void equals(double number) {
    for (int i = 0; i < N_fcts; i++) {
      fct_list[i]->equals(number);
    }
  };
  //===============================================
  // check whether states are equal
  //===============================================
  bool IsEqualTo(hydro_state * s) {
    bool fct_equal = true;
    bool state_equal = true;
    for (int i = 0; i < N_fcts; i++) {
       fct_equal = fct_list[i]->IsEqualTo(s->fct_list[i]);
       if (!fct_equal) {
	 cout << " HYDRO: Functions " << fct_list[i]->Name() 
	      << " in hydro_states "
	      << Name() << " and " << s->Name() << " are different!" << endl;
	 state_equal = false;
       }
    }
    if (state_equal == true)
      cout << " HYDRO: hydro_states " << Name() << " and " << s->Name()
	   << " are equal! " << endl; 
    return state_equal;
  };
  //===============================================
  // print list of all functions in hydro_state
  //===============================================
  void function_names() {
    cout << " HYDRO_STATE: List of all functions in hydro_state " << name << ": " << endl;
    for (int i = 0; i < N_fcts; i++)
      cout << "      " << fct_list[i]->Name() << endl;
  };
  //===============================================
  // return name of hydro_state
  //===============================================
  const char * Name() { return name; };
  //===============================================
  // Dump grid functions
  //===============================================
  void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
		 const char * suffix = "") {
    if (dump->time_to_dump(timestep)  || strcmp(suffix,"") ) {
      cout << " HYDRO_STATE: Dumping functions in hydro_state "
	   << name << " at time t = " << time  << endl;
      for (int i = 0; i < N_dump; i++) {
	dump->dump(time,prop_time,timestep,dump_list[i], suffix);
	dump->slice(time,prop_time,timestep,dump_list[i], suffix);
      }
    }
  };
  //===============================================
  // Find all functions to be dumped
  //===============================================
  int assemble_dump_list(const char * dump_list_file) {
    ifstream infile;
    infile.open(dump_list_file);
    if (!infile) {
      cerr << " HYDRO_STATE: can't oppen " << dump_list_file << " for input." << endl;
      return 0;
    }
    cout << " HYDRO_STATE: reading dump list from file " << dump_list_file << endl;
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
    if (N_dump > N_fcts) {
      cerr << " HYDRO: found too many grid functions in assemble_dump_list() " << endl;
      N_dump = N_fcts;
    }
    return N_dump;
  };
  //===============================================
  // Regrid
  //===============================================
  void Regrid(VecDoub r_new) {
    for (int i = 0; i < N_fcts; i++)
      fct_list[i]->Regrid(r_new);
  }
  //===============================================
  // Destructor
  //===============================================
  ~hydro_state() {
    delete fct_list;
    delete dump_list;
    cout << " HYDRO_STATE: ... closing hydro_state " << name << "... " << endl;
  };
};


#endif  /* HYDRO_STATE_H */
