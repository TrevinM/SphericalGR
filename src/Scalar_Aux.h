// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary functions for scalar
//
//================================================
#ifndef SCALAR_AUX_H
#define SCALAR_AUX_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class scalar_aux {
public:
  //
  // Note: will add functions to scalar_aux through derived classes
  //
  gf3d Omega, gradient_parts, pi_parts, potential_parts, K_parts, A2_parts, R_parts, X;
  //
  int N_fcts, N_dump, n_r, n_theta, n_phi;
  gf3d ** fct_list;   // list of all grid functions in scalar_aux
  gf3d ** dump_list;  // list of all grid functions to be dumped
  Grid *grid;
  dumper *dump;
  const char * name;
public:
  //===============================================
  // Constructor
  //===============================================
  scalar_aux(Grid *grid_i, dumper *dump_i, const char * name_i) :
    grid(grid_i), dump(dump_i), name(name_i)
  {
    N_fcts = 8;
    fct_list = new gf3d*[N_fcts];
    dump_list = new gf3d*[N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    //
    int gf_counter = 0;
    //
    fct_list[gf_counter] = Omega.setup(grid, 1, "Omega", gf_counter, -1, -1, 1);
    gf_counter++;

    fct_list[gf_counter] = gradient_parts.setup(grid, 1, "gradient_parts", gf_counter, +1, +1, +1);
    gf_counter++;

    fct_list[gf_counter] = pi_parts.setup(grid, 1, "pi_parts", gf_counter, +1, +1, +1);
    gf_counter++;

    fct_list[gf_counter] = potential_parts.setup(grid, 1, "potential_parts", gf_counter, +1, +1, +1);
    gf_counter++;

    fct_list[gf_counter] = K_parts.setup(grid, 1, "K_parts", gf_counter, +1, +1, +1);
    gf_counter++;

    fct_list[gf_counter] = A2_parts.setup(grid, 1, "A2_parts", gf_counter, +1, +1, +1);
    gf_counter++;

    fct_list[gf_counter] = R_parts.setup(grid, 1, "R_parts", gf_counter, +1, +1, +1);
    gf_counter++;

    fct_list[gf_counter] = X.setup(grid, 1, "X", gf_counter, +1, +1, +1);
    gf_counter++;

    n_r = Omega.dim1();
		n_theta = Omega.dim2();
		n_phi = Omega.dim3();
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN SCALAR_AUX!!! " << endl;
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
  void add(double factor, scalar_aux * rhs) {
    for (int i = 0; i < N_fcts; i++)
      (*fct_list)[i].add(factor, rhs->fct_list[i]);
  };
  //===============================================
  // addition 
  //===============================================
  void add(scalar_aux *scalar_aux1, double factor, scalar_aux * scalar_aux2) {
    for (int i = 0; i < N_fcts; i++)
      (*fct_list)[i].add(scalar_aux1->fct_list[i], factor, scalar_aux2->fct_list[i]);
  };
  //===============================================
  // addition
  //===============================================
  void equals(scalar_aux * rhs) {
    for (int i = 0; i < N_fcts; i++) {
      fct_list[i]->equals(rhs->fct_list[i]);
    }
  };
  //===============================================
  // print list of all functions in scalar_aux
  //===============================================
  void function_names() {
    cout << " SCALAR_AUX: List of all functions in scalar_aux " << name << ": " << endl;
    for (int i = 0; i < N_fcts; i++)
      cout << "      " << fct_list[i]->Name() << endl;
  };
  //===============================================
  // return name of scalar_aux
  //===============================================
  const char * Name() { return name; };
  //===============================================
  // Dump grid functions
  //===============================================
  void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
		 const char * suffix = "") {
    if (dump->time_to_dump(timestep)) {
      cout << " SCALAR_AUX: Dumping functions in scalar_aux "
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
      cerr << " SCALAR_AUX: can't oppen " << dump_list_file << " for input." << endl;
      return 0;
    }
    cout << " SCALAR_AUX: reading dump list from file " << dump_list_file << endl;
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
  // Compute X
  //===============================================

  void Compute_X(state * s, double a_friedmann) {
    for (int i = 0; i < n_r; i++)    
      for (int j = 0; j < n_theta; j++)
	      for (int k = 0; k < n_phi; k++) {
	        X[i][j][k] = exp(-2.0*s->phi(i,j,k)) / a_friedmann;
	  }
  }
  //===============================================
  // Destructor
  //===============================================
  ~scalar_aux() {
    delete fct_list;
    delete dump_list;
    cout << " SCALAR_AUX: ... closing scalar_aux " << name << "... " << endl;
  };
};


#endif  /* SCALAR_AUX_H */
