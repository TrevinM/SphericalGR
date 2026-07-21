// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary variables
//
//================================================
#ifndef ADM_SOURCE_TERMS_H
#define ADM_SOURCE_TERMS_H

#include "gridfunction.h"
#include "Grid.h"

class ADM_Source_Terms {
public:
  gf3d rho_ADM;
  //
  // Note: all lower indices, NOT rescaled!
  // 
  gf3d S_r, S_t, S_p;
  // 
  gf3d S_rr, S_rt, S_rp, S_tt, S_tp, S_pp, trace_S;
  //
  int N_g, N_r, N_t, N_p;  // ghost zones, *total* number grid points
  Grid *grid;
  dumper *dump;
  const char * name;
  int N_fcts, N_dump;
  gf3d ** fct_list;   // list of all grid functions in state
  gf3d ** dump_list;  // list of all grid functions to be dumped
public:
  //===============================================
  // Constructor
  //===============================================
  ADM_Source_Terms(Grid *grid_i, dumper *dump_i, const char * name_i) :
    grid(grid_i), dump(dump_i), name(name_i)
  {
    N_fcts = 11;
    fct_list = new gf3d*[N_fcts];
    dump_list = new gf3d*[N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    //
    //
    int gf_counter = 0;
    //
    fct_list[gf_counter] = rho_ADM.setup(grid, 1, "rho_ADM", gf_counter, +1, +1, +1);
    gf_counter++;
    //
    fct_list[gf_counter] = S_r.setup(grid, 2, "S_r", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = S_t.setup(grid, 2, "S_t", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = S_p.setup(grid, 2, "S_p", gf_counter,-1,-1,+1);
    gf_counter++;
    //
    // inverse metric 
    //
    fct_list[gf_counter] = S_rr.setup(grid, 1, "S_rr", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = S_rt.setup(grid, 1, "S_rt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = S_rp.setup(grid, 1, "S_rp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = S_tt.setup(grid, 1, "S_tt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = S_tp.setup(grid, 1, "S_tp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = S_pp.setup(grid, 1, "S_pp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = trace_S.setup(grid, 1, "trace_S", gf_counter);
    gf_counter++;
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN ADM_SOURCE_TERMS!!! " << endl;

    N_g = grid->N_ghosts();
    N_r = grid->N_r_tot();
    N_t = grid->N_theta_tot();
    N_p = grid->N_phi_tot();
    
  };
  //===============================================
  // print list of all functions in state
  //===============================================
  void function_names() {
    cout << " STATE: List of all functions in ADM_Source_Terms " << name << ": " << endl;
    for (int i = 0; i < N_fcts; i++)
      cout << "      " << fct_list[i]->Name() << endl;
  };
  //===============================================
  // return name of state
  //===============================================
  const char * Name() { return name; };
  //===============================================
  // Dump grid functions
  //===============================================
  void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
		 const char * suffix = "") {
    if (dump->time_to_dump(timestep) || strcmp(suffix,"") ) {
      cout << " ADM_SOURCE_TERMS: Dumping source terms " << name
	   << " at time t = " << time  << endl;
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
      cerr << " ADM_SOURCE_TERMS: can't open " << dump_list_file << " for input." << endl;
      return 0;
    }
    cout << " ADM_SOURCE_TERMS: reading dump list from file " << dump_list_file << endl;
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
  ~ADM_Source_Terms() {
    delete fct_list;
    delete dump_list;
    cout << " ADM_SOURCE_TERMS: ... closing adm_source_terms ... " << endl;
  };
  //===============================================
  // Updates
  //===============================================
};

#endif /* ADM_SOURCE_TERMS_H */
