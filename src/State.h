// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of system, i.e. collection
// of grid functions for dynamical variables
//
// Reference: Baumgarte, Montero, Cordero-Carrion & Mueller, PRD 87, 044026 (2013)  (BMCM)
//
//================================================
#ifndef STATE_H
#define STATE_H

#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"


class state {
public:
  // metric components (rescaled as in (BMCM.20))
  gf3d phi, h_rr, h_rt, h_rp, h_tt, h_tp, h_pp;
  // extrinsic curvature components (see (BMCM.21))
  gf3d K,   a_rr, a_rt, a_rp, a_tt, a_tp, a_pp;
  // connection functions (see (BMCM.22))
  gf3d lam_r, lam_t, lam_p;
  // Theta (for Z4)
  gf3d Theta;
  // lapse and shift (rescaled!), indices upstairs
  gf3d lapse, shift_r, shift_t, shift_p;
  gf3d B_r, B_t, B_p;
  // 
  int N_fcts, N_dump;
  gf3d ** fct_list;   // list of all grid functions in state
  gf3d ** dump_list;  // list of all grid functions to be dumped
  Grid *grid;
  dumper *dump;
  const char * name;
public:
  //===============================================
  // Constructor
  //===============================================
  state(Grid *grid_i, dumper *dump_i, const char * name_i) :
    grid(grid_i), dump(dump_i), name(name_i)
  {
    N_fcts = 25;
    fct_list = new gf3d*[N_fcts];
    dump_list = new gf3d*[N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    int gf_counter = 0;
    //
    // metric components, g_{ij} = exp(4 \phi) (1 + h_{ij}), RESCALED
    //
    fct_list[gf_counter] = phi.setup(grid,  1, "phi",  gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = h_rr.setup(grid, 1, "h_rr", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = h_rt.setup(grid, 1, "h_rt", gf_counter,-1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = h_rp.setup(grid, 1, "h_rp", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = h_tt.setup(grid, 1, "h_tt", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = h_tp.setup(grid, 1, "h_tp", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = h_pp.setup(grid, 1, "h_pp", gf_counter,+1,+1,+1);
    gf_counter++;
    //
    // extrinsic curvature components, A_{ij} = exp(4 \phi) a_{ij}, RESCALED
    //
    double K_alpha_speed = 1.0;
    fct_list[gf_counter] = K.setup(grid,    2, "K",    gf_counter,+1,+1,+1, K_alpha_speed);
    gf_counter++;
    fct_list[gf_counter] = a_rr.setup(grid, 1, "a_rr", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = a_rt.setup(grid, 1, "a_rt", gf_counter,-1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = a_rp.setup(grid, 1, "a_rp", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = a_tt.setup(grid, 1, "a_tt", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = a_tp.setup(grid, 1, "a_tp", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = a_pp.setup(grid, 1, "a_pp", gf_counter,+1,+1,+1);
    gf_counter++;
    //
    // connection functions, \Lambda^i, RESCALED
    //
    fct_list[gf_counter] = lam_r.setup(grid, 2, "lam_r", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = lam_t.setup(grid, 2, "lam_t", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = lam_p.setup(grid, 2, "lam_p", gf_counter,-1,-1,-1);
    gf_counter++;
    //
    // Theta function for Z4
    //
    fct_list[gf_counter] = Theta.setup(grid, 2, "Theta", gf_counter,+1,+1,+1);
    gf_counter++;
    //
    // lapse and shift
    // 
    double background = 1.0;
    fct_list[gf_counter] = lapse.setup(grid, 1, "lapse", gf_counter,+1,+1,+1,
				       K_alpha_speed, background);
    gf_counter++;
    fct_list[gf_counter] = shift_r.setup(grid, 2, "shift_r", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = shift_t.setup(grid, 2, "shift_t", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = shift_p.setup(grid, 2, "shift_p", gf_counter,-1,-1,+1);
    gf_counter++;
    fct_list[gf_counter] = B_r.setup(grid, 2, "B_r", gf_counter,-1,+1, -1);
    gf_counter++;
    fct_list[gf_counter] = B_t.setup(grid, 2, "B_t", gf_counter,+1,-1,+1);
    gf_counter++;
    fct_list[gf_counter] = B_p.setup(grid, 2, "B_p", gf_counter,-1,-1,+1);
    gf_counter++;
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN STATE!!! " << endl;
  };
  //===============================================
  // check properties
  //===============================================
  void check_props() {
    for (int i = 0; i < N_fcts; i++) {
      (*fct_list)[i].constants();
    }
      //      (*fct_list)[i].fill_ghosts();
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
  void add(double factor, state * rhs) {
    for (int i = 0; i < N_fcts; i++)
      (*fct_list)[i].add(factor, rhs->fct_list[i]);
  };
  //===============================================
  // addition 
  //===============================================
  void add(state *state1, double factor, state * state2) {
    for (int i = 0; i < N_fcts; i++)
      (*fct_list)[i].add(state1->fct_list[i], factor, state2->fct_list[i]);
  };
  //===============================================
  // equals
  //===============================================
  void equals(state * rhs) {
    for (int i = 0; i < N_fcts; i++) {
      fct_list[i]->equals(rhs->fct_list[i]);
    }
  };
  //===============================================
  // equals
  //===============================================
  void equals(double number) {
    for (int i = 0; i < N_fcts; i++) {
      fct_list[i]->equals(number);
    }
  };
  //===============================================
  // check whether states are equal
  //===============================================
  bool IsEqualTo(state * s) {
    bool fct_equal = true;
    bool state_equal = true;
    for (int i = 0; i < N_fcts; i++) {
       fct_equal = fct_list[i]->IsEqualTo(s->fct_list[i]);
       if (!fct_equal) {
	 cout << " STATE: Functions " << fct_list[i]->Name() << " in states "
	      << Name() << " and " << s->Name() << " are different!" << endl;
	 state_equal = false;
       }
    }
    if (state_equal == true)
      cout << " STATE: states " << Name() << " and " << s->Name()
	   << " are equal! " << endl; 
    return state_equal;
  };
  //===============================================
  // apply characteristic outer boundaries
  //===============================================
  void char_OB(state * last, double dt) {
    for (int i = 0; i < N_fcts; i++)
      fct_list[i]->fill_outerboundary(*last->fct_list[i], dt);
  }
  //===============================================
  // check whether functions are finite
  //===============================================
  bool FINITE() {
    bool fine = true;
    for (int i = 0; i < N_fcts; i++) {
      if (!fct_list[i]->FINITE() ) {
	cout << " STATE: function " << fct_list[i]->Name() << " in state "
	     << name << " is not finite!" << endl;
	fine = false;
      }
    }
    return fine;
  };
  //===============================================
  // print list of all functions in state
  //===============================================
  void function_names() {
    cout << " STATE: List of all functions in state " << name << ": " << endl;
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
      cout << " STATE: Dumping functions in state " << name
	   << " at time t = " << time  << endl;
      for (int i = 0; i < N_dump; i++) {
	dump->dump(time,prop_time,timestep,dump_list[i],suffix);
	dump->slice(time,prop_time,timestep,dump_list[i],suffix);
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
      cerr << " STATE: can't open " << dump_list_file << " for input." << endl;
      return 0;
    }
    cout << " STATE: reading dump list from file " << dump_list_file << endl;
    N_dump = 0;
    char fct_name[64];
    infile >> fct_name;
    while (!infile.eof()) {
      for (int i = 0; i < N_fcts; i++) {
	if (!strcmp(fct_name, fct_list[i]->Name())) {
	  cout << " ... found function name " << fct_name << endl;
	  dump_list[N_dump] = fct_list[i]->Address();
	  //	  cout << " function " << fct_name << " : " <<  dump_list[N_dump] << "  " << fct_list[i] << "  " << fct_list[i]->Address() << endl;

	  N_dump++;
	}
      }
      infile >> fct_name;
    }
    if (N_dump > N_fcts) {
      cerr << " STATE: found too many grid functions in assemble_dump_list() " << endl;
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
  ~state() {
    delete fct_list;
    delete dump_list;
    cout << " STATE: ... closing state " << name << "... " << endl;
  };
};

#endif  /* STATE_H */
