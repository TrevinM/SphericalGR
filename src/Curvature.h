// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary variables
//
//================================================
#ifndef CURVATURE_H
#define CURVATURE_H

#include "gridfunction.h"
#include "Grid.h"

class curvature {
public:
  gf3d det, det_init, trace_A;     // determinant of metric
  // covariant derivatives of conformal metric (with respect
  // to flat background metric)
  gf3d Dr_e_rr, Dr_e_rt, Dr_e_rp, Dr_e_tt, Dr_e_tp, Dr_e_pp;
  gf3d Dt_e_rr, Dt_e_rt, Dt_e_rp, Dt_e_tt, Dt_e_tp, Dt_e_pp;
  gf3d Dp_e_rr, Dp_e_rt, Dp_e_rp, Dp_e_tt, Dp_e_tp, Dp_e_pp;
  // relative connection symbols
  // \Delta \Gamma^i_{jk} \equiv \bar \Gamma^i_{jk} - Background \Gamma^i_{jk} 
  gf3d DG_r_rr, DG_r_rt, DG_r_rp, DG_r_tt, DG_r_tp, DG_r_pp;
  gf3d DG_t_rr, DG_t_rt, DG_t_rp, DG_t_tt, DG_t_tp, DG_t_pp;
  gf3d DG_p_rr, DG_p_rt, DG_p_rp, DG_p_tt, DG_p_tp, DG_p_pp;
  // trace of above connection functions
  gf3d DG_r, DG_t, DG_p;
  // inverse conformal metric
  gf3d gup_rr, gup_rt, gup_rp, gup_tt, gup_tp, gup_pp;
  // Ricci tensor
  gf3d R_rr, R_rt, R_rp, R_tt, R_tp, R_pp, trace_R;  
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
  curvature(Grid *grid_i, dumper *dump_i, const char * name_i) :
    grid(grid_i), dump(dump_i), name(name_i)
  {
    N_fcts = 55;
    fct_list = new gf3d*[N_fcts];
    dump_list = new gf3d*[N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    //
    //
    int gf_counter = 0;
    //
    fct_list[gf_counter] = det.setup(grid, 1, "det", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = det_init.setup(grid, 1, "det_init", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = trace_A.setup(grid, 1, "trace_A", gf_counter, +1, +1, +1);
    gf_counter++;
    //
    // Covariant derivatives of conformal metric, not rescaled!
    // Don't need parities since fill_ghosts() not called for these functions
    // Also, fall_off irrelevant since outer_boundary never called for these functions
    //
    fct_list[gf_counter] = Dr_e_rr.setup(grid, 1, "Dr_e_rr", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dr_e_rt.setup(grid, 1, "Dr_e_rt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dr_e_rp.setup(grid, 1, "Dr_e_rp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dr_e_tt.setup(grid, 1, "Dr_e_tt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dr_e_tp.setup(grid, 1, "Dr_e_tp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dr_e_pp.setup(grid, 1, "Dr_e_pp", gf_counter);
    gf_counter++;

    fct_list[gf_counter] = Dt_e_rr.setup(grid, 1, "Dt_e_rr", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dt_e_rt.setup(grid, 1, "Dt_e_rt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dt_e_rp.setup(grid, 1, "Dt_e_rp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dt_e_tt.setup(grid, 1, "Dt_e_tt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dt_e_tp.setup(grid, 1, "Dt_e_tp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dt_e_pp.setup(grid, 1, "Dt_e_pp", gf_counter);
    gf_counter++;

    fct_list[gf_counter] = Dp_e_rr.setup(grid, 1, "Dp_e_rr", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dp_e_rt.setup(grid, 1, "Dp_e_rt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dp_e_rp.setup(grid, 1, "Dp_e_rp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dp_e_tt.setup(grid, 1, "Dp_e_tt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dp_e_tp.setup(grid, 1, "Dp_e_tp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = Dp_e_pp.setup(grid, 1, "Dp_e_pp", gf_counter);
    gf_counter++;
    //
    // relative connection symbols
    // \Delta \Gamma^i_{jk} \equiv \bar \Gamma^i_{jk} - Background \Gamma^i_{jk} 
    // 
    fct_list[gf_counter] = DG_r_rr.setup(grid, 1, "DG_r_rr", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_r_rt.setup(grid, 1, "DG_r_rt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_r_rp.setup(grid, 1, "DG_r_rp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_r_tt.setup(grid, 1, "DG_r_tt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_r_tp.setup(grid, 1, "DG_r_tp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_r_pp.setup(grid, 1, "DG_r_pp", gf_counter);
    gf_counter++;

    fct_list[gf_counter] = DG_t_rr.setup(grid, 1, "DG_t_rr", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_t_rt.setup(grid, 1, "DG_t_rt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_t_rp.setup(grid, 1, "DG_t_rp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_t_tt.setup(grid, 1, "DG_t_tt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_t_tp.setup(grid, 1, "DG_t_tp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_t_pp.setup(grid, 1, "DG_t_pp", gf_counter);
    gf_counter++;

    fct_list[gf_counter] = DG_p_rr.setup(grid, 1, "DG_p_rr", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_p_rt.setup(grid, 1, "DG_p_rt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_p_rp.setup(grid, 1, "DG_p_rp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_p_tt.setup(grid, 1, "DG_p_tt", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_p_tp.setup(grid, 1, "DG_p_tp", gf_counter);
    gf_counter++;
    fct_list[gf_counter] = DG_p_pp.setup(grid, 1, "DG_p_pp", gf_counter);
    gf_counter++;
    //
    // trace of above connection functions (should be equal to lambdas, as long 
    // as constraint holds...)
    //
    fct_list[gf_counter] = DG_r.setup(grid, 2, "DG_r", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = DG_t.setup(grid, 2, "DG_t", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = DG_p.setup(grid, 2, "DG_p", gf_counter,-1,-1,+1);
    gf_counter++;
    //
    // inverse metric 
    //
    fct_list[gf_counter] = gup_rr.setup(grid, 1, "gup_rr", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = gup_rt.setup(grid, 1, "gup_rt", gf_counter, -1, -1, -1);
    gf_counter++;
    fct_list[gf_counter] = gup_rp.setup(grid, 1, "gup_rp", gf_counter, +1, -1, +1);
    gf_counter++;
    fct_list[gf_counter] = gup_tt.setup(grid, 1, "gup_tt", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = gup_tp.setup(grid, 1, "gup_tp", gf_counter, -1, +1, -1);
    gf_counter++;
    fct_list[gf_counter] = gup_pp.setup(grid, 1, "gup_pp", gf_counter, +1, +1, +1);
    gf_counter++;
    //
    // Ricci tensor (indices down)
    //
    const int third_order_fo = 3;
    fct_list[gf_counter] = R_rr.setup(grid, third_order_fo, "R_rr", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = R_rt.setup(grid, third_order_fo, "R_rt", gf_counter, -1, -1, -1);
    gf_counter++;
    fct_list[gf_counter] = R_rp.setup(grid, third_order_fo, "R_rp", gf_counter, +1, -1, +1);
    gf_counter++;
    fct_list[gf_counter] = R_tt.setup(grid, third_order_fo, "R_tt", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = R_tp.setup(grid, third_order_fo, "R_tp", gf_counter, -1, +1, -1);
    gf_counter++;
    fct_list[gf_counter] = R_pp.setup(grid, third_order_fo, "R_pp", gf_counter, +1, +1, +1);
    gf_counter++;
    fct_list[gf_counter] = trace_R.setup(grid, third_order_fo, "trace_R", gf_counter, +1, +1, +1);
    gf_counter++;
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN CURVATURE!!! " << endl;

    N_g = grid->N_ghosts();
    N_r = grid->N_r_tot();
    N_t = grid->N_theta_tot();
    N_p = grid->N_phi_tot();
    
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
  // check whether functions are finite
  //===============================================
  bool FINITE() {
    bool fine = true;
    for (int i = 0; i < N_fcts; i++) {
      if (!fct_list[i]->FINITE() ) {
	cout << " CURVATURE: function " << fct_list[i]->Name() << " in "
	     << name << " is not finite!" << endl;
	fine = false;
      }
    }
    return fine;
  };
  //===============================================
  // Dump grid functions
  //===============================================
  void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
		 const char * suffix = "") {
    if (dump->time_to_dump(timestep) || strcmp(suffix,"") ) {
      cout << " CURVATURE: Dumping curvature functions " << name << " at time t = " << time  << endl;
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
      cerr << " CURVATURE: can't oppen " << dump_list_file << " for input." << endl;
      return 0;
    }
    cout << " CURVATURE: reading dump list from file " << dump_list_file << endl;
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
      cerr << " CURVATURE: found too many grid functions in assemble_dump_list() " << endl;
      N_dump = N_fcts;
    }
    return N_dump;
  };
  //===============================================
  // Destructor
  //===============================================
  ~curvature() {
    delete fct_list;
    delete dump_list;
    cout << " CURVATURE: ... closing curvature ... " << endl;
  };
  //===============================================
  // Updates
  //===============================================
  void Compute_Metric_Derivatives(state *c);
  void Compute_Inverse_Metric(state *c);
  int Compute_Determinant(state *c);
  void Compute_Connection();
  void Compute_Ricci(state *c);
  void Compute_Trace(state *c);
  void Update(state * c) {
    Compute_Metric_Derivatives(c);
    Compute_Inverse_Metric(c);
    Compute_Connection();
    Compute_Ricci(c);
  };
  void InitDet() {
    det_init.equals(&det);
  }
};

#endif /* CURVATURE_H */
