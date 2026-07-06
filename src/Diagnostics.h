// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing routines to compute diagnostics
//
//================================================
#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "gridfunction.h"
#include "Grid.h"
#include "ADM_Source_Terms.h"
#include "Fluxes.h"
#include "HorizonFinder.h"
#include "Monitor.h"

class diagnostics {
public:
  gf3d Ham;                     // Hamiltonian constraint violations
  gf3d Mom_r, Mom_t, Mom_p;     // components of momentum constraint violations (rescaled, indices down)
  gf3d CFC_r, CFC_t, CFC_p;     // violation of connection constraints (rescaled, indices down)
  gf3d R_prop, s_axial;         // proper distance of from origin, computed on current slice, along lines 
                                // of constant angle, and axial radius 
  gf3d E_rr, E_rt, E_rp, E_tt, E_tp, E_pp;  // electric and magnetic parts of Weyl tensor
  gf3d B_rr, B_rt, B_rp, B_tt, B_tp, B_pp;  // *not* rescaled
  gf3d I_Re, I_Im, J_Re, J_Im;  // curvature invariants, I_Re = R_{abcd} R^{abcd} / 16 
  gf3d zeta;  // invariant in axisymmetry; see Ledvinka & Khirnov 2021
  // for Schwarzschild I_Re = 3 M^2/R^6 where R is areal radius 
  // (see more examples in Burko et.al., PRD 73, 024002 (2006) (Note: "64" in C2 should be "96")
  gf3d psi4_Re, psi4_Im;        // computed from E_ij and B_ij; see eq. (8.6.17) in Alcubierre's book  
  double I_Re_max, I_Re_max_max, I_Im_max, I_Im_max_max;
  double J_Re_max, J_Re_max_max, J_Im_max, J_Im_max_max;
  double I_Re_int, J_Re_int, I_Im_int, J_Im_int;
  double zeta_max;
  bool new_max_last_timestep;
  int N_g, N_r, N_t, N_p;  // ghost zones, *total* number of grid points
  double r_max;
  Grid *grid;
  dumper *dump;
  Monitor *monitor;
  Cosmology * cosmology;
  HorizonFinder * horizonfinder;
  const char * name;
  int N_fcts, N_dump;
  gf3d ** fct_list;   // list of all grid functions in state
  gf3d ** dump_list;  // list of all grid functions to be dumped
  double PI;
public:
  //===============================================
  // Constructor
  //===============================================
  diagnostics(Grid *grid_i, dumper *dump_i, Monitor * monitor_i, 
	      Cosmology * cosmology_i, HorizonFinder * horizonfinder_i,
	      const char * name_i) :
    grid(grid_i), dump(dump_i), monitor(monitor_i), cosmology(cosmology_i), 
    horizonfinder(horizonfinder_i), name(name_i), new_max_last_timestep(false)
  {
    N_fcts = 28;
    fct_list = new gf3d*[N_fcts];
    dump_list = new gf3d*[N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    //
    //
    int gf_counter = 0;
    //
    // Hamiltonian constraint
    //
    fct_list[gf_counter] = Ham.setup(grid, 1, "Ham", gf_counter, +1, +1, +1);
    gf_counter++;
    //
    // momentum constraint
    //
    fct_list[gf_counter] = Mom_r.setup(grid, 2, "Mom_r", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = Mom_t.setup(grid, 2, "Mom_t", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = Mom_p.setup(grid, 2, "Mom_p", gf_counter,-1,-1,+1);
    gf_counter++;
    //
    // connection constraint
    //
    fct_list[gf_counter] = CFC_r.setup(grid, 2, "CFC_r", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = CFC_t.setup(grid, 2, "CFC_t", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = CFC_p.setup(grid, 2, "CFC_p", gf_counter,-1,-1,+1);
    gf_counter++;
    //
    // proper radius and axial radius
    //
    fct_list[gf_counter] = R_prop.setup(grid, -1, "R_prop", gf_counter,-1,+1,+1);
    gf_counter++;
    // NOTE: this is xi_phi / sin(theta) - therefore ax_par = +1, not -1
    fct_list[gf_counter] = s_axial.setup(grid, -1, "s_axial", gf_counter,-1,+1,+1);
    gf_counter++;
    //
    // Electric and magnetic parts of Weyl tensor; curvature invariant
    //
    fct_list[gf_counter] = E_rr.setup(grid, 1, "E_rr", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = E_rt.setup(grid, 1, "E_rt", gf_counter,-1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = E_rp.setup(grid, 1, "E_rp", gf_counter,+1,-1,+1);
    gf_counter++;
    fct_list[gf_counter] = E_tt.setup(grid, 1, "E_tt", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = E_tp.setup(grid, 1, "E_tp", gf_counter,-1,+1,-1);
    gf_counter++;
    fct_list[gf_counter] = E_pp.setup(grid, 1, "E_pp", gf_counter,+1,+1,+1);
    gf_counter++;
    // NOTE: B_ij is pseudo tensor (involves epsilon_{ijk}) - so symmetry
    // across equatorial plane is opposite that for normal rank-2 tensors...
    fct_list[gf_counter] = B_rr.setup(grid, 1, "B_rr", gf_counter,+1,+1,-1);
    gf_counter++;
    fct_list[gf_counter] = B_rt.setup(grid, 1, "B_rt", gf_counter,-1,-1,+1);
    gf_counter++;
    fct_list[gf_counter] = B_rp.setup(grid, 1, "B_rp", gf_counter,+1,-1,-1);
    gf_counter++;
    fct_list[gf_counter] = B_tt.setup(grid, 1, "B_tt", gf_counter,+1,+1,-1);
    gf_counter++;
    fct_list[gf_counter] = B_tp.setup(grid, 1, "B_tp", gf_counter,-1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = B_pp.setup(grid, 1, "B_pp", gf_counter,+1,+1,-1);
    gf_counter++;
    fct_list[gf_counter] = I_Re.setup(grid, 1, "I_Re", gf_counter,+1,+1,+1);
    gf_counter++;
    // Imaginary parts are odd power in B_{ij}, so inherit "pseudo"-symmetry
    // from B_{ij}
    fct_list[gf_counter] = I_Im.setup(grid, 1, "I_Im", gf_counter,+1,+1,-1);
    gf_counter++;
    fct_list[gf_counter] = J_Re.setup(grid, 1, "J_Re", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = J_Im.setup(grid, 1, "J_Im", gf_counter,+1,+1,-1);
    gf_counter++;
    fct_list[gf_counter] = psi4_Re.setup(grid, 1, "psi4_Re", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = psi4_Im.setup(grid, 1, "psi4_Im", gf_counter,+1,+1,+1);
    gf_counter++;
    fct_list[gf_counter] = zeta.setup(grid, 1, "zeta", gf_counter,+1,+1,+1);
    gf_counter++;
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN DIAGNOSTICS!!! " << endl;
    //
    N_g = grid->N_ghosts();
    N_r = grid->N_r_tot();
    N_t = grid->N_theta_tot();
    N_p = grid->N_phi_tot();
    r_max = grid->r_max();
    I_Re_max_max = I_Im_max_max = J_Re_max_max = J_Im_max_max = 0.0;
    PI = acos(-1.0);
    // let horizon finder know about electric and magnetic parts 
    // of Weyl tensor:
    horizonfinder->AssignEB(&E_rr, &E_rt, &E_rp, &E_tt, &E_tp, &E_pp,
			    &B_rr, &B_rt, &B_rp, &B_tt, &B_tp, &B_pp);

  };
  //===============================================
  // print list of all functions in state
  //===============================================
  void function_names() {
    cout << " STATE: List of all functions in diagnostics " << name << ": " << endl;
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
      cout << " DIAGNOSTICS: Dumping diagnostic functions " << name << " at time t = " << time  << endl;
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
      cerr << " DIAGNOSTICS: can't oppen " << dump_list_file << " for input." << endl;
      return 0;
    }
    cout << " DIAGNOSTICS: reading dump list from file " << dump_list_file << endl;
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
      cerr << " DIAGNOSTIC: found too many grid functions in assemble_dump_list() " << endl;
      N_dump = N_fcts;
    }
    return N_dump;
  };
  //===============================================
  // Destructor
  //===============================================
  ~diagnostics() {
    delete fct_list;
    delete dump_list;
    cout << " DIAGNOSTICS: ... closing diagnostics ... " << endl;
  };
  //===============================================
  // Evaluate Diagnostics
  //===============================================
  double Hamiltonian(state *s, curvature *curve, auxiliary *aux,
		     ADM_Source_Terms * adm,
		     bool excise = false, double frac_grid = 0.9,
		     double time = 0.0);
  // Note: Momentum constraint needs state t for temporary storage - use carefully!!
  double MomentumConstraint(state * s, curvature * curve, state * t, 
			    ADM_Source_Terms * adm, 
			    double & Mom_r_norm, double & Mom_t_norm, double & Mom_p_norm,
			    bool excise = false, double frac_grid = 0.9, 
			    double time = 0.0);
  double ConnectionFunctionConstraint(state * s, curvature *c,
				      double & CFC_r_norm, double & CFC_t_norm, double & CFC_p_norm);
  bool FindHorizon(int timestep, double t, double tau_c,
		   state * s, curvature * c,
		   ADM_Source_Terms * adm, Fluxes * fluxes, auxiliary * aux, 
		   double adm_mass, double mom_guess = 0, bool force = false);
  double ADM_Mass_Surface(state * s, curvature * c, auxiliary * aux, 
			  int i, double time = 0.0);
  double Angular_Momentum(state * s, curvature * c, auxiliary * aux,
			  int i, double time = 0.0);
  double Linear_Momentum(state * s, curvature * c, auxiliary * aux,
			  int i, double time = 0.0);
  void Compute_Proper_Radius(state *s, curvature *c, auxiliary *aux, int sigma);
  bool CurvatureInvariant(state *s, curvature *c, ADM_Source_Terms *adm,
			  double t, double tau_c, int step);
  void Note_Invariants(int timestep, double time, double prop_time, 
		       Monitor *monitor) {
    monitor->note_invariants(timestep, time, prop_time,
			     I_Re(0.0,N_g,N_g), I_Re_max, I_Re_int,
			     I_Im(0.0,N_g,N_g), I_Im_max, I_Im_int,
			     J_Re(0.0,N_g,N_g), J_Re_max, J_Re_int,
			     J_Im(0.0,N_g,N_g), J_Im_max, J_Im_int,
			     zeta(0.0,N_g,N_g), zeta_max);
  };
};

#endif /* DIAGNOSTICS_H */
