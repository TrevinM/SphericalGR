// Tell emacs that this is -*-c++-*- mode
//
#include "Scalar_State.h"
#include "Scalar_Aux.h"
#include "Monitor.h"
#include "Potential.h"

//================================================
//
// Scalar Field
//
//================================================
class ScalarFieldWithPotential : public Matter {
public:
  double rho_center, rho_c_max, drhoddr;    // diagnostics 
  double sf_center_new, sf_center_old, sf_center_previous;
  scalar_state * last, * derivs, * inter, * updates;
  scalar_aux * aux;
  Monitor * monitor;
  ofstream monitorfile;
  double eta_KO;   // Kreiss-Oliger coefficient
  Potential* potential;
public:
  ScalarFieldWithPotential(Grid * grid_i, dumper * dump_i, InData * indata_i, 
	      int cowling_i, Cosmology * cosmology_i, 
	      Monitor *monitor_i, double eta_KO_i, Potential * potential_i) :
  Matter(grid_i, dump_i, indata_i, cowling_i, cosmology_i),
  monitor(monitor_i), eta_KO(eta_KO_i), potential(potential_i)
  {
    cout << " SCALARFIELDWITHPOTENTIAL: setting up scalar field... " << endl;
    // 
    // create states for dynamical variables
    // 
    cout << " SCALARFIELDWITHPOTENTIAL: setting up states... " << endl;
    last = new scalar_state(grid, dump, "scalar_last");
    derivs = new scalar_state(grid, dump, "scalar_derivs");
    inter = new scalar_state(grid, dump, "scalar_inter");
    updates = new scalar_state(grid, dump, "scalar_updates");
    // create dump list for state last:
    last->function_names();
    last->assemble_dump_list("Dump_List");
    //
    // create ADM sources
    // 
    adm_sources = new ADM_Source_Terms(grid,dump,"scalar_sources");
    adm_sources->assemble_dump_list("Dump_List");
    //
    // create auxiliary variables
    // 
    aux = new scalar_aux(grid, dump, "scalar_aux");
    aux->assemble_dump_list("Dump_List");
     //
    // finally create a monitor file...
    //
    ostringstream monfilename;
    monfilename << monitor->Filestem() << "_" << N_r - 2*N_g << "_" 
		<< N_t - 2*N_g << ".scalar_mon" << ends;
    monitorfile.open(monfilename.str().c_str());
    monitorfile.setf(ios::right);
    time_t clocktime;
    struct tm * currenttime;
    time(&clocktime);
    currenttime=localtime(&clocktime);
    monitorfile << "# File created on " << asctime(currenttime);  
    monitorfile << "# ScalarFieldWithPotential evolution " << endl;
    monitorfile << "# " << setw(18) << "time" 
		<< setw(18) << "pr time (r=0)" 
		<< setw(18) << "scal_fld_c" 
		<< setw(18) << "rho_ADM_c" 
		<< setw(18) << "rho_ADM_c_max" 
		<< endl;
    monitorfile << "#=======================================================================================================================================================================" << endl;
   //
    rho_center = rho_c_max = drhoddr = 0.0;
    sf_center_new = sf_center_old = sf_center_previous = 0.0;
  };
  ~ScalarFieldWithPotential() {
    delete last;
    delete derivs;
    delete inter;
    delete updates;
    delete adm_sources;
    delete aux;
    monitorfile.close();
    cout << " SCALARFIELDWITHPOTENTIAL: destructing derived class ScalarFieldWithPotential " << endl;
  };
  const char * Name() { return "ScalarFieldWithPotential's equations"; }
  //================================================
  // initialize
  //================================================
  void Initialize(state *s, curvature *c, diagnostics *d);
  //===============================================
  // Compute RHS sides for matter equations
  //===============================================
  void Compute_RHS(state *s, curvature *c, double time = 0.0);
  void dot_sf(scalar_state *m, state *s, double time = 0.0);
  void dot_pi(scalar_state *m, state *s, curvature *c, double time = 0.0);
  //===============================================
  // Start and finish RK steps
  //===============================================
  void Start_RK() {
    inter->equals(last);
    updates->equals(last);
  };
  void Finish_RK(state *s, curvature *c, double dt) {
    last->equals(updates);
    last->fill_ghosts();
    inter->equals(last);
    ADM_Sources(s, c);
  }
  //===============================================
  // Update matter state (adds dt * derivs to update)
  //===============================================
  void Update(double dt) {
    updates->add(dt, derivs);
  };
  //===============================================
  // Compute intermediate state (computes inter = last + dt * derivs )
  //===============================================
  void Compute_inter(double dt) {
    inter->add(last, dt, derivs);
    inter->fill_ghosts();
  };
  //===============================================
  // Regridding etc
  //===============================================
  int Regrid (VecDoub r_new) {
    last->Regrid(r_new);
    return 0;
  };
  double RegridCriterion() {
    const double tiny = 1.e-12;
    const double dfctddr = fabs(last->sf(N_g,N_g,N_g));
    const double scale = sqrt(1.0 / (dfctddr + tiny));
    const double grid_scale = grid->delta_r(N_g);
    return grid_scale / scale;
  };
  gf3d * MatterField() { return &last->sf; };
  //================================================
  // Note
  //================================================
  void Note(int time_step, double time, double tau_c) {
    if (rho_center > rho_c_max) rho_c_max = rho_center;
    if (time_step % monitor->Note_Step() == 0) {
      monitorfile.setf(ios::left);
      monitorfile << setw(18) << time
		  << setw(18) << tau_c
		  << setprecision(10) << setw(18) << rho_center
		  << setprecision(10) << setw(18) << rho_c_max
		  << endl;
    }
  };
  //================================================
  // dump grid functions
  //================================================
  void dump_fcts(double time, double prop_time, int timestep,
		 const char * suffix = "") 
  {
    if (dump->time_to_dump(timestep) || strcmp(suffix,"") ) {
      cout << " SCALARFIELDWITHPOTENTIAL: Dumping matter functions at time t = " 
	   << time  << endl;
      last->dump_fcts(time, prop_time, timestep, suffix);
      adm_sources->dump_fcts(time, prop_time, timestep, suffix);
      aux->dump_fcts(time, prop_time, timestep, suffix);
    }
  }
  //================================================
  // Diagnostics
  //================================================
  double Compute_Diagnostics(state * s, curvature * c) {
    for (int i = 0; i < N_r; i++) {
      const double rl = grid->r(i);
      for (int j = 0; j < N_t; j++) {
	// const double sinthetal = grid->sintheta(j);
	for (int k = 0; k < N_p; k++) {
	  double psi_l = exp(s->phi(i,j,k));
	  double h_tt_l = s->h_tt(i,j,k);
	  double r_areal = rl*psi_l*psi_l*sqrt(1.0 + h_tt_l);
	  aux->Omega[i][j][k] = 4.0*PI*r_areal*r_areal*
	    adm_sources->rho_ADM(i,j,k);
	}
      }
    }

    rho_components(s, c);
    Hamiltonian_components(s, c);

    return aux->Omega.max();
  };
  //================================================
  // ADM sources
  //================================================
  void ADM_Sources(state * s, curvature *c);

  void rho_components(state * s, curvature * c);

  void Hamiltonian_components(state *s, curvature *curve);

};
