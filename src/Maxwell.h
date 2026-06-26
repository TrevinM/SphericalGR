// Tell emacs that this is -*-c++-*- mode

#include "Maxwell_State.h"
#include "Maxwell_Aux.h"
#include "Monitor.h"
#include "CheckPoint.h"

//================================================
//
// Maxwell
//
//================================================
class Maxwell : public Matter {
public:
  double rho_center, rho_c_max, drhoddr;    // diagnostics 
  maxwell_state * last, * derivs, * inter, * updates;
  maxwell_aux * aux;
  Monitor * monitor;
  CheckPoint * checkpoint;
  ofstream monitorfile;
  double eta_KO;   // Kreiss-Oliger coefficient
  int char_OB;    // decides how Sommerfeld BCs are implemented
                  // 0: derivatives; 1: characteristic interpolation
public:
  Maxwell(Grid * grid_i, dumper * dump_i, InData * indata_i, 
	  int cowling_i, int char_OB_i, Cosmology * cosmology_i, 
	  Monitor *monitor_i, double eta_KO_i, CheckPoint * checkpoint_i) :
    Matter(grid_i, dump_i, indata_i, cowling_i, cosmology_i),
    char_OB(char_OB_i), monitor(monitor_i), eta_KO(eta_KO_i),
    checkpoint(checkpoint_i)
  {
    cout << " MAXWELL: setting up maxwell... " << endl;
    // 
    // create states for dynamical variables
    // 
    cout << " MAXWELL: setting up states... " << endl;
    last = new maxwell_state(grid, dump, "maxwell_last");
    derivs = new maxwell_state(grid, dump, "maxwell_derivs");
    inter = new maxwell_state(grid, dump, "maxwell_inter");
    updates = new maxwell_state(grid, dump, "maxwell_updates");
    // create dump list for state last:
    last->assemble_dump_list("Dump_List");
    // let checkpointer know about last
    checkpoint->CollectDynVariables(last);
    //
    // create ADM sources
    // 
    adm_sources = new ADM_Source_Terms(grid,dump,"maxwell_souces");
    adm_sources->assemble_dump_list("Dump_List");
    //
    // create auxiliary variables
    // 
    aux = new maxwell_aux(grid, dump, "maxwell_aux");
    aux->assemble_dump_list("Dump_List");
     //
    // finally create a monitor file...
    //
    ostringstream monfilename;
    monfilename << monitor->Filestem() << "_" << N_r - 2*N_g << "_" 
     		<< N_t - 2*N_g << ".maxwell_mon" << ends;
    monitorfile.open(monfilename.str().c_str());
    monitorfile.setf(ios::right);
    time_t clocktime;
    struct tm * currenttime;
    time(&clocktime);
    currenttime=localtime(&clocktime);
    monitorfile << "# File created on " << asctime(currenttime);  
    monitorfile << "# Maxwell evolution " << endl;
    monitorfile << "# " << setw(16) << "time" 
		<< setw(18) << "pr time (r=0)" 
		<< setw(18) << "rho_ADM_c" 
		<< setw(18) << "rho_ADM_c_max" 
		<< endl;
    monitorfile << "#=======================================================================================================================================================================" << endl;
   //
    rho_center = rho_c_max = drhoddr = 0.0;
  };
  ~Maxwell() {
    delete last;
    delete derivs;
    delete inter;
    delete updates;
    delete adm_sources;
    delete aux;
    monitorfile.close();
    cout << " MAXWELL: destructing derived class Maxwell " << endl;
  };
  const char * Name() { return "Maxwell's equations"; }
  //================================================
  // initialize
  //================================================
  void Initialize(state *s, curvature *c, diagnostics *d);
  //===============================================
  // Compute RHS sides for matter equations
  //===============================================
  void Compute_RHS(state *s, curvature *c, double time = 0.0);
  void dot_a_p(maxwell_state *m, state *s, double time = 0.0);
  void dot_e_p(maxwell_state *m, state *s, curvature *c, double time = 0.0);
  //===============================================
  // Start and finish RK steps
  //===============================================
  void Start_RK() {
    inter->equals(last);
    updates->equals(last);
  };
  void Finish_RK(state *s, curvature *c, double dt) {
    updates->char_OB(last, dt);
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
    inter->char_OB(last, dt);
  };
  //===============================================
  // Regridding etc
  //===============================================
  int Regrid (VecDoub r_new) { last->Regrid(r_new); return 0;};
  double RegridCriterion() { 
    const Doub tiny = 1.e-12;
    // const Doub drhoddr = abs(a_p_o.ddr(N_g,N_g,N_g));
    Doub rho_l = fabs(rho_center);
    if (rho_l < 1.0) rho_l = 1.0;
    const Doub scale = sqrt(rho_l/(fabs(drhoddr) + tiny));
    const Doub grid_scale = grid->delta_r(N_g);
    // if ( grid_scale / scale > 0.15 ) {
    //   cout << " rho_center = " << rho_center 
    // 	   << " drhoddr = " << drhoddr 
    // 	   << " scale = " << scale 
    // 	   << " grid_scale = " << grid->delta_r(N_g) 
    // 	   << endl;
    // }
    return grid_scale / scale;
  };
  gf3d * MatterField() { return &last->a_p; };
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
		 const char * suffix) 
  {
    if (dump->time_to_dump(timestep) || strcmp(suffix,"") ) {
      cout << " MAXWELL: Dumping matter functions at time t = " 
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
	const double sinthetal = grid->sintheta(j);
	for (int k = 0; k < N_p; k++) {
	  // recall that gup is *not* rescaled...
	  aux->A2[i][j][k] = exp(-4.0*s->phi(i,j,k))*c->gup_pp(i,j,k)*
	    last->a_p(i,j,k)*last->a_p(i,j,k)*rl*rl*sinthetal*sinthetal;
	  const double g_pp = exp(4.0 * s->phi(i,j,k)) 
	    * (1.0 + s->h_pp(i,j,k) ); 
	  aux->A_xi[i][j][k] = last->a_p(i,j,k) / sqrt(g_pp);
	}
      }
    }
    return aux->A2(0.0,N_g,N_g);
  };
  //================================================
  // ADM sources
  //================================================
  void ADM_Sources(state * s, curvature *c);
};




