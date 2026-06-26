// Tell emacs that this is -*-c++-*- mode

#include "RadHydro_State.h"
#include "RadHydro_Aux.h"
#include "Monitor.h"
#include "EOS.h"
#include "Slope_Limiter.h"
#include "Approximate_Riemann_Solver.h"
#include "Particles.h"
#include "CheckPoint.h"

//================================================
//
// RadHydro
//
//================================================
class RadHydro : public Matter {
public:
  radhydro_state * last, * derivs, * inter, * updates;
  radhydro_aux * aux;
  Monitor * monitor;
  ofstream monitorfile;
  //
  // EOS etc
  //
  EOS *eos;
  Slope *slope;
  Riemann_Solver *riemann_solver;
  //
  // radhydro parameters
  //
  int eos_type, slope_limiter_type, riemann_type;
  bool partial;
  double f_atm, f_thr;
  double rho_atm, rho_thr;
  //  double kappa, Gamma;
  //
  // some diagnostics
  //
  double rho_center, rho_c_max, drhoddr; 
  int recovery_error_counter;
  Particles * particles;
  gf3d ** fct_list;
  CheckPoint * checkpoint;
  //
  //
public:
  RadHydro(Grid * grid_i, dumper * dump_i, InData * indata_i, EOS *eos_i,
	  int cowling_i, Cosmology * cosmology_i, 
	   Monitor *monitor_i, CheckPoint * checkpoint_i) :
    Matter(grid_i, dump_i, indata_i, cowling_i, cosmology_i),
    monitor(monitor_i), eos(eos_i), recovery_error_counter(0),
    checkpoint(checkpoint_i)
  {
    cout << " RADHYDRO: setting up radhydro... " << endl;
    //
    // Read Input
    //
    int N_particles = 0;
    double R_max = 0.0;
    int read_error = ReadInput(N_particles, R_max);
    //
    // if read didn't work use default values...
    //
    if (read_error != 0) {
      eos_type = 2;
      slope_limiter_type = 1;
      riemann_type = 1;
      partial = false;
      f_atm = 1.e-12;
      f_thr = 1.e2;
    }
    //
    // set up slope limiter
    // 
    if (slope_limiter_type == minmod_lim) {
      slope = new minmod();
    } else if (slope_limiter_type == MC_lim) {
      slope = new MC();
    } else {
      cout << " RADHYDRO: Unknown slope limiter! " << endl;
    }
    //
    // set up Riemann Solver
    // 
    if (riemann_type == llf_solver) {
      riemann_solver = new LLF();
    } else if (riemann_type == hlle_solver) {
      riemann_solver = new HLLE();
    } else {
      cout << " RADHYDRO: Unknown Riemann Solver! " << endl;
    }
    //
    // echo out input stuff...
    //
    cout << " RADHYDRO: using " << eos->Name() << endl;
    cout << " RADHYDRO: using " << slope->Name() << " slope limiter and " 
	 << riemann_solver->Name() << " Riemann solver " << endl;
    if (partial) 
      cout << " RADHYDRO: using 'partial' approach " << endl;
    else 
      cout << " RADHYDRO: using 'full' approach " << endl;
    cout << " RADHYDRO: using f_atm = " << f_atm << " and f_thr = " << f_thr << endl;
    // 
    // create states for dynamical variables
    // 
    cout << " RADHYDRO: setting up states... " << endl;
    last = new radhydro_state(grid, dump, "radhydro_last");
    derivs = new radhydro_state(grid, dump, "radhydro_derivs");
    inter = new radhydro_state(grid, dump, "radhydro_inter");
    updates = new radhydro_state(grid, dump, "radhydro_updates");
    // create dump list for state last:
    last->assemble_dump_list("Dump_List");
    // let checkpointer know about last
    checkpoint->CollectDynVariables(last);
    //
    // create ADM sources
    // 
    adm_sources = new ADM_Source_Terms(grid,dump,"radhydro_sources");
    adm_sources->assemble_dump_list("Dump_List");
    //
    // create auxiliary variables
    // 
    aux = new radhydro_aux(grid, dump, "radhydro_aux");
    aux->assemble_dump_list("Dump_List");
    //
    // create fluid particle tracers
    //
    if (N_particles > 0) {
      cout << " RADHYDRO: setting up fluid traces..." << endl;
      particles = new Particles(N_particles, R_max, grid);
      particles->add_to_dump_list("Particle_Fct_List", last);
      particles->add_to_dump_list("Particle_Fct_List", aux);
      particles->add_to_dump_list("Particle_Fct_List", adm_sources);
    } else
      particles = NULL;
    //
    // finally create a monitor file...
    //
    ostringstream monfilename;
    monfilename << monitor->Filestem() << "_" << N_r - 2*N_g << "_" 
     		<< N_t - 2*N_g << ".radhydro_mon" << ends;
    monitorfile.open(monfilename.str().c_str());
    monitorfile.setf(ios::right);
    time_t clocktime;
    struct tm * currenttime;
    time(&clocktime);
    currenttime=localtime(&clocktime);
    monitorfile << "# File created on " << asctime(currenttime);  
    monitorfile << "# Radiation hydro evolution " << endl;
    monitorfile << "# " << setw(16) << "time" 
		<< setw(18) << "pr time (r=0)" 
		<< setw(18) << "rho_ADM_c" 
		<< setw(18) << "rho_ADM_c_max" 
		<< setw(18) << "P_c" 
		<< setw(18) << "epsilon_c" 
		<< setw(18) << "D_c" 
		<< setw(18) << "tau_c" 
		<< setw(18) << "M_0"       
		<< endl;
    monitorfile << "#=======================================================================================================================================================================" << endl;
   //
    rho_center = rho_c_max = drhoddr = 0.0;
  };
  ~RadHydro() {
    delete last;
    delete derivs;
    delete inter;
    delete updates;
    delete adm_sources;
    delete aux;
    delete particles;
    monitorfile << "# found " << recovery_error_counter 
		<< " errors in recovery - bye! " << endl; 
    monitorfile.close();
    cout << " RADHYDRO: destructing derived class RadHydro " << endl;
  };
  const char * Name() { return "relativistic radiation hydrodynamics"; }
  //================================================
  // Initialize
  //================================================
  void Initialize(state *s, curvature *c, diagnostics *d);
  //===============================================
  // Start and finish RK steps
  //===============================================
  void Start_RK() {
    inter->equals(last);
    updates->equals(last);
    if (particles != NULL) particles->Start_RK();
  };
  void Finish_RK(state *s, curvature *c, double dt) {
    last->equals(updates);
    last->fill_ghosts();
    inter->equals(last);
    ADM_Sources(s, c);
    // CHECK : compute primitives?
    if (particles != NULL) particles->Finish_RK();
 }
  //===============================================
  // Update matter state (adds dt * derivs to update)
  //===============================================
  void Update(double dt) {
    updates->add(dt, derivs);
    if (particles != NULL) particles->Update(dt);
  };
  //===============================================
  // Compute intermediate state (computes inter = last + dt * derivs )
  //===============================================
  void Compute_inter(double dt) {
    inter->add(last, dt, derivs);
    inter->fill_ghosts();
    if (particles != NULL) particles->Compute_inter(dt);
  };
  //===============================================
  // Read Input
  //===============================================
  int ReadInput(int & N_particles, double & R_max);
  //================================================
  // Compute conserved variables from primitive variables 
  // (both global and local version)
  //================================================
  int Compute_Conserved_Vars(radhydro_state *m, state *s, curvature *c);
  int Compute_Conserved_Vars(double rho_0, double eps, double p,
			     double v_r, double v_t, double v_p,
			     double E, double F_0, 
			     double F_r, double F_t, double F_p,
			     double h_rr, double h_rt, double h_rp,
			     double h_tt, double h_tp, double h_pp,
			     double lapse, double phi, double det,
			     double r, double sintheta,
			     double & D, 
			     double & S_r, double & S_t, double & S_p, 
			     double & tau, double & W, double & v2,
			     double & tau_rad, double & S_rad_r,
			     double & S_rad_t, double & S_rad_p,
			     bool verbose = false);
  //===============================================
  // Recovery: compute primitive from conserved variables
  //===============================================
  int Recovery(radhydro_state *m, state *s, curvature *c);
  double f(double p, double tau, double D, double S2, double & dfdp);
  //===============================================
  // Reconstruction
  //===============================================
  int Reconstruct_r(gf3d & fct, gf3d & fct_L, gf3d & fct_R);
  int Reconstruct_t(gf3d & fct, gf3d & fct_L, gf3d & fct_R);
  int Reconstruct_p(gf3d & fct, gf3d & fct_L, gf3d & fct_R);
  int Construct_Fluxes(state *s, curvature *c);
  int Compute_Fluxes(int index, double r, double sintheta,
		     double lapse, double phi, double det, double D,
		     double v, double shift, 
		     double S_r, double S_t, double S_p, 
		     double p, double h, double W, double tau,
		     double E, double F_0, double F, double tau_rad,
		     double S_rad_r, double S_rad_t, double S_rad_p, 
		     double & f_D, 
		     double & f_S_r, double & f_S_t, double & f_S_p,
		     double & f_tau, double & f_tau_rad,
		     double & f_S_rad_r, double & f_S_rad_t, 
		     double & f_S_rad_p, bool verbose = false);
  //===============================================
  // Compute right-hand sides
  //===============================================
  int dot_conserved_vars(radhydro_state *m, state *s, curvature *c);
  //===============================================
  // Compute eigenvalues
  //===============================================
  int lambda(double lapse, double v, double v2, 
	     double shift, double gamma, double cs,
	     double & lambda_0, double & lambda_p, double & lambda_m);
  //===============================================
  // Compute RHS sides for matter equations
  //===============================================
  void Compute_RHS(state *s, curvature *c, double time = 0.0);
  //===============================================
  // Regridding etc
  //===============================================
  int Regrid (VecDoub r_new) { return 0;};
  double RegridCriterion() { return 0.0; };
  gf3d * MatterField() { return &last->D; };
  //================================================
  // Note
  //================================================
  void Note(int time_step, double time, double tau_c) {
    rho_center = adm_sources->rho_ADM(0.0,N_g,N_g);
    if (rho_center > rho_c_max) rho_c_max = rho_center;
    if (time_step % monitor->Note_Step() == 0) {
      const double rest_mass = Rest_Mass(last);
      monitorfile.setf(ios::left);
      monitorfile << setw(18) << time
		  << setw(18) << tau_c
		  << setprecision(10) << setw(18) << rho_center
		  << setprecision(10) << setw(18) << rho_c_max
		  << setprecision(10) << setw(18) << aux->p(0.0,N_g,N_g)
		  << setprecision(10) << setw(18) << aux->eps(0.0,N_g,N_g)
		  << setprecision(10) << setw(18) << last->D(0.0,N_g,N_g)
		  << setprecision(10) << setw(18) << last->tau(0.0,N_g,N_g)
		  << setprecision(10) << setw(18) << rest_mass
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
      cout << " RADHYDRO: Dumping matter functions at time t = " 
	   << time  << endl;
      last->dump_fcts(time, prop_time, timestep, suffix);
      adm_sources->dump_fcts(time, prop_time, timestep, suffix);
      aux->dump_fcts(time, prop_time, timestep, suffix);
      if (particles != NULL)
	particles->Dump(time, prop_time, timestep, suffix);
    }
  }
  //================================================
  // Diagnostics: compute magnitude of F^\alpha, and u_a F^a
  //================================================
  double Compute_Diagnostics(state * s, curvature * c) {
    for (int i = N_g; i < N_r - N_g; i++) {
      const double rl = grid->r(i);
      for (int j = N_g; j < N_t - N_g; j++) {
	const double sinthetal = grid->sintheta(j);
	const double rst = rl * sinthetal;
	for (int k = N_g; k < N_p - N_g; k++) {
	  vect f(aux->f_r(i,j,k),aux->f_t(i,j,k),aux->f_p(i,j,k));
	  vect v(aux->v_r(i,j,k),aux->v_t(i,j,k),aux->v_p(i,j,k));
	  vect beta(s->shift_r(i,j,k),s->shift_t(i,j,k),s->shift_p(i,j,k));
	  tensor g_conf(1.0 + s->h_rr(i,j,k),
			s->h_rt(i,j,k),
			s->h_rp(i,j,k),
			1.0 + s->h_tt(i,j,k),
			s->h_tp(i,j,k),
			1.0 + s->h_pp(i,j,k));
	  const double lapsel = s->lapse(i,j,k);
	  const double psi = exp(s->phi(i,j,k));
	  const double psi4 = psi*psi*psi*psi;
	  const double fl = aux->f(i,j,k);
	  const double F0 = fl / lapsel;
	  double F_dot_F = - fl * fl;
	  double u_dot_F = - fl;
	  for (int a = 0; a < 3; a++) {
	    for (int b = 0; b < 3; b++) {
	      F_dot_F += psi4 * g_conf[a][b] * f[a] * f[b];
	      u_dot_F += psi4 * g_conf[a][b] * f[a] * v[b];
	    }
	  }
	  aux->FF[i][j][k] = F_dot_F;
	  aux->udotF[i][j][k] = u_dot_F * aux->W(i,j,k);

	}
      }
    }
    for (int i = 0; i < N_r; i++) {
      const double rl = grid->r(i);
      for (int j = 0; j < N_t; j++) {
	//     	const double sinthetal = grid->sintheta(j);
     	for (int k = 0; k < N_p; k++) {
	  // const double e2p = exp(2.0 * s->phi(i,N_g,N_g));
	  // const double R = rl * e2p * sqrt(1.0 + s->h_pp(i,N_g,N_g));
	  const double u = - aux->W(i,N_g,N_g) * 
	    (aux->v_r(i,N_g,N_g) - s->shift_r(i,N_g,N_g) / s->lapse(i,N_g,N_g));
	  const double e6p = exp(6.0 * s->phi(i,N_g,N_g));
	  const double det_l = e6p * sqrt( c->det(i,N_g,N_g) ) * rl * rl;
	  const double M_0_dot = 4.0 * PI * det_l * aux->rho_0(i,N_g,N_g) * u * s->lapse(i,N_g,N_g);	  
	  //
	  // CAREFUL: accretion rate correct only in spherical symmetry!
	  //
	  aux->flux[i][j][k] = M_0_dot;
	  aux->sound_speed[i][j][k] = eos->sound_speed(aux->rho_0(i,j,k),
						       aux->eps(i,j,k));
	  aux->temperature[i][j][k] = 
	    eos->temperature(aux->rho_0(i,j,k),aux->eps(i,j,k));
     	}
      }
    }
    return aux->rho_0(0.0,N_g,N_g);
  };
  double Rest_Mass(radhydro_state *m);
  //================================================
  // Test
  //================================================
  int Test_Recovery(state *s, curvature *c);
  //================================================
  // ADM sources
  //================================================
  void ADM_Sources(state * s, curvature *c);
};




