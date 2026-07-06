// Tell emacs that this is -*-c++-*- mode

#include "Hydro_State.h"
#include "Hydro_Aux.h"
#include "CheckPoint.h"
#include "Monitor.h"
#include "EOS.h"
#include "Slope_Limiter.h"
#include "Approximate_Riemann_Solver.h"
#include "Particles.h"
#include "HorizonFinder.h"

//================================================
//
// Hydro
//
//================================================
class Hydro : public Matter {
public:
  hydro_state * last, * derivs, * inter, * updates;
  hydro_aux * aux;
  Monitor * monitor;
  ofstream monitorfile;
  CheckPoint * checkpoint;
  //
  // EOS etc
  //
  EOS *eos;
  Slope *slope;
  Riemann_Solver *riemann_solver;
  //
  // hydro parameters
  //
  int slope_limiter_type, riemann_type;
  bool partial;
  double f_atm, f_thr;
  double rho_atm, rho_thr;
  //
  // flags for *not* evolving hydro variables inside black holes
  // 
  bool fix_hydro;     // if true, will not evolve first n_fixed radial points
  int n_fixed;
  //
  // some diagnostics
  //
  double rho_center, rho_c_max, drhoddr; 
  int recovery_error_counter, cs_error;
  Particles * particles;
  gf3d ** fct_list;  // list of functions to be tracked by particles
  HorizonFinder * horizonfinder;
  //
  //
public:
  Hydro(Grid * grid_i, dumper * dump_i, InData * indata_i, EOS *eos_i,
	int cowling_i, Cosmology * cosmology_i, 
	Monitor *monitor_i, HorizonFinder *horizonfinder_i,
	CheckPoint *checkpoint_i) :
    Matter(grid_i, dump_i, indata_i, cowling_i, cosmology_i),
    monitor(monitor_i), eos(eos_i), recovery_error_counter(0),
    horizonfinder(horizonfinder_i), checkpoint(checkpoint_i), cs_error(0)
  {
    cout << " HYDRO: setting up hydro... " << endl;
    //
    // make sure EOS has been set!
    // 
    if (eos == NULL) {
      cerr << " HYDRO: No EOS has been provided -- problem!! " << endl;
    }
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
      slope_limiter_type = 1;
      riemann_type = 1;
      partial = false;
      f_atm = 1.e-12;
      f_thr = 1.e2;
      // kappa = indata->Polytropic_K();
      // Gamma = indata->Polytropic_Gamma();
    }
    //
    // set up slope limiter
    // 
    if (slope_limiter_type == minmod_lim) {
      slope = new minmod();
    } else if (slope_limiter_type == MC_lim) {
      slope = new MC();
    } else {
      cout << " HYDRO: Unknown slope limiter! " << endl;
    }
    //
    // set up Riemann Solver
    // 
    if (riemann_type == llf_solver) {
      riemann_solver = new LLF();
    } else if (riemann_type == hlle_solver) {
      riemann_solver = new HLLE();
    } else {
      cout << " HYDRO: Unknown Riemann Solver! " << endl;
    }
    //
    // echo out input stuff...
    //
    cout << " HYDRO: using " << eos->Name() << endl;
    cout << " HYDRO: using " << slope->Name() << " slope limiter and " 
	 << riemann_solver->Name() << " Riemann solver " << endl;
    if (partial) 
      cout << " HYDRO: using 'partial' approach " << endl;
    else 
      cout << " HYDRO: using 'full' approach " << endl;
    cout << " HYDRO: using f_atm = " << f_atm << " and f_thr = " << f_thr << endl;
    if (fix_hydro)
      cout << " HYDRO: fixing first " << n_fixed << " radial grid points" << endl;
    else
      cout << " HYDRO: evolving hydro for all grid points" << endl;
    // 
    // create states for dynamical variables
    // 
    cout << " HYDRO: setting up states... " << endl;
    last = new hydro_state(grid, dump, "hydro_last");
    derivs = new hydro_state(grid, dump, "hydro_derivs");
    inter = new hydro_state(grid, dump, "hydro_inter");
    updates = new hydro_state(grid, dump, "hydro_updates");
    // let checkpointer know about aux:
    checkpoint->CollectDynVariables(last);
    // create dump list for state last:
    last->assemble_dump_list("Dump_List");
    //
    // create ADM sources
    // 
    adm_sources = new ADM_Source_Terms(grid,dump,"hydro_souces");
    adm_sources->assemble_dump_list("Dump_List");
    //
    // create fluxes (needed to compute accretion rate onto black holes)
    // 
    fluxes = new Fluxes(grid,dump,"hydro_fluxes");
    fluxes->assemble_dump_list("Dump_List");
    //
    // create auxiliary variables
    // 
    aux = new hydro_aux(grid, dump, "hydro_aux");
    aux->assemble_dump_list("Dump_List");
    // let checkpointer know about aux:
    checkpoint->CollectDynVariables(aux);
    //
    // create fluid particle tracers
    //
    if (N_particles > 0) {
      cout << " HYDRO: setting up fluid traces..." << endl;
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
    monfilename << "output/" << monitor->Filestem() << "_" << N_r - 2*N_g << "_" 
     		<< N_t - 2*N_g << ".hydro_mon" << ends;
    monitorfile.open(monfilename.str().c_str());
    monitorfile.setf(ios::right);
    time_t clocktime;
    struct tm * currenttime;
    time(&clocktime);
    currenttime=localtime(&clocktime);
    monitorfile << "# File created on " << asctime(currenttime);  
    monitorfile << "# Hydro evolution with EOS " << eos->Name() << endl;
    if (fix_hydro)
      monitorfile << "# fixing first " << n_fixed << " radial grid points" << endl;
    else
      monitorfile << "# evolving hydro for all grid points" << endl;
    monitorfile << "# " << setw(18) << "time" 
		<< setw(18) << "pr time (r=0)" 
		<< setw(18) << "rho_ADM_max" 
		<< setw(18) << "rho_0_max" 
		<< setw(18) << "P_max" 
		<< setw(18) << "epsilon_max" 
		<< setw(18) << "D_max" 
		<< setw(18) << "tau_max" 
		<< setw(18) << "M_0 (outside BH)"       
		<< setw(18) << "Unbound M_0"       
		<< endl;
    monitorfile << "#=========================================================================================================================================================================================" << endl;
   //
    rho_center = rho_c_max = drhoddr = 0.0;
  };
  ~Hydro() {
    delete last;
    delete derivs;
    delete inter;
    delete updates;
    delete adm_sources;
    delete aux;
    delete particles;
    monitorfile << "# found " << recovery_error_counter 
		<< " errors in recovery - bye! " << endl; 
    if (cs_error > 0) {
      cout << "HYDRO: Found cs >= 1 " << cs_error << " times... " << endl;
    }
    monitorfile << "# found cs >= 1 " << cs_error << " times! " << endl;
    monitorfile.close();
    cout << " HYDRO: destructing derived class Hydro " << endl;
  };
  const char * Name() { return "relativistic hydrodynamics"; }
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
    aux->hydro_errors.equals(0.0);  // set errors to zero...
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
    DownstreamBoundaries(updates);
    if (particles != NULL) particles->Update(dt);
  };
  //===============================================
  // Compute intermediate state (computes inter = last + dt * derivs )
  //===============================================
  void Compute_inter(double dt) {
    inter->add(last, dt, derivs);
    DownstreamBoundaries(inter);
    inter->fill_ghosts();
    if (particles != NULL) particles->Compute_inter(dt);
  };
  //===============================================
  // Downstream boundaries for conserved quantities
  //===============================================
  void DownstreamBoundaries(hydro_state *m) {
    for (int j = N_g; j < N_t - N_g; j++) {
      double thetal = aux->rho_0.theta(j);
      if (thetal > PI/2.0) {
	for (int k = N_g; k < N_p - N_g; k++) {
	  for (int i = N_r - N_g; i < N_r; i++) {
	    double rl = aux->rho_0.r(i);
	    double sinthetal = aux->rho_0.sintheta(j);
	    double sinthetap = aux->rho_0.sintheta(j-1);
	    double costhetal = aux->rho_0.costheta(j);
	    double costhetap = aux->rho_0.costheta(j-1);
	    double r_int = rl * sinthetal / sinthetap;
	    m->D[i][j][k]   = m->D(r_int,j-1,k);
	    m->tau[i][j][k] = m->tau(r_int,j-1,k);
	    double Sr_p = m ->S_r(r_int,j-1,k);
	    double St_p = m ->S_t(r_int,j-1,k);
	    double Sz = costhetap * Sr_p + sinthetap * St_p;
	    double Ss = sinthetap * Sr_p + costhetap * St_p;
	    m->S_r[i][j][k] = - costhetal * Sz + sinthetal * Ss;
	    m->S_t[i][j][k] =   sinthetal * Sz + costhetal * Ss;
	  }
	}
      }
    }
  }
  //===============================================
  // Read Input
  //===============================================
  int ReadInput(int & N_particles, double & R_max);
  //================================================
  // Compute conserved variables from primitive variables 
  // (both global and local version)
  //================================================
  int Compute_Conserved_Vars(hydro_state *m, state *s, curvature *c);
  int Compute_Conserved_Vars(double rho_0, double eps, double p,
			     double v_r, double v_t, double v_p,
			     double h_rr, double h_rt, double h_rp,
			     double h_tt, double h_tp, double h_pp,
			     double phi, double det,
			     double r, double sintheta,
			     double & D, 
			     double & S_r, double & S_t, double & S_p, 
			     double & tau, double & W, double & v2);
  //===============================================
  // Recovery: compute primitive from conserved variables
  //===============================================
  int Recovery(hydro_state *m, state *s, curvature *c);
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
		     double p, double tau,
		     double & f_D, 
		     double & f_S_r, double & f_S_t, double & f_S_p,
		     double & f_tau);
  //===============================================
  // Compute right-hand sides
  //===============================================
  int dot_conserved_vars(hydro_state *m, state *s, curvature *c);
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
      const double unbound_mass = Unbound_Mass(last);
      monitorfile.setf(ios::left);
      monitorfile << setw(18) << time
		  << setw(18) << tau_c
		  << setprecision(10) << setw(18) << adm_sources->rho_ADM.max()
		  << setprecision(10) << setw(18) << aux->rho_0.max()
		  << setprecision(10) << setw(18) << aux->p.max()
		  << setprecision(10) << setw(18) << aux->eps.max()
		  << setprecision(10) << setw(18) << last->D.max()
		  << setprecision(10) << setw(18) << last->tau.max()
		  << setprecision(10) << setw(18) << rest_mass
		  << setprecision(10) << setw(18) << unbound_mass
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
      cout << " HYDRO: Dumping matter functions at time t = " 
	   << time  << endl;
      last->dump_fcts(time, prop_time, timestep, suffix);
      adm_sources->dump_fcts(time, prop_time, timestep, suffix);
      aux->dump_fcts(time, prop_time, timestep, suffix);
      if (particles != NULL)
	particles->Dump(time, prop_time, timestep, suffix);
      fluxes->dump_fcts(time, prop_time, timestep, suffix);
    }
  }
  //================================================
  // Diagnostics
  //================================================
  double Compute_Diagnostics(state * s, curvature * c) {
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
	  //
	  // compute u_t (lower index!)
	  //
	  const double det_3D = exp(6.0*s->phi[i][j][k])*sqrt(c->det[i][j][k]);
	  const double Wl = aux->W[i][j][k];
	  const double rho_0l = aux->rho_0[i][j][k];
	  const double epsl = aux->eps[i][j][k];
	  const double hl = eos->h(rho_0l,epsl);
	  // compute u_i (lower index...)
	  const double ur = last->S_r[i][j][k] / (det_3D * Wl * rho_0l * hl);
	  const double ut = last->S_t[i][j][k] / (det_3D * Wl * rho_0l * hl);
	  const double up = last->S_p[i][j][k] / (det_3D * Wl * rho_0l * hl);
	  const double beta_dot_u = s->shift_r(i,j,k) * ur +
	    s->shift_t(i,j,k) * ut + s->shift_p(i,j,k) * up;
	  // get the following from
	  //     u^t = g^{ta} u_a = (- u_t + beta^i u_i ) / alpha^2
	  aux->ut_low[i][j][k] = beta_dot_u - s->lapse(i,j,k)*Wl; 
	}
      }
    }
    return aux->rho_0(0.0,N_g,N_g);
  };
  double Rest_Mass(hydro_state *m);
  double Unbound_Mass(hydro_state *m);
  //================================================
  // Test
  //================================================
  int Test_Recovery(state *s, curvature *c);
  //================================================
  // ADM sources
  //================================================
  void ADM_Sources(state * s, curvature *c);
  //================================================
  // Read grid function from checkpoint file
  //================================================
  bool read_gf(gf3d & gf, int timestep);
};




