// Tell emacs that this if is -*-c++-*- mode
//
//================================================
//================================================
//================================================
// Contains three layers of classes handling photon tracers:
//
// Photon: handles individual photons
//
// PhotonSwarm: handles all photons emitted at one instant
//
// Photons: handles all photons
//================================================
//================================================
//================================================
//
#ifndef PHOTONS_H
#define PHOTONS_H

#include "nr3.h"
#include "gridfunction.h"

#include "State.h"
#include "Curvature.h"
#include "Auxiliary.h"

//================================================
//================================================
//================================================
// Class for one individual photon -- 
//     a swarm of photons emitted at one time handled by 
//     PhotonSwarm below
//
// Use p^a = q^a + q n^a 
// with q = - n_a p^a = alpha p^t
//================================================
//================================================
//================================================
class Photon {
private:
  int number;        // photon number
  double theta_0;      // angle that photon was launched into
  double T_launch;   // self-similar time that photon was launched at
  double t_launch;   // coordinate time that photon was launched at
  double tau_launch; // proper time that photon was launched at
  double tau_star;
  //====================================================
  double r_ph, theta_ph;      // current coordinate radius and angle
  double q;                   // index up: d t / d lambda
  double q_r, q_t;        // index down
  double lambda;     // affine parameter
  double lambda_old; // affine parameter at last step
  double r_ph_old, theta_ph_old;
  //====================================================
  double fct1_val, fct2_val;
  double fct1_val_old, fct2_val_old; 
  //  double fct_max;    // maximum of function monitored during flight
  //  ofstream outfile;
  int N_g;
  bool photon_left_grid;
public:
  //================================================
  //================================================
  // Constructor -- doesn't do anything...
  //================================================
  //================================================
  Photon() : number(0), theta_0(0.0), T_launch(0.0), t_launch(0.0), 
	     tau_launch(0.0), r_ph(-1.0), theta_ph(0.0), q(0.0), q_r(0.0), q_t(0.0),
	     lambda(0.0), fct1_val(0.0), fct2_val(0.0), N_g(3), photon_left_grid(false)
  { 
    //    cout << " In Default constructor " << endl;
  }
  //================================================
  //================================================
  // Destructor
  //================================================
  //================================================
  ~Photon() {
  }
  //================================================
  //================================================
  // Launch photon
  //================================================
  //================================================
  int Launch(int n, double theta_0_i, 
	     double t, double tau, double tau_star_i,
	     gf3d * fct1, gf3d * fct2, state * s, curvature *curv) {
    tau_star = tau_star_i;
    number = n;
    theta_0 = theta_0_i;
    cout << " PHOTON: Launching photon # " << number
	 << " in direction theta = " << theta_0 << endl;
    T_launch = T(tau);
    t_launch = t;
    tau_launch = tau;
    //================================================
    // initialize photon parameters
    //================================================
    r_ph = 0.0;
    theta_ph = theta_0;
    lambda = T_launch;
    double lapse_l = s->lapse(r_ph,theta_ph,N_g);
    double phi_l = s->phi(r_ph,theta_ph,N_g);
    double gup_rr_l = curv->gup_rr(r_ph,theta_ph,N_g);
    q = (tau_star - tau);
    q_r = exp(2.0*phi_l)* q / sqrt(gup_rr_l);
    //================================================
    if (fct1 == NULL)
      fct1_val = 0.0;
    else
      fct1_val   = (*fct1)(r_ph,theta_ph,N_g);
    if (fct2 == NULL)
      fct2_val = 0.0;
    else
      fct2_val   = (*fct2)(r_ph,theta_ph,N_g);
    // fct_max = fct_val;
    return number;
  }
  //================================================
  //================================================
  // Advance photon
  //================================================
  //================================================
  bool Advance(state * s, curvature * c, auxiliary * aux, 
	       gf3d * fct1, gf3d * fct2, double dt, const double r_out) {
    if (!photon_left_grid) {
      //================================================
      // Remember old values...
      //================================================
      r_ph_old = r_ph;
      theta_ph_old = theta_ph;
      lambda_old = lambda;
      fct1_val_old = fct1_val;
      fct2_val_old = fct2_val;
      //================================================
      // assemble local functions
      //================================================
      const double lapse_l = s->lapse(r_ph,theta_ph,N_g);
      const double shift_r_l = s->shift_r(r_ph,theta_ph,N_g);  // upper component
      const double phi_l   = s->phi(r_ph,theta_ph,N_g);
      const double gup_rr_l = c->gup_rr(r_ph,theta_ph,N_g);
      const double gup_rt_l = c->gup_rt(r_ph,theta_ph,N_g);
      const double gup_tt_l = c->gup_tt(r_ph,theta_ph,N_g);
      const double em4phi = exp(-4.0*phi_l);
      const double gup_rr_phys  = em4phi * gup_rr_l;
      const double gup_rt_phys  = em4phi * gup_rt_l;
      const double gup_tt_phys  = em4phi * gup_tt_l;
      const double dlapse_dr_l = aux->dlapse_dr(r_ph,theta_ph,N_g);
      const double dlapse_dt_l = aux->dlapse_dt(r_ph,theta_ph,N_g);
      const double dphi_dr_l = aux->dphi_dr(r_ph,theta_ph,N_g);
      const double dphi_dt_l = aux->dphi_dt(r_ph,theta_ph,N_g);
      const double dshift_r_dr_l = aux->dshift_r_dr(r_ph,theta_ph,N_g);
      const double dshift_r_dt_l = aux->dshift_r_dt(r_ph,theta_ph,N_g);
      //
      double shift_t_l = 0.0;
      double dshift_t_dr_l = 0.0;
      double dshift_t_dt_l = 0.0;
      double qup_r = gup_rr_phys * q_r;
      double qup_t = 0.0;
      if (r_ph > 0.0) {
	shift_t_l = s->shift_t(r_ph,theta_ph,N_g) / r_ph;  // (physical) upper component
	dshift_t_dr_l = aux->dshift_t_dr(r_ph,theta_ph,N_g) / r_ph - shift_t_l / r_ph ;
	dshift_t_dt_l = aux->dshift_t_dt(r_ph,theta_ph,N_g) / r_ph;
	qup_r = gup_rr_phys * q_r + gup_rt_phys * q_t;
	qup_t = gup_rt_phys * q_r + gup_tt_phys * q_t;       
      }	
      //================================================
      // update photon parameters; see eqs. (7.6) in "Numerical Relativity".
      // Note: divide equations by q / lapse to get d / dt instead of d / d\lambda
      //================================================
      const double dlambda = lapse_l * dt / q;
      const double term1_r = q * ( dlapse_dr_l - 2.0 * lapse_l * dphi_dr_l ); 
      const double term1_t = q * ( dlapse_dt_l - 2.0 * lapse_l * dphi_dt_l );
      double term2_r, term2_t, term3_r, term3_t;
      if (r_ph > 0.0) {
	term2_r = q_r * dshift_r_dr_l + q_t * dshift_t_dr_l;
	term2_t = q_r * dshift_r_dt_l + q_t * dshift_t_dt_l;
	term3_r = ( ( c->DG_r_rr(r_ph,theta_ph,N_g)              ) * qup_r * q_r + 
		    ( c->DG_t_rr(r_ph,theta_ph,N_g)              ) * qup_r * q_t + 
		    ( c->DG_r_rt(r_ph,theta_ph,N_g)              ) * qup_t * q_r + 
		    ( c->DG_t_rt(r_ph,theta_ph,N_g) + 1.0 / r_ph ) * qup_t * q_t ) * lapse_l / q;
	term3_t = ( ( c->DG_r_rt(r_ph,theta_ph,N_g)              ) * qup_r * q_r + 
		    ( c->DG_t_rt(r_ph,theta_ph,N_g) + 1.0 / r_ph ) * qup_r * q_t + 
		    ( c->DG_r_tt(r_ph,theta_ph,N_g) - r_ph       ) * qup_t * q_r + 
		    ( c->DG_t_tt(r_ph,theta_ph,N_g) + 1.0 / r_ph ) * qup_t * q_t ) * lapse_l / q;
      } else {
	term2_r = q_r * dshift_r_dr_l; 
	term2_t = q_r * dshift_r_dt_l;
	term3_r = c->DG_r_rr(r_ph,theta_ph,N_g) * qup_r * q_r * lapse_l / q;
	term3_t = c->DG_r_rt(r_ph,theta_ph,N_g) * qup_r * q_r * lapse_l / q;
      }	
      const double dq_r = ( - term1_r + term2_r + term3_r ) * dt;
      const double dq_t = ( - term1_t + term2_t + term3_t ) * dt;
      const double dr =     ( lapse_l * ( gup_rr_phys * q_r + gup_rt_phys * q_t ) / q - shift_r_l ) * dt;
      const double dtheta = ( lapse_l * ( gup_rt_phys * q_r + gup_tt_phys * q_t ) / q - shift_t_l ) * dt;
      lambda += dlambda;
      r_ph += dr;
      //
      // CHECK: hack: keep photon at constant angle!!  Valid only for symmetries!
      // theta_ph += dtheta;
      q_r += dq_r;
      q_t += dq_t;
      q = sqrt( q_r * qup_r + q_t * qup_t );
      if (q < 0.0) cout << " q = " << q << ", q_r = " << q_r << ", gup_rr_phys = " 
			  << gup_rr_phys << ", lapse_l = " << lapse_l << endl;
      //================================================
      // find function value
      //================================================
      if (fct1 == NULL)
	fct1_val = 0.0;
      else
	fct1_val   = (*fct1)(r_ph,theta_ph,N_g);
      if (fct2 == NULL)
	fct2_val = 0.0;
      else
	fct2_val   = (*fct2)(r_ph,theta_ph,N_g);
      // if (fct_val > fct_max) fct_max = fct_val;
      //================================================
      // check whether photon left grid
      //================================================
      photon_left_grid = (r_ph > r_out);
    } 
    return photon_left_grid;
  }
  //================================================
  //================================================
  // Interpolate to new datapoint
  //================================================
  //================================================
  bool ReachedLambda(const double Lambda, double & r, double & theta, double & fct1, double & fct2) {
    bool update = ((lambda_old <= Lambda) && (lambda > Lambda));
    if (update) {
      const double frac = (Lambda - lambda_old) / (lambda - lambda_old);
      r     = r_ph_old +     frac * (r_ph -     r_ph_old);
      theta = theta_ph_old + frac * (theta_ph - theta_ph_old);
      fct1   = fct1_val_old +  frac * (fct1_val -  fct1_val_old);
      fct2   = fct2_val_old +  frac * (fct2_val -  fct2_val_old);
    }
    return update;
  }
  //================================================
  //================================================
  // Gather info...
  //================================================
  //================================================
  bool GatherInfo(double & r, double & theta, double & lam, double & fct1, double & fct2,
		  double & qr, double & qtheta, double & qt) { 
    r = r_ph;
    theta = theta_ph;
    lam = lambda;
    fct1 = fct1_val;
    fct2 = fct2_val;
    qr = q_r;
    qtheta = q_t;
    qt = q;
    return photon_left_grid;
  }
private:
  //================================================
  //================================================
  // Compute self-similar time T
  //================================================
  //================================================
  inline double T(double tau) { return - log(fabs(tau_star - tau)); }
};

//================================================
//================================================
//================================================
// Class that handles swarm of photons 
//    emitted at one instant of time
//================================================
//================================================
//================================================
class PhotonSwarm {
private:
  int number;        // photon number
  double T_launch;   // self-similar time that photons were launched at
  double t_launch;   // coordinate time that photons were launched at
  double tau_launch; // proper time that photons were launched at
  double tau_star;
  int N;             // number of data points that will be collected 
  int n_ax, n_eq;    // number of data points collected so far
  int n_min;         // min(n_ax,n_eq)
  double PI;
  double delta_lambda; // record results at intervals of delta_lambda
  double *lambda;
  double *r_ax, *theta_ax;
  double *r_eq, *theta_eq;
  double *fct1_ax, *fct2_ax;
  double *fct1_eq, *fct2_eq;
  Photon *PhotonAx;
  Photon *PhotonEq;
  bool monitor;
  ofstream outfile;
  double next_lambda_ax;
  double next_lambda_eq;
  bool photon_left_grid_ax;
  bool photon_left_grid_eq;
  bool printed_message;
  //================================================
  //================================================
  // Constructor - just calls constructor for individual photons
  //================================================
  //================================================
public:
  PhotonSwarm() : number(0) 
  { 
    PhotonAx = new Photon();
    PhotonEq = new Photon();
    PI = acos(-1.0);
    photon_left_grid_ax = false;
    photon_left_grid_eq = false;
    printed_message = false;
  };
  //================================================
  //================================================
  // Destructor
  //================================================
  //================================================
  ~PhotonSwarm() {
    //    cout << " PHOTONSWARM: deleting swarm # " << number << " ... ";
    if (!(lambda == NULL)) delete lambda;
    if (!(r_ax == NULL)) delete r_ax;
    if (!(r_eq == NULL)) delete r_eq;
    if (!(theta_ax == NULL)) delete theta_ax;
    if (!(theta_eq == NULL)) delete theta_eq;
    if (!(fct1_ax == NULL)) delete fct1_ax;
    if (!(fct1_eq == NULL)) delete fct1_eq;
    if (!(fct2_ax == NULL)) delete fct2_ax;
    if (!(fct2_eq == NULL)) delete fct2_eq;
    delete PhotonAx;
    delete PhotonEq;
    if (outfile) outfile.close();
    //    cout << " done! " << endl;
  }
  //================================================
  //================================================
  // Launch swarm
  //================================================
  //================================================
  void Launch(int num, int N_i, double Delta_T,
	     double t, double tau, double tau_star_i, bool mon,
	      gf3d * fct1, gf3d * fct2, state * s, curvature * curv) {
    cout << " PHOTONSWARM: Launching PhotonSwarm # " << num << endl; 
    number = num;
    N = N_i;
    tau_star = tau_star_i;
    T_launch = T(tau);
    t_launch = t;
    tau_launch = tau;
    delta_lambda = Delta_T;
    n_ax = n_eq = 0;
    next_lambda_ax = next_lambda_eq = T_launch;
    monitor = mon;
    //================================================
    // allocate memory for data collection
    //================================================
    lambda = new double[N];
    r_ax = new double[N];
    r_eq = new double[N];
    theta_ax = new double[N];
    theta_eq = new double[N];
    fct1_ax = new double[N];
    fct1_eq = new double[N];
    fct2_ax = new double[N];
    fct2_eq = new double[N];
    for (int i = 0; i < N; i++) 
      lambda[i] = T_launch + double(i+1)*delta_lambda;
    //================================================
    // Launch photons
    //================================================
    const double theta1 = 0.0;
    const double theta2 = PI / 2.0;
    PhotonAx->Launch(number, theta1, t, tau, tau_star, fct1, fct2, s, curv);
    PhotonEq->Launch(number, theta2, t, tau, tau_star, fct1, fct2, s, curv);
    if (monitor) {
      //================================================
      // create output file
      //================================================
      ostringstream filename;
      filename << "PhotonSwarm_" << number << ".mon" << ends;
      outfile.open(filename.str().c_str());
      ios init(NULL);
      init.copyfmt(outfile);
      //================================================
      // write header
      //================================================
      outfile << "# Data for photon swarm " << number << " launched at coordinate time " << t
    	      << ", proper time " << tau << endl;
      outfile << "# and self-similar time " << T_launch << endl;
      outfile << "# " << setw(48) << "|" << " photon emitted in direction theta = " << setw(8) << theta2 
	      << setw(59) << "|" << " photon emitted in direction theta = " << setw(8) << theta1 << endl; 
      outfile << "# " << setw(14) << "coord time"
    	      << setw(16) << "proper time"
    	      << setw(16) << "self-sim time"
	      << " |" 
    	      << setw(14) << "lambda"
    	      << setw(16) << "coord radius"
    	      << setw(16) << "coord theta"
    	      << setw(16) << "q_r"
    	      << setw(16) << "q"
    	      << setw(24) << fct1->Name()
    	      << setw(24) << fct2->Name()
	      << " |" 
    	      << setw(14) << "lambda"
    	      << setw(16) << "coord radius"
    	      << setw(16) << "coord theta"
    	      << setw(16) << "q_r"
    	      << setw(16) << "q"
    	      << setw(24) << fct1->Name()
    	      << setw(24) << fct2->Name()
    	      << endl;
      //      outfile << "#==================================================================================================================================================================" << endl;
      outfile << "#" << setfill ('=') << setw(255) << "=" << endl;
      outfile.copyfmt(init);
    }
  }
  //================================================
  //================================================
  // Advance photon swarm
  //================================================
  //================================================
  void Advance(state * s, curvature * curv, auxiliary * aux, 
	       gf3d * fct1, gf3d * fct2,  double dt, const double r_out) {
    photon_left_grid_ax = PhotonAx->Advance(s, curv, aux, fct1, fct2, dt, r_out);
    photon_left_grid_eq = PhotonEq->Advance(s, curv, aux, fct1, fct2, dt, r_out);
    //================================================
    // check whether we can update data points
    //================================================
    double r, theta, f1, f2;
    if (PhotonAx->ReachedLambda(lambda[n_ax],r,theta,f1, f2) && n_ax < N) {
      r_ax[n_ax] = r;
      fct1_ax[n_ax] = f1;
      fct2_ax[n_ax] = f2;
      theta_ax[n_ax] = theta;
      n_ax += 1;
    }
    if (PhotonEq->ReachedLambda(lambda[n_eq],r,theta,f1, f2) && n_eq < N) {
      r_eq[n_eq] = r;
      fct1_eq[n_eq] = f1;
      fct2_eq[n_eq] = f2;
      theta_ax[n_ax] = theta;
      n_eq += 1;
    }
  }
  //================================================
  //================================================
  // Find maximum difference
  //================================================
  //================================================
  double GatherSwarmInfo(double & t, double & tau, double & T, 
			 double & maxdiff, double & r_max, bool & maxlast) {
    t = t_launch;
    tau = tau_launch;
    T = T_launch;
    n_min = n_eq;
    if (n_ax < n_min) n_min = n_ax;
    maxdiff = 0.0;
    maxlast = false;
    r_max = -1.0;
    for (int i = 0; i < n_min; i++) {
      if (fabs(fct1_ax[i] - fct1_eq[i]) > fabs(maxdiff)) {
	maxdiff = fct1_ax[i] - fct1_eq[i];
	r_max = r_eq[i];
	if (i == n_min-1) maxlast = true;
      }
    }
    return maxdiff;
  }
  //================================================
  //================================================
  // Monitor
  //================================================
  //================================================
  void Monitor(const double t, const double tau) { 
    if (monitor) {
      if (!(photon_left_grid_ax && photon_left_grid_eq)) {
      double lambda_ax = 0.0;
      double lambda_eq = 0.0;
      double radius_ax = 0.0; 
      double radius_eq = 0.0;
      double angle_ax = 0.0;
      double angle_eq = 0.0;
      double fct1_val_ax = 0.0; 
      double fct1_val_eq = 0.0;
      double fct2_val_ax = 0.0; 
      double fct2_val_eq = 0.0;
      double qr_ax = 0.0;
      double qr_eq = 0.0;
      double qangle_ax = 0.0;
      double qangle_eq = 0.0;
      double q_ax = 0.0;
      double q_eq = 0.0;
      PhotonEq->GatherInfo(radius_eq,angle_eq,lambda_eq,fct1_val_eq,fct2_val_eq,qr_eq,qangle_eq,q_eq);
      PhotonAx->GatherInfo(radius_ax,angle_ax,lambda_ax,fct1_val_ax,fct2_val_ax,qr_ax,qangle_ax,q_ax);
      outfile << setw(16) << setprecision(8) << t
    	      << setw(16) << setprecision(8) << tau
    	      << setw(16) << setprecision(8) << T(tau)
    	      << setw(16) << setprecision(8) << lambda_eq
    	      << setw(16) << setprecision(8) << radius_eq
    	      << setw(16) << setprecision(8) << angle_eq
    	      << setw(16) << setprecision(8) << qr_eq
    	      << setw(16) << setprecision(8) << q_eq
    	      << setw(24) << setprecision(16) << fct1_val_eq
    	      << setw(24) << setprecision(16) << fct2_val_eq
    	      << setw(16) << setprecision(8) << lambda_ax
    	      << setw(16) << setprecision(8) << radius_ax
    	      << setw(16) << setprecision(8) << angle_ax
    	      << setw(16) << setprecision(8) << qr_ax
    	      << setw(16) << setprecision(8) << q_ax
    	      << setw(24) << setprecision(16) << fct1_val_ax
    	      << setw(24) << setprecision(16) << fct2_val_ax
    	      << endl;
      } else {
	if (!printed_message) {
	  outfile << "# The photons left the grid! " << endl;
	  printed_message = true;
	}
      }
    }
  }
  //================================================
  //================================================
  // Compute self-similar time T
  //================================================
  //================================================
private:
  inline double T(double tau) { return - log(fabs(tau_star - tau)); }
};

//================================================
//================================================
//================================================
// Class that handles all photons
//================================================
//================================================
//================================================
class Photons {
private:
  int N;                // number of photons in each direction
  int N_monitored;      // number of photons that will write monitor files
  int N_skip;           // number of photons between monitored ones
  double tau_star;      // assumed proper time of accumulation event
  double T_launch_min;
  double T_launch_max;
  PhotonSwarm **PhotonSwarmList;
  double Delta_T;       // times T between launches
  double T_next_launch; // you get it...
  int N_launched;       // number of photons launched already, in each direction
  double PI;
  char * filestem;      // stem of monitor output file
  ostringstream monitorfilename;
public:
  //================================================
  //================================================
  // Constructor 
  //================================================
  //================================================
  Photons(char * filename, Grid *grid) :
    N_launched(0),
    filestem(filename) {
    Read_Input(N, N_monitored, tau_star, T_launch_min, T_launch_max);  
    PhotonSwarmList = new PhotonSwarm*[N];
    for (int i = 0; i < N; i++) {
      PhotonSwarmList[i] = new PhotonSwarm();
    }
    Delta_T = (T_launch_max - T_launch_min)/double(N - 1);
    T_next_launch = T_launch_min;
    N_skip = N / N_monitored;
    cout << " PHOTONS: Will monitor every " << N_skip << "th photon " << endl;
    PI = acos(-1.0);
    int N_r = grid->N_r_int();
    int N_theta = grid->N_theta_int();
    monitorfilename << filestem << "_" << N_r << "_" << N_theta 
		    << ".photon_mon" << ends;
  }
  //================================================
  //================================================
  // Destructor
  //================================================
  //================================================
  ~Photons() {
    cout << " PHOTONS: Destructing Photons... ";
    for (int i = 0; i < N; i++) {
      delete PhotonSwarmList[i];
    }
    delete PhotonSwarmList;
    cout << " done! " << endl;
  }
  //================================================
  //================================================
  // Update photons
  //================================================
  //================================================
  void Update_Photons(const double t, const double tau, const double dt,
		      state *s, curvature *c, auxiliary *aux,
		      gf3d * fct1, gf3d * fct2, const double r_out) {
    //================================================
    // advance photons that have been launched already
    //================================================
    for (int i = 0; i < N_launched; i++) {
      PhotonSwarmList[i]->Advance(s, c, aux, fct1, fct2, dt, r_out);
    }
    //================================================
    // time to launch new photon swarm?
    //================================================
    if (T(tau) >= T_next_launch && N_launched < N ) {
      // decide whether new photon should be monitored:
      bool mon = (N_launched % N_skip == 0);
      PhotonSwarmList[N_launched]->Launch(N_launched, N, Delta_T, t, tau, 
					  tau_star, mon, fct1, fct2, s, c);
      N_launched++;
      T_next_launch += Delta_T;
    }
  }
  //================================================
  //================================================
  // Monitor photons
  //================================================
  //================================================
  void Monitor(const double t, const double tau) {
    //================================================
    // Monitor individual monitored photons
    //================================================
    for (int i = 0; i < N_launched; i++) {
      PhotonSwarmList[i]->Monitor(t, tau);
    }      
    //================================================
    // Collect data from all photons: keep overwriting same file
    //================================================
    ofstream outfile;
    outfile.open(monitorfilename.str().c_str());
    outfile << "# Photon data at time " << t << " proper time " << tau << " and self-similar time " << T(tau) << endl;
    outfile << "# " << setw(14) << "Photon Number" 
	    << setw(16) << "t_launch" 
	    << setw(16) << "tau_launch" 
	    << setw(16) << "T_launch"
	    << setw(24) << "Max Diff" 
	    << setw(16) << "r_eq at max"
	    << endl;
    outfile << "#===============================================================================================================" << endl;
    for (int i = 0; i < N_launched; i++) {
      double t_launch;
      double tau_launch;
      double T_launch;
      double max_diff = 0.0;
      bool max_last = false;
      double r_max = 0.0;
      PhotonSwarmList[i]->GatherSwarmInfo(t_launch, tau_launch, 
					  T_launch, max_diff, r_max, max_last);
      outfile << setw(16) << setprecision(8) << i 
	      << setw(16) << setprecision(8) << t_launch
	      << setw(16) << setprecision(8) << tau_launch
	      << setw(16) << setprecision(8) << T_launch
	      << setw(24) << setprecision(16) << max_diff
	      << setw(16) << setprecision(8) << r_max;
      if (max_last) outfile << "  Maximum on last data point!";
      outfile << endl;
    }      
    outfile.close();
  }    
private:
  //================================================
  //================================================
  // Compute self-similar time T
  //================================================
  //================================================
  inline double T(double tau) { return - log(fabs(tau_star - tau)); }
  //================================================
  //================================================
  // Read Input
  //================================================
  //================================================
  int Read_Input(int & n, int & n_monitor, 
		 double & tau_star, double & T_min, double & T_max) {
    int error = 0;
    ifstream infile;
    infile.open("Photon_Input");
    if (!infile) {
      cerr << " PHOTONS: Can't open file Photon_Input for Input. " << endl;
      error = 1;
      return error;
    } 
    int buf_size = 500;
    char buf[buf_size], c;
    infile.get(buf,buf_size,'='); infile.get(c); infile >> n;
    infile.get(buf,buf_size,'='); infile.get(c); infile >> n_monitor;
    infile.get(buf,buf_size,'='); infile.get(c); infile >> tau_star;
    infile.get(buf,buf_size,'='); infile.get(c); infile >> T_min;
    infile.get(buf,buf_size,'='); infile.get(c); infile >> T_max;
    infile.close();
    //================================================
    // echo out input
    //================================================
    cout << " PHOTONS: will track " << n << " photons, with " << n_monitor 
	 << " monitored." << endl;
    cout << " PHOTONS: will assume tau_star = " << tau_star << endl;
    cout << " PHOTONS: will emit photons between self-similar times "
	 << T_min << " and " << T_max << endl;
    return error;
  }
};
  

#endif  /* PHOTONS_H */
