// Tell emacs that this is -*-c++-*- mode
//================================================
// Classes for initial data
//================================================
#ifndef INDATA_H
#define INDATA_H

// #define _DEBUG_INDATA_
#include "nr3.h"
#include "gridfunction.h"
// #include "EllSolver3D.h" 
#include "FlatEllSolver3D2ndOrder.h" 
#include "VecLaplaceSolver.h"
#include "tensors.h"
#include "Cosmology.h"
#include "Grid.h"
#include "EOS.h"
#include "State.h"
#include "Curvature.h"
#include "Maxwell_State.h"

enum { linwave, schwarzschild, flat, tov, trumpet, rns, brill, 
       kerr, bowenyork, kerrschild, ken_trumpet, brill_lindquist, 
       shock, evanscoleman, rad_tov, os, rotperfectfluid, choptuik, 
       em_wave, radhydroshocktest, bondi, tov_bh, sms_tov, gauge_wave,
       wind_tunnel, chkpt, shibatawave };

//
//================================================
// Base class - doesn't do much...
//================================================
//
class InData {
protected:
  int indata_type;
  bool analytical;
  int N_g;   // number of ghost zones - hard-code to 2
  double Gamma; // polytropic exponent
  double Kappa; // polytropic constant, or P = Kappa rho for perfect fluid  
  double tau_star, xi;  // used for diagnostics in Scalar field evolution
  double PI;
  bool read_from_file;  // true if initial data are read from file
  double r_max;         // if so, need to set r_max and gridsize
  int nr, ntheta;     // to certain fixed values - routine FixGrid does that
  bool initialized;
  bool rescale_metric;
  Cosmology *cosmology;
  Grid *grid;
public:
  // Constructor
  InData(Grid * grid_i, Cosmology *cosmology_i) : 
    cosmology(cosmology_i), grid(grid_i) { 
    N_g = grid->N_ghosts(); 
    PI = acos(-1.0);
    initialized = false;
    // relevant only for matter:
    Gamma = 2.0; 
    Kappa = 0.0; 
    // relevant only for initial data that are read from file:
    read_from_file = false;
    r_max = 0.0;
    nr = 0;
    ntheta = 0;
  };
  // Destructor
  ~InData() {};
  //================================================
  // Virtual functions that will be overwritten by derived
  // classes below
  //================================================
  int Type() { return indata_type; }; 
  bool RescaleMetric(bool rescale) { return rescale_metric = rescale; }
  virtual double h_rr_analytical(double r, double theta, double phi, double t) = 0;
  virtual double h_rt_analytical(double r, double theta, double phi, double t) = 0;
  virtual double h_rp_analytical(double r, double theta, double phi, double t) = 0;
  virtual double h_tt_analytical(double r, double theta, double phi, double t) = 0;
  virtual double h_tp_analytical(double r, double theta, double phi, double t) = 0;
  virtual double h_pp_analytical(double r, double theta, double phi, double t) = 0;
  virtual double phi_analytical(double r, double theta, double phi, double t) = 0;
  virtual double lam_r_analytical(double r, double theta, double phi, double t, bool & done) = 0;
  virtual double lam_t_analytical(double r, double theta, double phi, double t, bool & done) = 0;
  virtual double lam_p_analytical(double r, double theta, double phi, double t, bool & done) = 0;
  virtual double a_rr_analytical(double r, double theta, double phi, double t) = 0;
  virtual double a_rt_analytical(double r, double theta, double phi, double t) = 0;
  virtual double a_rp_analytical(double r, double theta, double phi, double t) = 0;
  virtual double a_tt_analytical(double r, double theta, double phi, double t) = 0;
  virtual double a_tp_analytical(double r, double theta, double phi, double t) = 0;
  virtual double a_pp_analytical(double r, double theta, double phi, double t) = 0;
  virtual double K_analytical(double r, double theta, double phi, double t) = 0;
  virtual double lapse_analytical(double r, double theta, double phi, double t) = 0;
  virtual double shift_r_analytical(double r, double theta, double phi, double t) = 0;
  virtual double shift_t_analytical(double r, double theta, double phi, double t) = 0;
  virtual double shift_p_analytical(double r, double theta, double phi, double t) = 0;
  virtual double rho_0_analytical(double r, double theta, double phi, double t) = 0;
  virtual double P_analytical(double r, double theta, double phi, double t) = 0;
  virtual double v_r_analytical(double r, double theta, double phi, double t) = 0;
  virtual double v_t_analytical(double r, double theta, double phi, double t) = 0;
  virtual double v_p_analytical(double r, double theta, double phi, double t) = 0;
  virtual double sf_analytical(double r, double theta, double phi, double t) = 0;
  virtual double pi_analytical(double r, double theta, double phi, double t) = 0;
  virtual double e_p_analytical(double r, double theta, double phi, double t) = 0;
  virtual double a_p_analytical(double r, double theta, double phi, double t) = 0;
  virtual double E_analytical(double r, double theta, double phi, double t) = 0;
  virtual double F_0_analytical(double r, double theta, double phi, double t) = 0;
  virtual double F_r_analytical(double r, double theta, double phi, double t) = 0;
  virtual double F_t_analytical(double r, double theta, double phi, double t) = 0;
  virtual double F_p_analytical(double r, double theta, double phi, double t) = 0;
  virtual double Theta_analytical(double r, double theta, double phi, double t) = 0;
  virtual double B_r_analytical(double r, double theta, double phi, double t) = 0;
  virtual double B_t_analytical(double r, double theta, double phi, double t) = 0;
  virtual double B_p_analytical(double r, double theta, double phi, double t) = 0;
  virtual bool Initialize(gf3d & fct) = 0;
  void SelfSim(double & tau_star_r, double & xi_r) { 
    tau_star_r = tau_star; xi_r = xi; 
  }
  //================================================
  // Metric
  //================================================  
  bool Initialize_Metric(state * c) {
    if (!initialized) initialized = Initialize(c->h_rr);
    double t = 0.0;
    int N_r = c->h_rr.dim1();
    int N_theta = c->h_rr.dim2();
    int N_phi = c->h_rr.dim3();
    for (int i = N_g; i < N_r; i++) { // include outer boundary grid points    
      double rl = c->h_rr.r(i);
      for (int j = N_g; j < N_theta - N_g; j++) {
	double tl = c->h_rr.theta(j);
	for (int k = N_g; k < N_phi - N_g; k++) { 
	  double pl = c->h_rr.phi(k);
	  c->h_rr[i][j][k] = h_rr_analytical(rl,tl,pl,t);
	  c->h_rt[i][j][k] = h_rt_analytical(rl,tl,pl,t);
	  c->h_rp[i][j][k] = h_rp_analytical(rl,tl,pl,t);
	  c->h_tt[i][j][k] = h_tt_analytical(rl,tl,pl,t);
	  c->h_tp[i][j][k] = h_tp_analytical(rl,tl,pl,t);
	  c->h_pp[i][j][k] = h_pp_analytical(rl,tl,pl,t);
	  c->phi[i][j][k]  = phi_analytical(rl,tl,pl,t);
	  c->Theta[i][j][k] = Theta_analytical(rl,tl,pl,t);
	}
      }
    }
    c->h_rr.fill_ghosts();
    c->h_rt.fill_ghosts();
    c->h_rp.fill_ghosts();
    c->h_tt.fill_ghosts();
    c->h_tp.fill_ghosts();
    c->h_pp.fill_ghosts();
    c->phi.fill_ghosts();
    c->Theta.fill_ghosts();
    //
    // rescale metric if desired...
    //
    if (rescale_metric) {
      cout << " Rescaling initial metric so that det = 1 " << endl;
      Rescale_Metric(c->h_rr,c->h_rt,c->h_rp,c->h_tt,c->h_tp,c->h_pp,c->phi);
    }
    return initialized;
  };
  //================================================
  // Rescale metric so that determinant is r^4 sin^2 theta
  //================================================
  void Rescale_Metric(gf3d & h_rr, gf3d & h_rt, gf3d & h_rp, 
		      gf3d & h_tt, gf3d & h_tp, gf3d & h_pp,
		      gf3d & phi) {
    int N_r = h_rr.dim1();
    int N_theta = h_rr.dim2();
    int N_phi = h_rr.dim3();
    for (int i = 0; i < N_r; i++)    
      for (int j = 0; j < N_theta; j++)
	for (int k = 0; k < N_phi; k++) {
	  const double det = ( (1.0 + h_rr(i,j,k) ) * ( 1.0 + h_tt(i,j,k) ) * ( 1.0 + h_pp(i,j,k) ) + 
				 h_rt(i,j,k) * h_tp(i,j,k) * h_rp(i,j,k) + 
				 h_rp(i,j,k) * h_rt(i,j,k) * h_tp(i,j,k) - 
				 ( 1.0 + h_rr(i,j,k) ) * h_tp(i,j,k) * h_tp(i,j,k) -
				 h_rt(i,j,k) * h_rt(i,j,k) * ( 1.0 + h_pp(i,j,k) ) -
				 h_rp(i,j,k) * ( 1.0 + h_tt(i,j,k) ) * h_rp(i,j,k) );
	  const double factor = pow(1.0/det,1.0/3.0);
	  h_rr[i][j][k] = factor*(1.0 + h_rr(i,j,k)) - 1.0;
	  h_rt[i][j][k] *= factor;
	  h_rp[i][j][k] *= factor;
	  h_tt[i][j][k] = factor*(1.0 + h_tt(i,j,k)) - 1.0;
	  h_tp[i][j][k] *= factor;
	  h_pp[i][j][k] = factor*(1.0 + h_pp(i,j,k)) - 1.0;
	  phi[i][j][k] -= log(factor)/4.0;
	}
  };  
  //================================================
  // Extrinsic Curvature - interface for rotating equilibria
  //================================================
  // void Initialize_ExCurvature(gf3d & a_rr, gf3d & a_rt, gf3d & a_rp, 
  // 			     gf3d & a_tt, gf3d & a_tp, gf3d & a_pp,
  // 			     gf3d & K,
  // 			     gf3d & h_pp, gf3d & lapse, gf3d & shift_p) {
  //   cout << " Initializing Extrinsic Curvature for rotating equilibria " << endl;
  //   int N_r = a_rr.dim1();
  //   int N_theta = a_rr.dim2();
  //   int N_phi = a_rr.dim3();
  //   for (int i = N_g; i < N_r; i++) {
  //     for (int j = N_g; j < N_theta - N_g; j++)
  // 	for (int k = N_g; k < N_phi - N_g; k++) {    
  // 	  double rl = a_rr.r(i);
  // 	  double tl = a_rr.theta(j);
  // 	  double pl = a_rr.phi(k);
  // 	  double cottheta = cos(tl)/sin(tl);
  // 	  double factor = (h_pp(i,j,k) + 1.0) / (2.0 * lapse(i,j,k));
	  
  // 	  a_rr[i][j][k] = 0.0;
  // 	  a_rt[i][j][k] = 0.0;
  // 	  a_tt[i][j][k] = 0.0;
  // 	  a_pp[i][j][k] = 0.0;
  // 	  K[i][j][k]  = 0.0;
  // 	  a_tp[i][j][k] = factor * ( shift_p.dtheta(i,j,k) - shift_p(i,j,k)*cottheta ) / rl;
  // 	  if (i >= N_r - N_g) {  // use one-sided derivative (upwind) at outer boundary
  // 	    a_rp[i][j][k] = factor * ( shift_p.dr_OS(i,j,k,-1.0) - shift_p(i,j,k)/rl );
  // 	  } else { // ... centered derivative everywhere else
  // 	    a_rp[i][j][k] = factor * ( shift_p.dr(i,j,k) - shift_p(i,j,k)/rl );
  // 	  }
  // 	}
  //   }
  //   a_rr.fill_ghosts();
  //   a_rt.fill_ghosts();
  //   a_rp.fill_ghosts();
  //   a_tt.fill_ghosts();
  //   a_tp.fill_ghosts();
  //   a_pp.fill_ghosts();
  //   K.fill_ghosts();
  // }
  //================================================
  // Extrinsic Curvature - general interface
  //================================================
  void Initialize_ExCurvature(state *c) {
    //    cout << " Initializing Extrinsic Curvature " << endl;
    double t = 0.0;
    int N_r = c->a_rr.dim1();
    int N_theta = c->a_rr.dim2();
    int N_phi = c->a_rr.dim3();
    for (int i = N_g; i < N_r; i++) {     // include outer boundary grid points    
      double rl = c->a_rr.r(i);
      for (int j = N_g; j < N_theta - N_g; j++) {
	double tl = c->a_rr.theta(j);
	for (int k = N_g; k < N_phi - N_g; k++) {    
	  double pl = c->a_rr.phi(k);
	  c->a_rr[i][j][k] = a_rr_analytical(rl,tl,pl,t);
	  c->a_rt[i][j][k] = a_rt_analytical(rl,tl,pl,t);
	  c->a_rp[i][j][k] = a_rp_analytical(rl,tl,pl,t);
	  c->a_tt[i][j][k] = a_tt_analytical(rl,tl,pl,t);
	  c->a_tp[i][j][k] = a_tp_analytical(rl,tl,pl,t);
	  c->a_pp[i][j][k] = a_pp_analytical(rl,tl,pl,t);
	  c->K[i][j][k]  = K_analytical(rl,tl,pl,t);
	}
      }
    }
    //
    // if metric was rescaled, need to multiply A_ij with conformal factor...
    //
    if (rescale_metric) {
      cout << " Since we rescaled metric, will multiply A_ij with psi^{-4}..." << endl;
      for (int i = N_g; i < N_r; i++) {    
	double rl = c->a_rr.r(i);
	for (int j = N_g; j < N_theta - N_g; j++) {
	  double tl = c->a_rr.theta(j);
	  for (int k = N_g; k < N_phi - N_g; k++) {
	    // CHECK THIS LINE!!!
	    double pl = c->a_rr.phi(k);
	    double psim4 = exp(-4.0* ( c->phi(i,j,k) - phi_analytical(rl,tl,pl,t) ));
	    c->a_rr[i][j][k] *= psim4;
	    c->a_rt[i][j][k] *= psim4;
	    c->a_rp[i][j][k] *= psim4;
	    c->a_tt[i][j][k] *= psim4;
	    c->a_tp[i][j][k] *= psim4;
	    c->a_pp[i][j][k] *= psim4;
	  }
	}
      }
    }
    c->a_rr.fill_ghosts();
    c->a_rt.fill_ghosts();
    c->a_rp.fill_ghosts();
    c->a_tt.fill_ghosts();
    c->a_tp.fill_ghosts();
    c->a_pp.fill_ghosts();
    c->K.fill_ghosts();
  };
  //================================================
  // Lapse
  //================================================
  void Initialize_Lapse(state * c) {
    const double t = 0.0;
    int N_r = c->lapse.dim1();
    int N_theta = c->lapse.dim2();
    int N_phi = c->lapse.dim3();
    for (int i = N_g; i < N_r; i++)       // include outer boundary grid points
      for (int j = N_g; j < N_theta - N_g; j++) 
	for (int k = N_g; k < N_phi - N_g; k++) {
	  double rl = c->lapse.r(i);
	  double tl = c->lapse.theta(j);
	  double pl = c->lapse.phi(k);
	  c->lapse[i][j][k] = lapse_analytical(rl,tl,pl,t);
	}
    c->lapse.fill_ghosts();
  };
  //================================================
  // Shift
  //================================================
  void Initialize_Shift(state *c) {
    const double t = 0.0;
    int N_r = c->shift_r.dim1();
    int N_theta = c->shift_r.dim2();
    int N_phi = c->shift_r.dim3();
    for (int i = N_g; i < N_r; i++)    
      for (int j = N_g; j < N_theta - N_g; j++)
	for (int k = N_g; k < N_phi - N_g; k++) {
	  double rl = c->shift_r.r(i);
	  double tl = c->shift_r.theta(j);
	  double pl = c->shift_r.phi(k);
	  c->shift_r[i][j][k] = shift_r_analytical(rl,tl,pl,t);
	  c->shift_t[i][j][k] = shift_t_analytical(rl,tl,pl,t);
	  c->shift_p[i][j][k] = shift_p_analytical(rl,tl,pl,t);
	  c->B_r[i][j][k] = B_r_analytical(rl,tl,pl,t);
	  c->B_t[i][j][k] = B_t_analytical(rl,tl,pl,t);
	  c->B_p[i][j][k] = B_p_analytical(rl,tl,pl,t);
	}
    c->shift_r.fill_ghosts();
    c->shift_t.fill_ghosts();
    c->shift_p.fill_ghosts();
    c->B_r.fill_ghosts();
    c->B_t.fill_ghosts();
    c->B_p.fill_ghosts();
  };
  //================================================
  // Matter
  //================================================
  void Initialize_Matter(gf3d & rho_0, gf3d & v_r, gf3d & v_t, gf3d & v_p) {
    const double t = 0.0;
    int N_r = rho_0.dim1();
    int N_theta = rho_0.dim2();
    int N_phi = rho_0.dim3();
    for (int i = N_g; i < N_r; i++)    
      for (int j = N_g; j < N_theta - N_g; j++)
	for (int k = N_g; k < N_phi - N_g; k++) {    
	  double rl = rho_0.r(i);
	  double tl = rho_0.theta(j);
	  double pl = rho_0.phi(k);
	  rho_0[i][j][k] = rho_0_analytical(rl,tl,pl,t);
	  v_r[i][j][k] = v_r_analytical(rl,tl,pl,t);
	  v_t[i][j][k] = v_t_analytical(rl,tl,pl,t);
	  v_p[i][j][k] = v_p_analytical(rl,tl,pl,t);
	}
    rho_0.fill_ghosts();
    v_r.fill_ghosts();
    v_t.fill_ghosts();
    v_p.fill_ghosts();
  };
  //================================================
  // Radiation Hydro
  //================================================
  void Initialize_RadHydro(gf3d & rho_0, gf3d & p, gf3d & v_r, gf3d & v_t, gf3d & v_p,
			   gf3d & E, gf3d & F_0,
			   gf3d & F_r, gf3d & F_t, gf3d & F_p) {
    const double t = 0.0;
    int N_r = rho_0.dim1();
    int N_theta = rho_0.dim2();
    int N_phi = rho_0.dim3();
    for (int i = N_g; i < N_r; i++) {   
      double rl = rho_0.r(i);
      for (int j = N_g; j < N_theta - N_g; j++) {
	double tl = rho_0.theta(j);
	for (int k = N_g; k < N_phi - N_g; k++) {    
	  double pl = rho_0.phi(k);
	  rho_0[i][j][k] = rho_0_analytical(rl,tl,pl,t);
	  p[i][j][k]     = P_analytical(rl,tl,pl,t);
	  v_r[i][j][k]   = v_r_analytical(rl,tl,pl,t);
	  v_t[i][j][k]   = v_t_analytical(rl,tl,pl,t);
	  v_p[i][j][k]   = v_p_analytical(rl,tl,pl,t);
	  E[i][j][k]     = E_analytical(rl,tl,pl,t);
	  F_0[i][j][k]   = F_0_analytical(rl,tl,pl,t);
	  F_r[i][j][k]   = F_r_analytical(rl,tl,pl,t);
	  F_t[i][j][k]   = F_t_analytical(rl,tl,pl,t);
	  F_p[i][j][k]   = F_p_analytical(rl,tl,pl,t);
	}
      }
    }
    rho_0.fill_ghosts();
    v_r.fill_ghosts();
    v_t.fill_ghosts();
    v_p.fill_ghosts();
    E.fill_ghosts();
    F_0.fill_ghosts();
    F_r.fill_ghosts();
    F_t.fill_ghosts();
    F_p.fill_ghosts();
  };
  //================================================
  // Scalar Field
  //================================================
  void Initialize_Scalar_Field(gf3d & sf, gf3d & pi) {
    const double t = 0.0;
    int N_r = sf.dim1();
    int N_theta = sf.dim2();
    int N_phi = sf.dim3();
    for (int i = N_g; i < N_r; i++)    
      for (int j = N_g; j < N_theta - N_g; j++)
	for (int k = N_g; k < N_phi - N_g; k++) {    
	  double rl = sf.r(i);
	  double tl = sf.theta(j);
	  double pl = sf.phi(k);
	  sf[i][j][k] = sf_analytical(rl,tl,pl,t);
	  pi[i][j][k] = pi_analytical(rl,tl,pl,t);
	}
    sf.fill_ghosts();
    pi.fill_ghosts();
  };
  //================================================
  // Initialize EM fields
  //================================================
  void Initialize_Maxwell(gf3d & e_p, gf3d & a_p) {
    const double t = 0.0;
    int N_r = e_p.dim1();
    int N_t = e_p.dim2();
    int N_p = e_p.dim3();
    for (int i = N_g; i < N_r; i++) {   
      double rl = e_p.r(i);
      for (int j = N_g; j < N_t - N_g; j++) {
	double tl = e_p.theta(j);
	for (int k = N_g; k < N_p - N_g; k++) { 
	  double pl = e_p.phi(k);
	  e_p[i][j][k] = e_p_analytical(rl,tl,pl,t);
	  a_p[i][j][k] = a_p_analytical(rl,tl,pl,t);
	}
      }
    }
    e_p.fill_ghosts();
    a_p.fill_ghosts();
  };
  //================================================
  // Conformal Connection 
  //
  // use this version if we cannot compute metric analytically for 
  // any coordinate value; in that case, initialize from trace of
  // connection (needed for RNS initial data)
  //
  // NOTE: this version assumes that Compute_Connection() has already been 
  //       called.
  //================================================
  void Initialize_ConfConn(state *s, curvature *c) {
    //    cout << " INDATA: Initializing Conformal Connection Functions... " << endl; 
    double t = 0.0;
    int N_r = s->lam_r.dim1();
    int N_theta = s->lam_r.dim2();
    int N_phi = s->lam_r.dim3();
    bool done;
    for (int i = N_g; i < N_r; i++) {    
      const double rl = s->lam_r.r(i);
      for (int j = N_g; j < N_theta - N_g; j++) {
	const double tl = s->lam_r.theta(j);
	for (int k = N_g; k < N_phi - N_g; k++) {
	  const double pl = s->lam_r.phi(k);
	  s->lam_r[i][j][k] = lam_r_analytical(rl,tl,pl,t,done);
	  if (!done) s->lam_r[i][j][k] = c->DG_r(i,j,k);
	  s->lam_t[i][j][k] = lam_t_analytical(rl,tl,pl,t,done);
	  if (!done) s->lam_t[i][j][k] = c->DG_t(i,j,k);
	  s->lam_p[i][j][k] = lam_p_analytical(rl,tl,pl,t,done);
	  if (!done) s->lam_p[i][j][k] = c->DG_p(i,j,k);
	}
      }
    }
    s->lam_r.fill_ghosts();
    s->lam_t.fill_ghosts();
    s->lam_p.fill_ghosts();
  };
  // void Initialize_ConfConn(gf3d & DG_r, gf3d & DG_t, gf3d & DG_p,
  // 			  gf3d & lam_r, gf3d & lam_t, gf3d & lam_p) {
  //   //    cout << " INDATA: Initializing Conformal Connection Functions... " << endl; 
  //   int N_r = lam_r.dim1();
  //   int N_theta = lam_r.dim2();
  //   int N_phi = lam_r.dim3();
  //   for (int i = N_g; i < N_r; i++)    
  //     for (int j = N_g; j < N_theta - N_g; j++)
  // 	for (int k = N_g; k < N_phi - N_g; k++) {
  // 	  lam_r[i][j][k] = DG_r(i,j,k);
  // 	  lam_t[i][j][k] = DG_t(i,j,k);
  // 	  lam_p[i][j][k] = DG_p(i,j,k);
  // 	}
  //   lam_r.fill_ghosts();
  //   lam_t.fill_ghosts();
  //   lam_p.fill_ghosts();
  // };
  //================================================
  // Conformal Connection 
  //
  // otherwise compute "analytically" from metric everywhere
  //================================================
  void Initialize_ConfConn(state * cc) {
    double t = 0.0;
    int N_r = cc->lam_r.dim1();
    int N_theta = cc->lam_r.dim2();
    int N_phi = cc->lam_r.dim3();
    const double dr_min = 1.0e-6;
    const double dt_min = 1.0e-6;
    const double dp_min = 1.0e-6;
    double dr = dr_min;
    double dt = dt_min;
    double dp = dp_min;
    if (dr > cc->lam_r.r(N_g)/2.0)     dr = cc->lam_r.r(N_g)/2.0;
    if (dt > cc->lam_r.theta(N_g)/2.0) dt = cc->lam_r.theta(N_g)/2.0;
    if (dp > cc->lam_r.phi(N_g)/2.0)   dp = cc->lam_r.phi(N_g)/2.0;
    for (int i = N_g; i < N_r; i++)    
      for (int j = N_g; j < N_theta - N_g; j++)
	for (int k = N_g; k < N_phi - N_g; k++) {
	  double rl = cc->lam_r.r(i);
	  double tl = cc->lam_r.theta(j);
	  double pl = cc->lam_r.phi(k);
	  double st = cc->lam_r.sintheta(j);
	  double rl2 = rl*rl;
	  double st2 = st*st;
	  double ct = cc->lam_r.costheta(j);
	  //
	  // get h_{ij} at current location
	  //
	  const double h_rr = h_rr_analytical(rl,tl,pl,t);
	  const double h_rt = h_rt_analytical(rl,tl,pl,t);
	  const double h_rp = h_rp_analytical(rl,tl,pl,t);
	  const double h_tt = h_tt_analytical(rl,tl,pl,t);
	  const double h_tp = h_tp_analytical(rl,tl,pl,t);
	  const double h_pp = h_pp_analytical(rl,tl,pl,t);
	  //
	  // compute partial derivatives of h_ij
	  //
	  // const double d_r_h_rr = (h_rr_analytical(rl+dr,tl,pl,t) - h_rr_analytical(rl-dr,tl,pl,t))/(2.0*dr);
	  // const double d_r_h_rt = (h_rt_analytical(rl+dr,tl,pl,t) - h_rt_analytical(rl-dr,tl,pl,t))/(2.0*dr);
	  // const double d_r_h_rp = (h_rp_analytical(rl+dr,tl,pl,t) - h_rp_analytical(rl-dr,tl,pl,t))/(2.0*dr);
	  // const double d_r_h_tt = (h_tt_analytical(rl+dr,tl,pl,t) - h_tt_analytical(rl-dr,tl,pl,t))/(2.0*dr);
	  // const double d_r_h_tp = (h_tp_analytical(rl+dr,tl,pl,t) - h_tp_analytical(rl-dr,tl,pl,t))/(2.0*dr);
	  // const double d_r_h_pp = (h_pp_analytical(rl+dr,tl,pl,t) - h_pp_analytical(rl-dr,tl,pl,t))/(2.0*dr);

	  const double d_r_h_rr = (- h_rr_analytical(rl+2.0*dr,tl,pl,t)  + h_rr_analytical(rl-2.0*dr,tl,pl,t) 
				   + 8.0*(h_rr_analytical(rl+dr,tl,pl,t) - h_rr_analytical(rl-dr,tl,pl,t)) )/(12.0*dr);
	  const double d_r_h_rt = (- h_rt_analytical(rl+2.0*dr,tl,pl,t)  + h_rt_analytical(rl-2.0*dr,tl,pl,t) 
				   + 8.0*(h_rt_analytical(rl+dr,tl,pl,t) - h_rt_analytical(rl-dr,tl,pl,t)) )/(12.0*dr);
	  const double d_r_h_rp = (- h_rp_analytical(rl+2.0*dr,tl,pl,t)  + h_rp_analytical(rl-2.0*dr,tl,pl,t) 
				   + 8.0*(h_rp_analytical(rl+dr,tl,pl,t) - h_rp_analytical(rl-dr,tl,pl,t)) )/(12.0*dr);
	  const double d_r_h_tt = (- h_tt_analytical(rl+2.0*dr,tl,pl,t)  + h_tt_analytical(rl-2.0*dr,tl,pl,t) 
				   + 8.0*(h_tt_analytical(rl+dr,tl,pl,t) - h_tt_analytical(rl-dr,tl,pl,t)) )/(12.0*dr);
	  const double d_r_h_tp = (- h_tp_analytical(rl+2.0*dr,tl,pl,t)  + h_tp_analytical(rl-2.0*dr,tl,pl,t) 
				   + 8.0*(h_tp_analytical(rl+dr,tl,pl,t) - h_tp_analytical(rl-dr,tl,pl,t)) )/(12.0*dr);
	  const double d_r_h_pp = (- h_pp_analytical(rl+2.0*dr,tl,pl,t)  + h_pp_analytical(rl-2.0*dr,tl,pl,t) 
				   + 8.0*(h_pp_analytical(rl+dr,tl,pl,t) - h_pp_analytical(rl-dr,tl,pl,t)) )/(12.0*dr);

	  // const double d_t_h_rr = (h_rr_analytical(rl,tl+dt,pl,t) - h_rr_analytical(rl,tl-dt,pl,t))/(2.0*dt);
	  // const double d_t_h_rt = (h_rt_analytical(rl,tl+dt,pl,t) - h_rt_analytical(rl,tl-dt,pl,t))/(2.0*dt);
	  // const double d_t_h_rp = (h_rp_analytical(rl,tl+dt,pl,t) - h_rp_analytical(rl,tl-dt,pl,t))/(2.0*dt);
	  // const double d_t_h_tt = (h_tt_analytical(rl,tl+dt,pl,t) - h_tt_analytical(rl,tl-dt,pl,t))/(2.0*dt);
	  // const double d_t_h_tp = (h_tp_analytical(rl,tl+dt,pl,t) - h_tp_analytical(rl,tl-dt,pl,t))/(2.0*dt);
	  // const double d_t_h_pp = (h_pp_analytical(rl,tl+dt,pl,t) - h_pp_analytical(rl,tl-dt,pl,t))/(2.0*dt);

	  // CHECK!!!
	  // if (i == 4 && ( j == 5 || j == 6 ) && k == 2 ) {
	  //   //	    cout << " Lam_t(" << j << ") = " << setprecision(16) << cc->lam_t(i,j,k) << endl;
	  //   cout << " d_t_h_rr(" << j << ") = " << setprecision(16) << h_pp_analytical(rl,tl+dt,pl,t)  << "  " <<  setprecision(16) << h_pp_analytical(rl,tl-dt,pl,t) << "  " << (h_pp_analytical(rl,tl-dt,pl,t) - h_pp_analytical(rl,tl+dt,pl,t))/(2.0) << endl;
	  // }

	  const double d_t_h_rr = (- h_rr_analytical(rl,tl+2.0*dt,pl,t)  + h_rr_analytical(rl,tl-2.0*dt,pl,t) 
				   + 8.0*(h_rr_analytical(rl,tl+dt,pl,t) - h_rr_analytical(rl,tl-dt,pl,t)) )/(12.0*dt);
	  const double d_t_h_rt = (- h_rt_analytical(rl,tl+2.0*dt,pl,t)  + h_rt_analytical(rl,tl-2.0*dt,pl,t) 
				   + 8.0*(h_rt_analytical(rl,tl+dt,pl,t) - h_rt_analytical(rl,tl-dt,pl,t)) )/(12.0*dt);
	  const double d_t_h_rp = (- h_rp_analytical(rl,tl+2.0*dt,pl,t)  + h_rp_analytical(rl,tl-2.0*dt,pl,t) 
				   + 8.0*(h_rp_analytical(rl,tl+dt,pl,t) - h_rp_analytical(rl,tl-dt,pl,t)) )/(12.0*dt);
	  const double d_t_h_tt = (- h_tt_analytical(rl,tl+2.0*dt,pl,t)  + h_tt_analytical(rl,tl-2.0*dt,pl,t) 
				   + 8.0*(h_tt_analytical(rl,tl+dt,pl,t) - h_tt_analytical(rl,tl-dt,pl,t)) )/(12.0*dt);
	  const double d_t_h_tp = (- h_tp_analytical(rl,tl+2.0*dt,pl,t)  + h_tp_analytical(rl,tl-2.0*dt,pl,t) 
				   + 8.0*(h_tp_analytical(rl,tl+dt,pl,t) - h_tp_analytical(rl,tl-dt,pl,t)) )/(12.0*dt);
	  const double d_t_h_pp = (- h_pp_analytical(rl,tl+2.0*dt,pl,t)  + h_pp_analytical(rl,tl-2.0*dt,pl,t) 
				   + 8.0*(h_pp_analytical(rl,tl+dt,pl,t) - h_pp_analytical(rl,tl-dt,pl,t)) )/(12.0*dt);

	  // const double d_p_h_rr = (h_rr_analytical(rl,tl,pl+dp,t) - h_rr_analytical(rl,tl,pl-dp,t))/(2.0*dp);
	  // const double d_p_h_rt = (h_rt_analytical(rl,tl,pl+dp,t) - h_rt_analytical(rl,tl,pl-dp,t))/(2.0*dp);
	  // const double d_p_h_rp = (h_rp_analytical(rl,tl,pl+dp,t) - h_rp_analytical(rl,tl,pl-dp,t))/(2.0*dp);
	  // const double d_p_h_tt = (h_tt_analytical(rl,tl,pl+dp,t) - h_tt_analytical(rl,tl,pl-dp,t))/(2.0*dp);
	  // const double d_p_h_tp = (h_tp_analytical(rl,tl,pl+dp,t) - h_tp_analytical(rl,tl,pl-dp,t))/(2.0*dp);
	  // const double d_p_h_pp = (h_pp_analytical(rl,tl,pl+dp,t) - h_pp_analytical(rl,tl,pl-dp,t))/(2.0*dp);

	  const double d_p_h_rr = (- h_rr_analytical(rl,tl,pl+2.0*dp,t)  + h_rr_analytical(rl,tl,pl-2.0*dp,t) 
				   + 8.0*(h_rr_analytical(rl,tl,pl+dp,t) - h_rr_analytical(rl,tl,pl-dp,t)) )/(12.0*dp);
	  const double d_p_h_rt = (- h_rt_analytical(rl,tl,pl+2.0*dp,t)  + h_rt_analytical(rl,tl,pl-2.0*dp,t) 
				   + 8.0*(h_rt_analytical(rl,tl,pl+dp,t) - h_rt_analytical(rl,tl,pl-dp,t)) )/(12.0*dp);
	  const double d_p_h_rp = (- h_rp_analytical(rl,tl,pl+2.0*dp,t)  + h_rp_analytical(rl,tl,pl-2.0*dp,t) 
				   + 8.0*(h_rp_analytical(rl,tl,pl+dp,t) - h_rp_analytical(rl,tl,pl-dp,t)) )/(12.0*dp);
	  const double d_p_h_tt = (- h_tt_analytical(rl,tl,pl+2.0*dp,t)  + h_tt_analytical(rl,tl,pl-2.0*dp,t) 
				   + 8.0*(h_tt_analytical(rl,tl,pl+dp,t) - h_tt_analytical(rl,tl,pl-dp,t)) )/(12.0*dp);
	  const double d_p_h_tp = (- h_tp_analytical(rl,tl,pl+2.0*dp,t)  + h_tp_analytical(rl,tl,pl-2.0*dp,t) 
				   + 8.0*(h_tp_analytical(rl,tl,pl+dp,t) - h_tp_analytical(rl,tl,pl-dp,t)) )/(12.0*dp);
	  const double d_p_h_pp = (- h_pp_analytical(rl,tl,pl+2.0*dp,t)  + h_pp_analytical(rl,tl,pl-2.0*dp,t) 
				   + 8.0*(h_pp_analytical(rl,tl,pl+dp,t) - h_pp_analytical(rl,tl,pl-dp,t)) )/(12.0*dp);

	  //
	  // now compute partial derivatives of metric \bar \gamma_{ij}
	  //
	  double d_r_g_rr =             d_r_h_rr;
	  double d_r_g_rt = rl        * d_r_h_rt +                         h_rt ;
	  double d_r_g_rp = rl  * st  * d_r_h_rp +            st  *        h_rp ;
	  double d_r_g_tt = rl2       * d_r_h_tt + 2.0 * rl       * (1.0 + h_tt);
	  double d_r_g_tp = rl2 * st  * d_r_h_tp + 2.0 * rl * st  *        h_tp ;
	  double d_r_g_pp = rl2 * st2 * d_r_h_pp + 2.0 * rl * st2 * (1.0 + h_pp);

	  double d_t_g_rr =             d_t_h_rr;
	  double d_t_g_rt = rl        * d_t_h_rt;
	  double d_t_g_rp = rl  * st  * d_t_h_rp +       rl  * ct      *        h_rp ;
	  double d_t_g_tt = rl2       * d_t_h_tt;
	  double d_t_g_tp = rl2 * st  * d_t_h_tp +       rl2 * ct      *        h_tp ;
	  double d_t_g_pp = rl2 * st2 * d_t_h_pp + 2.0 * rl2 * st * ct * (1.0 + h_pp);

	  double d_p_g_rr =             d_p_h_rr;
	  double d_p_g_rt = rl        * d_p_h_rt;
	  double d_p_g_rp = rl  * st  * d_p_h_rp;
	  double d_p_g_tt = rl2       * d_p_h_tt;
	  double d_p_g_tp = rl2 * st  * d_p_h_tp;
	  double d_p_g_pp = rl2 * st2 * d_p_h_pp;
	  // if (i == 2 && ( j == 4 || j == 7 ) && k == 2 ) {
	  //   //	    cout << " Lam_t(" << j << ") = " << setprecision(16) << cc->lam_t(i,j,k) << endl;
	  //   cout << " d_t_h_rr(" << j << ") = " << setprecision(16) << d_t_h_rr << "  " <<  setprecision(16) << d_t_h_tt << "  "<< setprecision(16) << d_t_h_pp << endl;
	  // }
	  //
	  // now store things in tensors...
	  //
	  rank3tens D_g(d_r_g_rr,d_r_g_rt,d_r_g_rp,d_r_g_tt,d_r_g_tp,d_r_g_pp,
			d_t_g_rr,d_t_g_rt,d_t_g_rp,d_t_g_tt,d_t_g_tp,d_t_g_pp,
			d_p_g_rr,d_p_g_rt,d_p_g_rp,d_p_g_tt,d_p_g_tp,d_p_g_pp);
	  tensor g(1.0 + h_rr, rl*h_rt, rl*st*h_rp, rl2*(1.0 + h_tt), rl2*st*h_tp, rl2*st2*(1.0 + h_pp));
	  tensor gup;
	  gup = g.inverse();
	  cc->lam_r[i][j][k] = 0.0;
	  cc->lam_t[i][j][k] = 0.0;
	  cc->lam_p[i][j][k] = 0.0;
	  for (int b = 0; b < 3; b++)
	    for (int c = 0; c < 3; c++)
	      for (int d = 0; d < 3; d++) {
		cc->lam_r[i][j][k] += 0.5 * gup[b][c] * gup[0][d] * ( D_g[b][c][d] + D_g[c][b][d] - D_g[d][b][c] );
		cc->lam_t[i][j][k] += 0.5 * gup[b][c] * gup[1][d] * ( D_g[b][c][d] + D_g[c][b][d] - D_g[d][b][c] );
		cc->lam_p[i][j][k] += 0.5 * gup[b][c] * gup[2][d] * ( D_g[b][c][d] + D_g[c][b][d] - D_g[d][b][c] );
	      }
	  //
	  // subtract flat connection to obtain Lambda^i
	  //
	  cc->lam_r[i][j][k] -= gup[1][1] * (- rl ) + gup[2][2] * (- rl * st2 );
	  cc->lam_t[i][j][k] -= gup[2][2] * (-st*ct ) + 2.0 * gup[0][1] / rl;
	  cc->lam_p[i][j][k] -= 2.0 * gup[0][2] / rl + 2.0*gup[1][2] * ct/st;
	  //
	  // Now rescale to find lambda as used in code
	  //
	  cc->lam_t[i][j][k] *= rl;
	  cc->lam_p[i][j][k] *= rl*st;
	}
    //
    // CHECK this hack!!
    // 
    //    cc->lam_r.equals(0.0);
    //    cc->lam_t.equals(0.0);
    //    cc->lam_p.equals(0.0);
    cc->lam_r.fill_ghosts();
    cc->lam_t.fill_ghosts();
    cc->lam_p.fill_ghosts();
  };
//================================================
// some simple utilities...
//================================================
  bool Analytical_Solution_Exists() { return analytical; }
  virtual string Name() = 0; 
  //================================================
  // only relevant for matter:
  //================================================
  double Polytropic_Gamma() { return Gamma; }
  double Polytropic_K() { return Kappa; }
//================================================
// adjust grid for indata that are read from file
//================================================
  int Adjust_Grid(double & r_out, int & N_r, int & N_theta) {
    if (read_from_file) {     // don't do anything unless data read from file
      r_out = r_max;
      N_r = nr;
      if (ntheta > 0) {        // change ntheta only if needed (set to value > 0 in derived class)
	N_theta = ntheta;
	cout << " NOTE: adjusted grid to r_max = " << r_out << ", N_r = " << N_r 
	     << ", N_theta = " << N_theta << endl;
      } else {
	cout << " NOTE: adjusted grid to r_max = " << r_out << ", N_r = " << N_r << endl;
      }
    }
    return 0;
  };
//================================================
// and a transformation from cartesian to spherical polar coordinates 
// for a vector with upstairs indices
//================================================
vect Cartesian_to_Spherical_upper(vect V, Doub r, Doub theta, Doub phi) {
  const Doub st = sin(theta);
  const Doub ct = cos(theta);
  const Doub sp = sin(phi);
  const Doub cp = cos(phi);  
  //
  // compute transformation matrix \partial x^i / partial x^j'
  // 
  tensor Lambda(st*cp,      st*sp,     ct,
                ct*cp/r,    ct*sp/r,   -st/r,
                -sp/(r*st), cp/(r*st), 0.0);
  vect V_sc;
  for (int i=0; i<3; i++) {
    V_sc[i] = 0.0;
    for (int k=0; k<3; k++)
          V_sc[i] += Lambda[i][k] * V[k];
    }
  return V_sc;
};

//================================================
// and a transformation from cartesian to spherical polar coordinates 
// for rank-2 tensor with down-stairs indices
//================================================
  tensor Cartesian_to_Spherical(tensor g, double r, double theta, double phi) {
    Doub st = sin(theta);
    Doub ct = cos(theta);
    Doub sp = sin(phi);
    Doub cp = cos(phi);
    //
    // compute transformation matrix \partial x^i / partial x^j'
    // 
    tensor Lambda(st*cp,    st*sp,   ct,
		  r*ct*cp,  r*ct*sp, -r*st,
		  -r*st*sp, r*st*cp, 0.0);
    tensor g_sc;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++) {
	g_sc[i][j] = 0.0;
	for (int k=0; k<3; k++)
	  for (int l=0; l<3; l++)
	    g_sc[i][j] += Lambda[i][k] * Lambda[j][l] * g[k][l];
      }
    return g_sc;
  };
};
//
//=======================================================
// Now include files for individual types of initial data
//=======================================================
//
#include "INDATA/LinWave.h"
#include "INDATA/Schwarzschild.h"
#include "INDATA/Flat.h"
#include "INDATA/TOV.h"
#include "INDATA/Trumpet.h"
#include "INDATA/RNS.h"
#include "INDATA/Brill.h"
#include "INDATA/Kerr.h"
#include "INDATA/BowenYork.h"
#include "INDATA/KerrSchild.h"
#include "INDATA/KenTrumpet.h"
#include "INDATA/BrillLindquist.h"
#include "INDATA/Shock.h"
#include "INDATA/EvansColeman.h"
#include "INDATA/Rad_TOV.h"
#include "INDATA/OS.h"
#include "INDATA/RotPerfectFluid.h"
#include "INDATA/Choptuik.h"
#include "INDATA/EMWave.h"
#include "INDATA/RadHydroShockTest.h"
#include "INDATA/Bondi.h"
#include "INDATA/TOV_BH.h"
#include "INDATA/SMS_TOV.h"
#include "INDATA/GaugeWave.h"
#include "INDATA/WindTunnel.h"
#include "INDATA/ReadFromCheckPoint.h"
#include "INDATA/ShibataWave.h"

#endif  /* INDATA_H */
