// Tell emacs that this is -*-c++-*- mode
//================================================
// Classes for initial data
//================================================
#ifndef INDATA_H
#define INDATA_H

// #define _DEBUG_INDATA_
#include "nr3.h"
#include "gridfunction.h"
#include "FlatEllSolver3D2ndOrder.h" 
#include "VecLaplacian.h"    
#include "VecLaplaceSolver.h"  // old version that solves for phi comp only
#include "tensors.h"
#include "Cosmology.h"
#include "Grid.h"
#include "EOS.h"
#include "State.h"
#include "Curvature.h"
#include "Maxwell_State.h"

enum {
    linwave, schwarzschild, flat, tov, trumpet, rns, brill,
    kerr, bowenyork, kerrschild, ken_trumpet, brill_lindquist,
    shock, evanscoleman, rad_tov, os, rotperfectfluid, choptuik,
    em_wave, radhydroshocktest, bondi, tov_bh, sms_tov, gauge_wave,
    wind_tunnel, chkpt, shibatawave, disk
};

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
    Cosmology* cosmology;
    Grid* grid;
public:
    // Constructor
    InData(Grid* grid_i, Cosmology* cosmology_i);

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
    virtual double lam_r_analytical(double r, double theta, double phi, double t, bool& done) = 0;
    virtual double lam_t_analytical(double r, double theta, double phi, double t, bool& done) = 0;
    virtual double lam_p_analytical(double r, double theta, double phi, double t, bool& done) = 0;
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
    virtual double a_r_analytical(double r, double theta, double phi, double t) = 0;
    virtual double a_t_analytical(double r, double theta, double phi, double t) = 0;
    virtual double a_p_analytical(double r, double theta, double phi, double t) = 0;
    virtual double as_r_analytical(double r, double theta, double phi, double t) = 0;
    virtual double as_t_analytical(double r, double theta, double phi, double t) = 0;
    virtual double as_p_analytical(double r, double theta, double phi, double t) = 0;
    virtual double E_analytical(double r, double theta, double phi, double t) = 0;
    virtual double F_0_analytical(double r, double theta, double phi, double t) = 0;
    virtual double F_r_analytical(double r, double theta, double phi, double t) = 0;
    virtual double F_t_analytical(double r, double theta, double phi, double t) = 0;
    virtual double F_p_analytical(double r, double theta, double phi, double t) = 0;
    virtual double Theta_analytical(double r, double theta, double phi, double t) = 0;
    virtual double B_r_analytical(double r, double theta, double phi, double t) = 0;
    virtual double B_t_analytical(double r, double theta, double phi, double t) = 0;
    virtual double B_p_analytical(double r, double theta, double phi, double t) = 0;
    virtual bool Initialize(gf3d& fct) = 0;

    void SelfSim(double& tau_star_r, double& xi_r);

    //================================================
    // Metric
    //================================================  
    bool Initialize_Metric(state* c);

    //================================================
    // Rescale metric so that determinant is r^4 sin^2 theta
    //================================================
    void Rescale_Metric(gf3d& h_rr, gf3d& h_rt, gf3d& h_rp,
        gf3d& h_tt, gf3d& h_tp, gf3d& h_pp,
        gf3d& phi);

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
    void Initialize_ExCurvature(state* c);

    //================================================
    // Lapse
    //================================================
    void Initialize_Lapse(state* c);

    //================================================
    // Shift
    //================================================
    void Initialize_Shift(state* c);

    //================================================
    // Matter
    //================================================
    void Initialize_Matter(gf3d& rho_0, gf3d& v_r, gf3d& v_t, gf3d& v_p);

    //================================================
    // Radiation Hydro
    //================================================
    void Initialize_RadHydro(gf3d& rho_0, gf3d& p, gf3d& v_r, gf3d& v_t, gf3d& v_p,
        gf3d& E, gf3d& F_0,
        gf3d& F_r, gf3d& F_t, gf3d& F_p);

    //================================================
    // Scalar Field
    //================================================
    void Initialize_Scalar_Field(gf3d& sf, gf3d& pi);

    //================================================
    // Initialize EM fields
    //================================================
    void Initialize_Maxwell(gf3d& e_p, gf3d& a_p);

    //================================================
    // Initialize Dual EM fields
    //================================================
    void Initialize_DualMaxwell(gf3d& a_r, gf3d& a_t, gf3d& a_p,
        gf3d& as_r, gf3d& as_t, gf3d& as_p);

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
    void Initialize_ConfConn(state* s, curvature* c);

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
    void Initialize_ConfConn(state* cc);

    //================================================
    // some simple utilities...
    //================================================
    bool Analytical_Solution_Exists();
    virtual string Name() = 0;

    //================================================
    // only relevant for matter:
    //================================================
    double Polytropic_Gamma() { return Gamma; }
    double Polytropic_K() { return Kappa; }

    //================================================
    // adjust grid for indata that are read from file
    //================================================
    int Adjust_Grid(double& r_out, int& N_r, int& N_theta);

    //================================================
    // and a transformation from cartesian to spherical polar coordinates 
    // for a vector with upstairs indices
    //================================================
    vect Cartesian_to_Spherical_upper(vect V, Doub r, Doub theta, Doub phi);

    //================================================
    // and a transformation from cartesian to spherical polar coordinates 
    // for rank-2 tensor with down-stairs indices
    //================================================
    tensor Cartesian_to_Spherical(tensor g, double r, double theta, double phi);
};


#endif  /* INDATA_H */
