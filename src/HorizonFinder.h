// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Horizon finder
//
//================================================

#ifndef HORIZONFINDER_H
#define HORIZONFINDER_H

#include "EllSolver2D.h"

#include "nr3.h"
#include "surfacefunction.h"
#include "gridfunction.h"
// #include "Monitor.h"
#include "Grid.h"
#include "State.h"
#include "Curvature.h"
#include "Auxiliary.h"
#include "Fluxes.h"
#include <ctime>

enum { proper, unitsphere };

//
//================================================
// 
//================================================
//
class HorizonFinder {
private:
    // function h in level surface function tau = r - h(theta,phi)
    gf2d h, h_old;
    // projected metric
    gf2d gup_rr, gup_rt, gup_rp, gup_tt, gup_tp, gup_pp;
    // induced metric
    gf2d m_rr, m_rt, m_rp, m_tt, m_tp, m_pp;
    // conformal Christoffel symbols
    gf2d Gamma_r_rr, Gamma_r_rt, Gamma_r_rp, Gamma_r_tt, Gamma_r_tp, Gamma_r_pp;
    gf2d Gamma_t_rr, Gamma_t_rt, Gamma_t_rp, Gamma_t_tt, Gamma_t_tp, Gamma_t_pp;
    gf2d Gamma_p_rr, Gamma_p_rt, Gamma_p_rp, Gamma_p_tt, Gamma_p_tp, Gamma_p_pp;
    // extrinsic curvature
    gf2d A_rr, A_rt, A_rp, A_tt, A_tp, A_pp, K;
    // conformal exponent and derivatives
    gf2d phi_c, phi_r, phi_t, phi_p;
    // quantity needed for angular and linear momentum integrals
    gf2d integrand;
    // 
    gf2d lambda, expansion, res, rhs_grid, factor;
    // 
    // pointers to 3D grid functions
    //
    gf3d* gup_rr_3d, * gup_rt_3d, * gup_rp_3d;
    gf3d* gup_tt_3d, * gup_tp_3d, * gup_pp_3d;
    gf3d* DG_r_rr_3d, * DG_r_rt_3d, * DG_r_rp_3d, * DG_r_tt_3d, * DG_r_tp_3d, * DG_r_pp_3d;
    gf3d* DG_t_rr_3d, * DG_t_rt_3d, * DG_t_rp_3d, * DG_t_tt_3d, * DG_t_tp_3d, * DG_t_pp_3d;
    gf3d* DG_p_rr_3d, * DG_p_rt_3d, * DG_p_rp_3d, * DG_p_tt_3d, * DG_p_tp_3d, * DG_p_pp_3d;
    gf3d* phi_3d, * phi_r_3d, * phi_t_3d, * phi_p_3d;
    gf3d* a_rr_3d, * a_rt_3d, * a_rp_3d;
    gf3d* a_tt_3d, * a_tp_3d, * a_pp_3d;
    gf3d* K_3d;
    //
    // functions needed for vorticity and tencidity
    //
    bool EB_assigned;
    gf3d* E_rr_3d, * E_rt_3d, * E_rp_3d, * E_tt_3d, * E_tp_3d, * E_pp_3d;
    gf3d* B_rr_3d, * B_rt_3d, * B_rp_3d, * B_tt_3d, * B_tp_3d, * B_pp_3d;
    gf2d E_rr, E_rt, E_rp, E_tt, E_tp, E_pp;
    gf2d B_rr, B_rt, B_rp, B_tt, B_tp, B_pp;
    gf2d tendicity, vorticity;
    //
    int N_theta, N_phi;   // *include* ghost zones
    int N_g;
    bool foundhorizon;
    int timestep;
    Doub phystime;
#ifndef NoEllSolver
    EllSolver2D* ellsolver;
#endif
    //
    // parameters for solver (read in from file "AH_finder_Input")
    // 
    double tol_exp;    // expansion tolerance
    double tol_ell;    // tolerance for elliptic solver
    int max_iter_ell;  // maximum number of iterations in elliptic solver
    int max_iter_exp;  // maximum number of iterations for expansion
    double eta;      // parameter for expansion equation: \eta h 
    // added to both sides
    int find_steps;   // timesteps or times after which horizon should be found
    double find_times;
    int next_step, last_step;
    double next_time;
    double mass_guess;
    double a_friedmann;      // current value of cosmological expansion coefficient
    bool use_time_step_criterion;  // use step counter to decide whether to 
    // search for horizon; otherwise use time
//  char file_stem[64];
    ofstream monitorfile;
    ofstream surfacefile;
    double PI;
    double r_min, r_max;
    Grid* grid;
public:
    //================================================
    // Constructor
    //================================================
    HorizonFinder(Grid* grid_i, state* s, curvature* c, auxiliary* aux,
        const char* file_stem);
    //================================================
    // Destructor
    //================================================
    ~HorizonFinder() {
#ifndef NoEllSolver
        delete ellsolver;
#endif
        monitorfile.close();
        surfacefile.close();
    };
    //================================================
    // Read input
    //================================================
    int ReadInput();
    //================================================
    // time to find horizon?
    //================================================
    bool TimeToFindHorizon(int timestep, Doub time);
    //================================================
    // search for horizon (return true if horizon found)
    //================================================
    bool FindHorizon(int timestep, Doub time,
        state* s, Fluxes* fluxes,
        Doub adm_mass, Doub mom_guess, Doub a);
    //================================================
    // Horizon diagnostics
    //================================================
    double HorizonMass();
    double Spin();
    double LinearMomentum();
    double AccretionRate(gf3d& lapse, Fluxes* fluxes);
    //================================================
    // Adjust parameters
    //================================================
    double SetFindTimes(double time) { return find_times = time; };
    double SetFindSteps(int step) {
        find_steps = step;
        use_time_step_criterion = (find_steps > 0);
        return find_steps;
    };
    double SetMassGuess(double mass) { return mass_guess = mass; }
    void SetRMaxMin() {
        r_min = (*phi_t_3d).r(N_g);
        const int N_r = phi_t_3d->dim1();
        r_max = (*phi_t_3d).r(N_r - N_g);
    }
private:
    //================================================
    // Compute expansion
    //================================================
    double Expansion();
    //================================================
    // Compute expansion for surface h, return local exansion and normal
    //================================================
public:
    void Expansion(gf2d& h, gf2d& exp_at_h, gf2d& nr, gf2d& nt, gf2d& np);
    //================================================
    // Initial guess for horizon location
    // (use approximate result from Dennison et.al., PRD 74, 064016 (2006), eq. (24))
    //================================================
private:
    void InitialGuess() { InitialGuess(mass_guess); };
    void InitialGuess(double M_init, double P_z = 0);
    //================================================
    // Project 3D functions onto 2D surfaces
    //================================================
    void Project(gf2d& h);
    //================================================
    // Dump horizon
    //================================================
    void DumpHorizon();
    void Note(double time, double mass, double AngMom, double LinMom,
        double eqcir, double polcir,
        double h_eq, double h_pole,
        double tend_eq, double tend_pole, double accretion);

    //================================================
    // Surface integrals
    //================================================
    double SurfaceIntegral();
    double SurfaceIntegral(gf2d& function);
    double SurfaceIntegral2(gf2d& function, int type = proper);
    inline double SurfaceElement(int j, int k, int type = proper);
    double EquatorialCircumference();
    double PolarCircumference() { return PolarCircumference(N_g); };
    double PolarCircumference(int k);
    //================================================
    // Routines for tendicity and vorticity
    //================================================
public:
    bool AssignEB(gf3d* E_rr_p, gf3d* E_rt_p, gf3d* E_rp_p, gf3d* E_tt_p, gf3d* E_tp_p, gf3d* E_pp_p,
        gf3d* B_rr_p, gf3d* B_rt_p, gf3d* B_rp_p, gf3d* B_tt_p, gf3d* B_tp_p, gf3d* B_pp_p);
private:
    bool TendicityVorticity();
    //================================================
    // Utility to decide whether point is inside horizon
    //================================================
public:
    double HorizonLocation(int j, int k);
    bool InsideHorizon(double r, int j, int k);
};

#endif   /* HORIZONFINDER_H */

