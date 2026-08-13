// Tell emacs that this is -*-c++-*- mode
//====================================================
//
// Manages allocation of grid and regridding
//
// Uses uniform grids in x and y; r and theta are given as functions
// of x and y.  These functions are specified in r_fct and theta_fct.
//    
//====================================================
//
#include "nr3.h"
//
#ifndef GRID_H
#define GRID_H
//
class Grid {
public:
    int N_theta, N_phi, N_r, N_ghost;
private:
    double r_max_init, r_max_fin;
    int regrids, regrid_counter;
    int regrid_type, grid_type;
    double tau_star;    // for self-similar regridding
    double last_selfsim_ratio;
    bool selfsimregrid;
    double courant;
    double cutoff;
    double s_param, t_param, t_amp, theta_param;
    double r_focus, x_focus, r_tracker;
    bool tracking;
    double r_max_current, r_max_old, r_max_new;
    double r_max_factor;
    bool printed_warning;
    double PI;
    double dx, dy, dphi;
    VecDoub x_v, y_v;
    VecDoub r_v, theta_v, phi_v;
    VecDoub dxdr_v, dydtheta_v;
    VecDoub ddxdr_v, ddydtheta_v;
    VecDoub costheta_v, sintheta_v;
    double  eps;  // helps with conversion from coordinates to indices works
    int i_ind_sav;  // needed in "hunt" routine for i_ind
    //===================================================
    // Constructor
    //===================================================
public:
    Grid();
    ~Grid();

    //=================================================
    // Return grid specifications
    //=================================================
    int N_r_int() { return N_r; }
    int N_r_tot() { return N_r + 2 * N_ghost; }
    int N_theta_int() { return N_theta; }
    int N_theta_tot() { return N_theta + 2 * N_ghost; }
    int N_phi_int() { return N_phi; }
    int N_phi_tot() { return N_phi + 2 * N_ghost; }
    int N_ghosts() { return N_ghost; }
    double r_max() { return r_max_current; }
    double r_max_final() { return r_max_fin; }
    double r_parameter() { return s_param; }
    double theta_parameter() { return theta_param; }
    double courant_factor() { return courant; }
    double delta_x() { return dx; }
    double delta_y() { return dy; }
    double delta_phi() { return dphi; }
    double delta_r(int i) { return dx / dxdr_v[i]; }
    double delta_theta(int j) { return dy / dydtheta_v[j]; }
    double delta_phi(int k) { return dphi; }
    VecDoub* r() { return &r_v; }
    VecDoub* theta() { return &theta_v; }
    VecDoub* phi() { return &phi_v; }
    VecDoub* dxdr() { return &dxdr_v; }
    VecDoub* ddxdr() { return &ddxdr_v; }
    VecDoub* dydtheta() { return &dydtheta_v; }
    VecDoub* ddydtheta() { return &ddydtheta_v; }
    VecDoub* costheta() { return &costheta_v; }
    VecDoub* sintheta() { return &sintheta_v; }
    double r(int i) { return r_v[i]; }
    double theta(int j) { return theta_v[j]; }
    double phi(int k) { return phi_v[k]; }
    double dxdr(int i) { return dxdr_v[i]; }
    double ddxdr(int i) { return ddxdr_v[i]; }
    double sintheta(int j) { return sintheta_v[j]; }
    double costheta(int j) { return costheta_v[j]; }
    double dydtheta(int j) { return dydtheta_v[j]; }
    double ddydtheta(int j) { return ddydtheta_v[j]; }

    //=================================================
    // Set up grid
    //=================================================
    int Setup_Grid(VecDoub& r, VecDoub& r2, VecDoub& theta,
        VecDoub& sintheta, VecDoub& sin2theta,
        VecDoub& costheta, VecDoub& phi);

    //=================================================
    // regrid
    //=================================================
    bool TimeToRegrid(double criterion);
    int Regrid(VecDoub& r_new);
    double RegridCriterion();

    //=================================================
    // Setup radial grid
    //=================================================
    int Setup_Radial_Grid(VecDoub& r, VecDoub& r2);

    //=================================================
    // Setup angular grid
    //=================================================
    int Setup_Angular_Grid(VecDoub& theta, VecDoub& sintheta,
        VecDoub& sin2theta, VecDoub& costheta,
        VecDoub& phi);

    //================================================
    // Compute r as function of x -- specify function here!
    // Note:
    //    dx / d r = 1 / r'
    //    d^2 x/ d r^2 = - r'' / ( r' )^3
    //
    // Here: r = r_max sinh(r_param x) / sinh(r_param)
    //
    //================================================
    // double r_fct(Doub x, double & x_prime, double & x_dprime) {
    //   double r = 0.0;
    //   if (s_param > 0.) {
    //     const double sinhA = sinh(s_param);
    //     const double sinhAx = sinh(s_param * x);
    //     const double coshAx = cosh(s_param * x);
    //     r = r_max_current * sinhAx / sinhA;
    //     const double r_prime = r_max_current * s_param * coshAx / sinhA;
    //     const double r_dprime = s_param * s_param * r;
    //     x_prime = 1.0 / r_prime;
    //     x_dprime = - r_dprime / (r_prime * r_prime * r_prime); 
    //   } else {
    //     r = r_max_current * x;
    //     x_prime = 1.0 / r_max_current;
    //     x_dprime = 0.0;
    //   }
    //   return r;
    // }  
    //================================================
    // Compute r as function of x -- specify function here!
    // Note:
    //    dx / d r = 1 / r'
    //    d^2 x/ d r^2 = - r'' / ( r' )^3
    //
    // Here: r = r_max / (1 + A) ( A tanh(t_param x) / tanh(t_param)
    //                             + sinh(s_param x) / sinh(s_param) )
    //
    //================================================
    double r_fct_tanh(Doub x, Doub& x_prime, Doub& x_dprime);

    //================================================
    // Compute r as function of x -- specify function here!
    // Note:
    //    dx / d r = 1 / r'
    //    d^2 x/ d r^2 = - r'' / ( r' )^3
    //
    // Here: r = r_max (sinh(sx-sa)+sinh(sa))/(sinh(s-sa)+sinh(sa))
    //
    //  s_param -> flatness
    //  t_param -> r location of flatness 
    //  t_amp   -> x location of flatness (used if t_param=0)
    //
    //================================================
    double r_fct_offset(Doub x, Doub& x_prime, Doub& x_dprime);

    double r_fct(Doub x, Doub& x_prime, Doub& x_dprime);

    //================================================
    // Compute theta as function of y -- specify function here!
    // Note:
    //    dy / d theta = 1 / theta'
    //    d^2 y/ d theta^2 = - theta'' / ( theta' )^3
    //================================================
    double theta_fct(Doub y, Doub& y_prime, Doub& y_dprime);
    
    //
    //==================================================
    // finally routines that find indices for coordinates, eg., return
    //    index i with r[i] <= rl < r[i+1]
    //==================================================
    //
public:
    // int i_ind(const double rl) {
    //   if (s_param > 0.0) {
    //     // first find x
    //     const double temp = sinh(s_param) * rl / r_max_current;
    //     const double xl = log(temp + sqrt(temp*temp + 1.0)) / s_param;
    //     const double x_ind = xl + (N_ghost - 0.5)*dx;
    //     return int(x_ind/dx + eps);
    //   } else {
    //     const double xl = rl / r_max_current;
    //     const double x_ind = xl + (N_ghost - 0.5)*dx;
    //     return int(x_ind/dx + eps);
    //   }
    // }

    int i_ind(const double rl);
    int j_ind(const double thetal);
    int k_ind(double phil);

};

#endif /* GRID_H */
