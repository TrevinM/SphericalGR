// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing types for grid functions 
//
// Based on Mat3DDoub in numerical recipes, but includes
// functionality for spherical coordinate systems.
//
// Uses Delta_r_i = r_i - r_{i-1}
//================================================
#ifndef GF_H
#define GF_H

#define TEST

#include "nr3.h"
#include "Grid.h"

class gf3d {
private:
    int nr, nt, np;    // NOTE: these *include* ghost zones
    int N_g;  // number of ghost zones - hardcoded to three
    Doub*** v;
    VecDoub* r_p, * dxdr_p, * ddxdr_p, * theta_p, * ct_p, * dydtheta_p, * ddydtheta_p, * phi_p;
    VecDoub* st_p;
    // NOTE: r and phi are functions of uniform grids in x and y; see Grid.h
    double delta_x, delta_y, delta_phi;
    // Avoid divisions in derivatives...
    double Inv12dx, Inv12dy, Inv12dphi, Inv12dxdx, Inv12dydy, Inv12dphidphi;
    double Inv6dx, Inv6dy, Inv6dphi;
    double Inv60dx, Inv60dy;
    double Inv180dxdx, Inv180dydy;
    double Inv840dx, Inv840dy;
    double Inv5040dxdx, Inv5040dydy;
    int gfn; // optional gridfunction number
    int center_par, axis_par, eq_par;   // parity, set to 1 by default
    const char* name;
    int fall_off;
    double wave_speed;  // characteristic speed at outer boundary
    double background;  // background for outer boundary
    double PI;
    Grid* grid;
    int order;    // order of interpolations for operators ()
    double KO_factor;    // see arXiv:2104.06978, eq. (7)
    int KO_order;     // 
public:
    //================================================
    // Constructor
    //================================================
    gf3d();

    //================================================
    // Set up grid functions 
    //================================================
    //
    // Version that allows to specify function name, number, parity, wave speed and background 
    // 
    gf3d* setup(Grid* grid_i, int fall_off_i, const char* name_i,
        int gfn_i, int center_i, int axis_i, int eq_i, double wave_speed_i, double background_i);
    //
    // Version that allows to specify function name, number and parity and wave speed
    // 
    gf3d* setup(Grid* grid_i, int fall_off_i, const char* name_i,
        int gfn_i, int center_i, int axis_i, int eq_i, double wave_speed_i);

    //
    // Version that allows to specify function name, number and parity
    // 
    gf3d* setup(Grid* grid_i, int fall_off_i, const char* name_i,
        int gfn_i, int center_i, int axis_i, int eq_i);

    //
    // Version that allows to specify both function name and number
    // 
    gf3d* setup(Grid* grid_i, int fall_off_i, const char* name_i, int gfn_i);

    //
    // Version that allows to specify function number
    // 
    gf3d* setup(Grid* grid_i, int fall_off_i, int gfn_i);

    //
    // Minimum version
    //
    gf3d* setup(Grid* grid_i, int fall_off_i);

    //================================================
    // access memory...
    //================================================
    inline Doub** operator[](const int i) { return v[i]; } //subscripting: pointer to row i
    inline const Doub* const* operator[](const int i) const { return v[i]; }
    inline double operator()(int i, int j, int k) { return v[i][j][k]; }

    inline int dim1() const { return nr; }
    inline int dim2() const { return nt; }
    inline int dim3() const { return np; }
    inline int N_ghosts() { return N_g; }
    inline double r(int i) { return (*r_p)[i]; }
    inline double dxdr(int i) { return (*dxdr_p)[i]; }
    inline double ddxdr(int i) { return (*ddxdr_p)[i]; }
    inline double theta(int j) { return (*theta_p)[j]; }
    inline double dydtheta(int j) { return (*dydtheta_p)[j]; }
    inline double ddydtheta(int j) { return (*ddydtheta_p)[j]; }
    inline double costheta(int j) { return (*ct_p)[j]; }
    inline double sintheta(int j) { return (*st_p)[j]; }
    inline double phi(int k) { return (*phi_p)[k]; }
    inline void constants() { cout << " constant in " << name << " : " << Inv12dy << endl; }
    inline double dx() const { return delta_x; }
    inline double dy() const { return delta_y; }
    inline double dphi() const { return delta_phi; }
    inline int GridFunctionNumber() { return gfn; }
    inline int Center_Parity() const { return center_par; }
    inline int Axis_Parity() const { return axis_par; }
    inline int Eq_Parity() const { return eq_par; }
    inline int Fall_Off() const { return fall_off; }
    inline double WaveSpeed() const { return wave_speed; }
    inline double Background() const { return background; }
    inline int NumberGhostZones() { return N_g; }
    const char* Name() { return name; }
    inline gf3d* Address() { return this; }

    //================================================
    // Fill inner ghost zones, assuming a spherical grid.
    //================================================
    int fill_ghosts();

    //================================================
    // Compute time derivatives at outer boundary (provide function values)
    //================================================
    void derivs_outerboundary(gf3d* fct);

    //================================================
    // Fill outer boundary points
    //================================================
    int fill_outerboundary(gf3d& fct_old, double dt, double update_background);
    int fill_outerboundary(gf3d& fct_old, double dt);

    //    
    //================================================
    //================================================
    // Derivative operators 
    //================================================
    //================================================
    //
    // flat Laplace operator
    //
    double Laplace(int i, int j, int k);

    // flat Laplace operator to second order
    //
    double Laplace_so(int i, int j, int k);

    //
    // first derivatives
    // 
    double dr(int i, int j, int k);
    double dtheta(int i, int j, int k);
    double dphi(int i, int j, int k);

    //
    // second derivatives
    //
    double ddr(int i, int j, int k);
    double ddtheta(int i, int j, int k);
    double ddphi(int i, int j, int k);

    //
    // mixed second derivatives
    //  
    double dthetadr(int i, int j, int k);
    double drdtheta(int i, int j, int k);
    double dphidr(int i, int j, int k);
    double drdphi(int i, int j, int k);
    double dphidtheta(int i, int j, int k);
    double dthetadphi(int i, int j, int k);
    
    //
    // second-order versions of second derivatives
    //
    double ddr_so(int i, int j, int k);
    double ddtheta_so(int i, int j, int k);
    double ddphi_so(int i, int j, int k);

    //
    // second-order versions of first derivatives
    //
    double dr_so(int i, int j, int k);
    double dtheta_so(int i, int j, int k);
    double dphi_so(int i, int j, int k);

    //
    // second-order versions of mixed derivatives
    //
    double dthetadr_so(int i, int j, int k);
    double drdtheta_so(int i, int j, int k);
    double dphidr_so(int i, int j, int k);
    double drdphi_so(int i, int j, int k);
    double dphidtheta_so(int i, int j, int k);
    double dthetadphi_so(int i, int j, int k);
    
    //
    // first derivatives in UPWIND differencing (NOTE: same name for routine, but extra argument for shift) 
    // 
    double dr(int i, int j, int k, double shift);
    double dtheta(int i, int j, int k, double shift);
    double dphi(int i, int j, int k, double shift);

    //
    // first derivatives in UPWIND differencing (NOTE: same name for routine, but extra argument for shift) 
    // Now third-order versions
    //
    double dr_to(int i, int j, int k, double shift);
    double dtheta_to(int i, int j, int k, double shift);
    double dphi_to(int i, int j, int k, double shift);

    //
    // One-sided version (for use at boundaries)
    //
    double dr_OS(int i, int j, int k, double shift);
    //
    // NOTE: these return derivatives * (dx)^3 !
    //
    double d4r(int i, int j, int k);
    double d4theta(int i, int j, int k);
    double d4phi(int i, int j, int k);
    double d4(int i, int j, int k) ;

    //
    // Derivatives for Kreiss-Oliger terms - for a *uniform* grid,
    // these return (Delta x)^p \partial_x^{p+1} f 
    //
    double D8r(int i, int j, int k);
    double D8theta(int i, int j, int k);
    double D8phi(int i, int j, int k);

    double D6r(int i, int j, int k);
    double D6theta(int i, int j, int k);
    double D6phi(int i, int j, int k);
    double KO(int i, int j, int k);

    //================================================
    // a hand-full of global derivative operators...
    //================================================
    // take radial derivative of function and store in derivs
    //================================================
    int dr(gf3d* derivs);
    int dtheta(gf3d* derivs);

    //================================================
    // utilities...
    //================================================
    void add(double factor, gf3d* rhs);
    void add(gf3d* term1, double factor, gf3d* term2);
    double equals(double value);
    void equals(gf3d* rhs);
    bool IsEqualTo(gf3d* rhs);
    double min();
    double min(int& i_min, int& j_min, int& k_min);
    double max();
    double max(int& i_max, int& j_max, int& k_max);
    double abs_max();
    double abs_max(int& i_max, int& j_max, int& k_max);
    double max(double r_min);

    //
    // Use polynomial interpolation to find maximum in direction theta or j, fixed value of k (for phi)
    // 
    double radialmax(double theta, int k, double& r_max);
    double radialmax(int j, int k, double& r_max);
    //
    // Use polynomial interpolation to find maximum in slice for fixed value of k (for phi)
    // 
    double slicemax(int k);
    double slicemax(int k, double& r_max, double& theta_max);
    //
    // generic routine that returns maximum of polynomial fit through three points
    // 
    double max(double x_0, double x_1, double x_2, double f_0, double f_1, double f_2, double& x_max);
    //
    // Use polynomial interpolation to find minimum in direction theta, fixed value of k (for phi)
    // 
    double radialmin(double theta, int k, double& r_min);
    bool FINITE();
    double center(int j, int k);
    //
    // takes norm only inside r_frac * r_out
    // 
    double L2_norm(double r_frac = 1.1);
    double L2_norm(gf3d& det);

    //================================================
    // proper integrals
    //================================================
    double PropInt(gf3d& phi, gf3d& det, double low_lim_r = 0.0);
    double PropIntMag(gf3d& phi, gf3d& det, double low_lim_r = 0.0);

    //================================================
    // Interpolator to arbitrary value of rl > 0
    //================================================
    double operator()(double rl, int j, int k, int ilo = -1);

    //================================================
    // Interpolator to arbitrary value of theta >= 0
    //================================================
    double operator()(int i, double theta, int k);

    //================================================
    // Interpolator to arbitrary value of r > 0 AND theta >= 0
    //================================================
    double operator()(double rl, double theta, int k,
        int loc_order = 0, int ilo = -1);

    //================================================
    // Interpolator to arbitrary value of phi >= 0
    //================================================
    //
    double operator()(int i, int j, double phi);

    //================================================
    // surface integral
    //================================================
    double surface_integral(int i);

    //================================================
    // Write check-point file
    // NOTE: This logic has to match that in INDATA/ReadFromCheckPoint.setup_gf
    //================================================
    bool checkpoint(int timestep);

    //================================================
    // Regrid
    //================================================
    int Regrid(VecDoub r_new);

    //================================================
    // Destructor
    //================================================
    ~gf3d();

};

#endif  /* GF_H */











//     //    
//     //================================================
//     //================================================
//     // OLD INLINE Derivative operators 
//     //================================================
//     //================================================

//     //
//     // flat Laplace operator
//     //
//     inline double Laplace(int i, int j, int k) {
//         return ddr(i, j, k) + 2.0 * dr(i, j, k) / r(i)
//             + (ddtheta(i, j, k) + (*ct_p)[j] * dtheta(i, j, k) / (*st_p)[j]
//                 + ddphi(i, j, k) / ((*st_p)[j] * (*st_p)[j])) / (r(i) * r(i));
//     }

//     //
//     // flat Laplace operator to second order
//     //
//     inline double Laplace_so(int i, int j, int k) {
//         return ddr_so(i, j, k) + 2.0 * dr_so(i, j, k) / r(i)
//             + (ddtheta_so(i, j, k) + (*ct_p)[j] * dtheta_so(i, j, k) / (*st_p)[j]
//                 + ddphi_so(i, j, k) / ((*st_p)[j] * (*st_p)[j])) / (r(i) * r(i));
//     }

//     //
//     // first derivatives
//     // 
//     inline double dr(int i, int j, int k) {
//     #ifdef EIGHTHORDER
//         return Inv840dx * (-3.0 * (v[i + 4][j][k] - v[i - 4][j][k])
//             + 32.0 * (v[i + 3][j][k] - v[i - 3][j][k])
//             - 168.0 * (v[i + 2][j][k] - v[i - 2][j][k])
//             + 672.0 * (v[i + 1][j][k] - v[i - 1][j][k])) * dxdr(i);
//     #elif SIXTHORDER
//         return Inv60dx * (v[i + 3][j][k] - v[i - 3][j][k] -
//             9.0 * (v[i + 2][j][k] - v[i - 2][j][k]) +
//             45.0 * (v[i + 1][j][k] - v[i - 1][j][k])) * dxdr(i);
//     #else /* fourth-order by default... */
//         return (v[i - 2][j][k] - 8.0 * (v[i - 1][j][k] - v[i + 1][j][k]) - v[i + 2][j][k]) * Inv12dx * dxdr(i);
//     #endif  /* ORDER... */
//     }

//     inline double dtheta(int i, int j, int k) {
//     #ifdef EIGHTHORDER
//         return Inv840dy * (-3.0 * (v[i][j + 4][k] - v[i][j - 4][k])
//             + 32.0 * (v[i][j + 3][k] - v[i][j - 3][k])
//             - 168.0 * (v[i][j + 2][k] - v[i][j - 2][k])
//             + 672.0 * (v[i][j + 1][k] - v[i][j - 1][k])) * dydtheta(j);
//     #elif SIXTHORDER
//         return Inv60dy * (v[i][j + 3][k] - v[i][j - 3][k]
//             - 9.0 * (v[i][j + 2][k] - v[i][j - 2][k])
//             + 45.0 * (v[i][j + 1][k] - v[i][j - 1][k])) * dydtheta(j);
//     #else
//         return (v[i][j - 2][k] - 8.0 * (v[i][j - 1][k] - v[i][j + 1][k]) - v[i][j + 2][k]) * Inv12dy * dydtheta(j);
//     #endif  /* ORDER... */
//     }

//     inline double dphi(int i, int j, int k) {
//     #ifdef AXISYMMETRY
//         return 0.0;
//     #else
//         return (v[i][j][k - 2] - 8.0 * (v[i][j][k - 1] - v[i][j][k + 1]) - v[i][j][k + 2]) * Inv12dphi;
//     #endif  /* AXISYMMETRY */
//     }

//     //
//     // second derivatives
//     //
//     inline double ddr(int i, int j, int k) {
//     #ifdef EIGHTHORDER
//         return Inv5040dxdx * (-9.0 * (v[i + 4][j][k] + v[i - 4][j][k])
//             + 128.0 * (v[i + 3][j][k] + v[i - 3][j][k])
//             - 1008.0 * (v[i + 2][j][k] + v[i - 2][j][k])
//             + 8064.0 * (v[i + 1][j][k] + v[i - 1][j][k])
//             - 14350.0 * v[i][j][k]) * dxdr(i) * dxdr(i) +
//             Inv840dx * (-3.0 * (v[i + 4][j][k] - v[i - 4][j][k])
//                 + 32.0 * (v[i + 3][j][k] - v[i - 3][j][k])
//                 - 168.0 * (v[i + 2][j][k] - v[i - 2][j][k])
//                 + 672.0 * (v[i + 1][j][k] - v[i - 1][j][k])) * ddxdr(i);

//     #elif SIXTHORDER
//         return Inv180dxdx * (2.0 * (v[i + 3][j][k] + v[i - 3][j][k]) -
//             27.0 * (v[i + 2][j][k] + v[i - 2][j][k]) +
//             270.0 * (v[i + 1][j][k] + v[i - 1][j][k]) -
//             490.0 * v[i][j][k]) * dxdr(i) * dxdr(i) +
//             Inv60dx * (v[i + 3][j][k] - v[i - 3][j][k] -
//                 9.0 * (v[i + 2][j][k] - v[i - 2][j][k]) +
//                 45.0 * (v[i + 1][j][k] - v[i - 1][j][k])) * ddxdr(i);
//     #else
//         return (-(v[i - 2][j][k] + v[i + 2][j][k]) - 30.0 * v[i][j][k] + 16.0 * (v[i + 1][j][k] + v[i - 1][j][k])) *
//             Inv12dxdx * (dxdr(i) * dxdr(i)) +
//             (v[i - 2][j][k] - 8.0 * (v[i - 1][j][k] - v[i + 1][j][k]) - v[i + 2][j][k]) * Inv12dx * ddxdr(i);
//     #endif   /* ORDER... */
//     }

//     inline double ddtheta(int i, int j, int k) {
//     #ifdef EIGHTHORDER
//         return Inv5040dydy * (-9.0 * (v[i][j + 4][k] + v[i][j - 4][k])
//             + 128.0 * (v[i][j + 3][k] + v[i][j - 3][k])
//             - 1008.0 * (v[i][j + 2][k] + v[i][j - 2][k])
//             + 8064.0 * (v[i][j + 1][k] + v[i][j - 1][k])
//             - 14350.0 * v[i][j][k]) * dydtheta(j) * dydtheta(j) +
//             Inv840dy * (-3.0 * (v[i][j + 4][k] - v[i][j - 4][k])
//                 + 32.0 * (v[i][j + 3][k] - v[i][j - 3][k])
//                 - 168.0 * (v[i][j + 2][k] - v[i][j - 2][k])
//                 + 672.0 * (v[i][j + 1][k] - v[i][j - 1][k])) * ddydtheta(j);
//     #elif SIXTHORDER
//         return Inv180dydy * (2.0 * (v[i][j + 3][k] + v[i][j - 3][k]) -
//             27.0 * (v[i][j + 2][k] + v[i][j - 2][k]) +
//             270.0 * (v[i][j + 1][k] + v[i][j - 1][k]) -
//             490.0 * v[i][j][k]) * dydtheta(j) * dydtheta(j) +
//             Inv60dy * (v[i][j + 3][k] - v[i][j - 3][k] -
//                 9.0 * (v[i][j + 2][k] - v[i][j - 2][k]) +
//                 45.0 * (v[i][j + 1][k] - v[i][j - 1][k])) * ddydtheta(j);
//     #else
//         return (-(v[i][j - 2][k] + v[i][j + 2][k]) - 30.0 * v[i][j][k] + 16.0 * (v[i][j + 1][k] + v[i][j - 1][k])) *
//             Inv12dydy * (dydtheta(j) * dydtheta(j)) +
//             (v[i][j - 2][k] - 8.0 * (v[i][j - 1][k] - v[i][j + 1][k]) - v[i][j + 2][k]) * Inv12dy * ddydtheta(j);
//     #endif  /* ORDER */
//     }

//     inline double ddphi(int i, int j, int k) {
//     #ifdef AXISYMMETRY
//         return 0.0;
//     #else
//         return (-(v[i][j][k - 2] + v[i][j][k + 2]) - 30.0 * v[i][j][k] + 16.0 * (v[i][j][k + 1] + v[i][j][k - 1])) *
//             Inv12dphidphi;
//     #endif  /* AXISYMMETRY */
//     }

//     //
//     // mixed second derivatives
//     //  
//     inline double dthetadr(int i, int j, int k) { return drdtheta(i, j, k); }
//     inline double drdtheta(int i, int j, int k) {
//     #ifdef EIGHTHORDER
//         return (
//             -3.0 * (-3.0 * (v[i + 4][j + 4][k] - v[i + 4][j - 4][k]) + 32.0 * (v[i + 4][j + 3][k] - v[i + 4][j - 3][k]) - 168.0 * (v[i + 4][j + 2][k] - v[i + 4][j - 2][k]) + 672.0 * (v[i + 4][j + 1][k] - v[i + 4][j - 1][k]))
//             + 32.0 * (-3.0 * (v[i + 3][j + 4][k] - v[i + 3][j - 4][k]) + 32.0 * (v[i + 3][j + 3][k] - v[i + 3][j - 3][k]) - 168.0 * (v[i + 3][j + 2][k] - v[i + 3][j - 2][k]) + 672.0 * (v[i + 3][j + 1][k] - v[i + 3][j - 1][k]))
//             - 168.0 * (-3.0 * (v[i + 2][j + 4][k] - v[i + 2][j - 4][k]) + 32.0 * (v[i + 2][j + 3][k] - v[i + 2][j - 3][k]) - 168.0 * (v[i + 2][j + 2][k] - v[i + 2][j - 2][k]) + 672.0 * (v[i + 2][j + 1][k] - v[i + 2][j - 1][k]))
//             + 672.0 * (-3.0 * (v[i + 1][j + 4][k] - v[i + 1][j - 4][k]) + 32.0 * (v[i + 1][j + 3][k] - v[i + 1][j - 3][k]) - 168.0 * (v[i + 1][j + 2][k] - v[i + 1][j - 2][k]) + 672.0 * (v[i + 1][j + 1][k] - v[i + 1][j - 1][k]))
//             - 672.0 * (-3.0 * (v[i - 1][j + 4][k] - v[i - 1][j - 4][k]) + 32.0 * (v[i - 1][j + 3][k] - v[i - 1][j - 3][k]) - 168.0 * (v[i - 1][j + 2][k] - v[i - 1][j - 2][k]) + 672.0 * (v[i - 1][j + 1][k] - v[i - 1][j - 1][k]))
//             + 168.0 * (-3.0 * (v[i - 2][j + 4][k] - v[i - 2][j - 4][k]) + 32.0 * (v[i - 2][j + 3][k] - v[i - 2][j - 3][k]) - 168.0 * (v[i - 2][j + 2][k] - v[i - 2][j - 2][k]) + 672.0 * (v[i - 2][j + 1][k] - v[i - 2][j - 1][k]))
//             - 32.0 * (-3.0 * (v[i - 3][j + 4][k] - v[i - 3][j - 4][k]) + 32.0 * (v[i - 3][j + 3][k] - v[i - 3][j - 3][k]) - 168.0 * (v[i - 3][j + 2][k] - v[i - 3][j - 2][k]) + 672.0 * (v[i - 3][j + 1][k] - v[i - 3][j - 1][k]))
//             + 3.0 * (-3.0 * (v[i - 4][j + 4][k] - v[i - 4][j - 4][k]) + 32.0 * (v[i - 4][j + 3][k] - v[i - 4][j - 3][k]) - 168.0 * (v[i - 4][j + 2][k] - v[i - 4][j - 2][k]) + 672.0 * (v[i - 4][j + 1][k] - v[i - 4][j - 1][k]))
//             ) *
//             Inv840dx * Inv840dy * dxdr(i) * dydtheta(j);
//     #elif SIXTHORDER
//         return (
//             (v[i + 3][j + 3][k] - v[i + 3][j - 3][k] - 9.0 * (v[i + 3][j + 2][k] - v[i + 3][j - 2][k]) + 45.0 * (v[i + 3][j + 1][k] - v[i + 3][j - 1][k]))
//             - 9.0 * (v[i + 2][j + 3][k] - v[i + 2][j - 3][k] - 9.0 * (v[i + 2][j + 2][k] - v[i + 2][j - 2][k]) + 45.0 * (v[i + 2][j + 1][k] - v[i + 2][j - 1][k]))
//             + 45.0 * (v[i + 1][j + 3][k] - v[i + 1][j - 3][k] - 9.0 * (v[i + 1][j + 2][k] - v[i + 1][j - 2][k]) + 45.0 * (v[i + 1][j + 1][k] - v[i + 1][j - 1][k]))
//             - 45.0 * (v[i - 1][j + 3][k] - v[i - 1][j - 3][k] - 9.0 * (v[i - 1][j + 2][k] - v[i - 1][j - 2][k]) + 45.0 * (v[i - 1][j + 1][k] - v[i - 1][j - 1][k]))
//             + 9.0 * (v[i - 2][j + 3][k] - v[i - 2][j - 3][k] - 9.0 * (v[i - 2][j + 2][k] - v[i - 2][j - 2][k]) + 45.0 * (v[i - 2][j + 1][k] - v[i - 2][j - 1][k]))
//             - (v[i - 3][j + 3][k] - v[i - 3][j - 3][k] - 9.0 * (v[i - 3][j + 2][k] - v[i - 3][j - 2][k]) + 45.0 * (v[i - 3][j + 1][k] - v[i - 3][j - 1][k]))) *
//             Inv60dx * Inv60dy * dxdr(i) * dydtheta(j);
//     #else
//         return (
//             (v[i - 2][j - 2][k] - 8.0 * (v[i - 2][j - 1][k] - v[i - 2][j + 1][k]) - v[i - 2][j + 2][k])
//             - 8.0 * (v[i - 1][j - 2][k] - 8.0 * (v[i - 1][j - 1][k] - v[i - 1][j + 1][k]) - v[i - 1][j + 2][k])
//             + 8.0 * (v[i + 1][j - 2][k] - 8.0 * (v[i + 1][j - 1][k] - v[i + 1][j + 1][k]) - v[i + 1][j + 2][k])
//             - (v[i + 2][j - 2][k] - 8.0 * (v[i + 2][j - 1][k] - v[i + 2][j + 1][k]) - v[i + 2][j + 2][k])) *
//             Inv12dx * Inv12dy * dxdr(i) * dydtheta(j);
//     #endif
//     }

//     inline double dphidr(int i, int j, int k) { return drdphi(i, j, k); }
//     inline double drdphi(int i, int j, int k) {
//     #ifdef AXISYMMETRY
//         return 0.0;
//     #else
//         return ((v[i - 2][j][k - 2] - 8.0 * (v[i - 2][j][k - 1] - v[i - 2][j][k + 1]) - v[i - 2][j][k + 2])
//             - 8.0 * (v[i - 1][j][k - 2] - 8.0 * (v[i - 1][j][k - 1] - v[i - 1][j][k + 1]) - v[i - 1][j][k + 2])
//             + 8.0 * (v[i + 1][j][k - 2] - 8.0 * (v[i + 1][j][k - 1] - v[i + 1][j][k + 1]) - v[i + 1][j][k + 2])
//             - (v[i + 2][j][k - 2] - 8.0 * (v[i + 2][j][k - 1] - v[i + 2][j][k + 1]) - v[i + 2][j][k + 2])) *
//             Inv12dx * Inv12dphi * dxdr(i);
//     #endif  /* AXISYMMETRY */
//     }
//     inline double dphidtheta(int i, int j, int k) { return dthetadphi(i, j, k); }
//     inline double dthetadphi(int i, int j, int k) {
//     #ifdef AXISYMMETRY
//         return 0.0;
//     #else
//         return ((v[i][j - 2][k - 2] - 8.0 * (v[i][j - 2][k - 1] - v[i][j - 2][k + 1]) - v[i][j - 2][k + 2])
//             - 8.0 * (v[i][j - 1][k - 2] - 8.0 * (v[i][j - 1][k - 1] - v[i][j - 1][k + 1]) - v[i][j - 1][k + 2])
//             + 8.0 * (v[i][j + 1][k - 2] - 8.0 * (v[i][j + 1][k - 1] - v[i][j + 1][k + 1]) - v[i][j + 1][k + 2])
//             - (v[i][j + 2][k - 2] - 8.0 * (v[i][j + 2][k - 1] - v[i][j + 2][k + 1]) - v[i][j + 2][k + 2])) *
//             Inv12dy * Inv12dphi * dydtheta(j);
//     #endif  /* AXISYMMETRY */
//     }

//     //
//     // second-order versions of second derivatives
//     //
//     inline double ddr_so(int i, int j, int k) {
//         return (-2.0 * v[i][j][k] + (v[i + 1][j][k] + v[i - 1][j][k])) / (delta_x * delta_x)
//             * (dxdr(i) * dxdr(i)) +
//             (v[i + 1][j][k] - v[i - 1][j][k]) / (2.0 * delta_x) * ddxdr(i);
//     }

//     inline double ddtheta_so(int i, int j, int k) {
//         return (-2.0 * v[i][j][k] + (v[i][j + 1][k] + v[i][j - 1][k])) / (delta_y * delta_y)
//             * (dydtheta(j) * dydtheta(j)) +
//             (v[i][j + 1][k] - v[i][j - 1][k]) / (2.0 * delta_y) * ddydtheta(j);
//     }

//     inline double ddphi_so(int i, int j, int k) {
//     #ifdef AXISYMMETRY
//         return 0.0;
//     #else
//         return (-2.0 * v[i][j][k] + (v[i][j][k + 1] + v[i][j][k - 1])) / (delta_phi * delta_phi);
//     #endif  /* AXISYMMETRY */
//     }

//     //
//     // second-order versions of first derivatives
//     //
//     inline double dr_so(int i, int j, int k) {
//         return (v[i + 1][j][k] - v[i - 1][j][k]) / (2.0 * delta_x) * dxdr(i);
//     }
//     inline double dtheta_so(int i, int j, int k) {
//         return (v[i][j + 1][k] - v[i][j - 1][k]) / (2.0 * delta_y) * dydtheta(j);
//     }
//     inline double dphi_so(int i, int j, int k) {
//     #ifdef AXISYMMETRY
//         return 0.0;
//     #else
//         return (v[i][j][k + 1] - v[i][j][k - 1]) / (2.0 * delta_phi);
//     #endif
//     }

//     inline double dthetadr_so(int i, int j, int k) { return drdtheta_so(i, j, k); }
//     inline double drdtheta_so(int i, int j, int k) {
//         return (v[i + 1][j + 1][k] - v[i + 1][j - 1][k] - v[i - 1][j + 1][k] + v[i - 1][j - 1][k])
//             / (4.0 * delta_x * delta_y) * dxdr(i) * dydtheta(j);
//     }
//     inline double dphidr_so(int i, int j, int k) { return drdphi_so(i, j, k); }
//     inline double drdphi_so(int i, int j, int k) {
// #ifdef AXISYMMETRY
//         return 0.0;
// #else
//         return (v[i - 1][j][k - 1] - v[i - 1][j][k + 1] - v[i + 1][j][k - 1] + v[i + 1][j][k + 1])
//             / (4.0 * delta_x * delta_phi) * dxdr(i);
// #endif
//     }
//     inline double dphidtheta_so(int i, int j, int k) { return dthetadphi_so(i, j, k); }
//     inline double dthetadphi_so(int i, int j, int k) {
// #ifdef AXISYMMETRY
//         return 0.0;
// #else
//         return (v[i][j - 1][k - 1] - v[i][j - 1][k + 1] - v[i][j + 1][k - 1] + v[i][j + 1][k + 1])
//             / (4.0 * delta_y * delta_phi) * dydtheta(j);
// #endif
//     }
//     //
//     // first derivatives in UPWIND differencing (NOTE: same name for routine, but extra argument for shift) 
//     // 
//     inline double dr(int i, int j, int k, double shift) {
//         if (shift < 0.0) {
// #if defined SIXTHORDER || defined EIGHTHORDER
//             return (v[i - 4][j][k] - 8.0 * v[i - 3][j][k] + 30.0 * v[i - 2][j][k] - 80.0 * v[i - 1][j][k] + 35.0 * v[i][j][k] + 24.0 * v[i + 1][j][k] - 2.0 * v[i + 2][j][k])
//                 * Inv60dx * dxdr(i);
// #else
//             return   (-0.5 * v[i - 3][j][k] + 3.0 * v[i - 2][j][k] - 9.0 * v[i - 1][j][k] + 5.0 * v[i][j][k] + 1.5 * v[i + 1][j][k])
//                 * Inv6dx * dxdr(i);
// #endif  /* ORDER */
//         }
//         else {   // now shift > 0...
// #if defined SIXTHORDER || defined EIGHTHORDER
//             return -(v[i + 4][j][k] - 8.0 * v[i + 3][j][k] + 30.0 * v[i + 2][j][k] - 80.0 * v[i + 1][j][k] + 35.0 * v[i][j][k] + 24.0 * v[i - 1][j][k] - 2.0 * v[i - 2][j][k])
//                 * Inv60dx * dxdr(i);
// #else
//             return -(-0.5 * v[i + 3][j][k] + 3.0 * v[i + 2][j][k] - 9.0 * v[i + 1][j][k] + 5.0 * v[i][j][k] + 1.5 * v[i - 1][j][k])
//                 * Inv6dx * dxdr(i);
// #endif  /* ORDER */
//         }
//     }
//     inline double dtheta(int i, int j, int k, double shift) {
//         if (shift < 0.0) {
// #if defined SIXTHORDER || defined EIGHTHORDER
//             return (v[i][j - 4][k] - 8.0 * v[i][j - 3][k] + 30.0 * v[i][j - 2][k] - 80.0 * v[i][j - 1][k] + 35.0 * v[i][j][k] + 24.0 * v[i][j + 1][k] - 2.0 * v[i][j + 2][k])
//                 * Inv60dy * dydtheta(j);
// #else
//             return   (-0.5 * v[i][j - 3][k] + 3.0 * v[i][j - 2][k] - 9.0 * v[i][j - 1][k] + 5.0 * v[i][j][k] + 1.5 * v[i][j + 1][k])
//                 * Inv6dy * dydtheta(j);
// #endif  /* ORDER */
//         }
//         else { // now shift > 0...
// #if defined SIXTHORDER || defined EIGHTHORDER
//             return -(v[i][j + 4][k] - 8.0 * v[i][j + 3][k] + 30.0 * v[i][j + 2][k] - 80.0 * v[i][j + 1][k] + 35.0 * v[i][j][k] + 24.0 * v[i][j - 1][k] - 2.0 * v[i][j - 2][k])
//                 * Inv60dy * dydtheta(j);
// #else
//             return -(-0.5 * v[i][j + 3][k] + 3.0 * v[i][j + 2][k] - 9.0 * v[i][j + 1][k] + 5.0 * v[i][j][k] + 1.5 * v[i][j - 1][k])
//                 * Inv6dy * dydtheta(j);
// #endif  /* ORDER */
//         }
//     }
//     inline double dphi(int i, int j, int k, double shift) {
// #ifdef AXISYMMETRY
//         return 0.0;
// #else
//         if (shift < 0.0)
//             return   (-0.5 * v[i][j][k - 3] + 3.0 * v[i][j][k - 2] - 9.0 * v[i][j][k - 1] + 5.0 * v[i][j][k] + 1.5 * v[i][j][k + 1])
//             * Inv6dphi;
//         else
//             return -(-0.5 * v[i][j][k + 3] + 3.0 * v[i][j][k + 2] - 9.0 * v[i][j][k + 1] + 5.0 * v[i][j][k] + 1.5 * v[i][j][k - 1])
//             * Inv6dphi;
// #endif  /* AXISYMMETRY */
//     }
//     //
//     // first derivatives in UPWIND differencing (NOTE: same name for routine, but extra argument for shift) 
//     // Now third-order versions
//     //
//     inline double dr_to(int i, int j, int k, double shift) {
//         if (shift < 0.0)
//             return   (v[i - 2][j][k] - 6.0 * v[i - 1][j][k] + 3.0 * v[i][j][k] + 2.0 * v[i + 1][j][k])
//             * Inv6dx * dxdr(i);
//         else
//             return -(v[i + 2][j][k] - 6.0 * v[i + 1][j][k] + 3.0 * v[i][j][k] + 2.0 * v[i - 1][j][k])
//             * Inv6dx * dxdr(i);
//     }
//     inline double dtheta_to(int i, int j, int k, double shift) {
//         if (shift < 0.0)
//             return   (v[i][j - 2][k] - 6.0 * v[i][j - 1][k] + 3.0 * v[i][j][k] + 2.0 * v[i][j + 1][k])
//             * Inv6dy * dydtheta(j);
//         else
//             return -(v[i][j + 2][k] - 6.0 * v[i][j + 1][k] + 3.0 * v[i][j][k] + 2.0 * v[i][j - 1][k])
//             * Inv6dy * dydtheta(j);
//     }
//     inline double dphi_to(int i, int j, int k, double shift) {
// #ifdef AXISYMMETRY
//         return 0.0;
// #else
//         if (shift < 0.0)
//             return   (v[i][j][k - 2] - 6.0 * v[i][j][k - 1] + 3.0 * v[i][j][k] + 2.0 * v[i][j][k + 1])
//             * Inv6dphi;
//         else
//             return -(v[i][j][k + 2] - 6.0 * v[i][j][k + 1] + 3.0 * v[i][j][k] + 2.0 * v[i][j][k - 1])
//             * Inv6dphi;
// #endif  /* AXISYMMETRY */
//     }
//     //
//     // One-sided version (for use at boundaries)
//     //
//     inline double dr_OS(int i, int j, int k, double shift) {
//         // if (shift < 0.0) 
//         //   return ( Hm2 * v[i-2][j][k] + Hm1 * v[i-1][j][k] + H0 * v[i][j][k]) / Delta_r(i);
//         // else
//         //   return ( Gp2 * v[i+2][j][k] + Gp1 * v[i+1][j][k] + G0 * v[i][j][k]) / Delta_r(i);
//         if (shift < 0.0)
//             return (v[i - 2][j][k] - 4.0 * v[i - 1][j][k] + 3.0 * v[i][j][k]) / (2.0 * delta_x) * dxdr(i);
//         else
//             return -(v[i + 2][j][k] - 4.0 * v[i + 1][j][k] + 3.0 * v[i][j][k]) / (2.0 * delta_x) * dxdr(i);
//     }
//     //
//     // NOTE: these return derivatives * (dx)^3 !
//     //
//     inline double d4r(int i, int j, int k) {
//         return (v[i - 2][j][k] + v[i + 2][j][k] - 4.0 * (v[i - 1][j][k] + v[i + 1][j][k]) + 6.0 * v[i][j][k])
//             * dxdr(i) / delta_x;
//         // / (delta_r*delta_r*delta_r*delta_r);
//     }
//     inline double d4theta(int i, int j, int k) {
//         return (v[i][j - 2][k] + v[i][j + 2][k] - 4.0 * (v[i][j - 1][k] + v[i][j + 1][k]) + 6.0 * v[i][j][k])
//             * dydtheta(j) / delta_y;
//         // / (delta_theta*delta_theta*delta_theta*delta_theta);
//     }
//     inline double d4phi(int i, int j, int k) {
// #ifdef AXISYMMETRY
//         return 0.0;
// #else
//         return (v[i][j][k - 2] + v[i][j][k + 2] - 4.0 * (v[i][j][k - 1] + v[i][j][k + 1]) + 6.0 * v[i][j][k]) / delta_phi;
// #endif  /* AXISYMMETRY */
//     }
//     inline double d4(int i, int j, int k) {
//         return d4r(i, j, k) + d4theta(i, j, k) + d4phi(i, j, k);
//     }
//     //
//     // Derivatives for Kreiss-Oliger terms - for a *uniform* grid,
//     // these return (Delta x)^p \partial_x^{p+1} f 
//     //
//     inline double D8r(int i, int j, int k) {
//         const double zero = 70.0;
//         const double one = -56.0;
//         const double two = 28.0;
//         const double three = -8.0;
//         const double four = 1.0;
//         return (four * (v[i + 4][j][k] + v[i - 4][j][k]) +
//             three * (v[i + 3][j][k] + v[i - 3][j][k]) +
//             two * (v[i + 2][j][k] + v[i - 2][j][k]) +
//             one * (v[i + 1][j][k] + v[i - 1][j][k]) +
//             zero * v[i][j][k]) * dxdr(i) / delta_x;
//     }
//     inline double D8theta(int i, int j, int k) {
//         const double zero = 70.0;
//         const double one = -56.0;
//         const double two = 28.0;
//         const double three = -8.0;
//         const double four = 1.0;
//         return (four * (v[i][j + 4][k] + v[i][j - 4][k]) +
//             three * (v[i][j + 3][k] + v[i][j - 3][k]) +
//             two * (v[i][j + 2][k] + v[i][j - 2][k]) +
//             one * (v[i][j + 1][k] + v[i][j - 1][k]) +
//             zero * v[i][j][k]) * dydtheta(j) / delta_y;
//     }
//     inline double D8phi(int i, int j, int k) {
// #ifdef AXISYMMETRY
//         return 0.0;
// #else
//         const double zero = 70.0;
//         const double one = -56.0;
//         const double two = 28.0;
//         const double three = -8.0;
//         const double four = 1.0;
//         return (four * (v[i][j][k + 4] + v[i][j][k - 4]) +
//             three * (v[i][j][k + 3] + v[i][j][k - 3]) +
//             two * (v[i][j][k + 2] + v[i][j][k - 2]) +
//             one * (v[i][j][k + 1] + v[i][j][k - 1]) +
//             zero * v[i][j][k]) / delta_phi;
// #endif
//     }
//     //
//     inline double D6r(int i, int j, int k) {
//         const double zero = -20.0;
//         const double one = 15.0;
//         const double two = -6.0;
//         const double three = 1.0;
//         return (three * (v[i + 3][j][k] + v[i - 3][j][k]) +
//             two * (v[i + 2][j][k] + v[i - 2][j][k]) +
//             one * (v[i + 1][j][k] + v[i - 1][j][k]) +
//             zero * v[i][j][k]) * dxdr(i) / delta_x;
//     }
//     inline double D6theta(int i, int j, int k) {
//         const double zero = -20.0;
//         const double one = 15.0;
//         const double two = -6.0;
//         const double three = 1.0;
//         return (three * (v[i][j + 3][k] + v[i][j - 3][k]) +
//             two * (v[i][j + 2][k] + v[i][j - 2][k]) +
//             one * (v[i][j + 1][k] + v[i][j - 1][k]) +
//             zero * v[i][j][k]) * dydtheta(j) / delta_y;
//     }
//     inline double D6phi(int i, int j, int k) {
// #ifdef AXISYMMETRY
//         return 0.0;
// #else
//         const double zero = -20.0;
//         const double one = 15.0;
//         const double two = -6.0;
//         const double three = 1.0;
//         return (three * (v[i][j][k + 3] + v[i][j][k - 3]) +
//             two * (v[i][j][k + 2] + v[i][j][k - 2]) +
//             one * (v[i][j][k + 1] + v[i][j][k - 1]) +
//             zero * v[i][j][k]) / delta_phi;
// #endif  /* AXISYMMETRY */
//     }
//     inline double KO(int i, int j, int k) {
//         // see Appendix C in Babiuc et.al, arxiv:0709.3559 
//         // note wrong sign in (C7); see also eq. (7) in
//         // Bozzola & Paschalidis, arXiv:2104.06978
//         double KO_term = 0.0;
//         if (KO_order == 5)
//             KO_term = KO_factor * (D6r(i, j, k) + D6theta(i, j, k) + D6phi(i, j, k));
//         else if (KO_order == 7)
//             KO_term = KO_factor * (D8r(i, j, k) + D8theta(i, j, k) + D8phi(i, j, k));
//         else {
//             cerr << " Kreiss-Oliger of order " << KO_order << " not implemented! " << endl;
//             exit(1);
//         }
//         return KO_term;
//     }