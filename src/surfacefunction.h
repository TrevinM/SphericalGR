// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing types for surface functions (needed in horizon finder)
//
// Based on NRmatrix in numerical recipes, but includes
// functionality for spherical coordinate systems.
//
//================================================
#ifndef SF_H
#define SF_H

#include "nr3.h"
#include "Grid.h"

class gf2d {
private:
    int nt, np;   // include ghost zones
    int N_g;  // number of ghost zones - hardcoded to three
    Doub** v;
    VecDoub* theta_p, * dydtheta_p, * ddydtheta_p, * phi_p;
    VecDoub* st_p, * ct_p;
    Doub delta_y;
    Doub delta_phi;         // derivative operators
    int gfn;                // optional gridfunction number
    int axis_par, eq_par;   // parities, set to 1 by default
    const char* name;
    Grid* grid;
public:
    //================================================
    // Constructor
    //================================================
    gf2d() : nt(0), np(0), v(NULL), axis_par(1), eq_par(1)
        //	, r(NULL), theta(NULL), x(NULL), phi(NULL) 
    { // cout << " in default gridfunction constructor " << endl; 
    };
    //================================================
    // Set up grid functions 
    //================================================
    //
    //
    // Version that allows to specify function name, number and parity
    // 
    gf2d* setup(Grid* grid_i, const char* name_i,
        int gfn_i, int axis_i, int eq_i) {
        gf2d* address = setup(grid_i);
        axis_par = axis_i;
#ifdef EQSYMMETRY
        eq_par = eq_i;
#endif
        gfn = gfn_i;
        name = name_i;
        return address;
    }
    //
    // Version that allows to specify both function name and number
    // 
    gf2d* setup(Grid* grid_i, const char* name_i, int gfn_i) {
        gf2d* address = setup(grid_i);
        gfn = gfn_i;
        name = name_i;
        return address;
    }
    //
    // Version that allows to specify function number
    // 
    gf2d* setup(Grid* grid_i, int gfn_i) {
        gf2d* address = setup(grid_i);
        gfn = gfn_i;
        return address;
    }
    //
    // Minimum version
    //
    gf2d* setup(Grid* grid_i) {
        grid = grid_i;
        gfn = 0;
        N_g = grid->N_ghosts();
        nt = grid->N_theta_tot();
        np = grid->N_phi_tot();
        theta_p = grid->theta();
        dydtheta_p = grid->dydtheta();
        ddydtheta_p = grid->ddydtheta();
        phi_p = grid->phi();
        ct_p = grid->costheta();
        st_p = grid->sintheta();
        v = new Doub * [nt];
        v[0] = new Doub[nt * np];
        for (int j = 1; j < nt; j++) v[j] = v[j - 1] + np;
        delta_y = grid->delta_y();
        delta_phi = grid->delta_phi();
        return this;
    }
    // subscripting: pointer to row i
    inline Doub* operator[](const int i) { return v[i]; }
    inline Doub const* operator[](const int i) const { return v[i]; }
    inline Doub operator()(int j, int k) { return v[j][k]; }
    inline int dim2() const { return nt; }
    inline int dim3() const { return np; }
    //  inline VecDoub * theta_pointer() const { return theta_p;}
    //  inline VecDoub * x_pointer() const { return x_p;}
    //  inline VecDoub * phi_pointer() const { return phi_p;}
    inline Doub theta(int j) { return (*theta_p)[j]; }
    //  inline Doub x(int j) { return (*x_p)[j]; }
    inline Doub dydtheta(int j) { return (*dydtheta_p)[j]; }
    inline Doub ddydtheta(int j) { return (*ddydtheta_p)[j]; }
    inline Doub costheta(int j) { return (*ct_p)[j]; }
    inline Doub sintheta(int j) { return (*st_p)[j]; }
    inline Doub phi(int k) { return (*phi_p)[k]; }
    inline Doub dphi() const { return delta_phi; }
    //  inline Doub dtheta() const { return delta_theta}
    inline int GridFunctionNumber() { return gfn; }
    inline int Axis_Parity() const { return axis_par; }
    inline int Eq_Parity() const { return eq_par; }
    const char* Name() { return name; }
    //================================================
    // Fill inner ghost zones, assuming a spherical grid.
    //================================================
    void fill_ghosts() {
        //
        // both lower and upper theta:
        //
        for (int k = N_g; k < np - N_g; k++) {
#ifdef AXISYMMETRY
            const int K = k;
#else
            const int K = (k - N_g + (np - 2 * N_g) / 2) % (np - 2 * N_g) + N_g;
#endif   /* AXISYMMETRY */
            // do the upper two first, in case there are only two interior grid points...
            v[N_g - 2][k] = axis_par * v[N_g + 1][K];
            v[N_g - 1][k] = axis_par * v[N_g][K];
#ifdef EQSYMMETRY
            for (int j = 0; j < N_g; j++)
                v[nt - j - 1][k] = eq_par * v[nt - 2 * N_g + j][k];
#else
            for (int j = 0; j < N_g; j++)
                v[nt - j - 1][k] = axis_par * v[nt - 2 * N_g + j][k];
#endif /* EQSYMMETRY */
            // ... now fill in the lower ones
            for (int j = 0; j < N_g - 2; j++)
                v[j][k] = axis_par * v[2 * N_g - 1 - j][K];
        }

        //       v[0][k]    = axis_par * v[5][K];
        //       v[1][k]    = axis_par * v[4][K];
        //       v[2][k]    = axis_par * v[3][K];
        // #ifdef EQSYMMETRY
        //       v[nt-3][k] = eq_par * v[nt-4][k];
        //       v[nt-2][k] = eq_par * v[nt-5][k];
        //       v[nt-1][k] = eq_par * v[nt-6][k];
        // #else
        //       v[nt-3][k] = axis_par * v[nt-4][K];
        //       v[nt-2][k] = axis_par * v[nt-5][K];
        //       v[nt-1][k] = axis_par * v[nt-6][K];
        // #endif /* EQSYMMETRY */
            //
            // finally both lower and upper phi:
            //
        for (int j = 0; j < nt; j++) {
#ifdef AXISYMMETRY
            for (int k = 0; k < N_g; k++) {
                v[j][k] = v[j][N_g];
                v[j][np - 1 - k] = v[j][N_g];
            }
#else
            for (int k = 0; k < N_g; k++) {
                v[j][k] = v[j][np - 2 * N_g + k];
                v[j][np - 1 - k] = v[j][2 * N_g - 1 - k];
            }
#endif /* AXISYMMETRY */
            // #ifdef AXISYMMETRY
            //       v[j][0]    = v[j][N_g];
            //       v[j][1]    = v[j][N_g];
            //       v[j][2]    = v[j][N_g];
            //       v[j][np-1] = v[j][N_g];
            //       v[j][np-2] = v[j][N_g];
            //       v[j][np-3] = v[j][N_g];
            // #else
            //       v[j][0]    = v[j][np-6];
            //       v[j][1]    = v[j][np-5];
            //       v[j][2]    = v[j][np-4];
            //       v[j][np-1] = v[j][5];
            //       v[j][np-2] = v[j][4];
            //       v[j][np-3] = v[j][3];
            //  #endif /* AXISYMMETRY */
        }
    };
    //    
    //================================================
    // Derivative operators 
    // NOTE: this fourth-order implementation assumes equidistant grid!
    //================================================
    //
    // Laplace operator
    //
    inline Doub Laplace(int j, int k) {
        return ddtheta(j, k) + (*ct_p)[j] * dtheta(j, k) / (*st_p)[j]
            + ddphi(j, k) / ((*st_p)[j] * (*st_p)[j]);
    }
    //
    // first derivatives
    // 
    inline Doub dtheta(int j, int k) {
        return (v[j - 2][k] - 8.0 * (v[j - 1][k] - v[j + 1][k]) - v[j + 2][k])
            / (12.0 * delta_y) * dydtheta(j);
    }
    inline Doub dphi(int j, int k) {
#ifdef AXISYMMETRY
        return 0.0;
#else
        return (v[j][k - 2] - 8.0 * (v[j][k - 1] - v[j][k + 1]) - v[j][k + 2])
            / (12.0 * delta_phi);
#endif   /* AXISYMMETRY */
    }
    //
    // second derivatives
    // 
    inline Doub ddtheta(int j, int k) {
        return (-(v[j - 2][k] + v[j + 2][k]) - 30.0 * v[j][k] + 16.0 * (v[j + 1][k] + v[j - 1][k])) /
            (12.0 * delta_y * delta_y) * (dydtheta(j) * dydtheta(j)) +
            (v[j - 2][k] - 8.0 * (v[j - 1][k] - v[j + 1][k]) - v[j + 2][k]) / (12.0 * delta_y) * ddydtheta(j);
    }
    inline Doub ddphi(int j, int k) {
#ifdef AXISYMMETRY
        return 0.0;
#else
        return (-(v[j][k - 2] + v[j][k + 2]) - 30.0 * v[j][k] + 16.0 * (v[j][k + 1] + v[j][k - 1])) /
            (12.0 * delta_phi * delta_phi);
#endif   /* AXISYMMETRY */
    }
    //
    // mixed second derivatives
    //  
    inline Doub dphidtheta(int j, int k) { return dthetadphi(j, k); }
    inline Doub dthetadphi(int j, int k) {
#ifdef AXISYMMETRY
        return 0.0;
#else
        return ((v[j - 2][k - 2] - 8.0 * (v[j - 2][k - 1] - v[j - 2][k + 1]) - v[j - 2][k + 2])
            - 8.0 * (v[j - 1][k - 2] - 8.0 * (v[j - 1][k - 1] - v[j - 1][k + 1]) - v[j - 1][k + 2])
            + 8.0 * (v[j + 1][k - 2] - 8.0 * (v[j + 1][k - 1] - v[j + 1][k + 1]) - v[j + 1][k + 2])
            - (v[j + 2][k - 2] - 8.0 * (v[j + 2][k - 1] - v[j + 2][k + 1]) - v[j + 2][k + 2])) /
            (144.0 * delta_y * delta_phi) * dydtheta(j);
#endif   /* AXISYMMETRY */
    }
    //
    // first derivatives in UPWIND differencing (NOTE: same name for routine, but extra argument for direction)
    //
    // NOTE: If these routines are needed, implement 3rd-order version.
    // 
  //   inline Doub dtheta(int j, int k, double shift) {
  //     if (shift < 0.0) 
  //       return   (v[j-2][k] - 4.0*v[j-1][k] + 3.0*v[j][k])/(2.0*delta_theta);
  //     else
  //       return - (v[j+2][k] - 4.0*v[j+1][k] + 3.0*v[j][k])/(2.0*delta_theta);
  //   }
  //   inline Doub dphi(int j, int k, double shift) {
  // #ifdef AXISYMMETRY
  //     return 0.0;
  // #else
  //     if (shift < 0.0) 
  //       return   (v[j][k-2] - 4.0*v[j][k-1] + 3.0*v[j][k])/(2.0*delta_phi);
  //     else
  //       return - (v[j][k+2] - 4.0*v[j][k+1] + 3.0*v[j][k])/(2.0*delta_phi);
  // #endif   /* AXISYMMETRY */
  //  }
    //
    // fourth derivatives (for Kreiss-Oliger dissipation)
    //
    // NOTE: these return derivatives * (dx)^4 !
    //
    inline Doub d4theta(int j, int k) {
        return (v[j - 2][k] + v[j + 2][k] - 4.0 * (v[j - 1][k] + v[j + 1][k]) + 6.0 * v[j][k]);
        // / (delta_theta*delta_theta*delta_theta*delta_theta);
    }
    inline Doub d4phi(int j, int k) {
#ifdef AXISYMMETRY
        return 0.0;
#else
        return (v[j][k - 2] + v[j][k + 2] - 4.0 * (v[j][k - 1] + v[j][k + 1]) + 6.0 * v[j][k]);
        //  / (delta_phi*delta_phi*delta_phi*delta_phi);
#endif   /* AXISYMMETRY */
    }
    inline Doub d4(int j, int k) {
        return d4theta(j, k) + d4phi(j, k);
    }
    //================================================
    // several utilities
    //================================================
    double min() {
        double min = v[N_g][N_g];
        for (int j = N_g; j < nt - N_g; j++)
            for (int k = N_g; k < np - N_g; k++) {
                if (v[j][k] < min) min = v[j][k];
            }
        return min;
    };
    double max() {
        double max = v[N_g][N_g];
        for (int j = N_g; j < nt - N_g; j++)
            for (int k = N_g; k < np - N_g; k++) {
                if (v[j][k] > max) max = v[j][k];
            }
        return max;
    };
    bool FINITE() {
        bool fine = true;
        for (int j = N_g; j < nt - N_g; j++)
            for (int k = N_g; k < np - N_g; k++) {
                if (!isfinite(v[j][k])) fine = false;
            }
        return fine;
    };
    double L2_norm() {
        double norm = 0.0;
        double vol_int = 0.0;
        for (int j = N_g; j < nt - N_g; j++)
            for (int k = N_g; k < np - N_g; k++) {
                const double delta_theta = delta_y / dydtheta(j);
                const double SqrtJac = (*st_p)[j] * delta_theta * delta_phi;
                norm += SqrtJac * v[j][k] * v[j][k];
                vol_int += SqrtJac;
            }
        return sqrt(norm);
    }
    //================================================
    // Interpolator to arbitrary value of theta >= 0
    //================================================
    Doub operator()(double theta, int k) {
        int order = 4;
        // first find jlo
        int jlo = 0;
        while ((*theta_p)[jlo] < theta) jlo += 1;
        jlo -= 2;
        if (jlo < 0) jlo = 0;
        if (jlo > nt - order) jlo = nt - order; // take care of upper boundary
        //    cout << " jlo : " << jlo << endl;
        Doub* ta, * ya;
        ta = new Doub[order];      // allocate data arrays... 
        ya = new Doub[order];
        for (int j = 0; j < order; j++) {
            ta[j] = (*theta_p)[jlo + j];      // ...and fill them with data
            ya[j] = v[jlo + j][k];
        }
        //
        // Now follow numerical recipes routine Poly_interp
        // 
        Doub* c, * d;
        c = new Doub[order];      // allocate arrays that store differenes in tableaus
        d = new Doub[order];
        Doub dif = abs(theta - ta[0]);
        Doub dift;
        int ns = 0;
        for (int i = 0; i < order; i++) {
            if ((dift = abs(theta - ta[i])) < dif) {
                ns = i;
                dif = dift;
            }
            c[i] = ya[i];
            d[i] = ya[i];
        }
        Doub y = ya[ns--];
        Doub den, ho, hp, w, dy;
        for (int m = 1; m < order; m++) {
            for (int i = 0; i < order - m; i++) {
                ho = ta[i] - theta;
                hp = ta[i + m] - theta;
                w = c[i + 1] - d[i];
                den = ho - hp;
                den = w / den;
                d[i] = hp * den;
                c[i] = ho * den;
            }
            y += (dy = (2 * (ns + 1) < (order - m) ? c[ns + 1] : d[ns--]));
        }
        delete ta;
        delete ya;
        delete c;
        delete d;
        return y;
    };


    //================================================
    // Destructor
    //================================================
    ~gf2d() {
        if (v != NULL) {
            // cout << " deleting gf2d... " << endl; 
            delete[](v[0]);
            delete[](v);
        }
    }
};

#endif  /* SF_H */
