// Tell emacs that this is -*-c++-*- mode
//================================================
// Classes for (spatial) gauge conditions
//================================================
#ifndef GAUGE_H
#define GAUGE_H

#include "Grid.h"
#include "nr3.h"
#include "State.h"
#include "gridfunction.h"

enum { constantshift, gammadriver, advective_gammadriver, jena_driver, advective_jena_driver, covariant_advective_jena_driver, self_sim };
//
//================================================
// Base class - doesn't do much...
//================================================
//
class Gauge {
protected:
    int gauge_type;
    double eta_KO;
    Grid* grid;
    int N_g, N_r, N_t, N_p;
public:
    Gauge(Grid* grid_i) : grid(grid_i) {
        eta_KO = 0.0;
        N_g = grid->N_ghosts();
        N_r = grid->N_r_tot();
        N_t = grid->N_theta_tot();
        N_p = grid->N_phi_tot();
    };
    virtual ~Gauge() {};
    //================================================
    // Virtual functions that will be overwritten by derived
    // classes below
    //================================================
    virtual string Name() = 0;
    virtual void dot_shift(state* c, state* derivs) = 0;
    virtual void overwrite_shift(state* s, state* derivs, double tau_c) = 0;
    double Set_eta(double eta_in) {
        cout << " GAUGE: Using Kreiss-Oliger coefficient " << eta_in << " in Gauge.h " << endl;
        eta_KO = eta_in;
        return eta_KO;
    };
    int GaugeType() { return gauge_type; }
};

//
//================================================
// Derived class for constant shift
//================================================
//
class Constant_Shift : public Gauge {
public:
    Constant_Shift(Grid* grid_i) : Gauge(grid_i) {
        cout << " GAUGE: setting up constant shift..." << endl;
        gauge_type = constantshift;
    };
    ~Constant_Shift() {};
    string Name() {
        return "constant shift";
    };
    //================================================
    // Overwrite virtual functions
    //================================================
    void dot_shift(state* c, state* derivs) {
        derivs->shift_r.equals(0.0);
        derivs->shift_t.equals(0.0);
        derivs->shift_p.equals(0.0);
        derivs->B_r.equals(0.0);
        derivs->B_t.equals(0.0);
        derivs->B_p.equals(0.0);
    }
    void overwrite_shift(state* s, state* derivs, double tau_c) {};
};

//
//================================================
// Derived class for (non-advective) Gamma-driver shift, integrates
//
// CHECK: eta B^i term...
//
//  \partial_t B^i     = mu_S * \partial_t \Lambda^i - eta B^i
//  \partial_t \beta^i = B^i
//
//================================================
//
class Gamma_Driver : public Gauge {
private:
    double mu_S;
    double eta;    // coefficient for damping term
    ostringstream gauge_name;
public:
    Gamma_Driver(Grid* grid_i) : Gauge(grid_i) {
        cout << " GAUGE: setting up non-advective Gamma driver..." << endl;
        gauge_type = gammadriver;
        // default values for parameters
        eta = 2.0;
        mu_S = 0.75;
        // look for file "Shift_Driver_Input" for parameters
        ifstream infile;
        infile.open("Shift_Driver_Input");
        if (!infile)
            cout << " JENA_DRIVER: Can't open Shift_Driver_Input -- will use default values " << endl;
        else {
            char buf[100], c;
            infile.get(buf, 100, '='); infile.get(c); infile >> eta;
            infile.get(buf, 100, '='); infile.get(c); infile >> mu_S;
        }
        cout << " GAMMA_DRIVER: Using eta = " << eta << " and mu_S = " << mu_S << endl;
        if (infile) infile.close();
        gauge_name << "(non-advective) Gamma-driver with eta = " << eta << " and mu_S = " << mu_S << ends;
    };
    ~Gamma_Driver() {};
    string Name() {
        return gauge_name.str();
    };
    //================================================
    // Overwrite virtual functions
    //================================================
    void dot_shift(state* c, state* derivs) {
        for (int i = N_g; i < N_r - N_g; i++)
            for (int j = N_g; j < N_t - N_g; j++)
                for (int k = N_g; k < N_p - N_g; k++) {
                    derivs->shift_r[i][j][k] = c->B_r(i, j, k) + eta_KO * c->shift_r.KO(i, j, k);
                    derivs->shift_t[i][j][k] = c->B_t(i, j, k) + eta_KO * c->shift_t.KO(i, j, k);
                    derivs->shift_p[i][j][k] = c->B_p(i, j, k) + eta_KO * c->shift_p.KO(i, j, k);
                    //
                    derivs->B_r[i][j][k] = mu_S * derivs->lam_r(i, j, k) - eta * c->B_r(i, j, k)
                        + eta_KO * c->B_r.KO(i, j, k);
                    derivs->B_t[i][j][k] = mu_S * derivs->lam_t(i, j, k) - eta * c->B_t(i, j, k)
                        + eta_KO * c->B_t.KO(i, j, k);
                    derivs->B_p[i][j][k] = mu_S * derivs->lam_p(i, j, k) - eta * c->B_p(i, j, k)
                        + eta_KO * c->B_p.KO(i, j, k);
                }
    }
    void overwrite_shift(state* s, state* derivs, double tau_c) {};
};
//
//================================================
// Derived class for advective Gamma-driver shift, integrates
//
// CHECK damping term
//
//  \partial_t B^i     - beta^j \partial_j B^i     = mu_S * \partial_t \Lambda^i - eta B^i
//  \partial_t \beta^i - beta^j \partial_j \beta^i = B^i
//
//================================================
//
class Advective_Gamma_Driver : public Gauge {
private:
    double mu_S;
    double eta;    // coefficient for damping term
    ostringstream gauge_name;
public:
    Advective_Gamma_Driver(Grid* grid_i) : Gauge(grid_i) {
        cout << " GAUGE: setting up advective Gamma driver..." << endl;
        gauge_type = advective_gammadriver;
        eta = 2.0;
        mu_S = 0.75;
        // look for file "Shift_Driver_Input" for parameters
        ifstream infile;
        infile.open("Shift_Driver_Input");
        if (!infile)
            cout << " JENA_DRIVER: Can't open Shift_Driver_Input -- will use default values " << endl;
        else {
            char buf[100], c;
            infile.get(buf, 100, '='); infile.get(c); infile >> eta;
            infile.get(buf, 100, '='); infile.get(c); infile >> mu_S;
        }
        cout << " GAMMA_DRIVER: Using eta = " << eta << " and mu_S = " << mu_S << endl;
        if (infile) infile.close();
        gauge_name << "(advective) Gamma-driver with eta = " << eta << " and mu_S = " << mu_S;
    };
    ~Advective_Gamma_Driver() {};
    string Name() {
        return gauge_name.str();
    };
    //================================================
    // Overwrite virtual functions
    //================================================
    void dot_shift(state* c, state* derivs) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                const double ctl = grid->costheta(j);
                const double st2 = stl * stl;
                for (int k = N_g; k < N_p - N_g; k++) {
                    const double shift_r = c->shift_r(i, j, k);
                    const double shift_t = c->shift_t(i, j, k) / rl;
                    const double shift_p = c->shift_p(i, j, k) / (rl * stl);
                    //
                    // Note: compute advective terms in terms of physical quantities - rescale to code units in actual update...
                    //
                    double beta_d_beta_r =
                        shift_r * c->shift_r.dr(i, j, k, shift_r) +
                        shift_t * c->shift_r.dtheta(i, j, k, shift_t) +
                        shift_p * c->shift_r.dphi(i, j, k, shift_p);
                    double beta_d_beta_t =
                        shift_r * (c->shift_t.dr(i, j, k, shift_r) / rl - c->shift_t(i, j, k) / r2) +
                        shift_t * c->shift_t.dtheta(i, j, k, shift_t) / rl +
                        shift_p * c->shift_t.dphi(i, j, k, shift_p) / rl;
                    double beta_d_beta_p =
                        shift_r * (c->shift_p.dr(i, j, k, shift_r) / (rl * stl) - c->shift_p(i, j, k) / (r2 * stl)) +
                        shift_t * (c->shift_p.dtheta(i, j, k, shift_t) / (rl * stl) - c->shift_p(i, j, k) / (rl * st2) * ctl) +
                        shift_p * c->shift_p.dphi(i, j, k, shift_p) / (rl * stl);
                    //
                    double beta_d_B_r =
                        shift_r * c->B_r.dr(i, j, k, shift_r) +
                        shift_t * c->B_r.dtheta(i, j, k, shift_t) +
                        shift_p * c->B_r.dphi(i, j, k, shift_p);
                    double beta_d_B_t =
                        shift_r * (c->B_t.dr(i, j, k, shift_r) / rl - c->B_t(i, j, k) / r2) +
                        shift_t * c->B_t.dtheta(i, j, k, shift_t) / rl +
                        shift_p * c->B_t.dphi(i, j, k, shift_p) / rl;
                    double beta_d_B_p =
                        shift_r * (c->B_p.dr(i, j, k, shift_r) / (rl * stl) - c->B_p(i, j, k) / (r2 * stl)) +
                        shift_t * (c->B_p.dtheta(i, j, k, shift_t) / (rl * stl) - c->B_p(i, j, k) / (rl * st2) * ctl) +
                        shift_p * c->B_p.dphi(i, j, k, shift_p) / (rl * stl);
                    //
                        // rescale derivative terms
                    //
                    beta_d_beta_t *= rl;
                    beta_d_beta_p *= rl * stl;
                    beta_d_B_t *= rl;
                    beta_d_B_p *= rl * stl;
                    //
                    // now compute derivatives
                    // 
                    derivs->shift_r[i][j][k] = c->B_r(i, j, k) + beta_d_beta_r + eta_KO * c->shift_r.KO(i, j, k);
                    derivs->shift_t[i][j][k] = c->B_t(i, j, k) + beta_d_beta_t + eta_KO * c->shift_t.KO(i, j, k);
                    derivs->shift_p[i][j][k] = c->B_p(i, j, k) + beta_d_beta_p + eta_KO * c->shift_p.KO(i, j, k);
                    //
                    derivs->B_r[i][j][k] = mu_S * derivs->lam_r(i, j, k) + beta_d_B_r - eta * c->B_r(i, j, k)
                        + eta_KO * c->B_r.KO(i, j, k);
                    derivs->B_t[i][j][k] = mu_S * derivs->lam_t(i, j, k) + beta_d_B_t - eta * c->B_t(i, j, k)
                        + eta_KO * c->B_t.KO(i, j, k);
                    derivs->B_p[i][j][k] = mu_S * derivs->lam_p(i, j, k) + beta_d_B_p - eta * c->B_p(i, j, k)
                        + eta_KO * c->B_p.KO(i, j, k);
                }
            }
        }
    };
    void overwrite_shift(state* s, state* derivs, double tau_c) {};
};
//
//================================================
// Derived class for Jena version of (non-advective) Gamma-driver shift, integrates
//
//  \partial_t \beta^i = mu_S \tilde \Gamma^i - \eta \beta^i
//
// (see Thierfelder et.al., arXiv:1104.4751 (TBB) )
//================================================
//
class Jena_Gamma_Driver : public Gauge {
private:
    double eta;    // coefficient for damping term
    double mu_S;   // coefficient; see eq. (15) in TBB
    bool muequalsalpha2;
    ostringstream gauge_name;
public:
    Jena_Gamma_Driver(Grid* grid_i) : Gauge(grid_i) {
        cout << " GAUGE: setting up non-advective Jena Gamma driver..." << endl;
        gauge_type = jena_driver;
        // default values for parameters
        eta = 0.0;
        mu_S = 0.75;
        muequalsalpha2 = false;
        // look for file "Jena_Driver_Input" for parameters
        ifstream infile;
        infile.open("Shift_Driver_Input");
        if (!infile)
            cout << " JENA_DRIVER: Can't open Shift_Driver_Input -- will use default values " << endl;
        else {
            char buf[100], c;
            infile.get(buf, 100, '='); infile.get(c); infile >> eta;
            infile.get(buf, 100, '='); infile.get(c); infile >> mu_S;
            infile.get(buf, 100, '='); infile.get(c); infile >> muequalsalpha2;
        }
        if (muequalsalpha2)
            cout << " JENA_DRIVER: Using eta = " << eta << " and mu_S = alpha^2." << endl;
        else
            cout << " JENA_DRIVER: Using eta = " << eta << " and mu_S = " << mu_S << endl;
        if (infile)
            infile.close();
        if (muequalsalpha2)
            gauge_name << "(non-advective) Jena version of Gamma-driver with eta = " << eta << " and mu_S = alpha^2 ";
        else
            gauge_name << "(non-advective) Jena version of Gamma-driver with eta = " << eta << " and mu_S = " << mu_S;
    };
    ~Jena_Gamma_Driver() {};
    string Name() {
        return gauge_name.str();
    };
    //================================================
    // Overwrite virtual functions
    //================================================
    void dot_shift(state* c, state* derivs) {
        for (int i = N_g; i < N_r - N_g; i++)
            for (int j = N_g; j < N_t - N_g; j++)
                for (int k = N_g; k < N_p - N_g; k++) {
                    if (muequalsalpha2) mu_S = c->lapse(i, j, k) * c->lapse(i, j, k);
                    derivs->shift_r[i][j][k] = mu_S * c->lam_r(i, j, k) - eta * c->shift_r(i, j, k) + eta_KO * c->shift_r.KO(i, j, k);
                    derivs->shift_t[i][j][k] = mu_S * c->lam_t(i, j, k) - eta * c->shift_t(i, j, k) + eta_KO * c->shift_t.KO(i, j, k);
                    derivs->shift_p[i][j][k] = mu_S * c->lam_p(i, j, k) - eta * c->shift_p(i, j, k) + eta_KO * c->shift_p.KO(i, j, k);
                    derivs->B_r[i][j][k] = 0.0;
                    derivs->B_t[i][j][k] = 0.0;
                    derivs->B_p[i][j][k] = 0.0;
                };
    }
    void overwrite_shift(state* s, state* derivs, double tau_c) {};
};
//
//================================================
// Derived class for Jena version of advective Gamma-driver shift, integrates
//
//  \partial_t \beta^i - beta^j \partial_j \beta^i = mu_S \tilde \Gamma^i - \eta \beta^i
//
// (see Thierfelder et.al., arXiv:1104.4751 (TBB) )
//================================================
//
class Advective_Jena_Gamma_Driver : public Gauge {
private:
    double eta;    // coefficient for damping term
    double mu_S;   // coefficient; see eq. (15) in TBB
    bool muequalsalpha2;
    ostringstream gauge_name;
    bool fix;
    double z_fix;
    double sigma;   // damping lenth scale for fix-point term
public:
    Advective_Jena_Gamma_Driver(Grid* grid_i) : Gauge(grid_i) {
        cout << " GAUGE: setting up advective Jena Gamma driver..." << endl;
        gauge_type = advective_jena_driver;
        // default values for parameters
        eta = 0.0;
        mu_S = 0.75;
        muequalsalpha2 = false;
        // look for file "Shift_Driver_Input" for parameters
        ifstream infile;
        infile.open("Shift_Driver_Input");
        if (!infile)
            cout << " JENA_DRIVER: Can't open Shift_Driver_Input -- will use default values " << endl;
        else {
            char buf[100], c;
            infile.get(buf, 100, '='); infile.get(c); infile >> eta;
            infile.get(buf, 100, '='); infile.get(c); infile >> mu_S;
            infile.get(buf, 100, '='); infile.get(c); infile >> muequalsalpha2;
            infile.get(buf, 100, '='); infile.get(c); infile >> fix;
            infile.get(buf, 100, '='); infile.get(c); infile >> z_fix;
            infile.get(buf, 100, '='); infile.get(c); infile >> sigma;
        }
        if (muequalsalpha2)
            cout << " JENA_DRIVER: Using eta = " << eta << " and mu_S = alpha^2."
            << endl;
        else
            cout << " JENA_DRIVER: Using eta = " << eta << " and mu_S = "
            << mu_S << endl;
        if (fix)
            cout << " JENA_DRIVER: Keeping shift zero at fixed-point z_fix "
            << z_fix << " and damping length sigma = " << sigma << endl;
        if (infile)
            infile.close();
        if (muequalsalpha2)
            gauge_name << "advective Jena version of Gamma-driver with eta = " << eta << " and mu_S = alpha^2 ";
        else
            gauge_name << "advective Jena version of Gamma-driver with eta = " << eta << " and mu_S = " << mu_S;
        if (fix)
            gauge_name << ", using fixed-point at z_fix = " << z_fix << " with sigma = " << sigma;
    };
    ~Advective_Jena_Gamma_Driver() {};
    string Name() {
        return gauge_name.str();
    };
    //================================================
    // Overwrite virtual functions
    //================================================
    void dot_shift(state* c, state* derivs) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                const double ctl = grid->costheta(j);
                const double st2 = stl * stl;
                for (int k = N_g; k < N_p - N_g; k++) {
                    if (muequalsalpha2) mu_S = c->lapse(i, j, k) * c->lapse(i, j, k);
                    const double shift_r = c->shift_r(i, j, k);
                    const double shift_t = c->shift_t(i, j, k) / rl;
                    const double shift_p = c->shift_p(i, j, k) / (rl * stl);
                    //
                    // Note: compute advective terms in terms of physical quantities - rescale to code units in actual update...
                    //
                    double beta_d_beta_r =
                        shift_r * c->shift_r.dr(i, j, k, shift_r) +
                        shift_t * c->shift_r.dtheta(i, j, k, shift_t) +
                        shift_p * c->shift_r.dphi(i, j, k, shift_p);
                    double beta_d_beta_t =
                        shift_r * (c->shift_t.dr(i, j, k, shift_r) / rl - c->shift_t(i, j, k) / r2) +
                        shift_t * c->shift_t.dtheta(i, j, k, shift_t) / rl +
                        shift_p * c->shift_t.dphi(i, j, k, shift_p) / rl;
                    double beta_d_beta_p =
                        shift_r * (c->shift_p.dr(i, j, k, shift_r) / (rl * stl) - c->shift_p(i, j, k) / (r2 * stl)) +
                        shift_t * (c->shift_p.dtheta(i, j, k, shift_t) / (rl * stl) - c->shift_p(i, j, k) / (rl * st2) * ctl) +
                        shift_p * c->shift_p.dphi(i, j, k, shift_p) / (rl * stl);
                    //
                        // rescale derivative terms
                    //
                    beta_d_beta_t *= rl;
                    beta_d_beta_p *= rl * stl;
                    //
                    // now compute derivatives
                    //
                    derivs->shift_r[i][j][k] = mu_S * c->lam_r(i, j, k) - eta * c->shift_r(i, j, k) + beta_d_beta_r
                        + eta_KO * c->shift_r.KO(i, j, k);
                    derivs->shift_t[i][j][k] = mu_S * c->lam_t(i, j, k) - eta * c->shift_t(i, j, k) + beta_d_beta_t
                        + eta_KO * c->shift_t.KO(i, j, k);
                    derivs->shift_p[i][j][k] = mu_S * c->lam_p(i, j, k) - eta * c->shift_p(i, j, k) + beta_d_beta_p
                        + eta_KO * c->shift_p.KO(i, j, k);
                    derivs->B_r[i][j][k] = 0.0;
                    derivs->B_t[i][j][k] = 0.0;
                    derivs->B_p[i][j][k] = 0.0;
                }
            }
        }
        //
        // if using fixed-point, make sure that all time-derivatives vanish
        // at fixed-point (assumed to be on axis...)
        // 
        derivs->shift_r.fill_ghosts();
        int order = 2;  // just to try this out...
        if (fix) {
            // CHECK just a test...
            // const double dot_shift_z = -100.0;
            const double dot_shift_1 = derivs->shift_r(N_g - 2, 0.0, N_g);
            const double dot_shift_2 = derivs->shift_r(N_g - 1, 0.0, N_g);
            const double dot_shift_3 = derivs->shift_r(N_g, 0.0, N_g);
            const double dot_shift_4 = derivs->shift_r(N_g + 1, 0.0, N_g);
            const double dot_shift_5 = derivs->shift_r(N_g + 2, 0.0, N_g);
            const double dot_shift_6 = derivs->shift_r(N_g + 3, 0.0, N_g);
            const double dot_shift_7 = derivs->shift_r(N_g + 4, 0.0, N_g);
            const double dot_shift_z = derivs->shift_r(z_fix, 0.0, N_g);
            // HACK - CHECK!!!
            // const double dot_shift_z = 0.5 * (derivs->shift_r(0, 0.0, N_g) + derivs->shift_r(2*N_g-1, 0.0, N_g));
            // cout << " dot_shift_z = " << dot_shift_z
            // 	   << " -2 " << dot_shift_1  
            // 	   << " -1 " << dot_shift_2  
            // 	   << " N_g " << dot_shift_3  
            // 	   << " +1 " << dot_shift_4  
            // 	   << " +2 " << dot_shift_5  
            // 	   << " +3 " << dot_shift_6  
            // 	   << " +4 " << dot_shift_7  
            // 	   << endl;
            for (int i = N_g; i < N_r - N_g; i++) {
                const double rl = grid->r(i);
                for (int j = N_g; j < N_t - N_g; j++) {
                    const double stl = grid->sintheta(j);
                    const double ctl = grid->costheta(j);
                    const double zl = rl * ctl;
                    const double rhol = rl * stl;
                    const double dist2 = (zl - z_fix) * (zl - z_fix) + rhol * rhol;
                    const double factor = exp(-dist2 / (sigma * sigma));
                    for (int k = N_g; k < N_p - N_g; k++) {
                        derivs->shift_r[i][j][k] -= factor * dot_shift_z * ctl;
                        // CHECK sign!!
                        derivs->shift_t[i][j][k] += factor * dot_shift_z * stl;
                    }
                }
            }
        }
    }
    void overwrite_shift(state* s, state* derivs, double tau_c) {};
};
//
//================================================
// Derived class for Jena version of advective Gamma-driver shift, using covariant derivatives (with respect to 
// reference metric) integrates
//
//  \partial_t \beta^i - beta^j \hat D_j \beta^i = mu_S \tilde \Gamma^i - \eta \beta^i
//
// (see Thierfelder et.al., arXiv:1104.4751 (TBB) )
//================================================
//
class Covariant_Advective_Jena_Gamma_Driver : public Gauge {
private:
    double eta;    // coefficient for damping term
    double mu_S;   // coefficient; see eq. (15) in TBB
    bool muequalsalpha2;
    ostringstream gauge_name;
public:
    Covariant_Advective_Jena_Gamma_Driver(Grid* grid_i) : Gauge(grid_i) {
        cout << " GAUGE: setting up advective covariant Jena Gamma driver..." << endl;
        gauge_type = advective_jena_driver;
        // default values for parameters
        eta = 0.0;
        mu_S = 0.75;
        muequalsalpha2 = false;
        // look for file "Jena_Driver_Input" for parameters
        ifstream infile;
        infile.open("Jena_Driver_Input");
        if (!infile)
            cout << " JENA_DRIVER: Can't open Jena_Driver_Input -- will use default values " << endl;
        else {
            char buf[100], c;
            infile.get(buf, 100, '='); infile.get(c); infile >> eta;
            infile.get(buf, 100, '='); infile.get(c); infile >> mu_S;
            infile.get(buf, 100, '='); infile.get(c); infile >> muequalsalpha2;
        }
        if (muequalsalpha2)
            cout << " JENA_DRIVER: Using eta = " << eta << " and mu_S = alpha^2." << endl;
        else
            cout << " JENA_DRIVER: Using eta = " << eta << " and mu_S = " << mu_S << endl;
        if (infile)
            infile.close();
        if (muequalsalpha2)
            gauge_name << "covariant advective Jena version of Gamma-driver with eta = " << eta << " and mu_S = alpha^2 ";
        else
            gauge_name << "covariant advective Jena version of Gamma-driver with eta = " << eta << " and mu_S = " << mu_S;
    };
    ~Covariant_Advective_Jena_Gamma_Driver() {};
    string Name() {
        return gauge_name.str();
    };
    //================================================
    // Overwrite virtual functions
    //================================================
    void dot_shift(state* c, state* derivs) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                const double ctl = grid->costheta(j);
                const double st2 = stl * stl;
                for (int k = N_g; k < N_p - N_g; k++) {
                    if (muequalsalpha2) mu_S = c->lapse(i, j, k) * c->lapse(i, j, k);
                    const double shift_r = c->shift_r(i, j, k);
                    const double shift_t = c->shift_t(i, j, k) / rl;
                    const double shift_p = c->shift_p(i, j, k) / (rl * stl);
                    //
                    // Note: compute advective terms in terms of physical quantities - rescale to code units in actual update...
                    //
                    double beta_d_beta_r =
                        shift_r * c->shift_r.dr(i, j, k, shift_r) +
                        shift_t * c->shift_r.dtheta(i, j, k, shift_t) +
                        shift_p * c->shift_r.dphi(i, j, k, shift_p)
                        // add covariant derivative terms
                        - rl * shift_t * shift_t - rl * st2 * shift_p * shift_p;
                    double beta_d_beta_t =
                        shift_r * (c->shift_t.dr(i, j, k, shift_r) / rl - c->shift_t(i, j, k) / r2) +
                        shift_t * c->shift_t.dtheta(i, j, k, shift_t) / rl +
                        shift_p * c->shift_t.dphi(i, j, k, shift_p) / rl
                        // add covariant derivative terms
                        - stl * ctl * shift_p * shift_p + 2.0 * shift_r * shift_t / rl;
                    double beta_d_beta_p =
                        shift_r * (c->shift_p.dr(i, j, k, shift_r) / (rl * stl) - c->shift_p(i, j, k) / (r2 * stl)) +
                        shift_t * (c->shift_p.dtheta(i, j, k, shift_t) / (rl * stl) - c->shift_p(i, j, k) / (rl * st2) * ctl) +
                        shift_p * c->shift_p.dphi(i, j, k, shift_p) / (rl * stl)
                        // add covariant derivative terms
                        + 2.0 * shift_r * shift_p / rl + 2.0 * shift_t * shift_p * ctl / stl;
                    //
                        // rescale derivative terms
                    //
                    beta_d_beta_t *= rl;
                    beta_d_beta_p *= rl * stl;
                    //
                    // now compute derivatives
                    //
                    derivs->shift_r[i][j][k] = mu_S * c->lam_r(i, j, k) - eta * c->shift_r(i, j, k) + beta_d_beta_r
                        + eta_KO * c->shift_r.KO(i, j, k);
                    derivs->shift_t[i][j][k] = mu_S * c->lam_t(i, j, k) - eta * c->shift_t(i, j, k) + beta_d_beta_t
                        + eta_KO * c->shift_t.KO(i, j, k);
                    derivs->shift_p[i][j][k] = mu_S * c->lam_p(i, j, k) - eta * c->shift_p(i, j, k) + beta_d_beta_p
                        + eta_KO * c->shift_p.KO(i, j, k);
                    derivs->B_r[i][j][k] = 0.0;
                    derivs->B_t[i][j][k] = 0.0;
                    derivs->B_p[i][j][k] = 0.0;
                }
            }
        }
    }
    void overwrite_shift(state* s, state* derivs, double tau_c) {};
};

//
//================================================
// Derived class for self-similar shift: will supply
// zero time derivatives, but will overwrite shift after
// each time step.
//================================================
//
class Self_Sim_Shift : public Gauge {
private:
    double tau_acc;
    double beta_max;
    ostringstream gauge_name;

public:
    Self_Sim_Shift(Grid* grid_i) : Gauge(grid_i) {
        cout << " GAUGE: setting up self_similar shift..." << endl;
        gauge_type = self_sim;

        ifstream infile;
        infile.open("Self_Sim_Shift_Input");
        if (!infile)
            cout << " GAUGE: Can't open Self_Sim_Shift_Input -- will use default values " << endl;
        else {
            char buf[100], c;
            infile.get(buf, 100, '='); infile.get(c); infile >> tau_acc;
            infile.get(buf, 100, '='); infile.get(c); infile >> beta_max;
        }
        if (infile) infile.close();
        gauge_name << "Self-similar shift with tau_acc = " << tau_acc << " and beta_max = " << beta_max;


    };
    ~Self_Sim_Shift() {};
    string Name() {
        return gauge_name.str();
    };
    //================================================
    // Overwrite virtual functions
    //================================================
    void dot_shift(state* c, state* derivs) {
        derivs->shift_r.equals(0.0);
        derivs->shift_t.equals(0.0);
        derivs->shift_p.equals(0.0);
        derivs->B_r.equals(0.0);
        derivs->B_t.equals(0.0);
        derivs->B_p.equals(0.0);
    }
    //================================================
    // For self-similar shift, compute shift directly, rather than
    // integrating its time derivative 
    //================================================
    void overwrite_shift(state* s, state* derivs, double tau_c) {
        //
        // first: compute n^a nabla_a R having used current shift
        //
        derivs->h_rr.fill_ghosts();
        derivs->phi.fill_ghosts();
        for (int j = N_g; j < N_t - N_g; j++)
            for (int k = N_g; k < N_p - N_g; k++) {
                // for each ray integrate outwards: start with proper distance
                // of first interior grid point from origin
                double gamma_bar_rr = 1.0 + s->h_rr(0.0, j, k);
                double dRdr_mid = exp(2.0 * s->phi(0.0, j, k)) * sqrt(gamma_bar_rr);
                double R = 0.5 * dRdr_mid * grid->delta_r(N_g);
                // compute time derivative of (gamma_rr)^{1/2}
                // CHECK!!! 0.0 should be 0.5...
                double gam_dot = dRdr_mid * (2.0 * derivs->phi(0.0, j, k) +
                    0.5 * derivs->h_rr(0.0, j, k) / gamma_bar_rr);
                double R_dot = 0.5 * gam_dot * grid->delta_r(N_g);
                double dRdr = exp(2.0 * s->phi(N_g, j, k)) * sqrt(1.0 + s->h_rr(N_g, j, k));
                double alpha_n_nabla_R = R_dot - s->shift_r(N_g, j, k) * dRdr;
                // if (j == N_g && k == N_g) {
                //   cout << " tau_c = " << tau_c << endl;
                //   cout << " i = " << N_g << " R_dot = " << R_dot << " beta_term = "
                //        << s->shift_r(N_g,j,k)*dRdr
                //        << " alpha_n_nabla_R = " << alpha_n_nabla_R 
                //        << " dot h_rr = " << derivs->h_rr(0.0,j,k)
                //        << " shift_r = " << s->shift_r(0.0,j,k) << endl;
                // }
                // alpha_n_nabla_R = 0.0;
                // compute self-similar shift:
                double C = R / (tau_acc - tau_c);
                // NOTE: use central lapse in product with C...
                double shift_ss = -(s->lapse(0.0, j, k) * C + alpha_n_nabla_R) *
                    exp(-2.0 * s->phi(N_g, j, k)) / sqrt(1.0 + s->h_rr(N_g, j, k));
                // compute correction term:
                double R0 = 2 * (tau_acc - tau_c) * beta_max;
                double factor = R0 * R0 / (R0 * R0 + R * R);
                // and finally shift:
                s->shift_r[N_g][j][k] = shift_ss * factor;
                //
                // .. now go to all other points
                // 
                for (int i = N_g + 1; i < N_r; i++) {
                    const double r_mid = 0.5 * (grid->r(i) + grid->r(i - 1));
                    gamma_bar_rr = 1.0 + s->h_rr(r_mid, j, k);
                    dRdr_mid = exp(2.0 * s->phi(r_mid, j, k)) * sqrt(gamma_bar_rr);
                    R += dRdr_mid * grid->delta_r(i);
                    gam_dot = dRdr_mid * (2.0 * derivs->phi(r_mid, j, k) +
                        0.5 * derivs->h_rr(r_mid, j, k) / gamma_bar_rr);
                    R_dot += gam_dot * grid->delta_r(i);
                    dRdr = exp(2.0 * s->phi(i, j, k)) * sqrt(1.0 + s->h_rr(i, j, k));
                    alpha_n_nabla_R = (R_dot - s->shift_r(i, j, k) * dRdr);
                    if (i == N_g + 1 && j == N_g && k == N_g)
                        cout << " i = " << i << " R_dot = " << R_dot << " beta_term = "
                        << s->shift_r(i, j, k) * dRdr
                        << " alpha_n_nabla_R = " << alpha_n_nabla_R
                        << " dot h_rr = " << derivs->h_rr(r_mid, j, k) << endl;
                    // alpha_n_nabla_R = 0.0;
                    // compute self-similar shift:
                    C = R / (tau_acc - tau_c);
                    shift_ss = -(s->lapse(0.0, j, k) * C + alpha_n_nabla_R) *
                        exp(-2.0 * s->phi(i, j, k)) / sqrt(1.0 + s->h_rr(i, j, k));
                    // compute correction term:
                    factor = R0 * R0 / (R0 * R0 + R * R);
                    // and finally shift:
                    s->shift_r[i][j][k] = shift_ss * factor;
                };
            };
        s->shift_r.fill_ghosts();
    };
};

#endif /* GAUGE_H */
