//====================================================
//
// Bondi_Solution.h
//
//====================================================
#ifndef BONDI_SOLUTION_H
#define BONDI_SOLUTION_H

#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"
#include "SBN.h"

//
// structure needed for ODE integration...
// 

struct Bondi_derivatives {
    double Gamma;
    double kappa;
    double M;
    Bondi_derivatives(double M_in, double Gamma_in, double kappa_in) :
        M(M_in), Gamma(Gamma_in), kappa(kappa_in) {
    };
    void operator() (const double r, const VecDoub_I& y, VecDoub_O& dydr) {
        // unpack variables
        const double u = y[0];
        const double rho = y[1];
        const double u2 = u * u;
        const double temp1 = Gamma * kappa * pow(rho, Gamma - 1.0);
        const double a2 = temp1 / (1.0 + temp1 / (Gamma - 1.0));
        const double temp2 = 1.0 - 2.0 * M / r + u2;
        const double r2 = r * r;
        // see Eqs. (G.10)-(G.13) in "Compact Objects"
        const double D1 = (2.0 * temp2 * a2 / r - M / r2) / rho;
        const double D2 = (2.0 * u2 / r - M / r2) / u;
        const double D = (u2 - temp2 * a2) / (rho * u);
        if (D == 0.0) { // hack...
            dydr[0] = 0.0;
            dydr[1] = 0.0;
        } else {
            dydr[0] = D1 / D;
            dydr[1] = -D2 / D;
        }
    };
};

class Bondi_Solution {
protected:
    double PI;
    double R_hor;
    double x1_init, x2_init;
    double xmin;
    double xacc;
    double span;
    double Kappa;
    double Gamma;
    double M_dot, R_crit, M;
    double u_crit_l;
    double a_crit_l;
    double rho_0_crit_l;
    //================================================
    // Fluid parameters at critical radius
    //================================================
    double u_crit(void) {
        return sqrt(M / (2.0 * R_crit));
    }
    double a_crit(void) {
        return u_crit_l / sqrt(1.0 - 3.0 * u_crit_l * u_crit_l);
    }
    double rho_0_crit(void) {
        return M_dot / (4.0 * PI * R_crit * R_crit * u_crit_l);
    }

    //================================================
    // Radial component of fluid four-velocity computed in
    // normal Schwarzschild coords (used to compute u^r
    // in other coord systems, so it makes sense to define it
    // in the base class)
    //================================================
    double u_r_schw(const double r_in) {
        const double R = R_of_r(r_in);
        const double rho_0_l = rho_0(r_in);
        return -M_dot / (4.0 * PI * R * R * rho_0_l);
    }
public:
    //================================================
    // Constructor
    //================================================
    Bondi_Solution(const double M_dot_in, const double R_crit_in, const double M_in,
        const double Kappa_in, const double Gamma_in) :
        M_dot(M_dot_in), R_crit(R_crit_in), M(M_in), Kappa(Kappa_in), Gamma(Gamma_in) {
        //
             // Initialize variables, constants
             //
        PI = acos(-1.0);
        R_hor = 2.0 * M;
        x1_init = 0.0;
        x2_init = 1.0;
        span = 1.0e-6 * M;
        xmin = 0.0;
        xacc = 1.0e-12;
        //
        // Calculate values of fluid parameters at critical radius
        //
        u_crit_l = u_crit();
        a_crit_l = a_crit();
        rho_0_crit_l = rho_0_crit();
        // B_const_l = B_const();
        cout << " Kappa = " << Kappa << endl;
        cout << " Gamma = " << Gamma << endl;
        cout << " R_crit = " << R_crit << endl;
        cout << " u_crit = " << u_crit_l << endl;
        cout << " a_crit = " << a_crit_l << endl;
        cout << " rho_0_crit = " << rho_0_crit_l << endl;
    }
    //================================================
    // Destructor
    //================================================
    virtual ~Bondi_Solution() {}
    //================================================
    // Virtual functions (will be overwritten in derived classes)
    //================================================
    virtual double r_of_R(const double) = 0;
    virtual double R_of_r(const double) = 0;
    virtual double psi(const double) = 0;
    virtual double lapse(const double) = 0;
    virtual double shift_r(const double) = 0;
    virtual double shift_t(const double) = 0;
    virtual double shift_p(const double) = 0;
    virtual double rho_0(const double) = 0;
    virtual double u_0(const double) = 0;
    virtual double u_r(const double) = 0;
    virtual double u_t(const double) = 0;
    virtual double u_p(const double) = 0;
    virtual double h_rr(const double) = 0;
    virtual double h_rt(const double) = 0;
    virtual double h_rp(const double) = 0;
    virtual double h_tt(const double) = 0;
    virtual double h_tp(const double) = 0;
    virtual double h_pp(const double) = 0;
    virtual double a_rr(const double) = 0;
    virtual double a_rt(const double) = 0;
    virtual double a_rp(const double) = 0;
    virtual double a_tt(const double) = 0;
    virtual double a_tp(const double) = 0;
    virtual double a_pp(const double) = 0;
    virtual double K(const double) = 0;
    //================================================
    // Methods to calculate fluid parameters as functions of the
    // isotropic radius r (these are parameters that can be calculated
    // in the same way for all coordinate systems, and thus belong
    // in the base class)
    //================================================
    // phi = ln(psi)
    //================================================
    double phi(const double r_in) {
        return log(psi(r_in));
    }
    //================================================
    // Fluid rest-mass density
    //================================================
    //
    // USED TO BE: Intermediate variable alpha = rho_0^(1/3)
    // now: alpha = rho_0...
    //
    double alpha(const double r_in) {
        // find areal radius:
        const double R = R_of_r(r_in);
        // prepare for integration...
        VecDoub ystart(2);
        ystart[0] = -u_crit_l;
        ystart[1] = rho_0_crit_l;
        const Doub atol = 1.e-12, rtol = atol, h1 = 0.001, hmin = 0.0;
        Output out;
        Bondi_derivatives derivs(M, Gamma, Kappa);
        Odeint<StepperDopr5<Bondi_derivatives> > ode(ystart, R_crit, R,
            atol, rtol, h1, hmin, out,
            derivs);
        ode.integrate();
        return ystart[1];
    }
    //
    // General expression for rest-mass density rho_0
    //
    double rho_0_gen(const double r_in) {
        return alpha(r_in);
    }
    //================================================
    // General expression for time component of four-velocity
    // (computed in terms of conformal factor, lapse, and shift;
    // assumes shift_t = shift_p = 0)
    //================================================
    double u_0_gen(const double r_in) {
        const double psi_l = psi(r_in);
        const double psi4 = psi_l * psi_l * psi_l * psi_l;
        const double lapse_l = lapse(r_in);
        const double shift_r_l = shift_r(r_in);
        const double u_r_l = u_r(r_in);
        return (psi4 * shift_r_l * u_r_l +
            sqrt(psi4 * psi4 * shift_r_l * shift_r_l * u_r_l * u_r_l +
                (lapse_l * lapse_l - psi4 * shift_r_l * shift_r_l) *
                (psi4 * u_r_l * u_r_l + 1.0))) /
            (lapse_l * lapse_l - psi4 * shift_r_l * shift_r_l);
    }
    //================================================
    // Three-velocity of fluid
    //================================================
    double v_r(const double r_in) {
        return u_r(r_in) / u_0(r_in);
    }
    double v_t(const double r_in) {
        return u_t(r_in) / u_0(r_in);
    }
    double v_p(const double r_in) {
        return u_p(r_in) / u_0(r_in);
    }
    //================================================
    // Three-velocity of fluid as seen by a normal observer
    //================================================
    double v_norm_r(const double r_in) {
        double temp = (u_r(r_in) / u_0(r_in) + shift_r(r_in)) / lapse(r_in);
        if (!isfinite(temp)) {
            cout << " v_norm_r not finite..." << endl;
            cout << "  u_r = " << u_r(r_in)
                << "  u_0 = " << u_0(r_in)
                << "  lapse = " << lapse(r_in) << endl;
        }
        return (u_r(r_in) / u_0(r_in) + shift_r(r_in)) / lapse(r_in);
    }
    double v_norm_t(const double r_in) {
        return (u_t(r_in) / u_0(r_in) + shift_t(r_in)) / lapse(r_in);
    }
    double v_norm_p(const double r_in) {
        return (u_p(r_in) / u_0(r_in) + shift_p(r_in)) / lapse(r_in);
    }
    //================================================
    // Lorentz factor
    //================================================
    double W(const double r_in) {
        return lapse(r_in) * u_0(r_in);
    }
};

#endif /* BONDI_SOLUTION_H */
