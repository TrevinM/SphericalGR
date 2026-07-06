//====================================================
//
// Bondi_Solution_Iso.h
//
//====================================================
#ifndef BONDI_SOLUTION_ISO_H
#define BONDI_SOLUTION_ISO_H

#include "Bondi_Solution.h"

class Bondi_Solution_Iso : public Bondi_Solution {
private:
    bool use_artificial;
public:
    //================================================
    // Constructor
    //================================================
    Bondi_Solution_Iso(const double M_dot_in, const double R_crit_in, const double M_in,
		       const double Kappa_in, const double Gamma_in,
		       const bool use_artificial_in) :
    Bondi_Solution(M_dot_in, R_crit_in, M_in, Kappa_in, Gamma_in), use_artificial(use_artificial_in) {}
    //================================================
    // Destructor
    //================================================
    ~Bondi_Solution_Iso() {}
    //================================================
    // Methods to calculate fluid parameters as functions of
    // isotropic radius
    //================================================
    // Conversions between areal and isotropic radii
    //================================================
    double r_of_R(const double R_in) {
        return (R_in - M + sqrt(R_in * (R_in - 2.0 * M))) / 2.0;
    }
    double R_of_r(const double r_in) {
        const double psi_l = psi(r_in);
        return r_in * psi_l * psi_l;
    }
    //================================================
    // Conformal factor
    //================================================
    double psi(const double r_in) {
        return 1.0 + M / (2.0 * r_in);
    }
    //================================================
    // Lapse
    //================================================
    double lapse(const double r_in) {
        return (1.0 - M / (2.0 * r_in)) / (1.0 + M / (2.0 * r_in));
        double psi_l = psi(r_in);
        return 1.0 / (psi_l * psi_l);
    }
    //================================================
    // Shift
    //================================================
    double shift_r(const double r_in) {
        return 0.0;
    }
    double shift_t(const double r_in) {
        return 0.0;
    }
    double shift_p(const double r_in) {
        return 0.0;
    }
    //================================================
    // Deviation from conformally related spatial metric
    //================================================
    double h_rr(const double r_in) {
        return 0.0;
    }
    double h_rt(const double r_in) {
        return 0.0;
    }
    double h_rp(const double r_in) {
        return 0.0;
    }
    double h_tt(const double r_in) {
        return 0.0;
    }
    double h_tp(const double r_in) {
        return 0.0;
    }
    double h_pp(const double r_in) {
        return 0.0;
    }
    //================================================
    // Extrinsic curvature
    //================================================
    double a_rr(const double r_in) {
        return 0.0;
    }
    double a_rt(const double r_in) {
        return 0.0;
    }
    double a_rp(const double r_in) {
        return 0.0;
    }
    double a_tt(const double r_in) {
        return 0.0;
    }
    double a_tp(const double r_in) {
        return 0.0;
    }
    double a_pp(const double r_in) {
        return 0.0;
    }
    double K(const double r_in) {
        return 0.0;
    }
    //================================================
    // Rest-mass density of fluid (artificial approach for r < M
    // borrowed from Faber et al. 2007)
    //================================================
    double rho_0(const double r_in) {
        if (r_in < M && use_artificial == true) {
            const double eps = rho_0_gen(10.0 * M);
            const double a = rho_0_deriv(M);
            const double b = rho_0_gen(M);
            const double c = (b - eps) / 2.0 - a * M / 8.0;
            if (r_in < M / 2.0) {
                return c * (1.0 - cos(2.0 * PI * r_in / M)) + eps;
            } else {
                return a / M * r_in * r_in - a * r_in + b;
            }
        } else {
            return rho_0_gen(r_in);
        }
    }
    //================================================
    // Derivative of rest-mass density d(rho_0)/dr
    //================================================
    double rho_0_deriv(const double r_in) {
      cout << " Bondi_Solution_Iso NOT implemented... " << endl;
      exit(1);       
      /*        const double R = R_of_r(r_in);
        const double alpha_l = alpha(r_in);
        const double alpha6 = alpha_l * alpha_l * alpha_l * alpha_l * alpha_l * alpha_l;
        const double coef = 3.0 * alpha_l * alpha_l;
        const double term1 = 2.0 * M / (R * R);
        const double term2 = M_dot * M_dot / (4.0 * PI * PI * R * R * R * R * R * alpha6);
        const double term3 = 3.0 * M_dot * M_dot / (8.0 * PI * PI * R * R * R * R * alpha6 * alpha_l);
        const double term4 = 8.0 * Kappa * B_const_l * B_const_l / pow(1.0 + 4.0 * Kappa * alpha_l, 3);
        return coef * (term1 - term2) / (term3 - term4);
      */
    }
    //================================================
    // Four-velocity of fluid
    //================================================
    //
    // Time component
    //
    double u_0(const double r_in) {
        return u_0_gen(r_in);
    }
    //
    // General expression for radial component of four-velocity
    //
    double u_r_gen(const double r_in) {
        const double u_r_schw_l = u_r_schw(r_in);
        return u_r_schw_l / ((1.0 - M / (2.0 * r_in)) * (1.0 + M / (2.0 * r_in)));
    }
    //
    // Modify general expression for u^r inside r = M
    // (method borrowed from Faber et al. 2007)
    //
    double u_r(const double r_in) {
        if (r_in < M && use_artificial == true) {
            return u_r_gen(M) * r_in / M;
        } else {
            return u_r_gen(r_in);
        }
    }
    //
    // Theta and phi components
    //
    double u_t(const double r_in) {
        return 0.0;
    }
    double u_p(const double r_in) {
        return 0.0;
    }
};

#endif /* BONDI_SOLUTION_ISO_H */
