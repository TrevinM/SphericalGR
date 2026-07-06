//====================================================
//
// Bondi_Solution_Ken.h
//
//====================================================
#ifndef BONDI_SOLUTION_KEN_H
#define BONDI_SOLUTION_KEN_H

#include "Bondi_Solution.h"

class Bondi_Solution_Ken : public Bondi_Solution {
private:
    double R_0;
public:
    //================================================
    // Constructor
    //================================================
 Bondi_Solution_Ken(const double M_dot_in, const double R_crit_in, const double M_in,
		    const double Kappa_in, const double Gamma_in) :
    Bondi_Solution(M_dot_in, R_crit_in, M_in, Kappa_in, Gamma_in) {R_0 = M;}
    //================================================
    // Destructor
    //================================================
    ~Bondi_Solution_Ken() {}
    //================================================
    // Methods to calculate fluid parameters as functions of
    // isotropic radius
    //================================================
    // Function f1 that enters into spacetime metric;
    // defined here for convenience
    //================================================
    double f1(const double r_in) {
        return sqrt(2.0 * r_in * (M - R_0) + R_0 * (2.0 * M - R_0));
    }
    //================================================
    // Conversions between areal and isotropic radii
    //================================================
    double r_of_R(const double R_in) {
        return R_in - R_0;
    }
    double R_of_r(const double r_in) {
        return r_in + R_0;
    }
    //================================================
    // Conformal factor
    //================================================
    double psi(const double r_in) {
        return sqrt(1.0 + R_0 / r_in);
    }
    //================================================
    // Lapse
    //================================================
    double lapse(const double r_in) {
        return r_in / (r_in + R_0);
    }
    //================================================
    // Shift
    //================================================
    double shift_r(const double r_in) {
        const double f1_l = f1(r_in);
        return f1_l * r_in / ((r_in + R_0) * (r_in + R_0));
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
        double f1_l = f1(r_in);
        return -2.0 / 3.0 * ((3.0 * r_in + R_0) * (M - R_0) + 2.0 * M * R_0) /
               (f1_l * (r_in + R_0) * (r_in + R_0));
    }
    double a_rt(const double r_in) {
        return 0.0;
    }
    double a_rp(const double r_in) {
        return 0.0;
    }
    double a_tt(const double r_in) {
        double f1_l = f1(r_in);
        return (r_in * (M - R_0) + R_0 * (M - R_0 / 3.0)) /
               (f1_l * (r_in + R_0) * (r_in + R_0));
    }
    double a_tp(const double r_in) {
        return 0.0;
    }
    double a_pp(const double r_in) {
        double f1_l = f1(r_in);
        return (r_in * (M - R_0) + R_0 * (M - R_0 / 3.0)) /
               (f1_l * (r_in + R_0) * (r_in + R_0));
    }
    double K(const double r_in) {
        double f1_l = f1(r_in);
        return ((3.0 * r_in + 2.0 * R_0) * (M - R_0) + M * R_0) /
               (f1_l * (r_in + R_0) * (r_in + R_0));
    }
    //================================================
    // Rest-mass density of fluid
    //================================================
    double rho_0(const double r_in) {
        return rho_0_gen(r_in);
    }
    //================================================
    // Four-velocity of fluid
    //================================================
    double u_0(const double r_in) {
        const double R = R_of_r(r_in);
        const double f1_l = f1(r_in);
        const double u_r_l = u_r(r_in);
        if (fabs(R - R_hor) < span) {
            //
            // Near horizon, compute u_0 using expanded expression
            // to avoid numerical error introduced by 1/0 terms
            //
            cout << " Warning: within " << span << " of horizon" << endl;
            cout << " Computing u^0 using expansion..." << endl;
            return -(r_in + R_0) * (r_in + R_0) * u_r_l / (r_in * (r_in + f1_l)) - 1.0 / (2.0 * u_r_l);
        } else {
            //
            // Away from horizon, use general expression
            //
            return u_0_gen(r_in);
        }
    }
    double u_r(const double r_in) {
        return u_r_schw(r_in);
    }
    double u_t(const double r_in) {
        return 0.0;
    }
    double u_p(const double r_in) {
        return 0.0;
    }
};

#endif /* BONDI_SOLUTION_KEN_H */
