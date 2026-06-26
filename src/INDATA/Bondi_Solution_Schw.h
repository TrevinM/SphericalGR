//====================================================
//
// Bondi_Solution_Schw.h
//
//====================================================
#ifndef BONDI_SOLUTION_SCHW_H
#define BONDI_SOLUTION_SCHW_H

#include "Bondi_Solution.h"

class Bondi_Solution_Schw : public Bondi_Solution {
public:
    //================================================
    // Constructor
    //================================================
 Bondi_Solution_Schw(const double M_dot_in, const double R_crit_in,
		     const double M_in, const double Kappa_in, const double Gamma_in) :
  Bondi_Solution(M_dot_in, R_crit_in, M_in, Kappa_in, Gamma_in) {}
    //================================================
    // Destructor
    //================================================
    ~Bondi_Solution_Schw() {}
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
    // Conformal factor (really only applies to isotropic coords...)
    //================================================
    double psi(const double r_in) {
        return 1.0 + M / (2.0 * r_in);
    }
    //================================================
    // Lapse
    //================================================
    double lapse(const double r_in) {
        const double R = R_of_r(r_in);
        return sqrt(1.0 + 2.0 * M / R);
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
    // Rest-mass density of fluid
    //================================================
    double rho_0(const double r_in) {
        return rho_0_gen(r_in);
    }
    //================================================
    // Four-velocity of fluid
    //================================================
    double u_0(const double r_in) {
        return u_0_gen(r_in);
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

#endif /* _Bondi_Schw_H_ */
