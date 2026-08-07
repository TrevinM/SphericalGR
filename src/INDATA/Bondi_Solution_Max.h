//====================================================
//
// Bondi_Solution_Max.h
//
//====================================================
#ifndef BONDI_SOLUTION_MAX_H
#define BONDI_SOLUTION_MAX_H

#include "Bondi_Solution.h"
#include "SBN.h"

class Bondi_Solution_Max : public Bondi_Solution {
private:
    SBN* sbn;
    double C;
public:
    //================================================
    // Constructor
    //================================================
    Bondi_Solution_Max(const double M_dot_in, const double R_crit_in, const double M_in,
        const double Kappa_in, const double Gamma_in) :
        Bondi_Solution(M_dot_in, R_crit_in, M_in, Kappa_in, Gamma_in) {
        //
        // Create instance of SBN class
        //
        const int array_length = 1e6;
        const double R_max = 500.0;
        sbn = new SBN(R_max, array_length);
        //
        // Calculate constant C
        //
        C = 3.0 * sqrt(3.0) * M * M / 4.0;
    }
    //================================================
    // Destructor
    //================================================
    ~Bondi_Solution_Max() { delete sbn; }
    //================================================
    // Methods to calculate fluid parameters as functions of
    // isotropic radius
    //================================================
    // Conversions between areal and isotropic radii
    //================================================
    double r_of_R(const double R_in) {
        return (2.0 * R_in + M + sqrt(4.0 * R_in * R_in + 4.0 * M * R_in + 3.0 * M * M)) / 4.0 *
            pow((4.0 + 3.0 * sqrt(2.0)) * (2.0 * R_in - 3.0 * M) /
                (8.0 * R_in + 6.0 * M + 3.0 * sqrt(8.0 * R_in * R_in + 8.0 * M * R_in + 6.0 * M * M)),
                1.0 / sqrt(2.0));
    }
    double R_of_r(const double r_in) {
        return sbn->areal(r_in, M);
    }
    //================================================
    // Conformal factor
    //================================================
    double psi(const double r_in) {
        return sbn->conFactor(r_in, M);
    }
    //================================================
    // Lapse
    //================================================
    double lapse(const double r_in) {
        return sbn->lapse(r_in, M);
    }
    //================================================
    // Shift
    //================================================
    double shift_r(const double r_in) {
        return sbn->shift(r_in, M);
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
        const double R = R_of_r(r_in);
        return -2.0 * C / (R * R * R);
    }
    double a_rt(const double r_in) {
        return 0.0;
    }
    double a_rp(const double r_in) {
        return 0.0;
    }
    double a_tt(const double r_in) {
        const double R = R_of_r(r_in);
        return C / (R * R * R);
    }
    double a_tp(const double r_in) {
        return 0.0;
    }
    double a_pp(const double r_in) {
        const double R = R_of_r(r_in);
        return C / (R * R * R);
    }
    double K(const double r_in) {
        return 0.0;
    }
    //================================================
    // Function f used in velocity transformation
    //================================================
    double f(const double r_in) {
        const double R = R_of_r(r_in);
        return sqrt(1.0 - 2.0 * M / R + C * C / (R * R * R * R));
    }
    //================================================
    // Rest-mass density of fluid
    //================================================
    double rho_0(const double r_in) {
        // double rho_0 = rho_0_gen(r_in);
        // if (r_in < 0.001) cout << " rho at r = " << r_in << " = " << rho_0 << endl;
        return rho_0_gen(r_in);
    }
    //================================================
    // Four-velocity of fluid
    //================================================
    double u_0(const double r_in) {
        const double R = R_of_r(r_in);
        const double psi_l = psi(r_in);
        const double f_l = f(r_in);
        const double u_r_l = u_r(r_in);
        if (fabs(R - R_hor) < span) {
            //
            // Near horizon, compute u_0 using expanded expression
            // to avoid numerical error introduced by 1/0 terms
            //
            cout << " Warning: within " << span << " of horizon" << endl;
            cout << " Computing u^0 using expanded expression..." << endl;
            return -psi_l * psi_l * r_in * r_in / (2.0 * C * u_r_l) *
                (psi_l * psi_l * psi_l * psi_l * u_r_l * u_r_l + 1.0);
        } else {
            //
            // Away from horizon, use general expression
            //
            return u_0_gen(r_in);
        }
    }
    double u_r(const double r_in) {
        double psi_l = psi(r_in);
        double f_l = f(r_in);
        double u_r_schw_l = u_r_schw(r_in);
        return u_r_schw_l / (psi_l * psi_l * f_l);
    }
    double u_t(const double r_in) {
        return 0.0;
    }
    double u_p(const double r_in) {
        return 0.0;
    }
};

#endif /* BONDI_SOLUTION_MAX_H */
