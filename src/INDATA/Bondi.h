// Tell emacs that this is -*-c++-*- mode
//====================================================
// Bondi initial data, centered on origin
//====================================================
//
#include "Bondi_Solution.h"
#include "Bondi_Solution_Schw.h"
#include "Bondi_Solution_Iso.h"
#include "Bondi_Solution_Ken.h"
#include "Bondi_Solution_Max.h"

class Bondi : public InData {
private:
    int coord_type;
    double M_dot;
    double R_crit_factor, R_crit;
    double a_inf, v_inf, a_s, shield;
    double M;
    //    double x_C, y_C, z_C;
    Bondi_Solution* bondi_sol;
    bool all_clear;
    double PI;
public:
    //================================================
    // Constructor
    //================================================
    Bondi(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
        InData(grid_i, cosmology) {
        indata_type = bondi;
        analytical = true;
        all_clear = true;
        ifstream infile;
        infile.open(indata_input);
        if (!infile) {
            cerr << " BONDI: Can't open " << indata_input << " for input. This is bad." << endl;
            all_clear = false;
        } else {
            cout << " BONDI: Reading initial data parameters from file " << indata_input << endl;
        }
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> coord_type;
        infile.get(buf, 100, '='); infile.get(c); infile >> Gamma;
        infile.get(buf, 100, '='); infile.get(c); infile >> Kappa;
        infile.get(buf, 100, '='); infile.get(c); infile >> a_inf;
        infile.get(buf, 100, '='); infile.get(c); infile >> v_inf;
        infile.get(buf, 100, '='); infile.get(c); infile >> shield;
        infile.get(buf, 100, '='); infile.get(c); infile >> M;
        //        infile.get(buf, 100, '='); infile.get(c); infile >> x_C;
        //        infile.get(buf, 100, '='); infile.get(c); infile >> y_C;
        //        infile.get(buf, 100, '='); infile.get(c); infile >> z_C;
        string coord_str;
        if (coord_type == 1) coord_str = "Schwarzschild";
        else if (coord_type == 2) coord_str = "isotropic Schwarzschild";
        else if (coord_type == 3) coord_str = "Ken's trumpet";
        else if (coord_type == 4) coord_str = "maximal trumpet";
        R_crit = R_crit_factor * M;
        cout << " BONDI: Will set up Bondi flow in " << coord_str << " coords with" << endl;
        cout << " BONDI:   Gamma = " << Gamma << endl;
        cout << " BONDI:   Kappa = " << Kappa << endl;
        cout << " BONDI:   a_inf = " << a_inf << endl;
        cout << " BONDI:   v_inf = " << v_inf << endl;
        cout << " BONDI:  shield = " << shield << endl;
        cout << " BONDI: around a Schwarzschild black hole of mass M = " << M << endl;
        cout << "===================================================" << endl;
        PI = acos(-1.0);
    }
    //================================================
    // Destructor
    //================================================
    ~Bondi() { delete bondi_sol; }
    string Name() { return "Bondi initial data"; }
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) {
        //
        // Follow Richards, Baumgarte & Shapiro, MNRAS, 2021, to construct Bondi solution... 
        //
        // Eqs. (19)
        double a_inf2 = a_inf * a_inf;
        double A = (7.0 - 6.0 * Gamma) / 3.0;
        double B = (1.0 - Gamma) * (5.0 - 3.0 * Gamma) / 3.0;
        double C = a_inf2 * (2.0 * Gamma - 2.0 - a_inf2) / 3.0;
        // Eq. (24)
        double Q = (A * A - 3.0 * B) / 9.0;
        // Eq. (25)
        double R = (2.0 * A * A * A - 9.0 * A * B + 27.0 * C) / 54.0;
        // Eq. (29)
        double Theta = acos(R / (Q * sqrt(Q)));
        // Eq. (28) (third root...)
        double a_s2 = -2.0 * sqrt(Q) * cos((Theta - 2.0 * PI) / 3.0) - A / 3.0;
        a_s = sqrt(a_s2);
        //
        // now compute lambda_GR from Eq. (12)
        //
        double term1 = pow(a_s / a_inf, (5.0 - 3.0 * Gamma) / (Gamma - 1.0));
        double term2 = pow((Gamma - 1.0 - a_inf2) / (Gamma - 1.0 - a_s2), 1.0 / (Gamma - 1.0));
        double term3 = pow(1.0 + 3.0 * a_s2, 1.5) / 4.0;
        double lambda_GR = term1 * term2 * term3;
        cout << " BONDI: found a_s^2 = " << a_s2 << " and lambda_GR = " << lambda_GR << endl;
        //
        // now compute M_dot...
        //
        double term4 = a_inf2 / (1.0 - a_inf2 / (Gamma - 1.0));
        double rho_0_inf = pow(term4 / (Kappa * Gamma), 1.0 / (Gamma - 1.0));
        M_dot = 4.0 * PI * lambda_GR * M * M * rho_0_inf / (a_inf * a_inf2);
        cout << " BONDI: found M_dot = " << M_dot << ", rho_0_inf = " << rho_0_inf << endl;
        //
        // now compute R_crit from Shapiro & Teukolsky, 1983, App. G
        //
        double u_s2 = a_s2 / (1 + 3.0 * a_s2);
        double R_crit = M / (2.0 * u_s2);
        cout << " BONDI: found R_crit = " << R_crit << endl;
        //
        // finally set up Bondi classes:
        //
        //============================================
        // Choose derived class based on coord type flag
        // (Note: Bondi_Solution classes assume Gamma = 4/3)
        //============================================
        // For isotropic Schwarzschild: use artificial initial data inside r = M?
        const bool use_artificial = true;
        if (coord_type == 1) bondi_sol = new Bondi_Solution_Schw(M_dot, R_crit, M, Kappa, Gamma);
        else if (coord_type == 2) bondi_sol = new Bondi_Solution_Iso(M_dot, R_crit, M, Kappa, Gamma, use_artificial);
        else if (coord_type == 3) bondi_sol = new Bondi_Solution_Ken(M_dot, R_crit, M, Kappa, Gamma);
        else if (coord_type == 4) bondi_sol = new Bondi_Solution_Max(M_dot, R_crit, M, Kappa, Gamma);
        return all_clear;
    };
    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->h_rr(r);
    }
    double h_rt_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->h_rt(r);
    }
    double h_rp_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->h_rp(r);
    }
    double h_tt_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->h_tt(r);
    }
    double h_tp_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->h_tp(r);
    }
    double h_pp_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->h_pp(r);
    }
    double phi_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->phi(r);
    }
    //================================================
    // Analytical solution for connection coefficients
    //================================================
    double lam_r_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = false; return 0;
    }
    double lam_t_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = false; return 0;
    }
    double lam_p_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = false; return 0;
    }
    //================================================
    // Analytical solution for extrinsic curvature
    //================================================
    double a_rr_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->a_rr(r);
    }
    double a_rt_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->a_rt(r);
    }
    double a_rp_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->a_rp(r);
    }
    double a_tt_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->a_tt(r);
    }
    double a_tp_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->a_tp(r);
    }
    double a_pp_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->a_pp(r);
    }
    double K_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->K(r);
    }
    //================================================
    // Analytical solution for gauge
    //================================================
    double lapse_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->lapse(r);
    }
    double shift_r_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->shift_r(r);
    }
    double shift_t_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->shift_t(r);
    }
    double shift_p_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->shift_p(r);
    }
    //================================================
    // Analytical solution for auxiliary functions
    //================================================
    double Theta_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double B_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double B_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double B_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    //================================================
    // Analytical solution for hydro
    //================================================
    double rho_0_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->rho_0(r);
    }
    double P_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double v_r_analytical(double r, double theta, double phi, double t) {
        // a sanity check: make sure that v^2 < 1.0 ...
        const double rho = r / (shield * M);
        const double v_r = bondi_sol->v_norm_r(r) - tanh(rho) * v_inf * cos(theta);
        const double v_t = v_t_analytical(r, theta, phi, t);
        const double psi4 = exp(4.0 * phi_analytical(r, theta, phi, t));
        const double v2 = psi4 * (v_r * v_r + v_t * v_t); // conformal flatness...
        if (v2 >= 0.95) {
            cout << " BONDI: v2 >= 1.0 in initial data." << endl;
            cout << "     Either choose v_inf smaller, or variable shield larger"
                << endl;
            exit(1);
        }
        return v_r;
    }
    double v_t_analytical(double r, double theta, double phi, double t) {
        const double rho = r / (shield * M);
        return bondi_sol->v_norm_t(r) + tanh(rho) * v_inf * sin(theta);
    }
    double v_p_analytical(double r, double theta, double phi, double t) {
        return bondi_sol->v_norm_p(r);
    }
    //================================================
    // Analytical solution for scalar field
    //================================================
    double sf_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double pi_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for Maxwell OR for Dual Maxwell:
    // For Maxwell need e_p and a_p only, for Dual Maxwell all a_i and as_i
    //================================================
    double e_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double as_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double as_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double as_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for radiation
    //================================================
    double E_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double F_0_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double F_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double F_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double F_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
};
