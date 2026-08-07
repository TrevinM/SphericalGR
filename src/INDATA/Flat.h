// Tell emacs that this is -*-c++-*- mode
//================================================
// Flat initial data in funny coordinate system
//================================================

class Flat : public InData {
private:
    double A, B, C; // parameters for functions f(x), g(y), h(z)
public:
    //================================================
    // Constructor
    //================================================
    Flat(char* indata_input, Grid* grid_i, Cosmology* cosmology)
        : InData(grid_i, cosmology) {
        indata_type = flat;
        analytical = true;
        A = 0.1;  // parameters for functions f(x) etc.
        B = 0.3;
        C = 0.5;
        cout << " Will set up flat initial data, for testing purposes only. " << endl;
    };
    //================================================
    // Destructor
    //================================================
    ~Flat();
    string Name() { return "flat in funny coordinates"; };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) { return true; }
    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        // create flat metric in (non-trivial) cartesian coordinate system
        tensor g_cart(f(x), 0.0, 0.0, g(y), 0.0, h(z));
        // transform into spherical polar coordinates
        tensor g = Cartesian_to_Spherical(g_cart, r, theta, phi);
        // return (rr) component of h
        return g[0][0] - 1.0;
    };
    double h_rt_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        // create flat metric in (non-trivial) cartesian coordinate system
        tensor g_cart(f(x), 0.0, 0.0, g(y), 0.0, h(z));
        // transform into spherical polar coordinates
        tensor g = Cartesian_to_Spherical(g_cart, r, theta, phi);
        // return (rt) component of h
        return g[0][1] / r;
    };
    double h_rp_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        // create flat metric in (non-trivial) cartesian coordinate system
        tensor g_cart(f(x), 0.0, 0.0, g(y), 0.0, h(z));
        // transform into spherical polar coordinates
        tensor g = Cartesian_to_Spherical(g_cart, r, theta, phi);
        // return (rp) component of h
        return g[0][2] / (r * sin(theta));
    };
    double h_tt_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        // create flat metric in (non-trivial) cartesian coordinate system
        tensor g_cart(f(x), 0.0, 0.0, g(y), 0.0, h(z));
        // transform into spherical polar coordinates
        tensor g = Cartesian_to_Spherical(g_cart, r, theta, phi);
        // return (tt) component of h
        return g[1][1] / (r * r) - 1.0;
    };
    double h_tp_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        // create flat metric in (non-trivial) cartesian coordinate system
        tensor g_cart(f(x), 0.0, 0.0, g(y), 0.0, h(z));
        // transform into spherical polar coordinates
        tensor g = Cartesian_to_Spherical(g_cart, r, theta, phi);
        // return (tp) component of h
        return g[1][2] / (r * r * sin(theta));
    };
    double h_pp_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        // create flat metric in (non-trivial) cartesian coordinate system
        tensor g_cart(f(x), 0.0, 0.0, g(y), 0.0, h(z));
        // transform into spherical polar coordinates
        tensor g = Cartesian_to_Spherical(g_cart, r, theta, phi);
        // return (pp) component of h
        return g[2][2] / (r * r * sin(theta) * sin(theta)) - 1.0;
    };
    double phi_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
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
        return 0.0;
    };
    double a_rt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_rp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_tt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_tp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_pp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double K_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for gauge
    //================================================
    double lapse_analytical(double r, double theta, double phi, double t) {
        return 1.0;
    };
    double shift_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double shift_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double shift_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
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
        return 0.0;
    };
    double P_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double v_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double v_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double v_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
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
private:
    //================================================
    // Auxiliary functions f(x) = \gamma_{xx} etc.
    //================================================
    double f(double x) { return 1.0 + A * x * x + 0.2; }
    double g(double y) { return 1.0 + B * y * y + 0.5; }
    double h(double z) { return 1.0 + C * z * z + 0.1; }

};

//
