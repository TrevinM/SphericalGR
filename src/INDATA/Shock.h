// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Initial data for shock wave - either spherical or planar
//================================================
//
class Shock : public InData {
private:
    double rho_in, rho_out;
    double R_i, x_C, y_C, z_C;
    int shape;
public:
    //================================================
    // Constructor
    //================================================
    Shock(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
        InData(grid_i, cosmology) {
        indata_type = shock;
        analytical = false;
        ifstream infile;
        infile.open(indata_input);
        if (!infile)
            cerr << "Can't open " << indata_input
            << " for input. This is bad. " << endl;
        else
            cout << " Reading initial data parameters from file "
            << indata_input << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> rho_in;
        infile.get(buf, 100, '='); infile.get(c); infile >> rho_out;
        infile.get(buf, 100, '='); infile.get(c); infile >> shape;
        infile.get(buf, 100, '='); infile.get(c); infile >> R_i;
        infile.get(buf, 100, '='); infile.get(c); infile >> x_C;
        infile.get(buf, 100, '='); infile.get(c); infile >> y_C;
        infile.get(buf, 100, '='); infile.get(c); infile >> z_C;
        if (shape == 0) {
            cout << " Will set up spherical shock initial data for: " << endl;
            cout << "    inner density = " << rho_in << ", outer density = " << rho_out << endl;
            cout << "    initial radius = " << R_i << endl;
            cout << "    centered on x = " << x_C << ", y = " << y_C << ", z = " << z_C << endl;
        } else if (shape == 1) {
            cout << " Will set up planar shock initial data for: " << endl;
            cout << "    upper density = " << rho_out << ", lower density = " << rho_in << endl;
            cout << "    separated at z = " << z_C << endl;
        } else {
            cerr << " Don't know shape = " << shape << ". This is bad. " << endl;
        };
        cout << "===================================================" << endl;
    };
    //================================================
    // Destructor
    //================================================
    ~Shock();
    string Name() {
        if (shape == 0)
            return "spherical shock tube initial data";
        else
            return "planar shock tube initial data";
    };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) { return true; }
    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_rt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_rp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_tt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_tp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_pp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double phi_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for connection coefficients
    //================================================
    double lam_r_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = true; return 0;
    }
    double lam_t_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = true; return 0;
    }
    double lam_p_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = true; return 0;
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
        if (shape == 0) {
            // spherical shock wave
            const double x_l = r * sin(theta) * cos(phi);
            const double y_l = r * sin(theta) * sin(phi);
            const double z_l = r * cos(theta);
            const double r_l = sqrt((x_l - x_C) * (x_l - x_C) + (y_l - y_C) * (y_l - y_C) + (z_l - z_C) * (z_l - z_C));
            if (r_l < R_i)
                return rho_in;
            else
                return rho_out;
        } else if (shape == 1) {
            // planar shock wave
            const double z_l = r * cos(theta);
            if (z_l > z_C)
                return rho_out;
            else
                return rho_in;
        };
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
};
