// Tell emacs that this is -*-c++-*- mode
//
//================================================
// RNS initial data, read in from file output.dat_UniRot_S*
//================================================
//
class RNS : public InData {
private:
    int NR;
    int NTheta;
    MatDoub* rho_data;
    MatDoub* v_phi_up_data;
    MatDoub* lapse_data;
    MatDoub* beta_phi_up_data;
    MatDoub* w_data;
    MatDoub* aa_data;
    MatDoub* bb_data;
    MatDoub* eps_data;
    double* r_data;
    double* theta_data;
    double delta_r;
    double delta_theta;
public:
    //================================================
    // Constructor
    //================================================
    RNS(char* data_file, Grid* grid_i, Cosmology* cosmology)
        : InData(grid_i, cosmology) {
        indata_type = rns;
        analytical = true;
        Gamma = 2.0;
        Kappa = 2180866104449.0630;  // according to data file - will later rescale to K = 1
        double sqrtKappa = sqrt(Kappa);
        // parameters for grid
        read_from_file = true;
        //    r_max = 5748581.5/sqrtKappa;   // slightly larger than correct value, to make sure that divisions work...
        //    nr = 98;
        ifstream infile;
        infile.open(data_file);
        if (!infile)
            cerr << "Can't open " << data_file
            << " for input. This is bad. " << endl;
        else
            cout << " Reading RNS data from file "
            << data_file << endl;
        // first read NR and Ntheta
        infile >> NR;
        infile >> NTheta;
        cout << " Reading data for NR = " << NR << " and NTheta = "
            << NTheta << " gridpoints." << endl;
        nr = NR - 2;  // NR include outer ghost zones
        ntheta = 2 * NTheta;  // no equatorial symmetry in code...
        // if (!strcmp(data_file,"output.dat_UniRot_S1")) {
        //   NTheta = 32;
        //   ntheta = 64;
        // } else if (!strcmp(data_file,"output.dat_UniRot_S2")) {
        //   NTheta = 64;
        //   ntheta = 128;
        // } else 
        //   cerr << " I don't know this data file.  This is bad. " << endl;
        //
        // allocate arrays
        //
        r_data = new double[NR];
        theta_data = new double[NTheta];
        rho_data = new MatDoub(NR, NTheta);        // rest-mass density!
        v_phi_up_data = new MatDoub(NR, NTheta);
        lapse_data = new MatDoub(NR, NTheta);
        beta_phi_up_data = new MatDoub(NR, NTheta);
        w_data = new MatDoub(NR, NTheta);          // gamma factor
        aa_data = new MatDoub(NR, NTheta);
        bb_data = new MatDoub(NR, NTheta);
        eps_data = new MatDoub(NR, NTheta);
        //
        // read data into arrays
        //
        for (int i = 0; i < NR; i++) {
            for (int j = 0; j < NTheta; j++) {
                infile >> r_data[i];
                infile >> theta_data[j];
                infile >> (*rho_data)[i][j];
                infile >> (*v_phi_up_data)[i][j];
                infile >> (*lapse_data)[i][j];
                infile >> (*beta_phi_up_data)[i][j];
                infile >> (*w_data)[i][j];
                infile >> (*aa_data)[i][j];
                infile >> (*bb_data)[i][j];
                infile >> (*eps_data)[i][j];
            }
            // cout << " r(" << i << ") = " << r_data[i] << ", lapse(r) = " <<
            // 	(*lapse_data)[i][2] << endl;
        }
        //
        // r_max...
        //
        r_max = 0.5 * (r_data[NR - 2] + r_data[NR - 3]);
        //
        // now rescale to dimensionless units:
        //
        r_max /= sqrtKappa;
        for (int i = 0; i < NR; i++) {
            r_data[i] /= sqrtKappa;
            for (int j = 0; j < NTheta; j++) {
                (*rho_data)[i][j] *= Kappa;
            }
        }
        Kappa = 1.0;
        delta_r = r_data[1] - r_data[0];
        delta_theta = theta_data[1] - theta_data[0];
        // make deltas slightly smaller, so that we compute correct j's below...
        delta_r *= 0.999999999;
        delta_theta *= 0.999999999;
    };
    //================================================
    // Destructor
    //================================================
    ~RNS() {
        delete r_data; delete theta_data; delete rho_data; delete v_phi_up_data; delete lapse_data;
        delete beta_phi_up_data; delete w_data; delete aa_data; delete bb_data; delete eps_data;
    };
    string Name() { return "Rotating neutron star initial data"; };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) { return true; }
    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double r, double theta, double phi, double t) {
        int i = abs(r) / delta_r - 0.5;
        if (theta > PI / 2.0) theta = PI - theta;
        int j = abs(theta) / delta_theta - 0.5;
        return (*aa_data)[i][j] - 1.0;
    }
    double h_rt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_rp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_tt_analytical(double r, double theta, double phi, double t) {
        int i = abs(r) / delta_r - 0.5;
        if (theta > PI / 2.0) theta = PI - theta;
        int j = abs(theta) / delta_theta - 0.5;
        return (*aa_data)[i][j] - 1.0;
    }
    double h_tp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_pp_analytical(double r, double theta, double phi, double t) {
        int i = abs(r) / delta_r - 0.5;
        if (theta > PI / 2.0) theta = PI - theta;
        int j = abs(theta) / delta_theta - 0.5;
        return (*bb_data)[i][j] - 1.0;
    }
    double phi_analytical(double r, double theta, double phi, double t) {
        return 0.0;
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
    }  //================================================
    // Analytical solution for extrinsic curvature
    // NOTE: these should never be called; for rotating equilibria call interface
    // to InitializeCurvature that does *not* call these routines.  Only include the 
    // following here to make compiler happy...
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
        int i = abs(r) / delta_r - 0.5;
        if (theta > PI / 2.0) theta = PI - theta;
        int j = abs(theta) / delta_theta - 0.5;
        return (*lapse_data)[i][j];
    };
    double shift_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double shift_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double shift_p_analytical(double r, double theta, double phi, double t) {
        int i = abs(r) / delta_r - 0.5;
        if (theta > PI / 2.0) theta = PI - theta;
        int j = abs(theta) / delta_theta - 0.5;
        return -(*beta_phi_up_data)[i][j];
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
        int i = abs(r) / delta_r - 0.5;
        if (theta > PI / 2.0) theta = PI - theta;
        int j = abs(theta) / delta_theta - 0.5;
        return (*rho_data)[i][j];
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
        int i = abs(r) / delta_r - 0.5;
        if (theta > PI / 2.0) theta = PI - theta;
        int j = abs(theta) / delta_theta - 0.5;
        return (*v_phi_up_data)[i][j];
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
