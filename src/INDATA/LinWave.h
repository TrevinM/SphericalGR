// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Linear Wave Initial Data (Teukolsky wave)
// This version based on Rinne, 2008 arXiv:0809.1761
// (compare pages 320 ff in Numerical Relativity)
//================================================
//
class LinWave : public InData {
private:
    double Amp;
    double lambda;
    double r0;
    int l;
    int m;
    double r_exp;
    double alpha, delta_T, R_0;  // parameters for shift
    double expo;
    double A2_term0, A3_term1, A2_term2, A3_term3, A2_term4;
    double B2_term0, B3_term1, B2_term2, B3_term3, B2_term4;
    double C2_term0, C3_term1, C2_term2, C3_term3, C2_term4;
    double A2_term6, B2_term6, C2_term6;
    double A3_term5, B3_term5, C3_term5;
    double A4_term2, A4_term4, A4_term6, A4_term8;
    double B4_term2, B4_term4, B4_term6, B4_term8;
    double C4_term2, C4_term4, C4_term6, C4_term8;
    double delta_r; // for numerical derivatives...
    ostringstream indata_name;
public:
    //================================================
    // Constructor
    //================================================
    LinWave(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
        InData(grid_i, cosmology) {
        indata_type = linwave;
        analytical = true;
        ifstream infile;
        infile.open(indata_input);
        if (!infile)
            cerr << "Can't open " << indata_input
            << " for input. This is bad. " << endl;
        else
            cout << " Reading initial data parameters from file "
            << indata_input << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> Amp;
        infile.get(buf, 100, '='); infile.get(c); infile >> lambda;
        infile.get(buf, 100, '='); infile.get(c); infile >> r0;
        infile.get(buf, 100, '='); infile.get(c); infile >> l;
        infile.get(buf, 100, '='); infile.get(c); infile >> m;
        infile.get(buf, 100, '='); infile.get(c); infile >> r_exp;
        infile.get(buf, 100, '='); infile.get(c); infile >> alpha;
        infile.get(buf, 100, '='); infile.get(c); infile >> delta_T;
        infile.get(buf, 100, '='); infile.get(c); infile >> R_0;
        cout << " Will set up l = " << l << ", m = " << m << " Teukolsky wave with parameters " << endl;
        cout << "        Amplitude = " << Amp << endl;
        cout << "           lambda = " << lambda << endl;
        cout << "               r0 = " << r0 << endl;
        cout << "  Using Expansion inside r_exp = " << r_exp << endl;
        cout << "  Building shift from alpha = " << alpha << " delta_T = "
            << delta_T << " R_0 = " << R_0 << endl;
        cout << "===================================================" << endl;
        if (m != 0) {
            cout << " INDATA/LinWave.h currently implement for m = 0 only!!! "
                << endl;
            exit(0);
        }
        indata_name << "Teukolsky wave initial data for l = " << l << " with A = " << Amp << " and r0 = " << r0;
        //
        // Sanity check: make sure no equatorial symmetry is used for l = 3...
        //
#ifdef EQSYMMETRY
        if (l == 3) {
            cerr << " LINWAVE: equatorial symmetry is inconsistent with l = 3. "
                << endl;
            cerr << " LINWAVE: recompile without the -DEQSYMMETRY option. " << endl;
            exit(0);
        }
#endif
        //
        // compute expansion coefficients
        // 
        const double lam2 = lambda * lambda;
        const double lam4 = lam2 * lam2;
        const double lam6 = lam4 * lam2;
        const double lam8 = lam6 * lam2;
        const double lam10 = lam8 * lam2;
        const double lam12 = lam10 * lam2;
        const double lam14 = lam10 * lam4;
        const double lam16 = lam12 * lam4;
        const double lam18 = lam14 * lam4;
        const double lam22 = lam18 * lam4;
        const double lam26 = lam22 * lam4;
        const double lam30 = lam26 * lam4;
        const double r02 = r0 * r0;
        const double r04 = r02 * r02;
        const double r06 = r04 * r02;
        const double r08 = r06 * r02;
        const double r010 = r08 * r02;
        const double r012 = r010 * r02;
        const double r014 = r012 * r02;
        expo = Amp * exp(-(r02) / (lam2));
        //
        // Expansion coefficients for l = 2
        //
        A2_term0 = -64 * expo * (20 * r04 * lam2 - 60 * r02 * lam4 +
            15 * lam6) / (5 * lam10);
        A2_term2 = 64 * expo * (-56 * r06 * lam2 + 420 * r04 * lam4 -
            630 * r02 * lam6 + 105 * lam8) / (35 * lam14);
        A2_term4 = -16 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (315 * lam18);
        A2_term6 = 16 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4 -
            110880 * r06 * lam6 + 277200 * r04 * lam8 -
            207900 * r02 * lam10 + 20790 * lam12) / (10395 * lam22);
        //
        B2_term0 = -32 * expo * (20 * r04 * lam2 - 60 * r02 * lam4 +
            15 * lam6) / (5 * lam10);
        B2_term2 = 32 * expo * (-56 * r06 * lam2 + 420 * r04 * lam4 -
            630 * r02 * lam6 + 105 * lam8) / (21 * lam14);
        B2_term4 = -8.0 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (135 * lam18);
        B2_term6 = 8 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4 -
            110880 * r06 * lam6 + 277200 * r04 * lam8 -
            207900 * r02 * lam10 + 20790 * lam12) / (3465 * lam22);
        //
        C2_term0 = -32 * expo * (20 * r04 * lam2 - 60 * r02 * lam4 +
            15 * lam6) / (5 * lam10);
        C2_term2 = 352 * expo * (-56 * r06 * lam2 + 420 * r04 * lam4 -
            630 * r02 * lam6 + 105 * lam8) / (105 * lam14);
        C2_term4 = -184 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (945 * lam18);
        C2_term6 = 104 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4 -
            110880 * r06 * lam6 + 277200 * r04 * lam8 -
            207900 * r02 * lam10 + 20790 * lam12) / (10395 * lam22);
        //
        // Expansion coefficients for l = 3
        //
        A3_term1 = 64 * expo * (-112 * r06 * lam2 + 840 * r04 * lam4 - 1260 * r02 * lam6 +
            210 * lam8) / (7 * lam14);
        A3_term3 = -64 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (63 * lam18);
        A3_term5 = 32 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4
            - 110880 * r06 * lam6 + 277200 * r04 * lam8
            - 207900 * r02 * lam10 + 20790 * lam12) / (693 * lam22);
        //    
        B3_term1 = 64 * expo * (-112 * r06 * lam2 + 840 * r04 * lam4 - 1260 * r02 * lam6 +
            210 * lam8) / (21 * lam14);
        B3_term3 = -32 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (63 * lam18);
        B3_term5 = 32 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4
            - 110880 * r06 * lam6 + 277200 * r04 * lam8
            - 207900 * r02 * lam10 + 20790 * lam12) / (693 * lam22);
        //    
        C3_term1 = 32 * expo * (-112 * r06 * lam2 + 840 * r04 * lam4 - 1260 * r02 * lam6 +
            210 * lam8) / (21 * lam14);
        C3_term3 = -32 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (63 * lam18);
        C3_term5 = 464 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4
            - 110880 * r06 * lam6 + 277200 * r04 * lam8
            - 207900 * r02 * lam10 + 20790 * lam12) / (10395 * lam22);
        //
    // Expansion coefficients for l = 4
    //
        A4_term2 = -128 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (21 * lam18);
        A4_term4 = 128 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4 -
            110880 * r06 * lam6 + 277200 * r04 * lam8 -
            207900 * r02 * lam10 + 20790 * lam12) / (231 * lam22);
        A4_term6 = -64 * expo * (1664 * r012 * lam2 - 54912 * r010 * lam4 +
            617760 * r08 * lam6 - 2882880 * r06 * lam8 +
            5405400 * r04 * lam10 - 3243240 * r02 * lam12 +
            270270 * lam14) / (3003 * lam26);
        A4_term8 = 64 * expo * (-3840 * r014 * lam2 + 174720 * r012 * lam4 -
            2882880 * r010 * lam6 + 21621600 * r08 * lam8 -
            75675600 * r06 * lam10 + 113513400 * r04 * lam12
            - 56756700 * r02 * lam14
            + 4054050 * lam16) / (135135 * lam30);
        //
        B4_term2 = -32 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (21 * lam18);
        B4_term4 = 32 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4 -
            110880 * r06 * lam6 + 277200 * r04 * lam8 -
            207900 * r02 * lam10 + 20790 * lam12) / (165 * lam22);
        B4_term6 = -48 * expo * (1664 * r012 * lam2 - 54912 * r010 * lam4 +
            617760 * r08 * lam6 - 2882880 * r06 * lam8 +
            5405400 * r04 * lam10 - 3243240 * r02 * lam12 +
            270270 * lam14) / (5005 * lam26);
        B4_term8 = 16 * expo * (-3840 * r014 * lam2 + 174720 * r012 * lam4 -
            2882880 * r010 * lam6 + 21621600 * r08 * lam8 -
            75675600 * r06 * lam10 + 113513400 * r04 * lam12
            - 56756700 * r02 * lam14
            + 4054050 * lam16) / (61425 * lam30);
        //
        C4_term2 = -32 * expo * (288 * r08 * lam2 - 4032 * r06 * lam4 + 15120 * r04 * lam6 -
            15120 * r02 * lam8 + 1890 * lam10) / (63 * lam18);
        C4_term4 = 416 * expo * (-704 * r010 * lam2 + 15840 * r08 * lam4 -
            110880 * r06 * lam6 + 277200 * r04 * lam8 -
            207900 * r02 * lam10 + 20790 * lam12) / (3465 * lam22);
        C4_term6 = -1136. * expo * (1664. * r012 * lam2 - 54912. * r010 * lam4 +
            617760. * r08 * lam6 - 2882880. * r06 * lam8 +
            5405400. * r04 * lam10 - 3243240. * r02 * lam12 +
            270270. * lam14) / (135135. * lam26);
        C4_term8 = 592 * expo * (-3840 * r014 * lam2 + 174720 * r012 * lam4 -
            2882880 * r010 * lam6 + 21621600 * r08 * lam8 -
            75675600 * r06 * lam10 + 113513400 * r04 * lam12
            - 56756700 * r02 * lam14
            + 4054050 * lam16) / (2027025 * lam30);
        // recompute one term with better ordering just to try it out...
        // ... didn't make a difference in first 6 or so digits...
        C4_term6 = -1136. * expo * ((((((1664. * r02 * lam2 -
            54912. * lam4) * r02 +
            617760. * lam6) * r02 -
            2882880. * lam8) * r02 +
            5405400. * lam10) * r02 -
            3243240. * lam12) * r02 +
            270270. * lam14) / (135135. * lam26);
        // CHECK!!
        delta_r = 1.e-4;
    };
    //================================================
    // Destructor
    //================================================
    ~LinWave();
    string Name() { return indata_name.str(); };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) {
        return true;
    }
    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double r, double theta, double phi, double t) {
        return A(r, t) * Y(theta, phi);
    };
    double h_rt_analytical(double r, double theta, double phi, double t) {
        return B(r, t) * Y_t(theta, phi);
    };
    double h_rp_analytical(double r, double theta, double phi, double t) {
        return B(r, t) * Y_p(theta, phi);
    };
    double h_tt_analytical(double r, double theta, double phi, double t) {
        return C(r, t) * Y_tt(theta, phi) - 0.5 * A(r, t) * Y(theta, phi);
    };
    double h_tp_analytical(double r, double theta, double phi, double t) {
        return C(r, t) * Y_tp(theta, phi);
    };
    double h_pp_analytical(double r, double theta, double phi, double t) {
        return -C(r, t) * Y_tt(theta, phi) - 0.5 * A(r, t) * Y(theta, phi);
    };
    double phi_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for connection coefficients
    //================================================
    double lam_r_analytical(double r, double theta, double phi, double t,
        bool& done) {
        double lam_r, lam_t, lam_p;
        done = compute_lambdas(r, theta, phi, t, lam_r, lam_t, lam_p);
        return lam_r;
    }
    double lam_t_analytical(double r, double theta, double phi, double t,
        bool& done) {
        double lam_r, lam_t, lam_p;
        done = compute_lambdas(r, theta, phi, t, lam_r, lam_t, lam_p);
        return lam_t;
    }
    double lam_p_analytical(double r, double theta, double phi, double t,
        bool& done) {
        double lam_r, lam_t, lam_p;
        done = compute_lambdas(r, theta, phi, t, lam_r, lam_t, lam_p);
        return lam_p;
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
        double factor = R_0 * R_0 / (R_0 * R_0 + r * r);
        return -alpha * r * factor / delta_T;
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
    //================================================
    // Auxiliary functions for linear wave
    //================================================
private:
    //================================================
    // radial functions
    //================================================
    inline double A(double r, double t) {
        if (l == 2 && r0 == 0.0 && t == 0.0) {
            const double rl = r / lambda;
            return -192 * Amp * exp(-rl * rl);
        } else if (fabs(r) < r_exp && t == 0.0)
            return A_exp(r);
        else
            return A_not_exp(r, t);
    }
    inline double B(double r, double t) {
        if (l == 2 && r0 == 0.0 && t == 0.0) {
            const double rl = r / lambda;
            return 8 * Amp * exp(-rl * rl) * (8 * rl * rl - 12);
        } else if (fabs(r) < r_exp && t == 0.0)
            return B_exp(r);
        else
            return B_not_exp(r, t);
    }
    inline double C(double r, double t) {
        if (l == 2 && r0 == 0.0 && t == 0.0) {
            const double rl = r / lambda;
            return -4 * Amp * exp(-rl * rl) * (16 * rl * rl * rl * rl - 64 * rl * rl + 24);
        } else if (fabs(r) < r_exp && t == 0.0) {
            return C_exp(r);
        } else
            return C_not_exp(r, t);
    }
    //
    //
    //  
    inline double A_not_exp(double r, double t) {
        if (l == 2)
            return 24.0 * (r * (-r * F2(r, t) + 3.0 * F1(r, t)) - 3.0 * F(r, t)) / (r * r * r * r * r);
        else if (l == 3)
            return 120.0 * (r * (r * (-r * F3(r, t) + 6.0 * F2(r, t)) - 15.0 * F1(r, t)) + 15.0 * F(r, t)) / (r * r * r * r * r * r);
        else if (l == 4)
            return 360.0 * (r * (r * (r * (-r * F4(r, t) + 10.0 * F3(r, t)) - 45.0 * F2(r, t)) + 105.0 * F1(r, t)) - 105.0 * F(r, t)) / (r * r * r * r * r * r * r);
        else {
            cout << " LINWAVE: l = " << l << " not implemented! " << endl;
            exit(0);
        }
    };
    inline double B_not_exp(double r, double t) {
        if (l == 2)
            return 4.0 * (r * (r * (-r * F3(r, t) + 3.0 * F2(r, t)) - 6.0 * F1(r, t)) + 6.0 * F(r, t)) / (r * r * r * r * r);
        else if (l == 3)
            return 10.0 * (r * (r * (r * (-r * F4(r, t) + 6.0 * F3(r, t)) - 21.0 * F2(r, t)) + 45.0 * F1(r, t)) - 45.0 * F(r, t)) / (r * r * r * r * r * r);
        else if (l == 4)
            return 18.0 * (r * (r * (r * (r * (-r * F5(r, t) + 10.0 * F4(r, t)) - 55.0 * F3(r, t)) + 195.0 * F2(r, t)) - 420.0 * F1(r, t)) + 420.0 * F(r, t)) / (r * r * r * r * r * r * r);
        else {
            cout << " LINWAVE: l = " << l << " not implemented! " << endl;
            exit(0);
        }
    };
    inline double C_not_exp(double r, double t) {
        if (l == 2)
            return 2.0 * (r * (r * (r * (-r * F4(r, t) + 2.0 * F3(r, t)) - 3.0 * F2(r, t)) + 3.0 * F1(r, t)) - 3.0 * F(r, t)) / (r * r * r * r * r);
        else if (l == 3)
            return 2.0 * (r * (r * (r * (r * (-r * F5(r, t) + 5.0 * F4(r, t)) - 15.0 * F3(r, t)) + 30.0 * F2(r, t)) - 45.0 * F1(r, t)) + 45.0 * F(r, t)) / (r * r * r * r * r * r);
        else if (l == 4)
            return 2.0 * (r * (r * (r * (r * (r * (-r * F6(r, t) + 9.0 * F5(r, t)) - 45.0 * F4(r, t)) + 150.0 * F3(r, t)) - 360.0 * F2(r, t)) + 630.0 * F1(r, t)) - 630 * F(r, t)) / (r * r * r * r * r * r * r);
        else {
            cout << " LINWAVE: l = " << l << " not implemented! " << endl;
            exit(0);
        }
    };
    inline double F(double r, double t) {
        const double um = (r - t - r0);
        const double vm = (r + t - r0);
        const double uml = um / lambda;
        const double vml = vm / lambda;
        const double up = (r - t + r0);
        const double vp = (r + t + r0);
        const double upl = up / lambda;
        const double vpl = vp / lambda;
        const double uc = (r - t);
        const double vc = (r + t);
        // const double ucl = uc / lambda;
        // const double vcl = vc / lambda;
        const double pterm = (uc * exp(-upl * upl) +
            vc * exp(-vpl * vpl));
        const double mterm = (uc * exp(-uml * uml) +
            vc * exp(-vml * vml));
        return 0.5 * Amp * (pterm + mterm);
    }
    inline double F1(double r, double t) {
        const double uml = (r - t - r0) / lambda;
        const double vml = (r + t - r0) / lambda;
        const double upl = (r - t + r0) / lambda;
        const double vpl = (r + t + r0) / lambda;
        const double ucl = (r - t) / lambda;
        const double vcl = (r + t) / lambda;
        const double pterm = ((1.0 - 2.0 * ucl * upl) * exp(-upl * upl) +
            (1.0 - 2.0 * vcl * vpl) * exp(-vpl * vpl));
        const double mterm = ((1.0 - 2.0 * ucl * uml) * exp(-uml * uml) +
            (1.0 - 2.0 * vcl * vml) * exp(-vml * vml));
        return 0.5 * Amp * (pterm + mterm);
    }
    inline double F2(double r, double t) {
        const double uml = (r - t - r0) / lambda;
        const double vml = (r + t - r0) / lambda;
        const double upl = (r - t + r0) / lambda;
        const double vpl = (r + t + r0) / lambda;
        const double ucl = (r - t) / lambda;
        const double vcl = (r + t) / lambda;
        const double pterm = ((-4.0 * upl - 2.0 * ucl + 4.0 * ucl * upl * upl) * exp(-upl * upl) +
            (-4.0 * vpl - 2.0 * vcl + 4.0 * vcl * vpl * vpl) * exp(-vpl * vpl)) / lambda;
        const double mterm = ((-4.0 * uml - 2.0 * ucl + 4.0 * ucl * uml * uml) * exp(-uml * uml) +
            (-4.0 * vml - 2.0 * vcl + 4.0 * vcl * vml * vml) * exp(-vml * vml)) / lambda;
        return 0.5 * Amp * (pterm + mterm);
    }
    inline double F3(double r, double t) {
        const double uml = (r - t - r0) / lambda;
        const double vml = (r + t - r0) / lambda;
        const double upl = (r - t + r0) / lambda;
        const double vpl = (r + t + r0) / lambda;
        const double ucl = (r - t) / lambda;
        const double vcl = (r + t) / lambda;
        const double pterm = ((3.0 - 6.0 * upl * upl - 6.0 * ucl * upl + 4.0 * ucl * upl * upl * upl) * exp(-upl * upl) +
            (3.0 - 6.0 * vpl * vpl - 6.0 * vcl * vpl + 4.0 * vcl * vpl * vpl * vpl) * exp(-vpl * vpl)) / (lambda * lambda);
        const double mterm = ((3.0 - 6.0 * uml * uml - 6.0 * ucl * uml + 4.0 * ucl * uml * uml * uml) * exp(-uml * uml) +
            (3.0 - 6.0 * vml * vml - 6.0 * vcl * vml + 4.0 * vcl * vml * vml * vml) * exp(-vml * vml)) / (lambda * lambda);
        return -2.0 * 0.5 * Amp * (pterm + mterm);
    }
    inline double F4(double r, double t) {
        const double uml = (r - t - r0) / lambda;
        const double vml = (r + t - r0) / lambda;
        const double upl = (r - t + r0) / lambda;
        const double vpl = (r + t + r0) / lambda;
        const double ucl = (r - t) / lambda;
        const double vcl = (r + t) / lambda;
        const double pterm = ((12.0 * upl + 3.0 * ucl - 8.0 * upl * upl * upl - 12.0 * ucl * upl * upl + 4.0 * ucl * upl * upl * upl * upl) * exp(-upl * upl) +
            (12.0 * vpl + 3.0 * vcl - 8.0 * vpl * vpl * vpl - 12.0 * vcl * vpl * vpl + 4.0 * vcl * vpl * vpl * vpl * vpl) * exp(-vpl * vpl))
            / (lambda * lambda * lambda);
        const double mterm = ((12.0 * uml + 3.0 * ucl - 8.0 * uml * uml * uml - 12.0 * ucl * uml * uml + 4.0 * ucl * uml * uml * uml * uml) * exp(-uml * uml) +
            (12.0 * vml + 3.0 * vcl - 8.0 * vml * vml * vml - 12.0 * vcl * vml * vml + 4.0 * vcl * vml * vml * vml * vml) * exp(-vml * vml))
            / (lambda * lambda * lambda);
        return 4.0 * 0.5 * Amp * (pterm + mterm);
    }
    inline double F5(double r, double t) {
        const double uml = (r - t - r0) / lambda;
        const double vml = (r + t - r0) / lambda;
        const double upl = (r - t + r0) / lambda;
        const double vpl = (r + t + r0) / lambda;
        const double ucl = (r - t) / lambda;
        const double vcl = (r + t) / lambda;
        const double pterm = ((-32.0 * ucl * upl * upl * upl * upl * upl
            + 80.0 * upl * upl * upl * upl
            + 160.0 * ucl * upl * upl * upl
            - 120.0 * ucl * upl - 240.0 * upl * upl + 60.0) *
            exp(-upl * upl) +
            (-32.0 * vcl * vpl * vpl * vpl * vpl * vpl
                + 80.0 * vpl * vpl * vpl * vpl
                + 160.0 * vcl * vpl * vpl * vpl
                - 120.0 * vcl * vpl - 240.0 * vpl * vpl + 60.0) *
            exp(-vpl * vpl)) / (lambda * lambda * lambda * lambda);
        const double mterm = ((-32.0 * ucl * uml * uml * uml * uml * uml
            + 80.0 * uml * uml * uml * uml
            + 160.0 * ucl * uml * uml * uml
            - 120.0 * ucl * uml - 240.0 * uml * uml + 60.0) *
            exp(-uml * uml) +
            (-32.0 * vcl * vml * vml * vml * vml * vml
                + 80.0 * vml * vml * vml * vml
                + 160.0 * vcl * vml * vml * vml
                - 120.0 * vcl * vml - 240.0 * vml * vml + 60.0) *
            exp(-vml * vml)) / (lambda * lambda * lambda * lambda);
        return 0.5 * Amp * (pterm + mterm);
    }
    inline double F6(double r, double t) {
        const double uml = (r - t - r0) / lambda;
        const double vml = (r + t - r0) / lambda;
        const double upl = (r - t + r0) / lambda;
        const double vpl = (r + t + r0) / lambda;
        const double ucl = (r - t) / lambda;
        const double vcl = (r + t) / lambda;
        const double pterm = ((64.0 * ucl * upl * upl * upl * upl * upl * upl
            - 192.0 * upl * upl * upl * upl * upl
            - 480.0 * ucl * upl * upl * upl * upl
            + 720.0 * ucl * upl * upl
            + 960.0 * upl * upl * upl
            - 120.0 * ucl - 720.0 * upl) * exp(-upl * upl) +
            (64.0 * vcl * vpl * vpl * vpl * vpl * vpl * vpl
                - 192.0 * vpl * vpl * vpl * vpl * vpl
                - 480.0 * vcl * vpl * vpl * vpl * vpl
                + 720.0 * vcl * vpl * vpl
                + 960.0 * vpl * vpl * vpl
                - 120.0 * vcl - 720.0 * vpl) * exp(-vpl * vpl))
            / (lambda * lambda * lambda * lambda * lambda);
        const double mterm = ((64.0 * ucl * uml * uml * uml * uml * uml * uml
            - 192.0 * uml * uml * uml * uml * uml
            - 480.0 * ucl * uml * uml * uml * uml
            + 720.0 * ucl * uml * uml
            + 960.0 * uml * uml * uml
            - 120.0 * ucl - 720.0 * uml) * exp(-uml * uml) +
            (64.0 * vcl * vml * vml * vml * vml * vml * vml
                - 192.0 * vml * vml * vml * vml * vml
                - 480.0 * vcl * vml * vml * vml * vml
                + 720.0 * vcl * vml * vml
                + 960.0 * vml * vml * vml
                - 120.0 * vcl - 720.0 * vml) * exp(-vml * vml))
            / (lambda * lambda * lambda * lambda * lambda);
        return 0.5 * Amp * (pterm + mterm);
    }
    //================================================
    // angular functions
    //================================================
    inline double Y(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        const double st2 = stl * stl;
        const double st4 = st2 * st2;
        if (l == 2 && m == 0)
            return 2.0 - 3.0 * st2;
        else if (l == 2 && m == 2)
            return stl * stl * cos(2.0 * phi);
        else if (l == 3 && m == 0)
            return ctl * (2.0 - 5.0 * st2);
        else if (l == 4 && m == 0)
            return 35.0 * st4 - 40 * st2 + 8.0;
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double Y_t(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        const double st2 = stl * stl;
        if (l == 2 && m == 0)
            return -6.0 * stl * ctl;
        else if (l == 2 && m == 2)
            return 2.0 * stl * ctl * cos(2.0 * phi);
        else if (l == 3 && m == 0)
            return 3.0 * stl * (5.0 * st2 - 4.0);
        else if (l == 4 && m == 0)
            return 20.0 * ctl * stl * (7.0 * st2 - 4.0);
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double Y_p(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        if (l == 2 && m == 0)
            return 0.0;
        else if (l == 2 && m == 2)
            return -2.0 * stl * sin(2.0 * phi);
        else if (l == 3 && m == 0)
            return 0.0;
        else if (l == 4 && m == 0)
            return 0.0;
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double Y_tt(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        const double st2 = stl * stl;
        if (l == 2 && m == 0)
            return 3.0 * st2;
        else if (l == 2 && m == 2)
            return (2.0 - st2) * cos(2.0 * phi);
        else if (l == 3 && m == 0)
            return 15.0 * ctl * st2;
        else if (l == 4 && m == 0)
            return 30.0 * st2 * (6.0 - 7.0 * st2);
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double Y_tp(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        if (l == 2 && m == 0)
            return 0.0;
        else if (l == 2 && m == 2)
            return -2.0 * ctl * sin(2.0 * phi);
        else if (l == 3 && m == 0)
            return 0.0;
        else if (l == 4 && m == 0)
            return 0.0;
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    //======================================================
    // coefficients for expansion of A, B, and C about r = 0
    //======================================================
    inline double A_exp(double r) {
        if (l == 2)
            return A2_term0 + A2_term2 * r * r + A2_term4 * r * r * r * r +
            A2_term6 * r * r * r * r * r * r;
        else if (l == 3)
            return A3_term1 * r + A3_term3 * r * r * r + A3_term5 * r * r * r * r * r;
        else if (l == 4) {
            const double r2 = r * r;
            return r2 * (A4_term2 + r2 *
                (A4_term4 + r2 * (A4_term6 + r2 * A4_term8)));
        } else {
            cout << " LINWAVE: l = " << l << " not implemented! " << endl;
            exit(0);
            return 0;
        }
    };
    //======================================================
    inline double B_exp(double r) {
        if (l == 2)
            return B2_term0 + B2_term2 * r * r + B2_term4 * r * r * r * r +
            B2_term6 * r * r * r * r * r * r;
        else if (l == 3)
            return B3_term1 * r + B3_term3 * r * r * r + B3_term5 * r * r * r * r * r;
        else if (l == 4) {
            const double r2 = r * r;
            return r2 * (B4_term2 + r2 *
                (B4_term4 + r2 * (B4_term6 + r2 * B4_term8)));
        } else {
            cout << " LINWAVE: l = " << l << " not implemented! " << endl;
            exit(0);
            return 0;
        }
    };
    //======================================================
    inline double C_exp(double r) {
        if (l == 2)
            return C2_term0 + C2_term2 * r * r + C2_term4 * r * r * r * r +
            C2_term6 * r * r * r * r * r * r;
        else if (l == 3)
            return C3_term1 * r + C3_term3 * r * r * r + C3_term5 * r * r * r * r * r;
        else if (l == 4) {
            const double r2 = r * r;
            return r2 * (C4_term2 + r2 *
                (C4_term4 + r2 * (C4_term6 + r2 * C4_term8)));
        } else {
            cout << " LINWAVE: l = " << l << " not implemented! " << endl;
            exit(0);
            return 0;
        }
    };
    //======================================================
    // Derivatives of radial functions...
    //======================================================
    inline double dAdr(double r, double t) {
        if (l == 2 && r0 == 0.0 && t == 0.0) {
            const double rl = r / lambda;
            return 384 * Amp * rl / lambda * exp(-rl * rl);
        } else if (fabs(r) < r_exp && t == 0.0) {
            if (l == 2)
                return 2.0 * A2_term2 * r + 4.0 * A2_term4 * r * r * r + 6.0 * A2_term6 * r * r * r * r * r;
            else if (l == 3)
                return A3_term1 + 3.0 * A3_term3 * r * r + 5.0 * A3_term5 * r * r * r * r;
            else if (l == 4) {
                const double r2 = r * r;
                return r * (2.0 * A4_term2 +
                    r2 * (4.0 * A4_term4 +
                        r2 * (6.0 * A4_term6 +
                            r2 * 8.0 * A4_term8)));
            } else {
                cout << " LINWAVE: l = " << l << " not implemented! " << endl;
                exit(0);
                return 0;
            }
        } else
            return (1.0 * (A_not_exp(r + 3.0 * delta_r, t) - A_not_exp(r - 3.0 * delta_r, t)) -
                9.0 * (A_not_exp(r + 2.0 * delta_r, t) - A_not_exp(r - 2.0 * delta_r, t)) +
                45.0 * (A_not_exp(r + delta_r, t) - A_not_exp(r - delta_r, t))) /
            (60.0 * delta_r);
    };
    inline double dBdr(double r, double t) {
        if (l == 2 && r0 == 0.0 && t == 0.0) {
            const double rl = r / lambda;
            return -4 * Amp * rl / lambda * exp(-rl * rl) * (32 * rl * rl - 80);
        } else if (fabs(r) < r_exp && t == 0.0) {
            if (l == 2)
                return 2.0 * B2_term2 * r + 4.0 * B2_term4 * r * r * r + 6.0 * B2_term6 * r * r * r * r * r;
            else if (l == 3)
                return B3_term1 + 3.0 * B3_term3 * r * r + 5.0 * B3_term5 * r * r * r * r;
            else if (l == 4) {
                const double r2 = r * r;
                return r * (2.0 * B4_term2 +
                    r2 * (4.0 * B4_term4 +
                        r2 * (6.0 * B4_term6 +
                            r2 * 8.0 * B4_term8)));
            } else {
                cout << " LINWAVE: l = " << l << " not implemented! " << endl;
                exit(0);
                return 0;
            }
        } else
            return (1.0 * (B_not_exp(r + 3.0 * delta_r, t) - B_not_exp(r - 3.0 * delta_r, t)) -
                9.0 * (B_not_exp(r + 2.0 * delta_r, t) - B_not_exp(r - 2.0 * delta_r, t)) +
                45.0 * (B_not_exp(r + delta_r, t) - B_not_exp(r - delta_r, t))) /
            (60.0 * delta_r);

    };
    inline double dCdr(double r, double t) {
        if (l == 2 && r0 == 0.0 && t == 0.0) {
            const double rl = r / lambda;
            return 4 * Amp * rl / lambda * exp(-rl * rl) * (32 * rl * rl * rl * rl - 192 * rl * rl + 176);
        } else if (fabs(r) < r_exp && t == 0.0) {
            if (l == 2)
                return 2.0 * C2_term2 * r + 4.0 * C2_term4 * r * r * r + 6.0 * C2_term6 * r * r * r * r * r;
            else if (l == 3)
                return C3_term1 + 3.0 * C3_term3 * r * r + 5.0 * C3_term5 * r * r * r * r;
            else if (l == 4) {
                const double r2 = r * r;
                return r * (2.0 * C4_term2 +
                    r2 * (4.0 * C4_term4 +
                        r2 * (6.0 * C4_term6 +
                            r2 * 8.0 * C4_term8)));
            } else {
                cout << " LINWAVE: l = " << l << " not implemented! " << endl;
                exit(0);
                return 0;
            }
        } else
            return (1.0 * (C_not_exp(r + 3.0 * delta_r, t) - C_not_exp(r - 3.0 * delta_r, t)) -
                9.0 * (C_not_exp(r + 2.0 * delta_r, t) - C_not_exp(r - 2.0 * delta_r, t)) +
                45.0 * (C_not_exp(r + delta_r, t) - C_not_exp(r - delta_r, t))) /
            (60.0 * delta_r);
    };
    //================================================
    // derivatives of angular functions
    //================================================
    inline double dYdt(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        const double st2 = stl * stl;
        if (l == 2 && m == 0)
            return -6.0 * stl * ctl;
        else if (l == 2 && m == 2)
            return 2.0 * stl * ctl * cos(2.0 * phi);
        else if (l == 3 && m == 0)
            return -stl * (2.0 - 5.0 * stl * stl) - 10.0 * ctl * ctl * stl;
        else if (l == 4 && m == 0)
            return 20.0 * ctl * stl * (7.0 * st2 - 4.0);
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double dY_tdt(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        const double st2 = stl * stl;
        const double st4 = st2 * st2;
        if (l == 2 && m == 0)
            return -6.0 * ctl * ctl + 6.0 * stl * stl;
        else if (l == 2 && m == 2)
            return 2.0 * (ctl * ctl - stl * stl) * cos(2.0 * phi);
        else if (l == 3 && m == 0)
            return 3.0 * ctl * (5.0 * stl * stl - 4.0) + 30.0 * ctl * stl * stl;
        else if (l == 4 && m == 0)
            return -560.0 * st4 + 580 * st2 - 80.0;
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double dY_pdt(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        if (l == 2 && m == 0)
            return 0.0;
        else if (l == 2 && m == 2)
            return -ctl * sin(2.0 * phi);
        else if (l == 3 && m == 0)
            return 0.0;
        else if (l == 4 && m == 0)
            return 0.0;
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double dY_ttdt(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        const double st2 = stl * stl;
        if (l == 2 && m == 0)
            return 6.0 * stl * ctl;
        else if (l == 2 && m == 2)
            return -2.0 * stl * ctl * cos(2.0 * phi);
        else if (l == 3 && m == 0)
            return 30.0 * ctl * ctl * stl - 15.0 * stl * stl * stl;
        else if (l == 4 && m == 0)
            return 60.0 * stl * ctl * (6.0 - 14.0 * st2);
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    inline double dY_tpdt(double theta, double phi) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        if (l == 2 && m == 0)
            return 0.0;
        else if (l == 2 && m == 2)
            return 2.0 * stl * sin(2.0 * phi);
        else if (l == 3 && m == 0)
            return 0.0;
        else if (l == 4 && m == 0)
            return 0.0;
        else {
            cerr << " l = " << l << ", m = " << m << " has not been implemented! " << endl;
            return 0.0;
        }
    };
    //================================================
    // compute lambdas...
    //================================================
    bool compute_lambdas(double r, double theta, double phi, double t,
        double& lam_r, double& lam_t, double& lam_p) {
        double stl = sin(theta);
        double ctl = cos(theta);
        if (theta == grid->theta(grid->j_ind(theta))) {
            stl = grid->sintheta(grid->j_ind(theta));
            ctl = grid->costheta(grid->j_ind(theta));
        }
        const double r2 = r * r;
        const double st2 = stl * stl;
        //
        const double AA = A(r, t);
        const double BB = B(r, t);
        const double CC = C(r, t);
        //
        const double h_rr = AA * Y(theta, phi);
        const double h_rt = BB * Y_t(theta, phi);
        const double h_rp = BB * Y_p(theta, phi);
        const double h_tt = CC * Y_tt(theta, phi) - 0.5 * AA * Y(theta, phi);
        const double h_tp = CC * Y_tp(theta, phi);
        const double h_pp = -CC * Y_tt(theta, phi) - 0.5 * AA * Y(theta, phi);
        tensor g(1.0 + h_rr,
            r * h_rt,
            r * stl * h_rp,
            r2 * (1.0 + h_tt),
            r2 * stl * h_tp,
            r2 * st2 * (1.0 + h_pp));
        //    tensor g(1.0, 0.0, 0.0, r2, 0.0, r2*st2);
        tensor gup = g.inverse();
        // see eqs. (25) in Baumgarte et.al., PRD 87, 044026 (2013) 
        const double D_r_g_rr = dAdr(r, t) * Y(theta, phi);
        const double D_r_g_rt = r * dBdr(r, t) * Y_t(theta, phi);
        const double D_r_g_rp = r * stl * dBdr(r, t) * Y_p(theta, phi);
        const double D_r_g_tt = r2 * (dCdr(r, t) * Y_tt(theta, phi) - 0.5 *
            dAdr(r, t) * Y(theta, phi));
        const double D_r_g_tp = r2 * stl * dCdr(r, t) * Y_tp(theta, phi);
        const double D_r_g_pp = r2 * st2 * (-dCdr(r, t) * Y_tt(theta, phi) - 0.5 *
            dAdr(r, t) * Y(theta, phi));
        //
        const double D_t_g_rr = AA * dYdt(theta, phi) - 2.0 * h_rt;
        const double D_t_g_rt = r * (BB * dY_tdt(theta, phi) + h_rr - h_tt);
        const double D_t_g_rp = r * stl * (BB * dY_pdt(theta, phi) - h_tp);
        const double D_t_g_tt = r2 * (CC * dY_ttdt(theta, phi) - 0.5 *
            AA * dYdt(theta, phi) + 2.0 * h_rt);
        const double D_t_g_tp = r2 * stl * (CC * dY_tpdt(theta, phi) + h_rp);
        const double D_t_g_pp = r2 * st2 * (-CC * dY_ttdt(theta, phi) - 0.5 *
            AA * dYdt(theta, phi));
        // CHECK: assuming m = 0: no phi dependence, and h_rp = h_tp = 0.0;
        const double D_p_g_rr = 0.0;
        const double D_p_g_rt = 0.0;
        const double D_p_g_rp = r * stl * (stl * h_rr + ctl * h_rt - stl * h_pp);
        const double D_p_g_tt = 0.0;
        const double D_p_g_tp = r2 * stl * (stl * h_rt + ctl * h_tt - ctl * h_pp);
        const double D_p_g_pp = 0.0;
        //
        rank3tens D_g(D_r_g_rr, D_r_g_rt, D_r_g_rp, D_r_g_tt, D_r_g_tp, D_r_g_pp,
            D_t_g_rr, D_t_g_rt, D_t_g_rp, D_t_g_tt, D_t_g_tp, D_t_g_pp,
            D_p_g_rr, D_p_g_rt, D_p_g_rp, D_p_g_tt, D_p_g_tp, D_p_g_pp);
        //
        lam_r = 0.0;
        lam_t = 0.0;
        lam_p = 0.0;
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++) {
                    lam_r += gup[j][k] * gup[0][l] *
                        (D_g[j][k][l] + D_g[k][j][l] - D_g[l][j][k]);
                    lam_t += gup[j][k] * gup[1][l] *
                        (D_g[j][k][l] + D_g[k][j][l] - D_g[l][j][k]);
                    lam_p += gup[j][k] * gup[2][l] *
                        (D_g[j][k][l] + D_g[k][j][l] - D_g[l][j][k]);
                }
        lam_r *= 0.5;
        lam_t *= 0.5 * r;
        lam_p *= 0.5 * r * stl;
        return true;
    };
};

