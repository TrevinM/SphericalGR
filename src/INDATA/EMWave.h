// Tell emacs that this is -*-c++-*- mode
//
//================================================
// 
// Initial for electromagnetic wave
//
//================================================
//
class EMWave : public InData {
private:
#ifndef NoEllSolver
    FlatEllSolver3D* laplace;
#endif
    gf3d psi, rho, res, u, delta_psi;
    gf3d e_p, a_p;  // phi-components of E and A fields
    int n_r, n_theta, n_phi;
    int N_g;
    int n_psi_e;  // power of psi in E, i.e. E \propto psi^n
    int max_it;   // parameters for elliptic solver
    double tol;
    double A;   // parameter for initial data
    double r0;
    int flat;
    double PI;
    bool all_clear;
    ostringstream indata_name;
public:
    //================================================
    // Constructor
    //================================================
    EMWave(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
        InData(grid_i, cosmology) {
        N_g = grid->N_ghosts();
        indata_type = em_wave;
        all_clear = true;
        ifstream infile;
        infile.open(indata_input);
        if (!infile) {
            cerr << " INDATA Can't open " << indata_input
                << " for input. This is bad. " << endl;
            all_clear = false;
        } else
            cout << " INDATA: Reading initial data parameters from file "
            << indata_input << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> A;
        infile.get(buf, 100, '='); infile.get(c); infile >> r0;
        infile.get(buf, 100, '='); infile.get(c); infile >> flat;
        infile.get(buf, 100, '='); infile.get(c); infile >> max_it;
        infile.get(buf, 100, '='); infile.get(c); infile >> tol;
        infile.get(buf, 100, '='); infile.get(c); infile >> tau_star;
        infile.get(buf, 100, '='); infile.get(c); infile >> xi;
        if (flat == 0) {
            cout << " INDATA: Will set up E&M wave initial data with A = "
                << setprecision(16) << A << " and  r_0 = " << r0 << endl;
            cout << " INDATA: Will run elliptic solver with max_it = " << max_it
                << " and tol = " << tol << endl;
            cout << " INDATA: Will initialize lapse to precollapsed lapse: alpha = psi^{-2} " << endl;
            analytical = false;
        } else {
            analytical = true;
            bool changed_r_0 = false;
            if (r0 != 0.0) {
                changed_r_0 = true;
                r0 = 0.0;
            }
            cout << " INDATA: Will set up FLAT E&M wave initial data with A = "
                << setprecision(16) << A << " and  r_0 = " << r0 << endl;
            if (changed_r_0)
                cout << " INDATA: NOTE: changed r_0 to zero! " << endl;
        }
        cout << " INDATA: will analyze self-similarity with tau_star = "
            << tau_star << " and xi = " << xi << endl;
        cout << "===================================================" << endl;
        PI = acos(-1.0);
        n_psi_e = -6;
        if (flat == 0) {
            indata_name << "E&M wave initial data for A = "
                << setprecision(16) << A << " and r_0 = " << r0;
        } else {
            indata_name << "Flat E&M wave initial data for A = "
                << setprecision(16) << A << " and r_0 = " << r0;
        }
    };
    //================================================
    // Destructor
    //================================================
    ~EMWave() {};
    string Name() {
        return indata_name.str();
    };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) {
        if (!all_clear) return false;  // don't even try...
        cout << " INDATA: Initializing E&M wave initial data... " << endl;
        n_r = fct.dim1();
        n_theta = fct.dim2();
        n_phi = fct.dim3();
        //
        // set up grid functions
        //
        int gf_counter = 2000;
        psi.setup(grid, 1, "psi", gf_counter++, +1, +1, +1);
        delta_psi.setup(grid, 1, "delta_psi", gf_counter++, +1, +1, +1);
        rho.setup(grid, 1);
        res.setup(grid, 1);
        u.setup(grid, 1);
        e_p.setup(grid, 1, "e_p", gf_counter++, -1, -1, +1);
        a_p.setup(grid, 1, "a_p", gf_counter++, -1, -1, +1);
        //
        //================================================
        // compute source functions and initialize psi etc...
        //================================================
        //
        psi.equals(1.0);
        Compute_fields();
        //
        bool success = true;
        if (flat == 0)
            success = Solve_Hamiltonian();
        return success;
    }
    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double rl, double thetal, double phil, double t) {
        return 0.0;
    };
    double h_rt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_rp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_tt_analytical(double rl, double thetal, double phil, double t) {
        return 0.0;
    };
    double h_tp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double h_pp_analytical(double rl, double thetal, double phil, double t) {
        return 0.0;
    };
    double phi_analytical(double rl, double thetal, double phil, double tl) {
        // reality check:
        if (thetal == 0.0) cout << " ouch - found theta = 0 in phi_analytical " << endl;
        VecDoub* r = grid->r();
        VecDoub* theta = grid->theta();
        VecDoub* phi = grid->phi();
        if (rl != (*r)[grid->i_ind(rl)])
            cout << " i index not right - rl = " << rl << " i = " << grid->i_ind(rl) << " r(i) = "
            << (*r)[grid->i_ind(rl)] << endl;
        if (!((thetal == (*theta)[grid->j_ind(thetal)]) || ((*theta)[grid->j_ind(thetal)] == PI - thetal)))
            //      cout << " j index not right - thetal = " << thetal << " PI - thetal = " << PI - thetal <<  " j = " << grid->j_ind(thetal) << " theta(j) = " << (*theta)[grid->j_ind(thetal)] << endl;
            if (phil != (*phi)[grid->k_ind(phil)])
                cout << " k index not right - phil = " << phil << " k = " << grid->k_ind(phil) << " phi(k) = "
                << (*phi)[grid->k_ind(phil)] << endl;
        int i = grid->i_ind(rl);
        int j = grid->j_ind(thetal);
        int k = grid->k_ind(phil);
        double psil = psi(i, j, k);
        return log(psil);
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
    }  //================================================
    // Analytical solution for extrinsic curvature
    // NOTE: in initialize we compute \bar A_{ij} = psi^6 \tilde A_{ij},
    // now need BSSN rescaling of extrinsic curvature...
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
    double lapse_analytical(double rl, double thetal, double phil, double t) {
        int i = grid->i_ind(rl);
        int j = grid->j_ind(thetal);
        int k = grid->k_ind(phil);
        double psil = psi(i, j, k);
        // return 0.2 + 0.8*rl*rl/(1.0 + rl*rl);
        // return 1.0;
        return 1.0 / (psil * psil);
    };
    double shift_r_analytical(double r, double theta, double phi, double t) {
        // return 2.0 * r / (1.0 + r*r);
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
    double rho_0_analytical(double rl, double thetal, double phil, double tl) {
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
     // Analytical solution for Maxwell
     //================================================
     // returns correct function at t=0 only
    double e_p_analytical(double r, double theta, double phi, double t) {
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        return e_p(i, j, k);
    };
    //
    // provides flat-space analytical result at any time,
    // but only for r_0 = 0; see
    // eq. (16) in Knapp, Walker & Baumgarte, 2002
    //
    double a_p_analytical(double r, double theta, double phi, double t) {
        const double u = t + r;
        const double v = t - r;
        const double term1 = (exp(-v * v) - exp(-u * u)) / (r * r);
        const double term2 = (v * exp(-v * v) + u * exp(-u * u)) / r;
        return A * sin(theta) * (term1 - 2 * term2);
        // int i = grid->i_ind(r);
        // int j = grid->j_ind(theta);
        // int k = grid->k_ind(phi);
        // return a_p(i,j,k);
    };
    double a_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_t_analytical(double r, double theta, double phi, double t) {
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
    //====================================================
    // Solve Hamiltonian constraint
    //====================================================
    bool Solve_Hamiltonian() {
#ifndef NoEllSolver
        laplace = new FlatEllSolver3D(grid);
        //
        // Allocate and set up elliptic solver
        //
        int max_step = 50;
        Compute_rho();
        int step = 0;
        double residual = Hamiltonian_Residual();
        cout << " EM Wave: initial residual: " << residual << endl;
        while (residual > tol && step < max_step) {
            step++;
            Compute_u();
            laplace->SetupSolver(1.0, u);
            laplace->SetRHS(-1.0, res);
            //
            // Solve and get solution
            //
            int num_it;
            laplace->Solve(max_it, num_it, tol);
            // cout << " Solved Hamiltonian constraint to residual " << res_tri 
            //	   << " in " << num_it << " iterations. " << endl;
            laplace->GetSolution(delta_psi);
            update_psi();
            Compute_fields();
            Compute_rho();
            residual = Hamiltonian_Residual();
            cout << " EM Wave: residual after " << step << " steps: " << residual << endl;
        }
        delete laplace;
        bool success = false;
        if (residual < tol) success = true;
        return success;
#else
        cout << " Can't construct EM initial data without an Elliptic Solver!! " << endl;
        return false;
#endif
    }
    //
    // Compute fields (valid at t=0 only)
    //
    void Compute_fields() {
        for (int i = 0; i < n_r; i++) {
            const double rl = rho.r(i);
            for (int j = 0; j < n_theta; j++) {
                const double stl = rho.sintheta(j);
                for (int k = 0; k < n_phi; k++) {
                    const double psil = psi(i, j, k);
                    const double r_prop = rl;
                    const double rpr0 = r_prop + r0;
                    const double rmr0 = r_prop - r0;
                    // const double r_prop = psil * psil * rl;
                    // Note: upper rescaled component
                    const double psi_n = pow(psil, n_psi_e);
                    e_p[i][j][k] = -4.0 * A * r_prop * psi_n * stl *
                        (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
                    a_p[i][j][k] = 0.0;
                }
            }
        }
    };
    //
    // computes rho_ADM from EM wave; appears on right-hand side of Hamiltonian
    //
    void Compute_rho() {
        const double oo4p = 1.0 / (4.0 * PI);
        for (int i = N_g; i < n_r - N_g; i++) {
            const double rl = rho.r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < n_theta - N_g; j++) {
                const double stl = rho.sintheta(j);
                const double st2 = stl * stl;
                const double ctl = rho.costheta(j);
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double psil = psi(i, j, k);
                    const double psi4 = psil * psil * psil * psil;
                    const double E_phi = e_p(i, j, k);  // upper component
                    const double E2 = psi4 * E_phi * E_phi;
                    const double D_r_A_p = rl * stl * a_p.dr(i, j, k);
                    const double D_t_A_p = rl * stl * a_p.dtheta(i, j, k);
                    const double D_p_A_r = -stl * a_p(i, j, k);
                    const double D_p_A_t = -rl * ctl * a_p(i, j, k);
                    const double DA_rp = D_r_A_p - D_p_A_r;
                    const double DA_tp = D_t_A_p - D_p_A_t;
                    // const double DA_pr = - DA_rp;
                    // const double DA_pt = - DA_tp;
                    const double F2 = 2.0 * (DA_rp * DA_rp / (r2 * st2) +
                        DA_tp * DA_tp / (r2 * r2 * st2))
                        / (psi4 * psi4) - 2.0 * E2;
                    rho[i][j][k] = oo4p * (E2 + F2 / 4.0);
                }
            }
        }
    };
    //
    void Compute_u() {
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    // const double rl = rho.r(i);
                    const double psil = psi(i, j, k);
                    const double psi4 = psil * psil * psil * psil;
                    double factor = 0.0;
                    factor = 5 + 4 + 2 * n_psi_e;
                    u[i][j][k] = 2.0 * PI * psi4 * rho(i, j, k) * factor;
                }
    };
    void update_psi() {
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    psi[i][j][k] += delta_psi(i, j, k);
                }
    };
    double Hamiltonian_Residual() {
        for (int i = N_g; i < n_r - N_g; i++)
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double psil = psi(i, j, k);
                    const double psi5 = psil * psil * psil * psil * psil;
                    res[i][j][k] = psi.Laplace(i, j, k)
                        + 2.0 * PI * psi5 * rho(i, j, k);
                }
        return res.L2_norm();
    };
};


