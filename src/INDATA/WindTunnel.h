// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Black hole with momentum P in uniform fluid
//================================================
//

class WindTunnel : public InData {
private:
    EOS* eos;
    double rho_inf;
    double x_C, y_C, z_C;
    double m_BH;
    double P_BH;
    //  double n_scale;
    bool all_clear;
# ifndef NoEllSolver
    FlatEllSolver3D* ellsolver;
# endif
    gf3d psi_NS, psi, u, delta_u, f, res_gr, K;
    int n_r, n_theta, n_phi;
    int N_g;
    double K_back;
    double tol, tol_res, tol_tri;
    int max_it, max_step, max_iter;
    double eps;
    double PI;
public:
    //================================================
    // Constructor
    //================================================
    WindTunnel(char* indata_input, EOS* eos_i, Grid* grid_i, Cosmology* cosmology)
        : InData(grid_i, cosmology), eos(eos_i) {
        N_g = grid->N_ghosts();
        indata_type = wind_tunnel;
        analytical = false;
        all_clear = true;
        ifstream infile;
        infile.open(indata_input);
        if (!infile) {
            cerr << "Can't open " << indata_input
                << " for input. This is bad. " << endl;
            all_clear = false;
        } else
            cout << " INDATA: Reading initial data parameters from file "
            << indata_input << endl;
        char buf[100], c;
        double rho_c;
        infile.get(buf, 100, '='); infile.get(c); infile >> rho_inf;
        infile.get(buf, 100, '='); infile.get(c); infile >> m_BH;
        infile.get(buf, 100, '='); infile.get(c); infile >> P_BH;
        // infile.get(buf,100,'='); infile.get(c); infile >> n_scale; // scaling exponent of conformal factor
        infile.get(buf, 100, '='); infile.get(c); infile >> tol;
        infile.get(buf, 100, '='); infile.get(c); infile >> max_it;
        if (!eos) {
            cerr << " Can't set up wind tunnel data without an eos..." << endl;
            exit(1);
        }
        cout << " Will set up wind tunnel initial data with asymptotic density rho = " << rho_inf << endl;
        cout << "    using EOS " << eos->Name() << endl;
        cout << "    with BH of mass =  " << m_BH << " and momentum P_BH = " << P_BH << endl;
        cout << "    accurate within tolerance = " << tol << endl;
        cout << "    completing max iterations = " << max_it << endl;
        cout << "===================================================" << endl;
        int N_array = 10000;
        double dr_init = 1.e-6;
        PI = acos(-1.0);
    };
    //================================================
    // Destructory
    //================================================
    ~WindTunnel() {};
    string Name() { return "wind tunnel initial data"; };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) {
        if (!all_clear) return false;
        cout << " INDATA: Initializing wind tunnel with BH initial data... " << endl;
        n_r = fct.dim1();
        n_theta = fct.dim2();
        n_phi = fct.dim3();
        //
        // set up grid functions
        //
        int gf_counter = 2000;
        psi_NS.setup(grid, 1, "psi_NS", gf_counter++, +1, +1, +1);
        psi.setup(grid, 1, "psi", gf_counter++, +1, +1, +1);
        // rho_bar.setup(grid, 1, "rho_bar", gf_counter++, +1, +1, +1);
        u.setup(grid, 1, "u", gf_counter++, +1, +1, +1);
        delta_u.setup(grid, 1, "delta_u", gf_counter++, +1, +1, +1);
        f.setup(grid, 1, "function", gf_counter++, +1, +1, +1);
        res_gr.setup(grid, 1, "residual_grid", gf_counter++, +1, +1, +1);
        // Oo_alpha.setup(grid, 1, "one_over_alpha", gf_counter++, +1, +1, +1);
        u.equals(0.0); // initialize u
        K_back = 2.0 * sqrt(6.0 * PI * rho_inf);
        tol_res = tol; // residual tolerance
        tol_tri = tol; // trilinos tolerance 
        max_step = max_it; // iterations for updating u 
        max_iter = max_it; // iterations for trilinos
        //
        // initialize --> solve Hamiltonian constraint
        //
        double success = Solve_Hamiltonian();
        return success;
    };
    //
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_rt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_rp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_tt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_tp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_pp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double phi_analytical(double r, double theta, double phi, double t) {
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        return log(psi(i, j, k));
        /*
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
        return log(tov_sol->psi(r_C));
        */
    }
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
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        const double rl = grid->r(i);
        const double ctl = grid->costheta(j);
        const double psil = psi(i, j, k);
        const double psi6 = psil * psil * psil * psil * psil * psil;
        return 3.0 * P_BH * ctl / (rl * rl * psi6);
    };
    double a_rt_analytical(double r, double theta, double phi, double t) {
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        const double rl = grid->r(i);
        const double stl = grid->sintheta(j);
        const double psil = psi(i, j, k);
        const double psi6 = psil * psil * psil * psil * psil * psil;
        return -3.0 * P_BH * stl / (2.0 * rl * rl * psi6);
    };
    double a_rp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_tt_analytical(double r, double theta, double phi, double t) {
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        const double rl = grid->r(i);
        const double ctl = grid->costheta(j);
        const double psil = psi(i, j, k);
        const double psi6 = psil * psil * psil * psil * psil * psil;
        return -3.0 * P_BH * ctl / (2.0 * rl * rl * psi6);
    };
    double a_tp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_pp_analytical(double r, double theta, double phi, double t) {
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        const double rl = grid->r(i);
        // const double stl = grid->sintheta(j);
        // const double st2 = stl*stl;
        const double ctl = grid->costheta(j);
        const double psil = psi(i, j, k);
        const double psi6 = psil * psil * psil * psil * psil * psil;
        return -3.0 * P_BH * ctl / (2.0 * rl * rl * psi6);
    };
    double K_analytical(double r, double theta, double phi, double t) {
        return K_back;
    };
    //================================================
    // Analytical solution for gauge
    //================================================
    double lapse_analytical(double r, double theta, double phi, double t) {
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        return 1.0 / (psi(i, j, k) * psi(i, j, k));
        /*
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt( (x - x_C)*(x - x_C) + (y - y_C)*(y - y_C) + (z - z_C)*(z - z_C));
        return tov_sol->lapse(r_C);
        */
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
        return eos->rho_0_of_rho(rho_inf);
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
    //
    // double compute_rho_i(double rho_0) {
    //   //    double P = Kappa * pow(rho_0, Gamma); 
    //   //    double rho_i = P/(Gamma-1.0);
    //   double eps = eos->cold_eps(rho_0);
    //   return rho_0 * eps;
    // };
    //====================================================                                        
    // Solve Hamiltonian constraint                                                               
    //====================================================                                        
    bool Solve_Hamiltonian() {
#ifndef NoEllSolver
        ellsolver = new FlatEllSolver3D(grid);
        //                                                                         
        // Allocate and set up elliptic solver                                    
        //   
        int step = 0;
        compute_psi();
        double res_norm = Hamiltonian_Residual();
        cout << "Wind tunnel with BH: initial residual: " << res_norm << endl;
        while (res_norm > tol_res && step < max_step) {
            step++;
            compute_f(); // compute f
            ellsolver->SetupSolver(1.0, f);
            ellsolver->SetRHS(-1.0, res_gr);
            //                                                                      
            // Solve and get solution                                                
            //                                                                       
            int num_it;
            const double tol_tri = tol;
            ellsolver->Solve(max_iter, num_it, tol_tri);
            cout << " Solved Hamiltonian constraint to residual " << res_norm <<
                " in " << num_it << " iterations. " << endl;
            ellsolver->GetSolution(delta_u); //compute delta_u
            update_u(); //update u
            compute_psi(); //update psi
            res_norm = Hamiltonian_Residual(); //update residual
            cout << " INDATA: wind tunnel residual after " << step << " steps: " << res_norm << endl;
        }
        delete ellsolver;
        bool success = false;
        if (res_norm < tol) success = true;
        return success;
#else
        cout << " Can't construct wind tunnel with BH initial data without an Elliptic Solver!! " << endl;
        return false;
#endif
    };

    double Hamiltonian_Residual() {
        for (int i = N_g; i < n_r - N_g; i++)
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double rl = grid->r(i);
                    const double ctl = grid->costheta(j);
                    const double psi_l = psi(i, j, k);
                    const double psi5 = psi_l * psi_l * psi_l * psi_l * psi_l;
                    // const double psi5pn = pow(psi_l, 5+n_scale);   // 5pn: 5 plus n
                    const double psi7 = psi5 * psi_l * psi_l;
                    const double A2 = 9.0 * P_BH * P_BH / (2.0 * rl * rl * rl * rl) * (1.0 + 2.0 * ctl * ctl);
                    const double h_of_u = -A2 / psi7 / 8.0;  // Note: rho and K^2 terms cancel!
                    res_gr[i][j][k] = u.Laplace(i, j, k) - h_of_u;
                }
        return res_gr.L2_norm();
    };

    void compute_f() {
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    const double rl = grid->r(i);
                    const double ctl = grid->costheta(j);
                    const double psi_l = psi(i, j, k);
                    const double psi6 = psi_l * psi_l * psi_l * psi_l * psi_l * psi_l;
                    const double psi8 = psi6 * psi_l * psi_l;
                    const double A2 = 9.0 * P_BH * P_BH / (2.0 * rl * rl * rl * rl) * (1.0 + 2.0 * ctl * ctl);
                    // const double rho_bar_l = rho_bar(i, j, k);
                    f[i][j][k] = -7.0 * A2 / psi8 / 8.0;
                }
    };

    void compute_psi() {
        for (int i = 0; i < n_r; i++) {
            const double r_l = psi.r(i);
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    const double psi_NS_l = psi_NS[i][j][k];
                    const double u_l = u[i][j][k];
                    const double Oo_alpha_l = m_BH / (2.0 * r_l);
                    psi[i][j][k] = 1.0 + Oo_alpha_l + u_l;
                }
        }
    };

    void update_u() {
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    u[i][j][k] += delta_u(i, j, k);
                }
    };
};


