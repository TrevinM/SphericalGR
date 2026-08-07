// Tell emacs that this is -*-c++-*- mode
//
//================================================
// 
// Initial data for scalar field - see eq. (9) in Choptuik, Hirschmann, Liebling & Pretorius, 2003
//
// NOTE: will compute only scalar field analytically; all derivatives taken numerically
//
//================================================
//
class Choptuik : public InData {
private:
#ifndef NoEllSolver
    FlatEllSolver3D* laplace;
#endif
    gf3d psi, rho, res, u, delta_psi, tau;
    gf3d sf;  // scalar field - NOTE: psi is conformal factor...
    int n_r, n_theta, n_phi;
    int N_g;
    int max_it;   // parameters for elliptic solver
    double tol;
    double A, eps;   // parameter for initial data
    //  double r_0;          // ibid...
    double PI;
    ostringstream indata_name;
public:
    //================================================
    // Constructor
    //================================================
    Choptuik(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
        InData(grid_i, cosmology) {
        N_g = grid->N_ghosts();
        indata_type = choptuik;
        analytical = false;
        ifstream infile;
        infile.open(indata_input);
        if (!infile)
            cerr << " INDATA Can't open " << indata_input
            << " for input. This is bad. " << endl;
        else
            cout << " INDATA: Reading initial data parameters from file "
            << indata_input << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> A;
        infile.get(buf, 100, '='); infile.get(c); infile >> eps;
        infile.get(buf, 100, '='); infile.get(c); infile >> max_it;
        infile.get(buf, 100, '='); infile.get(c); infile >> tol;
        infile.get(buf, 100, '='); infile.get(c); infile >> tau_star;
        infile.get(buf, 100, '='); infile.get(c); infile >> xi;
        cout << " INDATA: Will set up scalar wave initial data " << endl;
        cout << " INDATA: ...  with A = " << setprecision(16) << A
            << " and eps = " << eps << endl;
        cout << " INDATA: Will run elliptic solver with max_it = " << max_it
            << " and tol = " << tol << endl;
        cout << " INDATA: Will initialize lapse to precollapsed lapse: alpha = psi^{-2} " << endl;
        cout << " INDATA: will analyze self-similarity with tau_star = "
            << tau_star << " and xi = " << xi << endl;
        cout << "===================================================" << endl;
        PI = acos(-1.0);
    };
    //================================================
    // Destructor
    //================================================
    ~Choptuik() {};
    string Name() {
        indata_name << "scalar wave initial data for A = " << setprecision(16) << A << " and eps = "
            << setprecision(16) << eps;
        return indata_name.str();
    };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) {
        cout << " SCALARFIELD: Initializing scalar field initial data... " << endl;
        n_r = fct.dim1();
        n_theta = fct.dim2();
        n_phi = fct.dim3();
        //
        // set up grid functions
        //
        // int gf_counter = 2000;
        psi.setup(grid, 1);
        delta_psi.setup(grid, 1);
        rho.setup(grid, 1);
        res.setup(grid, 1);
        u.setup(grid, 1);
        tau.setup(grid, 1);
        sf.setup(grid, 1);
        //
        //================================================
        // compute source functions and initialize psi etc...
        //================================================
        //
        Compute_SF();
        //
        psi.equals(1.0);
        Solve_Hamiltonian();
        return true;
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
        VecDoub* r = grid->r();
        VecDoub* theta = grid->theta();
        VecDoub* phi = grid->phi();
        if (rl != (*r)[grid->i_ind(rl)])
            cout << " i index not right - rl = " << rl << " i = " << grid->i_ind(rl) << " r(i) = "
            << (*r)[grid->i_ind(rl)] << endl;
        if (!((thetal == (*theta)[grid->j_ind(thetal)]) || ((*theta)[grid->j_ind(thetal)] == PI - thetal)))
            cout << " j index not right - thetal = " << thetal << " PI - thetal = " << PI - thetal << " j = " << grid->j_ind(thetal) << " theta(j) = "
            << (*theta)[grid->j_ind(thetal)] << endl;
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
    }
    //================================================
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
        //    return 1.0;
        return 1.0 / (psil * psil);
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
    // Initial data for scalar field
    //================================================
    double sf_analytical(double r, double theta, double phi, double t) {
        int i = grid->i_ind(r);
        int j = grid->j_ind(theta);
        int k = grid->k_ind(phi);
        return sf(i, j, k);
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
    //====================================================
    // Solve Hamiltonian constraint
    //====================================================
    double Solve_Hamiltonian() {
#ifndef NoEllSolver
        laplace = new FlatEllSolver3D(grid);
        //
        // Allocate and set up elliptic solver
        //
        int max_step = 20;
        Compute_tau();
        int step = 0;
        double residual = Hamiltonian_Residual();
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
            Compute_tau();
            residual = Hamiltonian_Residual();
            cout << " CHOPTUIK: residual after " << step << " steps: " << residual << endl;
        }
        delete laplace;
        return residual;
#else
        cout << " Can't construct scalar field initial data without an Elliptic Solver!! " << endl;
        return false;
#endif
    }
    //
    // Compute scalar field
    //
    void Compute_SF() {
        for (int i = 0; i < n_r; i++) {
            const double rl = rho.r(i);
            for (int j = 0; j < n_theta; j++) {
                const double stl = rho.sintheta(j);
                const double ctl = rho.costheta(j);
                for (int k = 0; k < n_phi; k++) {
                    // const double psil = psi(i,j,k);
                    // const double r_areal = rl*psil*psil;
                    // const double r_rel = r_areal / r_0;
                    sf[i][j][k] = A * exp(-rl * rl * (stl * stl + (1.0 - eps * eps) * ctl * ctl));
                    // if (j == 2 && k == 2) cout << " sf(" << i << ") = " << sf(i,j,k) << endl;
                }
            }
        }
    };
    //
    // computes rho_ADM = tau from scalar field; appears on right-hand side of Hamiltonian
    // NOTE: will assume that time-derivative is zero
    //
    void Compute_tau() {
        for (int i = N_g; i < n_r - N_g; i++)
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double rl = rho.r(i);
                    const double stl = rho.sintheta(j);
                    const double psil = psi(i, j, k);
                    const double psim4 = 1.0 / (psil * psil * psil * psil);
                    const double sf_r = sf.dr(i, j, k);
                    const double sf_t = sf.dtheta(i, j, k);
                    const double sf_p = sf.dphi(i, j, k);
                    tau[i][j][k] = 0.5 * psim4 * (sf_r * sf_r + (sf_t * sf_t + sf_p * sf_p / (stl * stl)) / (rl * rl));
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
                    u[i][j][k] = 2.0 * PI * psi4 * tau(i, j, k);
                }
    };
    void update_psi() {
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    psi[i][j][k] += delta_psi(i, j, k);
                    //	  if (j == 2 && k == 2) cout << " psi(" << i << ") = " << psi(i,j,k) << endl;
                }
    };
    double Hamiltonian_Residual() {
        for (int i = N_g; i < n_r - N_g; i++)
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double psil = psi(i, j, k);
                    const double psi5 = psil * psil * psil * psil * psil;
                    res[i][j][k] = psi.Laplace(i, j, k)
                        + 2.0 * PI * psi5 * tau(i, j, k);
                }
        return res.L2_norm();
    };
};


