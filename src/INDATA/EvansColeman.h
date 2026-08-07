// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Evans-Coleman initial data for radiation fluid, but 
// allows for non-spherical deformation (parameter eps)
//
// Reference: PRL 72, 1782 (1994)
//
// Note: with rho_0 = 0 and fluid at rest have rho_ADM = rho_i
//
//================================================
//

class EvansColeman : public InData {
private:
    double x_C, y_C, z_C;
    double M;
    double a, lambda, rho0;
    double eps;   // deformation parameter
    int ell;      // order of Legendre polynomial
#ifndef NoEllSolver
    FlatEllSolver3D* ellsolver;
#endif
    gf3d psi, rho, res, u, delta_psi;
    gf3d psilapse, delta_psilapse;
    int n_r, n_theta, n_phi;
    //  VecDoub *r, *Delta_r, *theta, *costheta, *x, *x_prime, *x_dprime, *phi; 
    //  double dtheta, dphi, c;
    int N_g;
    int max_it;   // parameters for elliptic solver
    int init_lapse;  // parameter for initial choice of lapse
    double tol;
    double tiny;  // small number, to make sure that conversion from coordinates to indices works
    double eta;   // parameter for initial data, see eq. (6)
    double r_0, r_c;   // ibid...
    double PI;
    ostringstream indata_name;
public:
    //================================================
    // Constructor
    //================================================
    EvansColeman(char* indata_input, Grid* grid_i, Cosmology* cosmology) :
        InData(grid_i, cosmology) {
        N_g = grid->N_ghosts();
        indata_type = evanscoleman;
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
        infile.get(buf, 100, '='); infile.get(c); infile >> eta;
        infile.get(buf, 100, '='); infile.get(c); infile >> r_0;
        infile.get(buf, 100, '='); infile.get(c); infile >> r_c;
        infile.get(buf, 100, '='); infile.get(c); infile >> eps;
        infile.get(buf, 100, '='); infile.get(c); infile >> ell;
        infile.get(buf, 100, '='); infile.get(c); infile >> max_it;
        infile.get(buf, 100, '='); infile.get(c); infile >> tol;
        infile.get(buf, 100, '='); infile.get(c); infile >> init_lapse;
        cout << " INDATA: Will set up Evans-Coleman Wave Initial data " << endl;
        cout << " INDATA: ...  with amplitude eta = " << setw(24) << setprecision(16) << eta
            << ", length-scale r_0 = " << r_0 << endl;
        cout << " INDATA: ...  centered on r_c = " << r_c << ", deformation parameter eps = " << eps
            << " and ell = " << ell << endl;
        cout << " INDATA: Will run elliptic solver with max_it = " << max_it
            << " and tol = " << tol << endl;
        if (init_lapse == 0)
            cout << " INDATA: Will initialize lapse to precollapsed lapse: alpha = psi^{-2} " << endl;
        else if (init_lapse == 1)
            cout << " INDATA: Will initialize lapse to maximal slicing: dot K = 0 " << endl;
        else {
            cout << " INDATA: Unknown choice for lapse : " << init_lapse << endl;
            init_lapse = 0;
        }
        cout << "===================================================" << endl;
        PI = acos(-1.0);
    };
    //================================================
    // Destructor
    //================================================
    ~EvansColeman() {};
    string Name() {
        if (init_lapse == 1)
            indata_name << "Evans Coleman initial data with eta = " << setprecision(16) << eta
            << ", initial lapse maximal";
        else
            indata_name << "Evans Coleman initial data with eta = " << setprecision(16) << eta
            << ", precollapsed lapse";
        return indata_name.str();
    };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct) {
        cout << " EVANSCOLEMAN: Initializing Evans-Coleman data... " << endl;
        tiny = 1.e-10;
        n_r = fct.dim1();
        n_theta = fct.dim2();
        n_phi = fct.dim3();
        // r = fct.r_pointer();
        // Delta_r = fct.Delta_r_pointer();
        // theta = fct.theta_pointer();
        // costheta = fct.costheta_pointer();
        // x = fct.x_pointer();
        // x_prime = fct.x_prime_pointer();
        // x_dprime = fct.x_dprime_pointer();
        // phi = fct.phi_pointer();
        // c = fct.log_factor();
        // //    dr = fct.r(1) - fct.r(0);
        // dtheta = fct.theta(1) - fct.theta(0);
        // dphi = fct.phi(1) - fct.phi(0);
        //
        // set up grid functions
        //
        psi.setup(grid, 1);
        delta_psi.setup(grid, 1);
        rho.setup(grid, 1);  // rho_ADM = rho_i
        res.setup(grid, 1);
        u.setup(grid, 1);
        psilapse.setup(grid, 1);
        delta_psilapse.setup(grid, 1);
        //
        //================================================
        // compute source functions and initialize psi
        //================================================
        //
        psi.equals(1.0);
#ifndef NoEllSolver
        int max_step = 100;
        Compute_rho();
        int step = 0;
        double residual = Residual();
        //    double residual = 0.0;
        while (residual > tol && step < max_step) {
            step++;
            //
            // Allocate and set up elliptic solver
            //
            ellsolver = new FlatEllSolver3D(grid);
            Compute_u();
            ellsolver->SetupSolver(1.0, u);
            ellsolver->SetRHS(-1.0, res);
            //
            // Solve and get solution
            //
            int num_it;
            ellsolver->Solve(max_it, num_it, tol);
            //	cout << " Solved Hamiltonian constraint to residual " << res_tri 
            //	     << " in " << num_it << " iterations. " << endl;
            ellsolver->GetSolution(delta_psi);
            update_psi();
            Compute_rho();
            residual = Residual();
            cout << " EVANSCOLEMAN: After " << step << " steps: residual = " << residual << endl;
            delete ellsolver;
        }
        //
        //================================================
        // now compute initial maximal lapse
        //
        // Note: eq. is linear, so could solve directly - except that 
        // we're using a second order solver.  Iterate to get higher order
        //
        // Also note: psilapse = alpha * psi - 1
        //================================================
        //
        if (init_lapse == 1) {
            psilapse.equals(0.0);
            residual = Lapse_Residual();
            cout << " EVANSCOLEMAN: initial maximal slicing residual: " << residual
                << endl;
            step = 0;
            for (int i = 0; i < n_r; i++)
                for (int j = 0; j < n_theta; j++)
                    for (int k = 0; k < n_phi; k++) {
                        const double psil = psi(i, j, k);
                        u[i][j][k] = 2.0 * PI * psil * psil * psil * psil * rho(i, j, k);
                    }
            while (residual > tol && step < max_step) {
                step++;
                ellsolver = new FlatEllSolver3D(grid);
                ellsolver->SetupSolver(-1.0, u);
                ellsolver->SetRHS(-1.0, res);
                int num_it;
                ellsolver->Solve(max_it, num_it, tol);
                ellsolver->GetSolution(delta_psilapse);
                for (int i = 0; i < n_r; i++)
                    for (int j = 0; j < n_theta; j++)
                        for (int k = 0; k < n_phi; k++)
                            psilapse[i][j][k] += delta_psilapse(i, j, k);
                residual = Lapse_Residual();
                cout << " EVANSCOLEMAN: after " << step << " iterations maximal slicing residual: " << residual
                    << endl;
                delete ellsolver;
            }
        }
        return true;
#else
        cout << " Can't construct Evans-Coleman initial data without an Elliptic Solver!! " << endl;
        return false;
#endif
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
            cout << " i index not right - rl = " << rl
            << " i = " << grid->i_ind(rl) << " r(i) = "
            << (*r)[grid->i_ind(rl)] << endl;
        if (!((thetal == (*theta)[grid->j_ind(thetal)]) ||
            ((*theta)[grid->j_ind(thetal)] == PI - thetal)))
            cout << " j index not right - thetal = " << thetal << " PI - thetal = "
            << PI - thetal << " j = " << grid->j_ind(thetal) << " theta(j) = "
            << (*theta)[grid->j_ind(thetal)] << endl;
        if (phil != (*phi)[grid->k_ind(phil)])
            cout << " k index not right - phil = " << phil << " k = "
            << grid->k_ind(phil) << " phi(k) = "
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
        if (init_lapse == 1)
            return (psilapse(i, j, k) + 1.0) / psil;
        else
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
    //
    // NOTE: we're assuming radiation fluid here; rho_0_analytical will return
    //       rho_i instead of rho_0
    //
    double rho_0_analytical(double rl, double thetal, double phil, double tl) {
        int i = grid->i_ind(rl);
        int j = grid->j_ind(thetal);
        int k = grid->k_ind(phil);
        return rho(i, j, k);
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
    // computes rho according to eq. (6)
    //
    void Compute_rho() {
        const double PIm32 = 1.0 / (PI * sqrt(PI));
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    const double theta_l = rho.theta(j);
                    const double costheta = cos(theta_l);
                    //	  const double P2 = (3.0*costheta*costheta - 1.0) / 2.0;
                    const double Pell = P_ell(costheta);
                    const double rl = rho.r(i);
                    const double r_areal_plus = rl * psi(i, j, k) * psi(i, j, k) + r_c;
                    const double r_areal_minus = rl * psi(i, j, k) * psi(i, j, k) - r_c;
                    const double r_0_areal = r_0;  // Interpret r_0 as areal radius!
                    const double factor_plus = exp(-r_areal_plus * r_areal_plus / (r_0_areal * r_0_areal));
                    const double factor_minus = exp(-r_areal_minus * r_areal_minus / (r_0_areal * r_0_areal));
                    const double factor = (factor_plus + factor_minus) / 2.0;
                    rho[i][j][k] = 0.5 * PIm32 * eta * factor / (r_0_areal * r_0_areal) *
                        (1.0 + eps * rl * rl / (1.0 + rl * rl) * Pell);
                }
    };
    void Compute_u() {
        const double PIm32 = 1.0 / (PI * sqrt(PI));
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    const double psil = psi(i, j, k);
                    const double psi4 = psil * psil * psil * psil;
                    const double theta_l = rho.theta(j);
                    const double costheta = cos(theta_l);
                    // const double P2 = (3.0*costheta*costheta - 1.0) / 2.0;
                    const double Pell = P_ell(costheta);
                    const double rl = rho.r(i);
                    const double r_areal_plus = rl * psil * psil + r_c;
                    const double r_areal_minus = rl * psil * psil - r_c;
                    const double r_0_areal = r_0;  // Interpret r_0 as areal radius!
                    const double factor_plus = exp(-r_areal_plus * r_areal_plus / (r_0_areal * r_0_areal));
                    const double factor_minus = exp(-r_areal_minus * r_areal_minus / (r_0_areal * r_0_areal));
                    // const double factor = (factor_plus + factor_minus)/2.0;
                    const double dfpdpsi = -4.0 * factor_plus * (psil * psil * rl + r_c) * psil * rl / (r_0 * r_0);
                    const double dfmdpsi = -4.0 * factor_minus * (psil * psil * rl - r_c) * psil * rl / (r_0 * r_0);
                    const double dfdpsi = (dfpdpsi + dfmdpsi) / 2.0;
                    const double drhodpsi = 0.5 * PIm32 * eta * dfdpsi / (r_0_areal * r_0_areal) *
                        (1.0 + eps * rl * rl / (1.0 + rl * rl) * Pell);
                    //	  u[i][j][k] = 10.0*PI*psi4*rho(i,j,k) - 8.0*PI*psi4*psi4*rho(i,j,k)*r_l*r_l/(r_0*r_0);
                    u[i][j][k] = 10.0 * PI * psi4 * rho(i, j, k) + 2.0 * psi4 * psil * PI * drhodpsi;
                    // u[i][j][k] = 0.0;
                }
    };
    void update_psi() {
        for (int i = 0; i < n_r; i++)
            for (int j = 0; j < n_theta; j++)
                for (int k = 0; k < n_phi; k++) {
                    psi[i][j][k] += delta_psi(i, j, k);
                }
    };
    double Residual() {
        for (int i = N_g; i < n_r - N_g; i++)
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double psil = psi(i, j, k);
                    const double psi5 = psil * psil * psil * psil * psil;
                    res[i][j][k] = psi.Laplace(i, j, k) + 2.0 * PI * psi5 * rho(i, j, k);
                }
        return res.L2_norm();
    };
    double Lapse_Residual() {
        for (int i = N_g; i < n_r - N_g; i++)
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double psil = psi(i, j, k);
                    const double psi4 = psil * psil * psil * psil;
                    res[i][j][k] = psilapse.Laplace(i, j, k) -
                        (psilapse(i, j, k) + 1.0) * 2.0 * PI * psi4 * rho(i, j, k);
                }
        return res.L2_norm();
    };
    //
    // Legendre polynomials
    // 
    inline double P_ell(double x) {
        if (ell == 2)
            return (3.0 * x * x - 1.0) / 2.0;
        else if (ell == 4)
            return (35.0 * x * x * x * x - 30.0 * x * x + 3.0) / 8.0;
        else
            cerr << " ell = " << ell << " not implemented in INDATA/EvansColeman.h " << endl;
        return 0.0;
    }
};


