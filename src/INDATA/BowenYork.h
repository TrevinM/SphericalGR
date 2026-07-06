// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Rotating Bowen-York initial data 
//================================================
//
class BowenYork : public InData {
private:
  double x_C, y_C, z_C;
  double M;
  double a, J; 
#ifndef NoEllSolver
  FlatEllSolver3D * ellsolver;
#endif
  gf3d psi0, u, h, A2, res, linterm;
  int n_r, n_theta, n_phi;
  int N_g;
  int max_it;   // parameters for elliptic solver
  double tol;
  double eps;   // small number, to make sure that conversion from coordinates to indices works
public:
  //================================================
  // Constructor
  //================================================
  BowenYork(char * indata_input, Grid * grid_i, Cosmology *cosmology) :
    InData(grid_i, cosmology) {
    N_g = grid->N_ghosts();
    indata_type = bowenyork;
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
    infile.get(buf,100,'='); infile.get(c); infile >> M;
    infile.get(buf,100,'='); infile.get(c); infile >> a;
    infile.get(buf,100,'='); infile.get(c); infile >> max_it;
    infile.get(buf,100,'='); infile.get(c); infile >> tol;
    cout << " Will set up Bowen-York black hole with mass M = " << M << endl;
    cout << "    and angular momentum parameter a = " << a << endl; 
    cout << "    Will run elliptic solver with max_it = " << max_it << " and tol = " << tol << endl;
    cout << "===================================================" << endl;
    J = a * M;
  };
  //================================================
  // Destructor
  //================================================
  ~BowenYork() {};
  string Name() { return "Bowen-York initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { 
    cout << " Initializing Bowen-York spinning black hole... " << endl;
    eps = 1.e-10;
    n_r = fct.dim1();
    n_theta = fct.dim2();
    n_phi = fct.dim3();
    //
    // set up grid functions
    //
    psi0.setup(grid, 1);
    h.setup(grid, 1);
    u.setup(grid, 1);
    A2.setup(grid, 1);
    res.setup(grid, 1);
    linterm.setup(grid, 1);
    //
    // compute source function A2 and initialize psi
    //
    Compute_A2(A2);
    Initialize_psi();
    //
    // compute initial residual
    //
    double residual = Residual();
#ifndef NoEllSolver
    //
    // if we have an elliptic solver, iterate...
    //
    cout << " Using puncture method for Bowen-York data - initial residual = " << residual << endl;

    int it = 0;
    while (residual > tol && it < max_it) {
      it++;
      //
      // Compute u
      //
      Compute_linterm();
      //
      // Allocate and set up elliptic solver
      //
      ellsolver = new FlatEllSolver3D(grid);
      ellsolver->SetupSolver(1.0,linterm);
      ellsolver->SetRHS(-1.0,res);
      //
      // Solve and get solution
      //
      int num_it;
      double res_tri = ellsolver->Solve(max_it,num_it,tol);
      cout << " residual after " << it << " loops: "  << res_tri << " in " << num_it << " iterations. " << endl;
      ellsolver->GetSolution(h);
      Update_u();
      residual = Residual();
      cout << " Our residual = " << residual << endl;
      delete ellsolver;
    }
    return true; 
#else
    cout << " Can't construct Brill waves without an Elliptic Solver!! " << endl;
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
  double phi_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    return log(psi0(i,j,k) + u(i,j,k));
  };
  //================================================
  // Analytical solution for connection coefficients
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true; return 0;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true; return 0;
  }
  double lam_p_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
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
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi0(i,j,k) + u(i,j,k);
    return 3.0*J*sin(theta)/(r*r*r*psil*psil*psil*psil*psil*psil);
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
    double psil = psi0(i,j,k) + u(i,j,k);
    return 1.0/(psil*psil);
    //    return 1.0;
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
  double psi_background(double r) { return 1.0 + M/(2.0*abs(r)); }
  void Update_u() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  u[i][j][k] += h(i,j,k);
	}
  };
  void Initialize_psi() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = h.r(i);
	  psi0[i][j][k] = psi_background(rl);
	  u[i][j][k] = 0.0;
	}
  };
  void Compute_A2(gf3d & A2) {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = u.r(i);
	  const double sintheta = sin(u.theta(j));
	  const double sin2theta = sintheta*sintheta;
	  A2[i][j][k] = 18.0*J*J*sin2theta/(rl*rl*rl*rl*rl*rl);
	}
  };
  void Compute_linterm() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double psil = psi0(i,j,k) + u(i,j,k);
	  linterm[i][j][k] = - 7.0*A2(i,j,k)/(8.0*psil*psil*psil*psil*psil*psil*psil*psil);
	}
  };
  double Residual() {
    for (int i = N_g; i < n_r-N_g; i++) 
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  const double psil = psi0(i,j,k) + u(i,j,k);
	  res[i][j][k] = u.Laplace(i,j,k) + A2(i,j,k)/(8.0*psil*psil*psil*psil*psil*psil*psil);
	}
    return res.L2_norm();
  };
};
