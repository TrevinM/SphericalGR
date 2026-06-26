// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Brill Wave initial data 
//
// Uses functions
//
//  q = A (rho/lambda)^2 exp(-((rho - rho_0)^2 - z^2)/lambda^2)
//    = A (r/lambda)^2 sin^2 theta exp( - ((r sin theta - rho_0)^2 - r^2 cos^2 theta)/lambda^2 )
//
// and 
//
// u = 1/4 * (\partial2 q / \partial \rho^2 + \partial^2 q / \partial z^2)
//
//================================================
//
class Brill : public InData {
private:
  double x_C, y_C, z_C;
  double M;
  double a, lambda, rho0; 
#ifndef NoEllSolver
  FlatEllSolver3D * ellsolver;
#endif
  gf3d psi, q, res, u, delta_psi;
  int n_r, n_theta, n_phi;
  int N_g;
  int max_it;   // parameters for elliptic solver
  double tol;
  double eps;   // small number, to make sure that conversion from coordinates to indices works
public:
  //================================================
  // Constructor
  //================================================
  Brill(char * indata_input, Grid * grid_i, Cosmology *cosmology) :
    InData(grid_i, cosmology) {			
    N_g = grid->N_ghosts();
    indata_type = brill;
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
    infile.get(buf,100,'='); infile.get(c); infile >> a;
    infile.get(buf,100,'='); infile.get(c); infile >> lambda;
    infile.get(buf,100,'='); infile.get(c); infile >> rho0;
    infile.get(buf,100,'='); infile.get(c); infile >> max_it;
    infile.get(buf,100,'='); infile.get(c); infile >> tol;
    cout << " Will set up Brill Wave Initial data " << endl;
    cout << "    with amplitude a = " << a << ", wavelength lambda = " << lambda 
	 << " and offset rho_0 = " << rho0 << endl;
    cout << "    Will run elliptic solver with max_it = " << max_it << " and tol = " << tol << endl;
    cout << "===================================================" << endl;
  };
  //================================================
  // Destructor
  //================================================
  ~Brill() {};
  string Name() { return "Brill wave initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { 
    cout << " Initializing Brill waves... " << endl;
    eps = 1.e-10;
    n_r = fct.dim1();
    n_theta = fct.dim2();
    n_phi = fct.dim3();
    //
    // set up grid functions
    //
    psi.setup(grid, 1);
    delta_psi.setup(grid, 1);
    u.setup(grid, 1);
    q.setup(grid, 1);
    res.setup(grid, 1);
    //
    // compute source functions and q
    //
    Compute_u(u);
    Compute_q(q);
    psi.equals(0.0);
#ifndef NoEllSolver
    //
    // need to iterate, since we're using 2nd order elliptic solver
    //
    double residual = Residual(psi,u,res);
    cout << " BRILL: Initial residual = " << residual << endl;
    int it = 0;
    int it_max = 100;
    while (residual > tol && it < it_max) {
      //
      // Allocate and set up elliptic solver
      //
      ellsolver = new FlatEllSolver3D(grid);
      ellsolver->SetupSolver(1.0,u);
      ellsolver->SetRHS(-1.0,res);
      //
      // Solve and get solution
      //
      int num_it;
      // double res_tri = 
      ellsolver->Solve(max_it,num_it,tol);
      //      cout << " Solved Brill equation to residual " << res_tri << " in " << num_it << " iterations. " << endl;
      ellsolver->GetSolution(delta_psi);
      Update_Psi(psi,delta_psi);
      //
      // Compute new residual
      //
      residual = Residual(psi,u,res);
      cout << " BRILL: Residual after " << it << " steps = " << residual << endl;
      it++;
      delete ellsolver;
    }
    //
    // if iteration converged: 
    //
    if (it < it_max)
      return true; 
    //
    // else...
    //
    else {
      cout << " BRILL: iteration did not converge! " << endl;
      return false;
    }
#else
    cout << " Can't construct Brill waves without an Elliptic Solver!! " << endl;
    return false;
#endif
  } 
  //================================================
  // Analytical solution for h_{ij}
  //================================================
  double h_rr_analytical(double rl, double thetal, double phil, double t) {
    // reality check:
    VecDoub *r = grid->r();
    VecDoub *theta = grid->theta();
    VecDoub *phi = grid->phi();
    if (rl != (*r)[grid->i_ind(rl)]) 
      cout << " i index not right - rl = " 			    
    	   << rl << " i = " << grid->i_ind(rl) << " r(i) = " 
    	   << (*r)[grid->i_ind(rl)] << endl;
    if (thetal != (*theta)[grid->j_ind(thetal)]) 
      cout << " j index not right - thetal = " 
    	   << thetal << " j = " << grid->j_ind(thetal) << " theta(j) = " 
    	   << (*theta)[grid->j_ind(thetal)] << endl;
    if (phil != (*phi)[grid->k_ind(phil)]) 
      cout << " k index not right - phil = " << phil 
    	   << " k = " << grid->k_ind(phil) << " phi(k) = " 
    	   << (*phi)[grid->k_ind(phil)] << endl;
    int i = grid->i_ind(rl);
    int j = grid->j_ind(thetal);
    int k = grid->k_ind(phil);
    // double psil = psi(i,j,k) + 1.0;
    // double psi4 = psil*psil*psil*psil;
    double e2q = exp(2.0 * q(i,j,k));
    // return psi4 * e2q - 1.0;
    return e2q - 1.0;
  };
  double h_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_tt_analytical(double rl, double thetal, double phil, double t) {
    int i = grid->i_ind(rl);
    int j = grid->j_ind(thetal);
    int k = grid->k_ind(phil);
    double e2q = exp(2.0 * q(i,j,k));
    // return psi4 * e2q - 1.0;
    return e2q - 1.0;
  };
  double h_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_pp_analytical(double rl, double thetal, double phil, double t) {
    return 0.0;
  };
  double phi_analytical(double rl, double thetal, double phil, double tl) {
    int i = grid->i_ind(rl);
    int j = grid->j_ind(thetal);
    int k = grid->k_ind(phil);
    double psil = psi(i,j,k) + 1.0;
    // return 0.0;
    return log(psil);
  };
  //================================================
  // Analytical solution for connection coefficients CHECK THESE!!!
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; 
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double em2q = exp(- 2.0 * q(i,j,k));
    return (1.0 - em2q) / r;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = false; 
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double em2q = exp(- 2.0 * q(i,j,k));
    return (1.0 - em2q) * cos(theta) / sin(theta);
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
    double psil = psi(i,j,k) + 1.0;
    // return 1.0/(psil*psil);
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
  void Compute_u(gf3d & u) {
    for (int i = 0; i < n_r-N_g; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = u.r(i);
	  // const double r2 = (rl/lambda)*(rl/lambda);
	  const double sintheta = u.sintheta(j);
	  const double costheta = u.costheta(j);
	  //	  const double sin2theta = sintheta*sintheta;
	  const double rhol = rl * sintheta / lambda;
	  const double zl = rl * costheta/ lambda;
	  const double delta_rho = rhol - rho0 / lambda;
	  // u[i][j][k] = 0.5 * a * ( 1.0 - (6.0*r2 - 2.0*r2*r2)*sin2theta) * exp(-r2);
	  u[i][j][k] = 0.5 * a * (1.0 - 2.0 * rhol * rhol - 4.0 * rhol * delta_rho + 
				  2.0 * rhol * rhol * ( delta_rho*delta_rho + zl*zl )) *
	    exp(- delta_rho*delta_rho - zl*zl) / (lambda*lambda);
	}
  };
  void Compute_q(gf3d & q) {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = u.r(i);
	  // const double r2 = (rl/lambda)*(rl/lambda);
	  const double sintheta = u.sintheta(j);
	  const double costheta = u.costheta(j);
	  //	  const double sin2theta = sintheta*sintheta;
	  //	  const double cos2theta = costheta*costheta;
	  const double rho_l = rl * sintheta / lambda;
	  const double z_l = rl * costheta / lambda;
	  const double delta_rho_l = rho_l - rho0 / lambda;
	  q[i][j][k] = a * rho_l * rho_l * exp(- delta_rho_l*delta_rho_l - z_l*z_l);
	}
  };
  void Update_Psi(gf3d & psi, gf3d & delta_psi) {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  psi[i][j][k] += delta_psi(i,j,k);
	}
  };
  double Residual(gf3d & psi, gf3d & u, gf3d & res) {
    for (int i = N_g; i < n_r-N_g; i++) 
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  res[i][j][k] = psi.Laplace(i,j,k) + u(i,j,k) * psi(i,j,k) + u(i,j,k);
	}
    return res.L2_norm();
  };
};

