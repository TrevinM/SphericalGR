// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Brill Wave initial data 
//
// Uses function
//
//  q = A (rho/lambda)^2 exp(-((rho - rho_0)^2 - z^2)/lambda^2)
//    = A (r/lambda)^2 sin^2 theta exp( - ((r sin theta - rho_0)^2 - r^2 cos^2 theta)/lambda^2 )
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
  gf3d psi, q, res, u;
  int n_r, n_theta, n_phi;
  VecDoub *r, *Delta_r, *theta, *x, *phi; 
  double dtheta, dphi, c;
  int N_g;
  int max_it;   // parameters for elliptic solver
  double tol;
  double eps;   // small number, to make sure that conversion from coordinates to indices works
public:
  //================================================
  // Constructor
  //================================================
  Brill(char * indata_input, Cosmology *cosmology) :
    InData(cosmology), N_g(2) {
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
  const char * Name() { return "Brill wave initial data"; };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { 
    cout << " Initializing Brill waves... " << endl;
    eps = 1.e-10;
    n_r = fct.dim1();
    n_theta = fct.dim2();
    n_phi = fct.dim3();
    r = fct.r_pointer();
    Delta_r = fct.Delta_r_pointer();
    theta = fct.theta_pointer();
    x = fct.x_pointer();
    phi = fct.phi_pointer();
    c = fct.log_factor();
    dtheta = fct.theta(1) - fct.theta(0);
    dphi = fct.phi(1) - fct.phi(0);
    //
    // set up grid functions
    //
    psi.setup(n_r, n_theta, n_phi, c, r, Delta_r, theta, x, phi, 1);
    u.setup(n_r, n_theta, n_phi, c, r, Delta_r, theta, x, phi, 1);
    q.setup(n_r, n_theta, n_phi, c, r, Delta_r, theta, x, phi, 1);
    res.setup(n_r, n_theta, n_phi, c, r, Delta_r, theta, x, phi, 1);
    //
    for (int i = 0; i < n_r; i++) {
      //      cout << " i = " << i << " r(i) = " << (*r)[i] << " ind = " << i_ind((*r)[i]) << endl;
    }
    //
    // compute source functions and q
    //
    Compute_u(u);
    Compute_q(q);
#ifndef NoEllSolver
    //
    // Allocate and set up elliptic solver
    //
    ellsolver = new FlatEllSolver3D(n_r, n_theta, n_phi, c,
				    r, Delta_r, theta, phi);
    ellsolver->SetupSolver(1.0,u);
    ellsolver->SetRHS(-1.0,u);
    //
    // Solve and get solution
    //
    int num_it;
    double res_tri = ellsolver->Solve(max_it,num_it,tol);
    cout << " Solved Brill equation to residual " << res_tri << " in " << num_it << " iterations. " << endl;
    ellsolver->GetSolution(psi);
    cout << " Our residual = " << Residual(psi,u,res) << endl;
    for (int i = N_g; i < n_r; i++) {
      //      for (int j = N_g; j < n_theta-N_g; j++) {
    	int k = 2;
	int j = 2;
	//	cout << i << "  " << (*r)[i] << "  " << " psi = " << setw(16) << psi(i,j,k) << " res = " << setw(16) << res(i,j,k) << endl;
    }
    delete ellsolver;
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
    // reality check:
    if (rl != (*r)[i_ind(rl)]) cout << " i index not right - rl = " 
				    << rl << " i = " << i_ind(rl) << " r(i) = " 
				    << (*r)[i_ind(rl)] << endl;
    if (thetal != (*theta)[j_ind(thetal)]) cout << " j index not right - thetal = " 
						<< thetal << " j = " << j_ind(thetal) << " theta(j) = " 
						<< (*theta)[j_ind(thetal)] << endl;
    if (phil != (*phi)[k_ind(phil)]) cout << " k index not right - phil = " << phil 
					  << " k = " << k_ind(phil) << " phi(k) = " 
					  << (*phi)[k_ind(phil)] << endl;
    int i = i_ind(rl);
    int j = j_ind(thetal);
    int k = k_ind(phil);
    double psil = psi(i,j,k) + 1.0;
    double psi4 = psil*psil*psil*psil;
    double e2q = exp(2.0 * q(i,j,k));
    return psi4 * e2q - 1.0;
  };
  double h_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_tt_analytical(double rl, double thetal, double phil, double t) {
    int i = i_ind(rl);
    int j = j_ind(thetal);
    int k = k_ind(phil);
    double psil = psi(i,j,k) + 1.0;
    double psi4 = psil*psil*psil*psil;
    double e2q = exp(2.0 * q(i,j,k));
    return psi4 * e2q - 1.0;
  };
  double h_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_pp_analytical(double rl, double thetal, double phil, double t) {
    int i = i_ind(rl);
    int j = j_ind(thetal);
    int k = k_ind(phil);
    double psil = psi(i,j,k) + 1.0;
    double psi4 = psil*psil*psil*psil;
    return psi4 - 1.0;
  };
  double phi_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
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
    // int i = i_ind(rl);
    // int j = j_ind(thetal);
    // int k = k_ind(phil);
    // double psil = psi(i,j,k) + 1.0;
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
  double rho_0_analytical(double r, double theta, double phi, double t) {
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
	    exp(- delta_rho*delta_rho - zl*zl);
	}
  };
  void Compute_q(gf3d & q) {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = u.r(i);
	  const double r2 = (rl/lambda)*(rl/lambda);
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
  double Residual(gf3d & psi, gf3d & u, gf3d & res) {
    for (int i = N_g; i < n_r-N_g; i++) 
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  res[i][j][k] = psi.Laplace_so(i,j,k) + u(i,j,k) * psi(i,j,k) + u(i,j,k);
	}
    return res.L2_norm();
  };
  inline int i_ind(double r) { return res.i_ind(r); }
  //
  // CHECK: hardcodying equatorial symmetry!!!
  //
  // inline int j_ind(double theta) { 
  //   if (theta > PI/2.0) theta = PI - theta;
  //   return int(theta/dtheta + 1.5 + eps); 
  // }
  inline int j_ind(double theta) { return int(theta/dtheta + 1.5 + eps); }
  inline int k_ind(double phi) { return int(phi/dphi + 1.5 + eps); }

};

