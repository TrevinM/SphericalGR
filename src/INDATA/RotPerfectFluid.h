// Tell emacs that this is -*-c++-*- mode
//
//================================================
// 
// Rotating Perfect Fluid
//
// Note on rescaling...  In this version all quantities are (geometrically)
// rescaled
//
// Note: we're using S^phi here; rest of code uses S_phi.  
//
//
// Note: for radiation fluid have tau = rho_ADM.  For *static* radiation
// fluid also have rho = rho_ADM = tau, but not for the rotating configurations
// considered here.  We still have S_ADM_i = S_i, though.  Will set up 
// initial data for given conserved fluid quantities S^i and tau, and then
// use a recovery step in the end to compute the primitive fluid variables.
//
//================================================
//
class RotPerfectFluid : public InData {
private:
#ifndef NoEllSolver
  FlatEllSolver3D * laplace;
  VecLaplace * veclaplace;
#endif
  gf3d psi, rho, res, u, delta_psi, tau, delta_w;
  gf3d w_phi, s_phi;   // both defined with index upstairs! 
  gf3d v_phi;          // index upstairs
  gf3d a_rp, a_tp;     // indices downstairs, rescaled!
  int n_r, n_theta, n_phi;
  //  VecDoub *r, *Delta_r, *theta, *costheta, *x, *x_prime, *x_dprime, *phi; 
  // double dtheta, dphi, c;
  int N_g;
  int max_it;   // parameters for elliptic solver
  double tol;
  double tiny;  // small number, to make sure that conversion from coordinates to indices works
  double eta, Omega;   // parameter for initial data
  double r_0;          // ibid...
  double PI;
  ostringstream indata_name;
public:
  //================================================
  // Constructor
  //================================================
  RotPerfectFluid(char * indata_input, Grid * grid_i, Cosmology * cosmology) :
    InData(grid_i, cosmology) {
    N_g = grid->N_ghosts();
    indata_type = rotperfectfluid;
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
    infile.get(buf,100,'='); infile.get(c); infile >> Kappa;
    infile.get(buf,100,'='); infile.get(c); infile >> eta;
    infile.get(buf,100,'='); infile.get(c); infile >> Omega;
    infile.get(buf,100,'='); infile.get(c); infile >> r_0;
    infile.get(buf,100,'='); infile.get(c); infile >> max_it;
    infile.get(buf,100,'='); infile.get(c); infile >> tol;
    cout << " INDATA: Will set up Rotating Perfect Fluid Initial data for Kappa = " << Kappa << endl;
    cout << " INDATA: ...  with eta = " << setprecision(16) << eta 
	 << ", Omega = " << Omega << " and r_0 = " << r_0 << endl;
    cout << " INDATA: Will run elliptic solver with max_it = " << max_it 
	 << " and tol = " << tol << endl;
    cout << " INDATA: Will initialize lapse to precollapsed lapse: alpha = psi^{-2} " << endl;
    cout << "===================================================" << endl;
    PI = acos(-1.0);
  };
  //================================================
  // Destructor
  //================================================
  ~RotPerfectFluid() {};
  string Name() { 
    indata_name << "rotating perfect fluid initial data for Kappa = " << Kappa << " with eta = " << setprecision(16) << eta 
		<< " and Omega = " << Omega;
    return indata_name.str();
  };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) { 
    cout << " ROTPERFECTFLUID: Initializing rotating radiation fluid data... " << endl;
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
    // dtheta = fct.theta(1) - fct.theta(0);
    // dphi = fct.phi(1) - fct.phi(0);
    //
    // set up grid functions
    //
    int gf_counter = 2000;
    psi.setup(grid, 1);
    delta_psi.setup(grid, 1);
    rho.setup(grid, 1); 
    res.setup(grid, 1);
    u.setup(grid, 1);
    tau.setup(grid, 1);
    w_phi.setup(grid, 3, "w_phi", gf_counter++, -1, -1, +1);
    delta_w.setup(grid, 3, "delta_w", gf_counter++, -1, -1, +1);
    s_phi.setup(grid, 1, "s_phi", gf_counter++, -1, -1, +1 );
    a_rp.setup(grid, 1, "a_rp",gf_counter++, +1, -1, +1 );
    a_tp.setup(grid, 1, "a_tp",gf_counter++, -1, +1, -1 );
    v_phi.setup(grid, 1, "v_phi", gf_counter++, -1, -1, +1);
    //
    //================================================
    // compute source functions and initialize psi etc...
    //================================================
    //
    psi.equals(1.0);
    int it_max = 20;
    int it = 0;
    // to get iteration started...
    double residual = 10.0 * tol;
#ifndef NoEllSolver
    while (it < it_max && residual > tol) {
      it++;
      Solve_Hamiltonian();
      Solve_momentum();
      Compute_A_ij();
      residual = Hamiltonian_Residual();
      cout << " ROTPERFECTFLUID: residual = " << residual << " after " << it 
	   << " loops. " << endl;
    }
    if (it < it_max) {
      cout << " ROTPERFECTFLUID: iteration converted to residual = " 
	   << residual << endl;
      int error = Recovery();
      cout << " ROTPERFECTFLUID: found " << error << " errors in recovery step."
	   << endl;
      return true;
    } else {
      cout << " ROTPERFECTFLUID: iteration did not converge!!!" << endl;
      return false;
    }
#else
    cout << " ROTPERFECTFLUID: Can't construct Rotating Perfect Fluid initial data without Elliptic Solvers!!" << endl;
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
    VecDoub * r = grid->r();
    VecDoub * theta = grid->theta();
    VecDoub * phi = grid->phi();
    if (rl != (*r)[grid->i_ind(rl)]) 
      cout << " i index not right - rl = " << rl << " i = " << grid->i_ind(rl) << " r(i) = " 
	   << (*r)[grid->i_ind(rl)] << endl;
    if ( !( (thetal == (*theta)[grid->j_ind(thetal)]) || ((*theta)[grid->j_ind(thetal)] == PI - thetal)) ) 
      cout << " j index not right - thetal = " << thetal << " PI - thetal = " << PI - thetal <<  " j = " << grid->j_ind(thetal) << " theta(j) = " 
	   << (*theta)[grid->j_ind(thetal)] << endl;
    if (phil != (*phi)[grid->k_ind(phil)]) 
      cout << " k index not right - phil = " << phil << " k = " << grid->k_ind(phil) << " phi(k) = " 
	   << (*phi)[grid->k_ind(phil)] << endl;
    int i = grid->i_ind(rl);
    int j = grid->j_ind(thetal);
    int k = grid->k_ind(phil);
    double psil = psi(i,j,k);
    return log(psil);
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
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k);
    const double psi6 = psil*psil*psil*psil*psil*psil;
    return a_rp(i,j,k)/psi6;
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k);
    const double psi6 = psil*psil*psil*psil*psil*psil;
    return a_tp(i,j,k)/psi6;
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
    double psil = psi(i,j,k);
    return 1.0/(psil*psil);
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
    return rho(i,j,k);
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
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    return v_phi(i,j,k);  
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
  //====================================================
  // Solve Hamiltonian constraint
  //====================================================
  void Solve_Hamiltonian() {
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
      laplace->SetupSolver(1.0,u);
      laplace->SetRHS(-1.0,res);
      //
      // Solve and get solution
      //
      int num_it;
      laplace->Solve(max_it,num_it,tol);
      //	cout << " Solved Hamiltonian constraint to residual " << res_tri 
      //	     << " in " << num_it << " iterations. " << endl;
      laplace->GetSolution(delta_psi);
      update_psi();
      Compute_tau();    
      residual = Hamiltonian_Residual();
      // cout << " ROTPERFECTFLUID: Hamiltonian constraint residual after " 
      // 	   << step << " steps: " << residual << endl;
    }
    delete laplace;
#else
    cout << " Can't construct rotating radiation fluid initial data without an Elliptic Solver!! " << endl;
#endif
  }    
  //====================================================
  // Solve momentum constraint
  //
  // Momentum constraint is linear in W^phi, but our solver is only 2nd order.
  // Therefore will solve iteratively anyway...
  //====================================================
  void Solve_momentum() {
#ifndef NoEllSolver
#ifdef VECLAPLACEINSTALLED
    //
    // Allocate and set up elliptic solver
    //
    veclaplace = new VecLaplace(grid);
    veclaplace->SetFallOff(3);
    int max_step = 20;
    Compute_S();    
    int step = 0;
    double residual = momentum_Residual();
    while (residual > tol && step < max_step) {
      step++;
      veclaplace->SetupSolver(0.0,u);
      veclaplace->SetRHS(-1.0,res);
      //
      // Solve and get solution
      //
      int num_it;
      double res_tri = veclaplace->Solve(max_it,num_it,tol);
      veclaplace->GetSolution(delta_w);
      update_W();
      residual = momentum_Residual();
      //      cout << " ROTPERFECTFLUID: momentum constraint residual after " 
      //	    << step << " steps : " << residual << endl;
    }
    delete veclaplace;
#endif
#else
    cout << " Can't construct rotating radiation fluid initial data without an Elliptic Solver!! " << endl;
#endif
  }    
  //====================================================
  // Compute *rescaled* extrinsic curvature from W^phi
  // NOTE: this is conformally rescaled according to CTT decomposition;
  // will translate to BSSN rescaling in analytical functions...
  //====================================================
  void Compute_A_ij() {
    for (int i = N_g; i < n_r - N_g; i++) 
      for (int j = N_g; j < n_theta - N_g; j++) 
	for (int k = N_g; k < n_phi - N_g; k++) {
	  const double rl = rho.r(i);
	  const double sintheta = rho.sintheta(j);
	  const double costheta = rho.costheta(j);
	  const double cottheta = costheta / sintheta;
	  a_rp[i][j][k] =   w_phi.dr(i,j,k) - w_phi(i,j,k)/rl;
	  a_tp[i][j][k] = ( w_phi.dtheta(i,j,k) - cottheta * w_phi(i,j,k) ) / rl;
	}
  }
  //====================================================
  // Recovery: compute primitive fluid variables
  //====================================================
  int Recovery() {

    int error = 0;
    for (int i = N_g; i < n_r; i++) 
      for (int j = N_g; j < n_theta - N_g; j++) 
	for (int k = N_g; k < n_phi - N_g; k++) {
	  //
	  // compute S2 and tau2 - recall - s^phi is rescaled
	  //
	  // const double rl = rho.r(i);
	  // const double r2 = rl*rl;
	  // const double sintheta = rho.sintheta(j);
	  // const double sin2theta = sintheta*sintheta;
	  const double psil = psi(i,j,k);
	  const double psi4 = psil*psil*psil*psil;
	  const double S_p = s_phi(i,j,k);
	  const double S2 = psi4*S_p*S_p;
	  const double taul = tau(i,j,k);
	  const double tau2 = taul*taul;
	  //
	  // compute discriminant and gamma factor
	  //
	  // if (i == 6 && j == 6 && k == 2) {
	  //   const double s_phi_low = rl*sintheta*psi4*S_p * psi4*psil*psil;
	  //   cout << " S_p in INDATA : " << s_phi_low << endl;
	  // }
	  const double opk = 1.0 + Kappa;
	  const double opk2 = opk * opk;
	  const double opk4 = opk2 * opk2;
	  const double temp1 = opk4*tau2*tau2 - 4.0*opk2*Kappa*tau2*S2;
	  if (temp1 < 0.0) {
	    cout << " ROTPERFECTFLUID: Trouble in Recovery: temp1 < 0.0 " << endl;
	    error++;
	  }
	  const double temp2 = opk2*tau2 - 2.0*Kappa*S2 + sqrt(temp1);
	  double v2 = 0.0;
	  if (temp2 > 0.0)
	    v2 = 2.0*S2/temp2;
	  else
	    v2 = 0.0;
	  if (v2 >= 1.0) {
	    cout << " ROTPERFECTFLUID: Trouble in Recovery: v2 >= 1.0 " << endl;
	    error++;
	  }
	  const double Wl = 1.0/sqrt( 1.0 - v2 ); 
	  //
	  // compute primitive variables, including *rescaled* version of 
	  // v^phi
	  //
	  const double pl = Kappa * taul/( opk * Wl*Wl - Kappa);
	  if (!isfinite(pl)) cout << v2 << "    "
				<< Wl << "    "
				<< taul << "    "
				<< endl;
	  rho[i][j][k] = pl/Kappa;
	  const double taup = taul + pl;
	  if (taup > 0.0) {
	    v_phi[i][j][k] = S_p / taup;
	  } else
	    v_phi[i][j][k] = 0.0;
	  //
	  if (v_phi(i,j,k) < 0.0) 
	    cout << " v_p < 0 " << S_p << "  " << taup << endl;
	  // check...
	  // if (i == 6 && j == 6 && k == 2) { 
	  //   const double S_p_up = Wl*Wl * 4.0/3.0 * rho[i][j][k] * v_phi[i][j][k] 
	  //     / (rl * sintheta);
	  //   const double S_p_down = rl*sintheta*psi4*S_p_up * psi4*psil*psil;
	  //   cout << " S_p_down = " << S_p_down << endl;
	  // }
	}
    return error;
  }
  //
  // computes rho_ADM = tau; appears on right-hand side of Hamiltonian
  //
  void Compute_tau() {
    const double PIm32 = 1.0/(PI*sqrt(PI));
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = rho.r(i);
	  const double psil = psi(i,j,k);
	  const double r_areal = rl*psil*psil;
	  const double r_rel = r_areal / r_0;
	  tau[i][j][k] = 0.5 * PIm32 * eta * exp( - r_rel * r_rel)/(r_0*r_0);
	}
  };
  //
  // computes s_phi; appears on right-hand side of momentum constraint (index upstairs)
  // NOTE: assumes that Compute_tau has been called already
  //
  // NOTE: v_tilde is nothing physical (and not the same as v_phi) - just 
  // used to compute some physically reasonable s_phi
  //
  void Compute_S() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = rho.r(i);
	  const double sintheta = rho.sintheta(j);
	  const double psil = psi(i,j,k);
	  const double r_areal = rl*psil*psil;
	  const double r_rel = r_areal / r_0;
	  const double denominator = 1.0 + r_rel*r_rel;
	  const double v_tilde_phi = Omega / denominator;
	  // don't forget to rescale...
	  s_phi[i][j][k] = (1.0 + Kappa) * tau(i,j,k) * v_tilde_phi * rl * sintheta;
	}
  };

  void Compute_u() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  const double rl = rho.r(i);
	  const double psil = psi(i,j,k);
	  const double psi4 = psil*psil*psil*psil;
	  const double psim8 = 1.0/(psi4*psi4);
	  const double a_rpl = a_rp(i,j,k);
	  const double a_tpl = a_tp(i,j,k);
	  const double a2 = 2.0*(a_rpl*a_rpl + a_tpl*a_tpl);
	  u[i][j][k] = 10.0*PI*psi4*tau(i,j,k) 
	    - 8.0*PI*psi4*psi4*tau(i,j,k)*rl*rl/(r_0*r_0)
	    - 7.0/8.0*psim8*a2;
	}
  };
  void update_psi() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  psi[i][j][k] += delta_psi(i,j,k);
	}
  };
  void update_W() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  w_phi[i][j][k] += delta_w(i,j,k);
	}
  };
  double Hamiltonian_Residual() {
    for (int i = N_g; i < n_r-N_g; i++) 
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  const double psil = psi(i,j,k);
	  const double psi5 = psil*psil*psil*psil*psil;
	  const double psim7 = 1.0/(psi5*psil*psil);
	  const double a_rpl = a_rp(i,j,k);
	  const double a_tpl = a_tp(i,j,k);
	  const double a2 = 2.0*(a_rpl*a_rpl + a_tpl*a_tpl);
	  res[i][j][k] = psi.Laplace(i,j,k) 
	    + 2.0*PI*psi5*tau(i,j,k)
	    + psim7*a2/8.0;
	}
    return res.L2_norm();
  };
  double momentum_Residual() {
    for (int i = N_g; i < n_r-N_g; i++) 
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  const double rl = rho.r(i);
	  const double r2 = rl*rl;
	  const double sintheta = rho.sintheta(j);
	  const double sin2theta = sintheta*sintheta;
	  const double costheta = rho.costheta(j);
	  const double cottheta = costheta / sintheta;
	  const double psil = psi(i,j,k);
	  const double psi5 = psil*psil*psil*psil*psil;
	  const double psi10 = psi5*psi5;
	  res[i][j][k] = w_phi.ddr(i,j,k)
	    + 2.0 * w_phi.dr(i,j,k)/rl
	    + w_phi.ddtheta(i,j,k)/r2
	    + cottheta * w_phi.dtheta(i,j,k)/r2
	    - w_phi(i,j,k)/(r2*sin2theta) 
	    - 8.0 * PI * psi10 * s_phi(i,j,k);
	}
    return res.L2_norm();
  };
};


