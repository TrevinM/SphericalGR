// Tell emacs that this is -*-c++-*- mode
//
//================================================
// Linear Wave Initial Data, based on Nakamura 1984
// (see also Shibata & Nakamura, PRD 52, 5428 (1995)):
//
// constructs analytical solution to momentum constraint, then solves
// Hamiltonian constraint assuming conformal flatness
//
// NOTE: we construct m = 0 waves here, unlike the m = +- 2 waves in the
// above references
// 
//================================================
//
class ShibataWave : public InData {
private:
  double C_amp;   
  double lambda;
  double r0;
  int l;
  int m;
  int n_r, n_theta, n_phi, N_g;
  int family;
  double alpha, delta_T, R_0;  // parameters for shift
  double tol;
#ifndef NoEllSolver
  FlatEllSolver3D * ellsolver;
#endif
  gf3d psi, A2, res, delta_psi, u;
  //
  ostringstream indata_name;
public:
  //================================================
  // Constructor
  //================================================
  ShibataWave(char * indata_input, Grid * grid_i, Cosmology *cosmology) :
    InData(grid_i, cosmology) {
    N_g = grid->N_ghosts();
    indata_type = linwave;
    analytical = false;
    ifstream infile;
    infile.open(indata_input);
    if (!infile)
      cerr << "Can't open " << indata_input 
	   << " for input. This is bad. " << endl;
    else 
      cout << " Reading initial data parameters from file " 
	   << indata_input << endl;
    char buf[200], c;
    infile.get(buf,200,'='); infile.get(c); infile >> C_amp;
    infile.get(buf,200,'='); infile.get(c); infile >> lambda;
    infile.get(buf,200,'='); infile.get(c); infile >> r0;
    infile.get(buf,200,'='); infile.get(c); infile >> l;
    infile.get(buf,200,'='); infile.get(c); infile >> m;
    infile.get(buf,200,'='); infile.get(c); infile >> family;
    infile.get(buf,200,'='); infile.get(c); infile >> alpha;
    infile.get(buf,200,'='); infile.get(c); infile >> delta_T;
    infile.get(buf,200,'='); infile.get(c); infile >> R_0;
    infile.get(buf,200,'='); infile.get(c); infile >> tol;
    if ( (l == 2 || l == 4) && (family == 1 || family == 2 ) ) {
      cout << " SHIBATA WAVE: Will set up family " << family
	   << " for l = " << l
	   << " with amplitude C = " << C_amp << endl;
    } else {
      cout << " SHIBATA WAVE: family " << family << " for  l = "
	   << l << " not implemented " << endl;
      exit(0);
    }
    cout << " Building shift from alpha = " << alpha << " delta_T = "
	 << delta_T << " R_0 = " << R_0 << endl;
    cout << " Iterating to tolerance tol = " << tol << endl;
    cout << "===================================================" << endl;
    indata_name << "l = " << l << " family " << family
		<< " Nakamura wave initial data with C = " << C_amp;
  };
  //================================================
  // Destructor
  //================================================
  ~ShibataWave();
  string Name() { return indata_name.str(); };
  //================================================
  // Initializer
  //================================================
  bool Initialize(gf3d & fct) {
    cout << " SHIBATA WAVE: initializing..." << endl;
    n_r = fct.dim1();
    n_theta = fct.dim2();
    n_phi = fct.dim3();
    //
    // set up grid functions
    // 
    psi.setup(grid, 1);
    delta_psi.setup(grid, 1);
    A2.setup(grid, 1);
    res.setup(grid, 1);
    u.setup(grid, 1);
    //
    // compute A2
    // 
    Compute_A2();
    //
    // initialize psi...
    //
    psi.equals(0.0);
    //
    // set up iteration:
    //
    double residual = Residual();
#ifndef NoEllSolver
    ellsolver = new FlatEllSolver3D(grid);
    int it = 0;
    int it_max = 100;
    while (residual > tol && it < it_max) {
      Compute_u();
      ellsolver->SetupSolver(-1.0, u);
      ellsolver->SetRHS(-1.0, res);
      int max_it = 100;
      int num_it;
      ellsolver->Solve(max_it, num_it, tol/100.0);
      ellsolver->GetSolution(delta_psi);
      Update_Psi();
      residual = Residual();
      it++;
      cout << " SHIBATA WAVE: Residual after " << it << " steps = " << residual << endl;
    }
    delete ellsolver;
#else
    cout << " SHIBATA WAVE: can't construct waves without elliptic solver!! " << endl;
    return false;
#endif
    return true; 
  } 
  //================================================
  // Analytical solution for h_{ij}
  //================================================
  double h_rr_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_rt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_rp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_tt_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_tp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double h_pp_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  double phi_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    double psil = psi(i,j,k) + 1.0;
    return log(psil);
  };
  //================================================
  // Analytical solution for connection coefficients
  //================================================
  double lam_r_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true;
    return 0.0;
  }
  double lam_t_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true;
    return 0.0;
  }
  double lam_p_analytical(double r, double theta, double phi, double t, 
			  bool & done) {
    done = true;
    return 0.0;
  }
  //================================================
  // Analytical solution for extrinsic curvature
  //================================================
  double a_rr_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k) + 1.0;
    const double psim6 = 1.0/(psil*psil*psil*psil*psil*psil);
    return psim6*A_rr(r,theta,phi);
  };
  double a_rt_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k) + 1.0;
    const double psim6 = 1.0/(psil*psil*psil*psil*psil*psil);
    return psim6*A_rt(r,theta,phi);
  };
  double a_rp_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k) + 1.0;
    const double psim6 = 1.0/(psil*psil*psil*psil*psil*psil);
    return psim6*A_rp(r,theta,phi);
  };
  double a_tt_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k) + 1.0;
    const double psim6 = 1.0/(psil*psil*psil*psil*psil*psil);
    return psim6*A_tt(r,theta,phi);
  };
  double a_tp_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k) + 1.0;
    const double psim6 = 1.0/(psil*psil*psil*psil*psil*psil);
    return psim6*A_tp(r,theta,phi);
  };
  double a_pp_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    const double psil = psi(i,j,k) + 1.0;
    const double psim6 = 1.0/(psil*psil*psil*psil*psil*psil);
    return psim6*A_pp(r,theta,phi);
  };
  double K_analytical(double r, double theta, double phi, double t) {
    return 0.0;
  };
  //================================================
  // Analytical solution for gauge
  //================================================
  double lapse_analytical(double r, double theta, double phi, double t) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = grid->k_ind(phi);
    double psil = psi(i,j,k) + 1.0;
    // return 1.0/(psil*psil);
    return 1.0;
  };
  double shift_r_analytical(double r, double theta, double phi, double t) {
    double factor = R_0*R_0 / (R_0*R_0 + r*r);
    return - alpha * r * factor / delta_T;
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
  // rescaled spherical polar components of \hat A_{ij} for m=0 version of model (A)
  // (compare (3.7) in Shibata & Nakamura)
  //================================================
  inline double A_rr(double r, double theta, double phi) {
    const double r2 = r*r;
    if (l == 2) {
      if (family == 1) {
	const double costheta = cos(theta);
	return C_amp * exp(-r2/2.0) * (1. - 3.0 * costheta * costheta);
      } else if (family == 2) {
	const double costheta = cos(theta);
	return C_amp * exp(-r2/2.0) * (5. - r2) *
	  (1. - 3.0 * costheta * costheta);
      }
    } else if (l == 4) {
      const double cos2theta = cos(2.*theta);
      const double cos4theta = cos(4.*theta);
      const double r4 = r2*r2;
      const double term1 = 70 - 19*r2 + r4;
      const double term2 = 9 + 20*cos2theta + 35*cos4theta;
      return C_amp * r2 * exp(-r2/2) * term1 * term2;
    } else {
      cout << " SHIBATAWAVE: l = " << l << " not implemented!" << endl;
      exit(0);
    }
    return 0.0;
  };
  inline double A_rt(double r, double theta, double phi) {
    const double r2 = r*r;
    if (l == 2) {
      if (family == 1) {
	return C_amp * exp(-r2/2.0) * (3.0 - r2) * sin(theta) * cos(theta);
      } else if (family == 2) {
	return C_amp * exp(-r2/2.0) * (15.0 - 10.*r2 + r2*r2) *
	  sin(theta) * cos(theta);
      }
    } else if (l == 4) {
      const double cos2theta = cos(2.*theta);
      const double sin2theta = sin(2.*theta);
      const double r4 = r2*r2;
      const double r6 = r4*r2;
      const double term1 = -350 + 203*r2 - 28*r4 + r6;
      const double term2 = (1 + 7 * cos2theta) * sin2theta;
      return 2 * r2 * C_amp * exp(-r2/2) * term1 * term2;
    } else {
      cout << " SHIBATAWAVE: l = " << l << " not implemented!" << endl;
      exit(0);
    }
    return 0.0;
  };
  inline double A_rp(double r, double theta, double phi) {
    return 0.0;
  }
  inline double A_tt(double r, double theta, double phi) {
    const double r2 = r*r;
    const double r4 = r2*r2;
    if (l == 2) {
      if (family == 1) {
	const double costheta = cos(theta);
	const double sintheta = sin(theta);
	return - C_amp * exp(-r2/2.0) *
	  (2 - 6.*costheta*costheta + (6 - 8*r2 + r4)*sintheta*sintheta) / 4.0;
      } else if (family == 2) {
      const double cos2theta = cos(2.*theta);
      const double r6 = r4*r2;
      return - C_amp * exp(-r2/2.0) *
	(20. - 60.*r2 + 17.*r4 - r6 -
	 (60. - 68.*r2 + 17.*r4 - r6) * cos2theta) / 8.0;
      }
    } else if (l == 4) { 
      const double cos2theta = cos(2.*theta);
      const double cos4theta = cos(4.*theta);
      const double r6 = r4*r2;
      const double r8 = r4*r4;
      const double term1 = 3*(420 - 1410*r2 + 436*r4 - 39*r6 + r8);
      const double term2 = 4*r2*(-1296 + 430*r2 - 39*r4 + r6)*cos2theta;
      const double term3 = -7*(2100 - 1866*r2 + 460*r4 - 39*r6 + r8)*cos4theta;
      return C_amp * r2 / 6 * exp(-r2/2) * (term1 + term2 + term3);
    } else {
      cout << " SHIBATAWAVE: l = " << l << " not implemented!" << endl;
      exit(0);
    }
    return 0.0;
  };
  inline double A_tp(double r, double theta, double phi) {
    return 0.0;
  }
  inline double A_pp(double r, double theta, double phi) {
    const double r2 = r*r;
    const double r4 = r2*r2;
    if (l == 2) {
      if (family == 1) {
	const double costheta = cos(theta);
	const double sintheta = sin(theta);
	return C_amp * exp(-r2/2.0) *
	  (- 2 + 6.*costheta*costheta +
	   (6 - 8*r2 + r4)*sintheta*sintheta) / 4.0;
      } else if (family == 2) {
      const double cos2theta = cos(2.*theta);
      const double r6 = r4*r2;
	return C_amp * exp(-r2/2.0) *
	  ( 40. - 64.*r2 + 17.*r4 - r6 +
	   r2*(56 - 17.*r2 + r4) * cos2theta ) / 8.0;      
      }
    } else if (l == 4) { 
      const double cos2theta = cos(2.*theta);
      const double cos4theta = cos(4.*theta);
      const double r6 = r4*r2;
      const double r8 = r4*r4;
      const double term1 = 3*(1680 - 1752*r2 + 454*r4 - 39*r6 + r8);
      const double term2 = 4*(2100 - 1866*r2 + 460*r4 - 39*r6 + r8)*cos2theta;
      const double term3 = - 7*r2*(-1296 + 430*r2 - 39*r4 + r6)*cos4theta;
      return - C_amp * r2 / 6 * exp(-r2/2) * (term1 + term2 + term3);
    } else {
      cout << " SHIBATAWAVE: l = " << l << " not implemented!" << endl;
      exit(0);
    }
    return 0.0;
  };
  //================================================
  // Compute A2
  //================================================
  void Compute_A2() {
    for (int i = 0; i < n_r-N_g; i++)
      for (int j = 0; j < n_theta; j++)
	for (int k = 0; k < n_phi; k++) {
	  const double rl = A2.r(i);
	  const double tl = A2.theta(j);
	  const double pl = A2.phi(k);
	  const double Arrl = A_rr(rl, tl, pl);
	  const double Artl = A_rt(rl, tl, pl);
	  const double Arpl = A_rp(rl, tl, pl);
	  const double Attl = A_tt(rl, tl, pl);
	  const double Atpl = A_tp(rl, tl, pl);
	  const double Appl = A_pp(rl, tl, pl);
	  A2[i][j][k] = Arrl*Arrl + 2.0*Artl*Artl +
	    2.0*Arpl*Arpl + Attl*Attl + 2.0*Atpl*Atpl + Appl*Appl;
	}
  }
  //================================================
  // Compute residual
  //================================================
  double Residual() { 
    for (int i = N_g; i < n_r-N_g; i++) 
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  const double psil = psi(i,j,k) + 1.0;
	  const double psim7 = 1.0/(psil*psil*psil*psil*psil*psil*psil);
	  res[i][j][k] = psi.Laplace(i,j,k) + psim7*A2(i,j,k)/8.0;
	}
    return res.L2_norm();
  }
  //================================================
  // Compute function u
  //================================================
  void Compute_u() { 
    for (int i = N_g; i < n_r-N_g; i++) 
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  const double psil = psi[i][j][k] + 1.0;
	  const double psim8 = 1.0/(psil*psil*psil*psil*psil*psil*psil*psil);
	  u[i][j][k] = 7.0*psim8*A2(i,j,k)/8.0;
	}
  }
  //================================================
  // Update psi
  //================================================
  void Update_Psi() {
    for (int i = 0; i < n_r; i++) 
      for (int j = 0; j < n_theta; j++) 
	for (int k = 0; k < n_phi; k++) {
	  psi[i][j][k] += delta_psi(i,j,k);
	}
  };
  
};

