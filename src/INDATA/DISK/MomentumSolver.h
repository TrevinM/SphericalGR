  //
  //================================================
  // functions for Momentum constraint
  //================================================
  //
  double Solve_Momentum(double tol_tri = 1.e-10, double tol_res = 1.e-8) {
#ifndef NoEllSolver
    // first update A_ij (which also updates aop6 and aop6_BH)
    Compute_Aij();
    // set up right-hand sides 
    Momentum_RHS();
    double res_norm = Momentum_Residual(grad_ln_ap_r, grad_ln_ap_t, grad_ln_ap_p);
    cout << " DISK - Momentum constraint: initial residual: "
	 << res_norm << endl;
    int step = 0;
    int max_step = 1;
    while (res_norm > tol_res && step < max_step) {
      step++;
      cout << " DISK: calling Setu..." << endl;
      veclaplacian->Setu(grad_ln_ap_r, grad_ln_ap_t, grad_ln_ap_p);
      cout << " DISK: setting up solver..." << endl;
      veclaplacian->SetupSolver();
      // cout << " DISK: calling SetRHS..." << endl;
      veclaplacian->SetRHS(RHS_r, RHS_t, RHS_p);
      // veclaplacian->SetRHS(-1.0, res_r, res_t, res_p);
      int max_it = 500;
      int num_it;
      double tol = 1.0e-8;
      // cout << " DISK: solving..." << endl;
      veclaplacian->Solve(max_it, num_it, tol);
      // cout << " DISK: extracting solution..." << endl;
      veclaplacian->GetSolution(beta_r, beta_t, beta_p);
      // res_norm = Momentum_Residual(grad_ln_ap_r, grad_ln_ap_t, grad_ln_ap_p);
      // cout << " DISK - Momentum residual right after solve = " << res_norm << endl;
      // beta_p.equals(0.0);
      // compute residual
      Compute_Aij();
      Momentum_RHS();
      res_norm = Momentum_Residual(grad_ln_ap_r, grad_ln_ap_t, grad_ln_ap_p);
      cout << " DISK - Momentum residual after " << step <<
	" steps = " << res_norm << endl;
    }
    return res_norm;
#else
    cout << " DISK: cannot solve momentum constraint without elliptic solver!"
	 << endl;
    return 0.0;   // just to make compiler happy...
#endif  /* NoEllSolver */
  }
  void Momentum_RHS() {
    for (int i = N_g; i < n_r - N_g; i++) {
      const double rl = psi_D.r(i);
      for (int j = N_g; j < n_theta - N_g; j++) {
	const double sinthetal = psi_D.sintheta(j);
	const double costhetal = psi_D.costheta(j);
	const double sin2theta = sinthetal*sinthetal;
	const double cottheta = costhetal / sinthetal;
	for (int k = N_g; k < n_phi - N_g; k++) {
	  const double aop6l = aop6(i,j,k);
	  const double aop6_BHl = aop6_BH(i,j,k);
	  const double psil = psi_D(i,j,k) + psi_BH(i,j,k);
	  const double psi10 = psil*psil*psil*psil*psil*psil*psil*psil*psil*psil;
	  const double d_ln_aop6_BHdr = dlapse_BHdr(i,j,k)/lapse_BH(i,j,k)
	    - 6.0*dpsi_BHdr(i,j,k)/psi_BH(i,j,k);
	  vect grad_ln_ap_D(grad_ln_ap_D_r(i,j,k), 
			    grad_ln_ap_D_t(i,j,k),
			    grad_ln_ap_D_p(i,j,k));
	  //vect d_ln_aop6_BH(ln_aop6_BH.dr(i,j,k),
	  //			 ln_aop6_BH.dtheta(i,j,k)/rl,
	  //			 ln_aop6_BH.dphi(i,j,k)/(rl*sinthetal));
	  vect d_K(K_BH.dr(i,j,k), 0.0, 0.0);  // FIX!!!
	  vect S(0.0,0.0,0.1*q(i,j,k)*rl*sinthetal);   // FIX both q and scaling!
	  tensor L_beta_D(L_beta_rr(i,j,k), L_beta_rt(i,j,k), L_beta_rp(i,j,k),
			  L_beta_tt(i,j,k), L_beta_tp(i,j,k), L_beta_pp(i,j,k));
	  tensor A_BH(A_rr_BH(i,j,k), 0.0, 0.0, A_tt_BH(i,j,k), 0.0, A_tt_BH(i,j,k));
	  vect RHS(0.0, 0.0, 0.0);
	  for (int a = 0; a < 3; a++) {
	    RHS[a] = 0.0;
	    for (int b = 0; b < 3; b++) {
	      RHS[a] += // L_beta_D[a][b]*d_ln_aop6[b]
		+ 2.0*aop6_BHl*A_BH[a][b]*grad_ln_ap_D[b]
		+ 4./3.*lapse_D(i,j,k)*d_K[a]
		+ 16.*PI*psi10*S[a];
	    }
	  }
	  RHS_r[i][j][k] = RHS[0];
	  RHS_t[i][j][k] = RHS[1];
	  RHS_p[i][j][k] = RHS[2];
	}
      }
    }
  }
  double Momentum_Residual(gf3d & u_r, gf3d & u_t, gf3d & u_p) {
    double res_r_norm2 = 0.0;
    double res_t_norm2 = 0.0;
    double res_p_norm2 = 0.0;
    for (int i = N_g; i < n_r - N_g; i++) {
      const double rl = psi_D.r(i);
      const double r2 = rl*rl;
      for (int j = N_g; j < n_theta - N_g; j++) {
	const double sintheta = psi_D.sintheta(j);
	const double costheta = psi_D.costheta(j);
	const double sin2theta = sintheta*sintheta;
	const double cottheta = costheta / sintheta;
	for (int k = N_g; k < n_phi - N_g; k++) {
	  // compute divergence of shift (in terms of rescaled components)
	  const double div = beta_r.dr_so(i,j,k) + 2.*beta_r(i,j,k)/rl
	    + beta_t.dtheta_so(i,j,k)/rl + cottheta*beta_t(i,j,k)/rl;
	  // compute D^i \beta^j (rescaled components!)
	  const double Drbetar = beta_r.dr_so(i,j,k);
	  const double Dtbetar = (beta_r.dtheta_so(i,j,k) - beta_t(i,j,k))/rl;
	  const double Dpbetar = - beta_p(i,j,k)/rl;
	  const double Drbetat = beta_t.dr_so(i,j,k);
	  const double Dtbetat = (beta_t.dtheta_so(i,j,k) + beta_r(i,j,k))/rl;
	  const double Dpbetat = -cottheta*beta_p(i,j,k)/rl;
	  const double Drbetap = beta_p.dr_so(i,j,k);
	  const double Dtbetap = beta_p.dtheta_so(i,j,k)/rl;
	  const double Dpbetap = (beta_r(i,j,k) + cottheta*beta_t(i,j,k))/rl;
	  // compute (L beta)^ij (still rescaled!)
	  const double L_b_rr = 2.0*Drbetar - 2./3.*div;
	  const double L_b_rt = Drbetat + Dtbetar;
	  const double L_b_rp = Drbetap + Dpbetar;
	  const double L_b_tt = 2.0*Dtbetat - 2./3.*div;
	  const double L_b_tp = Dtbetap + Dpbetat;
	  const double L_b_pp = 2.0*Dpbetap - 2./3.*div;
	  
	  res_r[i][j][k] = 4./3.*beta_r.ddr_so(i,j,k) + 8./3.*beta_r.dr_so(i,j,k)/rl
	    + beta_r.ddtheta_so(i,j,k)/r2 + cottheta*beta_r.dtheta_so(i,j,k)/r2 
	    - 8./3.*beta_r(i,j,k)/r2 + 1./3.*beta_t.drdtheta_so(i,j,k)/rl
	    - 7./3.*beta_t.dtheta_so(i,j,k)/r2 - 7./3.*cottheta*beta_t(i,j,k)/r2
	    + cottheta/3.*beta_t.dr_so(i,j,k)/rl
	    - (L_b_rr * u_r(i,j,k)  + L_b_rt * u_t(i,j,k)  + L_b_rp * u_p(i,j,k) )
	    - RHS_r(i,j,k);
	  res_r_norm2 += res_r(i,j,k)*res_r(i,j,k);
	  res_t[i][j][k] = beta_t.ddr_so(i,j,k) + 2.*beta_t.dr_so(i,j,k)/rl
	    + 4./3.*beta_t.ddtheta_so(i,j,k)/r2
	    + 4.*cottheta/3.*beta_t.dtheta_so(i,j,k)/r2
	    - 4./3.*beta_t(i,j,k)/(r2*sin2theta) + 8./3.*beta_r.dtheta_so(i,j,k)/r2
	    + 1./3.*beta_r.drdtheta_so(i,j,k)/rl
	    - (L_b_rt * u_r(i,j,k)  + L_b_tt * u_t(i,j,k)  + L_b_tp * u_p(i,j,k) )
	    - RHS_t(i,j,k);
	  res_t_norm2 += res_t(i,j,k)*res_t(i,j,k);
	  res_p[i][j][k] = beta_p.ddr_so(i,j,k) + 2.*beta_p.dr_so(i,j,k)/rl
	    + beta_p.ddtheta_so(i,j,k)/r2 + cottheta*beta_p.dtheta_so(i,j,k)/r2
	    - beta_p(i,j,k)/(r2*sin2theta)
	    - (L_b_rp * u_r(i,j,k)  + L_b_tp * u_t(i,j,k)  + L_b_pp * u_p(i,j,k) )
	    - RHS_p(i,j,k);
	  res_p_norm2 += res_p(i,j,k)*res_p(i,j,k);
	}
      }
    }
    return sqrt(res_r_norm2 + res_t_norm2 + res_p_norm2);
  }
