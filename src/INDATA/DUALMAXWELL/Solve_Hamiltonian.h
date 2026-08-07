// Tell emacs that this is -*-c++-*- mode
//====================================================
// Solve Hamiltonian constraint
//====================================================
bool Solve_Hamiltonian(double tol_tri = 1.e-10, double tol_res = 1.e-8,
		       bool verbose = true) {
#ifndef NoEllSolver
  laplace = new FlatEllSolver3D(grid);
  //
  // Allocate and set up elliptic solver
  //
  int max_step = 50;
  int step = 0;
  Compute_Fields();
  Compute_Sources();    
  double residual = Hamiltonian_Residual();
  if (verbose)
    cout << " DUALEMWAVE: initial Hamiltonian residual: " << residual << endl;
  while (residual > tol_res && step < max_step) {
    step++;
    Compute_u();
    laplace->SetupSolver(1.0,u);
    laplace->SetRHS(-1.0,res);
    //
    // Solve and get solution
    //
    int num_it;
    laplace->Solve(max_it,num_it,tol_tri);
    // if (verbose)
    //   cout << " DUALEMWAVE: Solved Hamiltonian constraint to residual "
    // 	   << res_tri << " in "
    // 	   << num_it << " iterations. " << endl;
    laplace->GetSolution(delta_psi);
    update_psi();
    Compute_Fields();
    Compute_Sources();
    residual = Hamiltonian_Residual();
    if (verbose)
      cout << " DUALEMWAVE: Hamiltonian residual after "
	   << step << " steps: " << residual << endl;
  }
  delete laplace;
  if (residual < tol_res) {
    cout << " DUALEMWAVE - Hamiltonian residual converged to " << residual
	 << " in " << step << " steps" << endl;
    return true;
  } else {
    cout << " DUALEMWAVE - Hamiltonian constraint did not converge! " << endl;
    cout << "              Residual " << residual << " after " << step
	 << " steps " << endl;
    return false;
  }
#else
  cout << " DUALEMWAVE: Can't construct EM initial data without an Elliptic Solver!! " << endl;
  return false;
#endif
}

//
void Compute_u() {
  for (int i = 0; i < n_r; i++) 
    for (int j = 0; j < n_theta; j++) 
      for (int k = 0; k < n_phi; k++) {
	// const double rl = rho.r(i);
	const double psil = psi(i,j,k);
	const double psi4 = psil*psil*psil*psil;
	const double psim8 = 1.0/pow(psil, 8);
	double factor = 0.0;
	factor = 5 + 4 + 2*n_psi;
	u[i][j][k] = 2.0*PI*psi4*rho(i,j,k) * factor
	  - 7.0/8.0*psim8*A2(i,j,k);
      }
};
void update_psi() {
  for (int i = 0; i < n_r; i++) 
    for (int j = 0; j < n_theta; j++) 
      for (int k = 0; k < n_phi; k++) {
	psi[i][j][k] += 0.6 * delta_psi(i,j,k);
      }
};
double Hamiltonian_Residual() {
  for (int i = N_g; i < n_r-N_g; i++) 
    for (int j = N_g; j < n_theta-N_g; j++) 
      for (int k = N_g; k < n_phi-N_g; k++) {
	const double psil = psi(i,j,k);
	const double psi5 = psil*psil*psil*psil*psil;
	const double psim7 = 1.0/pow(psil, 7);
	res[i][j][k] = psi.Laplace(i,j,k) 
	  + 2.0*PI*psi5*rho(i,j,k) + psim7*A2(i,j,k)/8.0;
      }
  return res.L2_norm();
};
