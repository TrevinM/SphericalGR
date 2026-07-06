// Tell emacs that this is -*-c++-*- mode
//====================================================
// Solve Hamiltonian constraint
//====================================================
bool Solve_Psi(double tol_tri = 1.e-10, double tol_res = 1.e-8,
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
  double residual = Hamiltonian_Psi_Residual();
  if (verbose)
    cout << " INFLATIONGW: initial Hamiltonian Psi residual: " << residual << endl;
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
    //   cout << " INFLATIONGW: Solved Hamiltonian constraint to residual "
    // 	   << res_tri << " in "
    // 	   << num_it << " iterations. " << endl;
    laplace->GetSolution(delta_psi);
    update_psi();
    //Compute_Fields();
    //Compute_Sources();
    residual = Hamiltonian_Psi_Residual();
    if (verbose)
      cout << " INFLATIONGW: Hamiltonian Psi residual after "
	   << step << " steps: " << residual << endl;
  }
  delete laplace;
  if (residual < tol_res) {
    cout << " INFLATIONGW - Hamiltonian Psi residual converged to " << residual
	 << " in " << step << " steps" << endl;
    return true;
  } else {
    cout << " INFLATIONGW - Hamiltonian Psi constraint did not converge! " << endl;
    cout << "              Residual " << residual << " after " << step
	 << " steps " << endl;
    return false;
  }
#else
  cout << " INFLATIONGW: Can't construct Inflation GW initial data without an Elliptic Solver!! " << endl;
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
  const double sfl = sf(i,j,k);
  const double Vl = potential->V(sfl);
  // const double Vl = V(i,j,k); 
	u[i][j][k] = 10. * PI * epsilon * Vl * psi4;
      }
};
void update_psi() {
  for (int i = 0; i < n_r; i++) 
    for (int j = 0; j < n_theta; j++) 
      for (int k = 0; k < n_phi; k++) {
	psi[i][j][k] += delta_psi(i,j,k);
      }
};
double Hamiltonian_Psi_Residual() {
  for (int i = N_g; i < n_r-N_g; i++) 
    for (int j = N_g; j < n_theta-N_g; j++) 
      for (int k = N_g; k < n_phi-N_g; k++) {
	const double psil = psi(i,j,k);
	const double psi5 = psil*psil*psil*psil*psil;
	const double psim7 = 1.0/pow(psil, 7);
	res[i][j][k] = psi.Laplace(i,j,k) 
	  + 2.0*PI*psi5*epsilon*V(i,j,k);
      }
  return res.L2_norm();
};
