// Tell emacs that this is -*-c++-*- mode
//====================================================
// Solve Hamiltonian constraint
//====================================================
void Solve_Psi(double tol_tri = 1.e-10, double tol_res = 1.e-8,
		       bool verbose = true) {	
	double psi0 = 1.0;
	double psi0_last = psi0;
	double psi_asym_error_init = Integrate(psi0) - 1.0;
	cout << " INFLATIONGW: Inital psi_asym_error = " << psi_asym_error_init << endl;
 	double psi_asym_error = psi_asym_error_init;
	int iter = 0;
	cout << " INFLATIONGW: Finding weak branch bounds for psi0" << endl;

	while (psi_asym_error * psi_asym_error_init > 0.0 && iter < max_it) {
		iter++;
		psi0_last = psi0;
		psi0 *= 1.01; //Change from hardcoded
		psi_asym_error = Integrate(psi0) - 1.0;
		cout << "INFLATIONGW: psi0 = " << psi0 << " gives asym_error = " << psi_asym_error << endl;
	}
	if (iter >= max_it) {
		cerr << " INFLATIONGW: Could not find bounds on psi0 for weak branch" << endl;
	}

	if (branch == 1) {
		cout << " INFLATIONGW: Finding strong branch" << endl;
		iter = 0;
		psi_asym_error_init = psi_asym_error;
		while (psi_asym_error * psi_asym_error_init > 0.0 && iter < max_it) {
			iter++;
			psi0_last = psi0;
			psi0 *= 1.01; //Change from hardcoded
			psi_asym_error = Integrate(psi0) - 1.0;
			cout << "INFLATIONGW: psi0 = " << psi0 << " gives asym_error = " << psi_asym_error << endl;
		}
		if (iter >= max_it) {
			cerr << " INFLATIONGW: Could not find bounds on psi0 for strong branch" << endl;
		}
	}

	double psi0_high = psi0;
	double psi0_low = psi0_last;
	double psi_asym_error_low = Integrate(psi0_low) - 1.0;
	double psi_asym_error_high = Integrate(psi0_high) - 1.0;
	cout << " INFLATIONGW: Narrowed psi0 to between " << psi0_low << " and " << psi0_high << endl;
	double psi0_mid = (psi0_high + psi0_low) / 2.0;
	while (abs(psi_asym_error) > tol_tri && psi0_high - psi0_low > 1.e-2*tol_tri)
	{
		psi0_mid = (psi0_high + psi0_low) / 2.0;
		double psi_asym_error_mid = Integrate(psi0_mid) - 1.0;
		if (psi_asym_error_mid * psi_asym_error_low < 0.0) {
			psi0_high = psi0_mid;
			psi_asym_error_high = psi_asym_error_mid;
		}
		else {
			psi0_low = psi0_mid;
			psi_asym_error_low = psi_asym_error_mid;
		}
		psi_asym_error = psi_asym_error_mid;
		cout << " INFLATIONGW: Narrowed psi0 now " << psi0_mid << endl;
		cout << " INFLATIONGW: psi_asym_error =  " << psi_asym_error << endl;
		//cout << " INFLATIONGWL Hamiltonian psi residual = " << Hamiltonian_Psi_Residual() << endl;
	}
	if (abs(psi_asym_error) > tol_tri) {
		cerr << " INFLATIONGW: ERROR: Root finding failed to converge with bound difference " << psi0_high - psi0_low << endl;
	}
	else {
		cout << " INFLATIONGW: Found psi0 = " << psi0_mid << endl;
	}
	Update_Psi();
}


void Update_Psi() {
	for (int i = N_g; i < n_r; i++) {
		const double rl = K.r(i);
		for (int j = N_g; j < n_theta - N_g; j++) {
			const double thetal = K.theta(j);
			for (int k = N_g; k < n_phi - N_g; k++) {
				psi[i][j][k] = psi_r[i];
			}
		}
	}
	psi.fill_ghosts();
}


double Integrate(double psi0) {
	pair<double, double> vars = {psi0, 0.0};
	double delta_r = (grid->delta_r(N_g - 1)) / 2.0;
	pair<double, double> k1 = Ham_RHS_0(vars, 0.0);
	pair<double, double> vars2 = {vars.first + k1.first * delta_r * 0.5, vars.second + k1.second * delta_r * 0.5};
	pair<double, double> k2 = Ham_RHS(vars2, delta_r * 0.5);
	pair<double, double> vars3 = {vars.first + k2.first * delta_r * 0.5, vars.second + k2.second * delta_r * 0.5};
	pair<double, double> k3 = Ham_RHS(vars3, delta_r * 0.5);
	pair<double, double> vars4 = {vars.first + k3.first * delta_r, vars.second + k3.second * delta_r};
	pair<double, double> k4 = Ham_RHS(vars4, delta_r);
	vars.first = vars.first + delta_r * (k1.first + 2. * k2.first + 2. * k3.first + k4.first) / 6.;
	vars.second = vars.second + delta_r * (k1.second + 2. * k2.second + 2. * k3.second + k4.second) / 6.;	

  	for (int i = N_g; i < n_r; i++) {
		const double rl = grid->r(i);
		delta_r = grid->delta_r(i);
		if (vars.first < 0.0) {
			//cerr << " INFLATIONGW: ERROR: found some negative psi... Bad..." << endl;
		} 
		psi_r[i] = vars.first;
		pair<double, double> k1 = Ham_RHS(vars, rl);
		pair<double, double> vars2 = {vars.first + k1.first * delta_r * 0.5, vars.second + k1.second * delta_r * 0.5};
		pair<double, double> k2 = Ham_RHS(vars2, rl + delta_r * 0.5);
		pair<double, double> vars3 = {vars.first + k2.first * delta_r * 0.5, vars.second + k2.second * delta_r * 0.5};
		pair<double, double> k3 = Ham_RHS(vars3, rl + delta_r * 0.5);
		pair<double, double> vars4 = {vars.first + k3.first * delta_r, vars.second + k3.second * delta_r};
		pair<double, double> k4 = Ham_RHS(vars4, rl + delta_r);
		vars.first = vars.first + delta_r * (k1.first + 2. * k2.first + 2. * k3.first + k4.first) / 6.;
		vars.second = vars.second + delta_r * (k1.second + 2. * k2.second + 2. * k3.second + k4.second) / 6.;
	}
	return vars.first;

}


pair<double, double> Ham_RHS(pair<double, double> vars, double r) {
	const double sf = compute_sf(r);
	const double V = potential->V(sf);
	const double psi = vars.first;
	const double dpsidr = vars.second;
	const double psi5 = psi * psi * psi * psi * psi;
	return pair<double, double> {dpsidr, -(2. / r) * dpsidr - 2. * PI * psi5 * V * epsilon};			
}

pair<double, double> Ham_RHS_0(pair<double, double> vars, double r) {
	const double sf = compute_sf(r);
	const double V = potential->V(sf);
	const double psi = vars.first;
	const double psi5 = psi * psi * psi * psi * psi;
	return pair<double, double> {0.0, -2.0 * PI * psi5 * epsilon * V / 3.0};			
}



double Hamiltonian_Psi_Residual() {
  for (int i = N_g; i < n_r-N_g; i++) 
    for (int j = N_g; j < n_theta-N_g; j++) 
      for (int k = N_g; k < n_phi-N_g; k++) {
	const double psil = psi(i,j,k);
	const double psi5 = psil*psil*psil*psil*psil;
	const double psim7 = 1.0/pow(psil, 7);
  const double sfl = sf(i,j,k);
  const double Vl = potential->V(sfl);
	res[i][j][k] = psi.Laplace(i,j,k) 
	  + 2.0*PI*psi5*epsilon*Vl;
      }
	dump(&res);
	dump(&psi);
  return res.L2_norm();
}


void Solve_K() {
	const double w = step_factor;
	for (int i = N_g; i < n_r; i++) {
		const double rl = K.r(i);
		for (int j = N_g; j < n_theta - N_g; j++) {
			const double thetal = K.theta(j);
			const double Knew = compute_K(rl, thetal);
			for (int k = N_g; k < n_phi - N_g; k++) {
				K[i][j][k] = w * Knew + (1.0 - w) * K[i][j][k];

			}
		}
	}
	//cout << " INFLATIONGW: K_error at r = " << K.r(10) << " and theta = " << K.theta(n_theta / 2) 
	//<< " is " << K(10, n_theta/2, 0) + 3.296513422807648 << endl;
	//dump_error_monitor(K(10, n_theta/2, 0) + 3.296513422807648);
	K.fill_ghosts();
}


double compute_K(double r, double theta) {
	int i = grid->i_ind(r);
	int j = grid->j_ind(theta);
	int k = N_g;
  const double sfl = sf(i,j,k);
	const double sf2 = sfl * sfl;
  const double Vl = potential->V(sfl);
  const double r2 = r*r;
  const double sigma4 = sigma * sigma * sigma * sigma;
  const double sigmam4 = 1.0 / sigma4;
  const double psim4 = pow(psi(i,j,k), -4);
  const double psim12 = pow(psi(i,j,k), -12);
 	const double A2l = A2(i,j,k);

	return -sqrt(24.0 * PI * (1.0 - epsilon) * Vl + 
  12.0 * PI * 4.0 * r2 * sigmam4 * psim4 * sf2 + 
  1.5 * psim12 * A2l);
}

double Hamiltonian_K_Residual() {
  for (int i = N_g; i < n_r-N_g; i++) 
    for (int j = N_g; j < n_theta-N_g; j++) 
      for (int k = N_g; k < n_phi-N_g; k++) {
		const double rl = psi.r(i);
		const double r2 = rl * rl;
		const double sigma4 = sigma * sigma * sigma * sigma;
		const double sigmam4 = 1.0 / sigma4;
		const double psil = psi(i,j,k);
		const double psim4 = pow(psil, -4);
		const double psim12 = psim4 * psim4 * psim4;
  		const double sfl = sf(i,j,k);
		const double sf2 = sfl * sfl;
  		const double Vl = potential->V(sfl);
		const double Kl = K(i,j,k);
		const double K2 = Kl * Kl;
		const double A2l = A2(i,j,k);
		res[i][j][k] = K2 - (24.0 * PI * (1.0 - epsilon) * Vl + 
  		12.0 * PI * 4.0 * r2 * sigmam4 * psim4 * sf2 + 
  		1.5 * psim12 * A2l);
      }
	  //cout << " INFLATIONGW: Hamiltonian K residual = " << res.L2_norm() << endl;
  	return res.L2_norm();
};