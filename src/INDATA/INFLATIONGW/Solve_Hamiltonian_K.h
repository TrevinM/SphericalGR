// Tell emacs that this is -*-c++-*- mode
//====================================================
// Algebraically Solve Hamiltonian constraint for K
//====================================================

void Solve_K() {
	for (int i = N_g; i < n_r - N_g; i++) {
		const double rl = K.r(i);
		for (int j = N_g; j < n_theta - N_g; j++) {
			const double thetal = K.theta(j);
			for (int k = N_g; k < n_phi - N_g; k++) {
				K[i][j][k] = compute_K(rl, thetal);

			}
		}
	}
}


double compute_K(double r, double theta) {
	int i = grid->i_ind(r);
	int j = grid->j_ind(theta);
	int k = N_g;
	const double phi2 = pow(sf(i, j, k), 2);
  	const double Vl = V(i,j,k);
  	const double r2 = r*r;
  	const double sigma4 = sigma * sigma * sigma * sigma;
  	const double sigmam4 = 1.0 / sigma4;
  	const double psim4 = pow(psi(i,j,k), -4);
  	const double psim12 = pow(psi(i,j,k), -12);
 	const double A2l = A2(i,j,k);

	return sqrt(24.0 * PI * (1.0 - epsilon) * Vl + 
  12.0 * PI * 4.0 * r2 * sigmam4 * psim4 * phi2 + 
  1.5 * psim12 * A2l);
}