// Tell emacs that this is -*-c++-*- mode

//===============================================================
// Compute fields...
//===============================================================
void Compute_Fields() {
	for (int i = N_g; i < n_r - N_g; i++) {
		const double rl = rho.r(i);
		for (int j = N_g; j < n_theta - N_g; j++) {
			const double thetal = rho.theta(j);
			for (int k = N_g; k < n_phi - N_g; k++) {
				sf[i][j][k] = compute_sf(rl);
			}
		}
	}
	sf.fill_ghosts();
	dump(&sf);
}

//===============================================================
// Compute fields: Inflaton
//===============================================================
double compute_sf(double r) {
	return phi_0 * exp(- r * r / (sigma * sigma));
}



//===============================================================
// compute rho_ADM and S^i for Inflation GW
// appear on right-hand sides of constraints
// S^i = 0 for 0 conjugate momentum
// Compute only once after solving for psi
//===============================================================
void Compute_Sources() {
	for (int i = N_g; i < n_r - N_g; i++) {
		const double rl = rho.r(i);
		const double r2 = rl * rl;
		for (int j = N_g; j < n_theta - N_g; j++) {
			for (int k = N_g; k < n_phi - N_g; k++) {
				const double phi = sf(i, j, k);
				const double phi2 = phi * phi;
				const double psil = psi(i, j, k);
				const double psi4 = psil * psil * psil * psil;
				const double psim4 = 1. / psi4;
				const double sigma4 = sigma * sigma * sigma * sigma;
				const double sigmam4 = 1. / sigma4;
				const double Vl = V(i, j, k); 
				//
				// compute energy density
				//
				rho[i][j][k] = (2. * r2 * sigmam4 * psim4 * phi2 + Vl);
				//
				// compute *rescaled* upstairs momentum densities 
				//
				s_r[i][j][k] = 0.0;
				s_t[i][j][k] = 0.0;
				s_p[i][j][k] = 0.0;
			}
		}
	}
	rho.fill_ghosts();
	s_r.fill_ghosts();
	s_t.fill_ghosts();
	s_p.fill_ghosts();
	dump(&rho);
	dump(&s_r);
	dump(&s_t);
	dump(&s_p);
};