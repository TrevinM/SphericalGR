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
				sf[i][j][k] = compute_sf(rl, thetal);
				V[i][j][k] = compute_V(rl, thetal);
			}
		}
	}
	sf.fill_ghosts();
	v.fill_ghosts();
	dump(&sf);
	dump(&V);
}


//Arbitrary function used in wave equation solution
double f_pm(double x) {
	return 1/2 * (exp(-pow(x+r0, 2)) + exp(-pow(x-r0, 2))) * x;
}

double f_pm_dx(double x) {
	return 1/2 * (exp(-pow(x+r0, 2)) + exp(-pow(x-r0, 2)))
		+ 1/2 * (-2 * (x+r0) * exp(-pow(x+r0, 2)) - 2 * (x-r0) * exp(-pow(x-r0, 2))) * x;
}

double f_pm_dx2(double x) {
	return -2 * (x+r0) * exp(-pow(x+r0, 2)) - 2 * (x-r0) * exp(-pow(x-r0, 2))
		+ x/2 * (-2 * exp(-pow(x+r0, 2)) - 2 * exp(-pow(x-r0, 2))
		+ 4 * pow(x+r0, 2) * exp(-pow(x+r0, 2)) + 4 * pow(x-r0, 2) * exp(-pow(x+r0, 2)));
}

double f_pm_dx3(double x) {
	return 	3/2 * (-2 * exp(-pow(x+r0, 2)) - 2 * exp(-pow(x-r0, 2))
		+ 4 * pow(x+r0, 2) * exp(-pow(x+r0, 2)) + 4 * pow(x-r0, 2) * exp(-pow(x+r0, 2)))
		+ x/2 * (12 * exp(-pow(x+r0, 2)) * (x+r0) + 12 * exp(-pow(x-r0, 2)) * (x-r0)
		- 8 * exp(-pow(x+r0, 2)) * pow(x+r0, 3) - 8 * exp(-pow(x-r0, 2)) * pow(x-r0, 3));
}

//Spherical harmonic solutions
double l1(double r, double theta) {
	return sin(theta) * ((f_pm(r) / pow(r, 2)) + f_pm_dx(r) / r);
}

double l2(double r, double theta) {
	return sin(theta) * cos(theta) 
		* (f_pm(r) / pow(r, 3) - f_pm_dx(r) / pow(r, 2) + f_pm_dx2(r) / (3 * r));
}

double l3(double r, double theta) {
	return (
		5 * pow(cos(theta), 2) * sin(theta)
		* (f_pm(r) / pow(r, 4) - f_pm_dx(r) / pow(r, 3)
		+ 2 * f_pm_dx2(r) / (5 * pow(r, 2)) - f_pm_dx3(r) / (15 * r))
	);
}

//===============================================================
// Compute fields: Inflaton
//===============================================================
double compute_sf(double r, double theta) {
	return phi_0 * exp(-pow(r, 2) * pow(sigma, -2));
}

double compute_V(double r, double theta) {
	int i = grid->i_ind(r);
	int j = grid->j_ind(theta);
	int k = N_g;
	const double phi2 = pow(sf(i, j, k), 2);

	return 0.5 * m * m * phi2;
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

//===============================================================
// compute psi for given epsilon 
// part of the hamiltonian constraint
// equation 9a. in Brady Baumgarte & Clough
//===============================================================

void Compute_Psi() {
	for (int i = N_g; i < n_r - N_g; i++) {
		for (int j = N_g; j < n_theta - N_g; j++) {
			for (int k = N_g; k < n_phi - N_g; k++) {
				const double psil = psi(i, j, k);
				const double psi5 = psil * psil * psil * psil * psil;
				const double Vl = V(i, j, k); 
				const double LHS = psi.Laplace();
				const double RHS = -2.0 * PI * psi5 * epsilon * Vl;
				cout << "Solving for Psi with tol = " << tol << " and max_it = " << max_it << endl;
				while (res > tol && step < max_it) {

				}
				//
				// compute psi
				//
				psi[i][j][k] = (2. * r2 * sigmam4 * psim4 * phi2 + potential);
			}
		}
	}
	psi.fill_ghosts();
	dump(&psi);
};

