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
				a_r[i][j][k] = compute_a_r(rl, thetal);
				a_t[i][j][k] = compute_a_t(rl, thetal);
				a_p[i][j][k] = compute_a_p(rl, thetal);
				as_r[i][j][k] = compute_as_r(rl, thetal);
				as_t[i][j][k] = compute_as_t(rl, thetal);
				as_p[i][j][k] = compute_as_p(rl, thetal);
			}
		}
	}
	a_r.fill_ghosts();
	a_t.fill_ghosts();
	a_p.fill_ghosts();
	as_r.fill_ghosts();
	as_t.fill_ghosts();
	as_p.fill_ghosts();
	dump(&a_r);
	dump(&a_t);
	dump(&a_p);
	dump(&as_r);
	dump(&as_t);
	dump(&as_p);
}


//Arbitrary function used in spherical harmonic solutions
double f_pm(double x) {
	const double rpr0 = x + r0;
	const double rmr0 = x - r0;
	return 0.5 * (exp(-rpr0*rpr0) + exp(-rmr0*rmr0))* x;
}

double f_pm_dx(double x) {
	const double rpr0 = x + r0;
	const double rmr0 = x - r0;
	return 0.5 * (exp(-rpr0*rpr0) + exp(-rmr0*rmr0)) * (1. - 2. * x*x);
}

double f_pm_dx2(double x) {
	const double rpr0 = x + r0;
	const double rmr0 = x - r0;
	return (exp(-rpr0*rpr0) + exp(-rmr0*rmr0)) * x * (-3. + 2. * x*x);
}

double f_pm_dx3(double x) {
	const double rpr0 = x + r0;
	const double rmr0 = x - r0;
	return (exp(-rpr0*rpr0) + exp(-rmr0*rmr0)) * (-3. + 12. * x*x - 4. * x*x*x*x);
}

double f_pm_dx4(double x) {
	const double rpr0 = x + r0;
	const double rmr0 = x - r0;
	return 2. * (exp(-rpr0*rpr0) + exp(-rmr0*rmr0)) * x * (15. - 20. * x*x + 4. * x*x*x*x);
}

//Spherical harmonic solutions
double l1(double r, double theta) {
	return sin(theta) * ((f_pm(r) / (r*r)) - f_pm_dx(r) / r);
}

double l2(double r, double theta) {
	return sin(theta) * cos(theta) 
		* (f_pm(r) / (r*r*r)
		- f_pm_dx(r) / (r*r)
		+ f_pm_dx2(r) / (3. * r));
}

double l3(double r, double theta) {
	return (5. * pow(cos(theta), 2) - 1.) * sin(theta)
		* (f_pm(r) / (r*r*r*r) 
		- f_pm_dx(r) / (r*r*r)
		+ 2.* f_pm_dx2(r) / (5. * r*r) 
		- f_pm_dx3(r) / (15. * r));
}

//===============================================================
// Compute fields: A
//===============================================================
double compute_a_p(double r, double theta) {
	int i = grid->i_ind(r);
	int j = grid->j_ind(theta);
	int k = N_g;
	const double psin = pow(psi(i, j, k), n_psi);

	return psin * (a1_amp * l1(r, theta) + a2_amp * l2(r, theta) + a3_amp * l3(r, theta));
}
double compute_a_r(double r, double theta) {
	return 0;
}
double compute_a_t(double r, double theta) {
	return 0;
}

//===============================================================
// Compute fields: *A
//===============================================================
double compute_as_p(double r, double theta) {
	int i = grid->i_ind(r);
	int j = grid->j_ind(theta);
	int k = N_g;
	const double psin = pow(psi(i, j, k), n_psi);

	return psin * (as1_amp * l1(r, theta) + as2_amp * l2(r, theta) + as3_amp * l3(r, theta));

}
double compute_as_r(double r, double theta) {
	return 0;
}
double compute_as_t(double r, double theta) {
	return 0;
}

//===============================================================
// compute rho_ADM and S^i for EM wave;
// appear on right-hand sides of constraints
//===============================================================
void Compute_Sources() {
	const double oo4p = 1.0 / (4.0 * PI);
	for (int i = N_g; i < n_r - N_g; i++) {
		const double rl = rho.r(i);
		const double r2 = rl * rl;
		for (int j = N_g; j < n_theta - N_g; j++) {
			const double stl = rho.sintheta(j);
			const double st2 = stl * stl;
			const double ctl = rho.costheta(j);
			for (int k = N_g; k < n_phi - N_g; k++) {
				const double psil = psi(i, j, k);
				const double psi2 = psil * psil;
				const double psi4 = psil * psil * psil * psil;
				const double psim4 = 1. / psi4;
				// const double psim8 = 1./(psi4*psi4);
				//
				// get B^i and E^i from curl of A and *A
				// 
				double b_r, b_t, b_p;   // scaled, indices upstairs
				curl(a_r, a_t, a_p, b_r, b_t, b_p, i, j, k);
				double e_r, e_t, e_p;   // scaled, indices upstairs
				curl(as_r, as_t, as_p, e_r, e_t, e_p, i, j, k);
				//
				// compute dot products (for initial data \bar gamma_ij = eta_ij)
				//
				const double b2 = psi4 * (b_r * b_r + b_t * b_t + b_p * b_p);
				const double e2 = psi4 * (e_r * e_r + e_t * e_t + e_p * e_p);
				//
				// compute energy density
				//
				rho[i][j][k] = oo4p * (e2 + b2) / 2.0;
				//
				// compute *rescaled* upstairs momentum densities
				//
				s_r[i][j][k] = psi2 * oo4p * (e_t * b_p - e_p * b_t);
				s_t[i][j][k] = psi2 * oo4p * (e_p * b_r - e_r * b_p);
				s_p[i][j][k] = psi2 * oo4p * (e_r * b_t - e_t * b_r);
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
//
