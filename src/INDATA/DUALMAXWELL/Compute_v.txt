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
//===============================================================
// Compute fields: A
//===============================================================
double compute_a_p(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psin = pow(psi(i, j, k), n_psi);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  return a_amp * pow(r, n_a) * sin(theta) * psin *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
}
double compute_a_r(double r, double theta) { // curl of V_A
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double ctt = cos(theta) / sin(theta);
  const double psil = psi(i, j, k);
  const double v_a_p = compute_v_a_p(r, theta);
  const double v_a_p_t = compute_v_a_p_t(r, theta); // \partial \theta (V_A)^phi
  return (v_a_p_t + ctt * v_a_p) / (r * psil * psil);
}
double compute_a_t(double r, double theta) { // curl of V_A
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psil = psi(i, j, k);
  const double v_a_p = compute_v_a_p(r, theta);
  const double v_a_p_r = compute_v_a_p_r(r, theta);  // \partial r (V_A)^phi
  return -(v_a_p_r + v_a_p / r) / (psil * psil);
}
//===============================================================
// Compute fields: *A
//===============================================================
double compute_as_p(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psin = pow(psi(i, j, k), n_psi);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  return as_amp * pow(r, n_as) * sin(theta) * psin *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
}
double compute_as_r(double r, double theta) { // curl of V_A*
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double ctt = cos(theta) / sin(theta);
  const double psil = psi(i, j, k);
  const double v_as_p = compute_v_as_p(r, theta);
  const double v_as_p_t = compute_v_as_p_t(r, theta);
  return (v_as_p_t + ctt * v_as_p) / (r * psil * psil);
}
double compute_as_t(double r, double theta) { // curl of V_A*
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psil = psi(i, j, k);
  const double v_as_p = compute_v_as_p(r, theta);
  const double v_as_p_r = compute_v_as_p_r(r, theta);
  return -(v_as_p_r + v_as_p / r) / (psil * psil);
}
//===============================================================
// Compute fields: V_A
//===============================================================
double compute_v_a_p(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psin = pow(psi(i, j, k), n_psi);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  return v_a_amp * pow(r, n_v_a) * sin(theta) * psin *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
}
double compute_v_a_p_r(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psil = psi(i, j, k);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  const double v_a_p = compute_v_a_p(r, theta);
  const double term1 = n_v_a * v_a_p / r;
  const double term2 = v_a_amp * sin(theta) * pow(r, n_v_a) *
    (-2. * rpr0 * exp(-rpr0 * rpr0) - 2. * rmr0 * exp(-rmr0 * rmr0));
  const double term3 = n_psi * v_a_p / psil * psi.dr(i, j, k);
  return term1 + term2 + term3;
}
double compute_v_a_p_t(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psin = pow(psi(i, j, k), n_psi);
  const double psil = psi(i, j, k);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  const double term1 = v_a_amp * pow(r, n_v_a) * cos(theta) * psin *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
  const double term2 = v_a_amp * pow(r, n_v_a) * sin(theta) *
    n_psi * psin / psil * psi.dtheta(i, j, k) *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
  return term1 + term2;
}
//===============================================================
// Compute fields: V_*A
//===============================================================
double compute_v_as_p(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psin = pow(psi(i, j, k), n_psi);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  return v_as_amp * pow(r, n_v_as) * sin(theta) * psin *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
}
double compute_v_as_p_r(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psil = psi(i, j, k);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  const double v_as_p = compute_v_as_p(r, theta);
  const double term1 = n_v_as * v_as_p / r;
  const double term2 = v_as_amp * sin(theta) * pow(r, n_v_as) *
    (-2. * rpr0 * exp(-rpr0 * rpr0) - 2. * rmr0 * exp(-rmr0 * rmr0));
  const double term3 = n_psi * v_as_p / psil * psi.dr(i, j, k);
  return term1 + term2 + term3;
}
double compute_v_as_p_t(double r, double theta) {
  int i = grid->i_ind(r);
  int j = grid->j_ind(theta);
  int k = N_g;
  const double psin = pow(psi(i, j, k), n_psi);
  const double psil = psi(i, j, k);
  const double rpr0 = r + r0;
  const double rmr0 = r - r0;
  const double term1 = v_as_amp * pow(r, n_v_as) * cos(theta) * psin *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
  const double term2 = v_as_amp * pow(r, n_v_as) * sin(theta) *
    n_psi * psin / psil * psi.dtheta(i, j, k) *
    (exp(-rpr0 * rpr0) + exp(-rmr0 * rmr0));
  return term1 + term2;
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
