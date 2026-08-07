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


double l1(double r, double theta) {
    double rmr0 = r - r0;
    double rpr0 = r + r0;

    double Em = exp(-rmr0 * rmr0);
    double Ep = exp(-rpr0 * rpr0);

    return 2. * sin(theta) * (Em * rmr0 + Ep * rpr0);
}


double l2(double r, double theta) {
    double rmr0 = r - r0;
    double rpr0 = r + r0;

    double Em = exp(-rmr0 * rmr0);
    double Ep = exp(-rpr0 * rpr0);

    return 1. / 3. * cos(theta) * sin(theta) * (
        Em * (-4. + (4. / r) * rmr0 + 8. * rmr0 * rmr0)
        + Ep * (-4. + (4. / r) * rpr0 + 8. * rpr0 * rpr0)
        );
}

double l3(double r, double theta) {
    double rmr0 = r - r0;
    double rpr0 = r + r0;

    double Em = exp(-rmr0 * rmr0);
    double Ep = exp(-rpr0 * rpr0);

    // Final simplified expression
    return (5. * cos(theta) * cos(theta) - 1.) * sin(theta) * (
        Em * (-2. / (5. * r) + (2. / (5. * r * r) - 12. / 15.) * rmr0 + (4. / (5. * r) * rmr0 * rmr0 + (8. / 15.) * rmr0 * rmr0 * rmr0))
        + Ep * (-2. / (5. * r) + (2. / (5. * r * r) - 12. / 15.) * rpr0 + (4. / (5. * r) * rpr0 * rpr0 + (8. / 15.) * rpr0 * rpr0 * rpr0))
        );
}

//This gives a potential in the AE_r component, that obeys curlAE = -dtAB for l1 (minus convention)
double A_r_l1m_dt(double r, double theta) {
    return -8. * exp(-r * r) * r * r * cos(theta);
}

//===============================================================
// Compute fields: A (B)
//===============================================================
double compute_a_p(double r, double theta) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = N_g;
    const double psin = pow(psi(i, j, k), n_psi);

    if (dual_sign == -1) {
        return 0;
    } else {
        return psin * (a1_amp * l1(r, theta) + a2_amp * l2(r, theta) + a3_amp * l3(r, theta));
    }
}
double compute_a_r(double r, double theta) {
    if (dual_sign == -1) {
        return -as1_amp * A_r_l1m_dt(r, theta);
    } else {
        return 0;
    }
}
double compute_a_t(double r, double theta) {
    return 0;
}

//===============================================================
// Compute fields: *A (E)
//===============================================================
double compute_as_p(double r, double theta) {
    int i = grid->i_ind(r);
    int j = grid->j_ind(theta);
    int k = N_g;
    const double psin = pow(psi(i, j, k), n_psi);

    return psin * (as1_amp * l1(r, theta) + as2_amp * l2(r, theta) + as3_amp * l3(r, theta));

}
double compute_as_r(double r, double theta) {
    if (dual_sign == -1) {
        return a1_amp * A_r_l1m_dt(r, theta);
    } else {
        return 0;
    }
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
