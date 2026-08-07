// Tell emacs that this is -*-c++-*- mode
//====================================================
// Solve Momentum constraint
//====================================================
bool Solve_Momentum(double tol_tri = 1.e-10, double tol_res = 1.e-8,
    bool verbose = true) {
#ifdef NoEllSolver
    cout << " DUALEMWAVE: Can't construct EM initial data without an Elliptic Solver!! " << endl;
    return false;
#else
    double res_norm = Momentum_Residual();
    if (verbose) cout << " DUALMAXWELL - Momentum constraint: initial residual: "
        << res_norm << endl;
    int step = 0;
    int max_step = 50;
    while (res_norm > tol_res && step < max_step) {
        step++;
        veclaplacian->SetupSolver();
        veclaplacian->SetRHS(-1.0, res_r, res_t, res_p);
        int max_it = 500;
        int num_it;
        double tol = 1.0e-8;
        veclaplacian->Solve(max_it, num_it, tol_tri);
        veclaplacian->GetSolution(del_W_r, del_W_t, del_W_p);
        update_W();
        res_norm = Momentum_Residual();
        if (verbose) cout << " DUALEMWAVE - Momentum residual after " << step
            << " steps = " << res_norm << endl;
    }
    if (res_norm < tol_res) {
        cout << " DUALEMWAVE - Momentum residual converged to " << res_norm
            << " in " << step << " steps" << endl;
        return true;
    } else {
        cout << " DUALEMWAVE - Momentum constraint did not converge! " << endl;
        cout << "              Residual " << res_norm << " after " << step
            << " steps " << endl;
        return false;
    }
#endif  /* NoEllSolver */
}
//====================================================
// Solve Momentum constraint
//====================================================
double Momentum_Residual() {
    double res_r_norm2 = 0.0;
    double res_t_norm2 = 0.0;
    double res_p_norm2 = 0.0;
    for (int i = N_g; i < n_r - N_g; i++) {
        const double rl = A2.r(i);
        const double r2 = rl * rl;
        for (int j = N_g; j < n_theta - N_g; j++) {
            const double sintheta = A2.sintheta(j);
            const double costheta = A2.costheta(j);
            const double sin2theta = sintheta * sintheta;
            const double cottheta = costheta / sintheta;
            for (int k = N_g; k < n_phi - N_g; k++) {
                const double psil = psi(i, j, k);
                const double psi10 = pow(psil, 10);
                const double RHS_r = 8.0 * PI * psi10 * s_r(i, j, k);
                const double RHS_t = 8.0 * PI * psi10 * s_t(i, j, k);
                const double RHS_p = 8.0 * PI * psi10 * s_p(i, j, k);
                //
                res_r[i][j][k] = 4. / 3. * W_r.ddr(i, j, k)
                    + 8. / 3. * W_r.dr(i, j, k) / rl
                    + W_r.ddtheta(i, j, k) / r2 + cottheta * W_r.dtheta(i, j, k) / r2
                    - 8. / 3. * W_r(i, j, k) / r2 + 1. / 3. * W_t.drdtheta(i, j, k) / rl
                    - 7. / 3. * W_t.dtheta(i, j, k) / r2 - 7. / 3. * cottheta * W_t(i, j, k) / r2
                    + cottheta / 3. * W_t.dr(i, j, k) / rl
                    - RHS_r;
                res_r_norm2 += res_r(i, j, k) * res_r(i, j, k);
                res_t[i][j][k] = W_t.ddr(i, j, k) + 2. * W_t.dr(i, j, k) / rl
                    + 4. / 3. * W_t.ddtheta(i, j, k) / r2
                    + 4. * cottheta / 3. * W_t.dtheta(i, j, k) / r2
                    - 4. / 3. * W_t(i, j, k) / (r2 * sin2theta) + 8. / 3. * W_r.dtheta(i, j, k) / r2
                    + 1. / 3. * W_r.drdtheta(i, j, k) / rl
                    - RHS_t;
                res_t_norm2 += res_t(i, j, k) * res_t(i, j, k);
                res_p[i][j][k] = W_p.ddr(i, j, k) + 2. * W_p.dr(i, j, k) / rl
                    + W_p.ddtheta(i, j, k) / r2 + cottheta * W_p.dtheta(i, j, k) / r2
                    - W_p(i, j, k) / (r2 * sin2theta)
                    - RHS_p;
                res_p_norm2 += res_p(i, j, k) * res_p(i, j, k);
            }
        }
    }
    return sqrt(res_r_norm2 + res_t_norm2 + res_p_norm2);
}


//====================================================
// Update W's
//====================================================
void update_W() {
    for (int i = 0; i < n_r; i++)
        for (int j = 0; j < n_theta; j++)
            for (int k = 0; k < n_phi; k++) {
                W_r[i][j][k] += del_W_r(i, j, k);
                W_t[i][j][k] += del_W_t(i, j, k);
                W_p[i][j][k] += del_W_p(i, j, k);
            }
}




