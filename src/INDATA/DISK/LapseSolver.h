//
//================================================
// functions for lapse solver
//================================================
//
bool Solve_Lapse(double tol_tri = 1.e-10, double tol_res = 1.e-8,
    bool verbose = false) {
    // first update A_ij (which also updates aop6 and aop6_BH)
    Compute_Aij();
#ifndef NoEllSolver
    ellsolver = new FlatEllSolver3D(grid);
    //                                                                         
    // Allocate and set up elliptic solver                                    
    //   
    int step = 0;
    compute_eta();
    double res_norm = Hamiltonian_Residual();
    int max_step = 100;
    // if (verbose)
    cout << " DISK - Hamiltonian constraint: initial residual: "
        << res_norm << endl;
    while (res_norm > tol_res && step < max_step) {
        step++;
        compute_f(); // compute f
        ellsolver->SetupSolver(1.0, f);
        ellsolver->SetRHS(-1.0, res);
        //                                                              
        // Solve and get solution                                       
        //                                                                       
        int num_it;
        int max_iter = 500;
        double tri_res = ellsolver->Solve(max_iter, num_it, tol_tri);
        if (verbose) cout << "   Solved Hamiltonian constraint to residual "
            << tri_res << " in " << num_it << " iterations. "
            << endl;
        ellsolver->GetSolution(delta); //compute delta
        compute_psi(); //update psi
        res_norm = Hamiltonian_Residual(); //update residual
        if (verbose) cout << "   Residual after " << step << " steps: "
            << res_norm << endl;
    }
    delete ellsolver;
    bool success = false;
    if (res_norm < tol_res) {
        success = true;
        cout << " DISK : Hamiltonian constraint converged to " << res_norm
            << " in " << step << " iterations " << endl;
    }
    return success;
#else
    cout << " DISK: cannot construct initial data without elliptic solver! "
        << endl;
    exit(1);
#endif /* NoEllSolver */    
}
double Hamiltonian_Residual(bool verbose = false) {
    int n_scale = 0;
    for (int i = N_g + 1; i < n_r - N_g; i++)
        for (int j = N_g; j < n_theta - N_g; j++)
            for (int k = N_g; k < n_phi - N_g; k++) {
                const double psil = psi_D(i, j, k) + psi_BH(i, j, k);
                const double psi5pn = pow(psil, 5 + n_scale); // 5pn: 5 + n
                const double psi5 = psil * psil * psil * psil * psil;
                const double psim7 = 1. / (psil * psil * psil * psil * psil * psil * psil);
                const double psiBHl = psi_BH(i, j, k);
                const double psiBHm7 = 1. / (psiBHl * psiBHl * psiBHl * psiBHl * psiBHl * psiBHl * psiBHl);
                const double psiBH5 = psiBHl * psiBHl * psiBHl * psiBHl * psiBHl;
                const double Kl = K_BH(i, j, k);
                const double rhol = q(i, j, k);  // FIX THIS!!!
                const double temp = 1. / 12. * (psi5 - psiBH5) * Kl * Kl
                    - 1. / 8. * (psim7 * A2(i, j, k) - psiBHm7 * A2_BH(i, j, k)) - 2 * PI * psi5pn * rhol;
                res[i][j][k] = psi_D.Laplace_so(i, j, k) - temp;
                if (i == N_g && j == N_g && k == N_g && verbose) {
                    cout << "    temp = " << temp << " res = " << res[i][j][k] << endl;
                    cout << "    A2_BH = " << A2_BH(i, j, k) << " A2 = " << A2[i][j][k] << endl;
                }
                if (!isfinite(res[i][j][k])) {
                    cout << "Oooops!  " << i << " " << j << " " << k << endl;
                    cout << "Laplace = " << psi_D.Laplace(i, j, k) << ", temp = " << temp << endl;
                    cout << "psi5pn = " << psi5pn << endl;
                    cout << "rho = " << rhol << " psim7 = " << psim7 << endl;
                    exit(0);
                }
            }
    return res.L2_norm();
}
void compute_f() {
    int n_scale = 0;
    for (int i = 0; i < n_r; i++)
        for (int j = 0; j < n_theta; j++)
            for (int k = 0; k < n_phi; k++) {
                const double psil = psi_D(i, j, k) + psi_BH(i, j, k);
                const double psi4 = psil * psil * psil * psil;
                const double psim8 = 1. / (psi4 * psi4);
                const double rho_bar_l = 0.0;
                const double Kl = K_BH(i, j, k);
                const double A2l = A2(i, j, k);
                const double rhol = q(i, j, k);  // FIX THIS!!!
                f[i][j][k] = -5. / 12. * psi4 * Kl * Kl - 7. / 8. * psim8 * A2l + 10. * PI * psi4 * rhol;
            }
};

void compute_eta() {
    for (int i = 0; i < n_r; i++)
        for (int j = 0; j < n_theta; j++)
            for (int k = 0; k < n_phi; k++) {
                const double psil = psi_D(i, j, k) + psi_BH(i, j, k);
                const double psi4 = psil * psil * psil * psil;
                const double psim8 = 1. / (psi4 * psi4);
                const double rho_bar_l = 0.0;
                const double Kl = K_BH(i, j, k);
                const double A2l = A2(i, j, k);
                const double rhol = q(i, j, k);  // FIX THIS!!!
                f[i][j][k] = -5. / 12. * psi4 * Kl * Kl - 7. / 8. * psim8 * A2l + 10. * PI * psi4 * rhol;
            }
};

void compute_lapse() {
    for (int i = 0; i < n_r; i++)
        for (int j = 0; j < n_theta; j++)
            for (int k = 0; k < n_phi; k++) {
                const double psil = psi_D(i, j, k) + psi_BH(i, j, k);
                lapse_D[i][j][k] = (eta_D(i, j, k) - psi_D(i, j, k) * lapse_BH(i, j, k)) / psil;
            }
};
