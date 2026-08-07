// Tell emacs that this is -*-c++-*- mode
//================================================
// Initial data for disk around black hole
//================================================
//
#include "SBN.h"

class Disk : public InData {
private:
    double M, J;     // mass and angular momentum of black hole
    int type;        // trumpet type
    double r_i, r_o; // inner and outer edge of disk
    double kappa, n; // polytropic parameters
    SBN* sbn;
    EOS* eos;
    int array_length;
    double Rmax;
    double C;
#ifndef NoEllSolver
    FlatEllSolver3D* ellsolver;
    VecLaplacian* veclaplacian;
#endif
    gf3d q;       // polytropic matter variable
    gf3d psi_D, lapse_D, beta_r, beta_t, beta_p; // rescaled
    gf3d eta_D, eta_BH;
    gf3d f, delta;
    gf3d aop6, aop6_BH, ln_aop6, ln_aop6_BH;   // alpha / psi^6
    gf3d res, res_r, res_t, res_p;
    gf3d RHS, RHS_r, RHS_t, RHS_p;
    gf3d temp1, temp2, temp3, temp4;
    gf3d psi_BH, lapse_BH, beta_BH, K_BH, dpsi_BHdr, dlapse_BHdr;
    gf3d A_rr, A_rt, A_rp, A_tt, A_tp, A_pp, A2;
    gf3d L_beta_rr, L_beta_rt, L_beta_rp, L_beta_tt, L_beta_tp, L_beta_pp, L_beta2;
    gf3d grad_ln_ap_r, grad_ln_ap_t, grad_ln_ap_p;
    gf3d grad_ln_ap_D_r, grad_ln_ap_D_t, grad_ln_ap_D_p;
    gf3d u_r, u_t, u_p;
    gf3d A_rr_BH, A_tt_BH, A2_BH, L_beta_BH2; // components of A^{ij} for BH
    int n_r, n_theta, n_phi, N_g;   // grid parameters
    double PI;
    double eps;
public:
    //================================================
    // Constructor
    //================================================
    Disk(char* indata_input, EOS* eos_i, Grid* grid_i, Cosmology* cosmology) :
        InData(grid_i, cosmology), eos(eos_i) {
        indata_type = disk;
        N_g = grid->N_ghosts();
        analytical = false;
        ifstream infile;
        infile.open(indata_input);
        if (!infile)
            cerr << " DISK: Can't open " << indata_input
            << " for input. This is bad. " << endl;
        else
            cout << " DISK: Reading initial data parameters from file "
            << indata_input << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> M;
        infile.get(buf, 100, '='); infile.get(c); infile >> J;
        infile.get(buf, 100, '='); infile.get(c); infile >> type;
        infile.get(buf, 100, '='); infile.get(c); infile >> r_i;
        infile.get(buf, 100, '='); infile.get(c); infile >> r_o;
        infile.get(buf, 100, '='); infile.get(c); infile >> kappa;
        infile.get(buf, 100, '='); infile.get(c); infile >> n;
        cout << " DISK: Will set up disk around trumpet black hole with M = "
            << M << " and J = " << J;
        if (type == 1) {
            cout << " in maximal trumpet slicing " << endl;
        } else if (type == 2) {
            cout << " in analytical trumpet slicing " << endl;
        } else {
            cout << " OOOOPS -- no such trumpet slicing " << endl;
            exit(0);
        }
        cout << " DISK: edges of disk at r_i = " << r_i
            << " and r_o = " << r_o << endl;
        cout << " DISK: Using kappa = " << kappa << " and n = " << n << endl;
        cout << "===================================================" << endl;
        PI = acos(-1.0);
    };
    //================================================
    // Destructor
    //================================================
    ~Disk() {};
    string Name() { return "disk initial data"; };
    //================================================
    // Methods for solvers
    //================================================
#include "DISK/HamiltonianSolver.h"
#include "DISK/MomentumSolver.h"
  //================================================
  // Initializer
  //================================================
    bool Initialize(gf3d& fct) {
        cout << " DISK: initializing disk initial data... " << endl;
        if (eos != NULL) {
            cout << " DISK: Double check: EOS name = " << eos->Name() << endl;
        }
        n_r = fct.dim1();
        n_theta = fct.dim2();
        n_phi = fct.dim3();
        //
        // set up grid functions
        // 
        int gf_counter = 2000;
        q.setup(grid, 1, "q", gf_counter++, +1, +1, +1);
        temp1.setup(grid, 1, "temp1", gf_counter++, +1, +1, +1);
        temp2.setup(grid, 1, "temp2", gf_counter++, +1, +1, +1);
        temp3.setup(grid, 1, "temp3", gf_counter++, +1, +1, +1);
        temp4.setup(grid, 1, "temp4", gf_counter++, +1, +1, +1);
        aop6.setup(grid, 1, "aop6", gf_counter++, +1, +1, +1);
        aop6_BH.setup(grid, 1, "aop6_BH", gf_counter++, +1, +1, +1);
        ln_aop6.setup(grid, 1, "ln_aop6", gf_counter++, +1, +1, +1);
        ln_aop6_BH.setup(grid, 1, "ln_aop6_BH", gf_counter++, +1, +1, +1);
        psi_D.setup(grid, 1, "psi_D", gf_counter++, +1, +1, +1);
        lapse_D.setup(grid, 1, "lapse_D", gf_counter++, +1, +1, +1);
        beta_r.setup(grid, 1, "beta_r", gf_counter++, -1, +1, +1);
        beta_t.setup(grid, 1, "beta_t", gf_counter++, +1, -1, -1);
        beta_p.setup(grid, 1, "beta_p", gf_counter++, -1, -1, +1);
        u_r.setup(grid, 1, "u_r", gf_counter++, -1, +1, +1);
        u_t.setup(grid, 1, "u_t", gf_counter++, +1, -1, -1);
        u_p.setup(grid, 1, "u_p", gf_counter++, -1, -1, +1);
        grad_ln_ap_r.setup(grid, 1, "grad_ln_ap_r", gf_counter++, -1, +1, +1);
        grad_ln_ap_t.setup(grid, 1, "grad_ln_ap_t", gf_counter++, +1, -1, -1);
        grad_ln_ap_p.setup(grid, 1, "grad_ln_ap_p", gf_counter++, -1, -1, +1);
        grad_ln_ap_D_r.setup(grid, 1, "grad_ln_ap_D_r", gf_counter++, -1, +1, +1);
        grad_ln_ap_D_t.setup(grid, 1, "grad_ln_ap_D_t", gf_counter++, +1, -1, -1);
        grad_ln_ap_D_p.setup(grid, 1, "grad_ln_ap_D_p", gf_counter++, -1, -1, +1);
        delta.setup(grid, 1, "delta", gf_counter++, +1, +1, +1);
        f.setup(grid, 1, "f", gf_counter++, +1, +1, +1);
        L_beta_rr.setup(grid, 1, "L_beta_rr", gf_counter++, +1, +1, +1);
        L_beta_rt.setup(grid, 1, "L_beta_rt", gf_counter++, -1, -1, -1);
        L_beta_rp.setup(grid, 1, "L_beta_rp", gf_counter++, +1, -1, +1);
        L_beta_tt.setup(grid, 1, "L_beta_tt", gf_counter++, +1, +1, +1);
        L_beta_tp.setup(grid, 1, "L_beta_tp", gf_counter++, -1, +1, -1);
        L_beta_pp.setup(grid, 1, "L_beta_pp", gf_counter++, +1, +1, +1);
        L_beta2.setup(grid, 1, "L_beta2", gf_counter++, +1, +1, +1);
        A_rr.setup(grid, 1, "A_rr", gf_counter++, +1, +1, +1);
        A_rt.setup(grid, 1, "A_rt", gf_counter++, -1, -1, -1);
        A_rp.setup(grid, 1, "A_rp", gf_counter++, +1, -1, +1);
        A_tt.setup(grid, 1, "A_tt", gf_counter++, +1, +1, +1);
        A_tp.setup(grid, 1, "A_tp", gf_counter++, -1, +1, -1);
        A_pp.setup(grid, 1, "A_pp", gf_counter++, +1, +1, +1);
        A2.setup(grid, 1, "A2", gf_counter++, +1, +1, +1);
        L_beta_BH2.setup(grid, 1, "L_beta_BH2", gf_counter++, +1, +1, +1);
        res.setup(grid, 1, "res", gf_counter++, +1, +1, +1);
        res_r.setup(grid, 1, "res_r", gf_counter++, -1, +1, +1);
        res_t.setup(grid, 1, "res_t", gf_counter++, +1, -1, -1);
        res_p.setup(grid, 1, "res_p", gf_counter++, -1, -1, +1);
        RHS.setup(grid, 1, "RHS", gf_counter++, +1, +1, +1);
        RHS_r.setup(grid, 1, "RHS_r", gf_counter++, -1, +1, +1);
        RHS_t.setup(grid, 1, "RHS_t", gf_counter++, +1, -1, -1);
        RHS_p.setup(grid, 1, "RHS_p", gf_counter++, -1, -1, +1);
        psi_BH.setup(grid, 1, "psi_BH", gf_counter++, +1, +1, +1);
        lapse_BH.setup(grid, 1, "lapse_BH", gf_counter++, +1, +1, +1);
        K_BH.setup(grid, 1, "K_BH", gf_counter++, +1, +1, +1);
        dpsi_BHdr.setup(grid, 1, "dpsi_BHdr", gf_counter++, -1, +1, +1); // r component!
        dlapse_BHdr.setup(grid, 1, "dlapse_BHdr", gf_counter++, -1, +1, +1); // r component!
        beta_BH.setup(grid, 1, "beta_BH", gf_counter++, -1, +1, +1); // r component!
        A_rr_BH.setup(grid, 1, "A_rr_BH", gf_counter++, +1, +1, +1);
        A_tt_BH.setup(grid, 1, "A_tt_BH", gf_counter++, +1, +1, +1);
        A2_BH.setup(grid, 1, "A2_BH", gf_counter++, +1, +1, +1);

#ifndef NoEllSolver 
        veclaplacian = new VecLaplacian(grid, true);
#endif
        //=============================================
        // set up black hole background and initial guess for density
        //=============================================
        BH_solution();
        Initialize_q(J);
        //=============================================
        // now solve constraints...
        //=============================================
        // TestShiftSolver();
        // exit(0);   // FIX THIS!!!


        Solve_Constraints();


        exit(0);  // FIX THIS!!!

        return true;
    }
    //================================================
    //================================================
    // functions for disk solution...
    //================================================
    //================================================
    //
    double Solve_Constraints(double tol_tri = 1.e-10, double tol_res = 1.e-8) {
        double residual = 0.0;
        int rounds = 15;
        int step = 0;
        while (step < rounds) {
            step++;
            Solve_Momentum();
            Solve_Hamiltonian();
            dump(&beta_r, step);
            dump(&beta_t, step);
            dump(&beta_p, step);
            dump(&res_r, step);
            dump(&res_t, step);
            dump(&res_p, step);
            dump(&RHS_r, step);
            dump(&RHS_t, step);
            dump(&RHS_p, step);
            dump(&psi_D, step);
            dump(&res, step);
        }
        return residual;
    }


    //================================================
    // Compute A_ij
    //================================================
    void Compute_Aij(bool verbose = 1) {
        if (verbose) cout << " DISK: computing A^ij" << endl;
        // first update aop6 including on outer boundaries 
        for (int i = N_g; i < n_r; i++)
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++) {
                    const double psil = psi_D(i, j, k) + psi_BH(i, j, k);
                    const double psi6 = psil * psil * psil * psil * psil * psil;
                    const double psi_BHl = psi_BH(i, j, k);
                    const double psi_BH6 = psi_BHl * psi_BHl * psi_BHl * psi_BHl * psi_BHl * psi_BHl;
                    aop6[i][j][k] = (lapse_D(i, j, k) + lapse_BH(i, j, k)) / psi6;
                    aop6_BH[i][j][k] = lapse_BH(i, j, k) / psi_BH6;
                    ln_aop6[i][j][k] = log(aop6(i, j, k));
                    ln_aop6_BH[i][j][k] = log(aop6_BH(i, j, k));
                }

        // now deal with rest
        for (int i = N_g; i < n_r - N_g; i++) {
            const double rl = psi_D.r(i);
            const double r2 = rl * rl;
            const double r3 = r2 * rl;
            for (int j = N_g; j < n_theta - N_g; j++) {
                const double sintheta = psi_D.sintheta(j);
                const double costheta = psi_D.costheta(j);
                const double sin2theta = sintheta * sintheta;
                const double cottheta = costheta / sintheta;
                for (int k = N_g; k < n_phi - N_g; k++) {
                    // compute divergence of shift (in terms of rescaled components)
                    const double div = beta_r.dr_so(i, j, k) + 2. * beta_r(i, j, k) / rl
                        + beta_t.dtheta_so(i, j, k) / rl + cottheta * beta_t(i, j, k) / rl;
                    // compute D^i \beta^j (rescaled components!)
                    const double Drbetar = beta_r.dr_so(i, j, k);
                    const double Dtbetar = (beta_r.dtheta_so(i, j, k) - beta_t(i, j, k)) / rl;
                    const double Dpbetar = -beta_p(i, j, k) / rl;
                    const double Drbetat = beta_t.dr_so(i, j, k);
                    const double Dtbetat = (beta_t.dtheta_so(i, j, k) + beta_r(i, j, k)) / rl;
                    const double Dpbetat = -cottheta * beta_p(i, j, k) / rl;
                    const double Drbetap = beta_p.dr_so(i, j, k);
                    const double Dtbetap = beta_p.dtheta_so(i, j, k) / rl;
                    const double Dpbetap = (beta_r(i, j, k) + cottheta * beta_t(i, j, k)) / rl;
                    // compute (L beta)^ij (still rescaled!)
                    L_beta_rr[i][j][k] = 2.0 * Drbetar - 2. / 3. * div;
                    L_beta_rt[i][j][k] = Drbetat + Dtbetar;
                    L_beta_rp[i][j][k] = Drbetap + Dpbetar;
                    L_beta_tt[i][j][k] = 2.0 * Dtbetat - 2. / 3. * div;
                    L_beta_tp[i][j][k] = Dtbetap + Dpbetat;
                    L_beta_pp[i][j][k] = 2.0 * Dpbetap - 2. / 3. * div;
                    if (i == N_g && j == N_g && k == N_g - 1)
                        cout << " reality check: "
                        << L_beta_rr(i, j, k) + L_beta_tt(i, j, k) + L_beta_pp(i, j, k) << endl;
                    // now compute A^ij (still rescaled)
                    const double factor1 = 1. / aop6(i, j, k);
                    const double factor2 = 1. / aop6_BH(i, j, k);
                    A_rr[i][j][k] = factor1 * L_beta_rr(i, j, k) / 2.0 + factor1 / factor2 * A_rr_BH(i, j, k);
                    A_rt[i][j][k] = factor1 * L_beta_rt(i, j, k) / 2.0;
                    A_rp[i][j][k] = factor1 * L_beta_rp(i, j, k) / 2.0;
                    A_tt[i][j][k] = factor1 * L_beta_tt(i, j, k) / 2.0 + factor1 / factor2 * A_tt_BH(i, j, k);
                    A_tp[i][j][k] = factor1 * L_beta_rr(i, j, k) / 2.0;
                    A_pp[i][j][k] = factor1 * L_beta_pp(i, j, k) / 2.0 + factor1 / factor2 * A_tt_BH(i, j, k);
                    if (i == N_g && j == N_g && k == N_g - 1)
                        cout << " reality check: "
                        << A_rr(i, j, k) + A_tt(i, j, k) + A_pp(i, j, k) << endl;
                    // finally compute A2 = A_{ij} A^{ij}
                    A2[i][j][k] = A_rr(i, j, k) * A_rr(i, j, k)
                        + 2.0 * A_rt(i, j, k) * A_rt(i, j, k)
                        + 2.0 * A_rp(i, j, k) * A_rp(i, j, k)
                        + A_tt(i, j, k) * A_tt(i, j, k)
                        + 2.0 * A_tp(i, j, k) * A_tp(i, j, k)
                        + A_pp(i, j, k) * A_pp(i, j, k);
                    L_beta2[i][j][k] = L_beta_rr(i, j, k) * L_beta_rr(i, j, k)
                        + 2.0 * L_beta_rt(i, j, k) * L_beta_rt(i, j, k)
                        + 2.0 * L_beta_rp(i, j, k) * L_beta_rp(i, j, k)
                        + L_beta_tt(i, j, k) * L_beta_tt(i, j, k)
                        + 2.0 * L_beta_tp(i, j, k) * L_beta_tp(i, j, k)
                        + L_beta_pp(i, j, k) * L_beta_pp(i, j, k);
                    //  compute components of gradient of alpha / psi^6
                    const double psil = psi_D(i, j, k) + psi_BH(i, j, k);
                    const double psi6 = psil * psil * psil * psil * psil * psil;
                    const double psi_BHl = psi_BH(i, j, k);
                    const double psi_BH6 = psi_BHl * psi_BHl * psi_BHl * psi_BHl * psi_BHl * psi_BHl;
                    const double psi_rat = psi_D(i, j, k) / psi_BH(i, j, k);
                    grad_ln_ap_D_r[i][j][k] = -6.0 * (psi_D.dr_so(i, j, k) - psi_rat * dpsi_BHdr(i, j, k)) / psil;
                    grad_ln_ap_D_t[i][j][k] = -6.0 * psi_D.dtheta_so(i, j, k) / psil;
                    grad_ln_ap_D_p[i][j][k] = 0.0;   // axisymmetry!
                    grad_ln_ap_r[i][j][k] = grad_ln_ap_D_r(i, j, k) + dlapse_BHdr(i, j, k) / lapse_BH(i, j, k) - 6.0 * dpsi_BHdr(i, j, k) / psi_BHl;
                    grad_ln_ap_t[i][j][k] = grad_ln_ap_D_t(i, j, k);
                    grad_ln_ap_p[i][j][k] = grad_ln_ap_D_p(i, j, k);
                }
            }
        }
        aop6.fill_ghosts();
        aop6_BH.fill_ghosts();
        ln_aop6.fill_ghosts();
        ln_aop6_BH.fill_ghosts();

        dump(&lapse_BH);
        dump(&psi_BH);
        dump(&grad_ln_ap_r);
        dump(&grad_ln_ap_D_r);
        dump(&temp2);
        dump(&A_rr_BH);
        dump(&A_tt_BH);
        dump(&aop6);
        dump(&aop6_BH);
    }

    void BH_solution() {
        if (type == 1) {
            // allocate SBN class
            array_length = 1000000;
            double r_max = grid->r_max();
            Rmax = 1.5 * r_max;
            sbn = new SBN(Rmax, array_length);
            // finally compute "magical" value of C
            C = 3.0 * sqrt(3.0) * M * M / 4.0;
            //
            for (int i = 0; i < n_r; i++) {
                double rl = psi_D.r(i);
                double psil = sbn->conFactor(abs(rl), M);
                double lapsel = sbn->lapse(abs(rl), M);
                for (int j = 0; j < n_theta; j++)
                    for (int k = 0; k < n_phi; k++) {
                        psi_BH[i][j][k] = psil;
                        lapse_BH[i][j][i] = lapsel;
                        K_BH[i][j][k] = 0.0;
                    }
            }
            delete sbn;
        } else if (type == 2) {
            for (int i = 0; i < n_r; i++) {
                // double rl = fabs(psi_D.r(i));
                double rl = psi_D.r(i);
                double psil = sqrt(1 + M / rl);
                double lapsel = rl / (rl + M);
                for (int j = 0; j < n_theta; j++)
                    for (int k = 0; k < n_phi; k++) {
                        const double Mpr = M + rl;
                        const double Mpr2 = Mpr * Mpr;
                        const double Mpr6 = Mpr2 * Mpr2 * Mpr2;
                        psi_BH[i][j][k] = psil;
                        lapse_BH[i][j][k] = lapsel;
                        A_rr_BH[i][j][k] = -4. * M * (M + rl) / (3. * rl * rl * rl);
                        A_tt_BH[i][j][k] = 2. * M * (M + rl) / (3. * rl * rl * rl);
                        K_BH[i][j][k] = M / ((rl + M) * (rl + M));
                        A2_BH[i][j][k] = 8. * M * M * (M + rl) * (M + rl) / (3. * rl * rl * rl * rl * rl * rl);
                        dpsi_BHdr[i][j][k] = -M / (2. * psil * rl * rl);
                        dlapse_BHdr[i][j][k] = -rl / ((rl + M) * (rl + M)) + 1. / (rl + M);
                        L_beta_BH2[i][j][k] = 32 * M * M * rl * rl / (3 * Mpr6);
                        eta_BH[i][j][k] = lapse_BH(i, j, k) * psi_BH(i, j, k);
                    }
            }
        }
        psi_BH.fill_ghosts();
        lapse_BH.fill_ghosts();
        A_rr_BH.fill_ghosts();
        A_tt_BH.fill_ghosts();
        K_BH.fill_ghosts();
        A2_BH.fill_ghosts();
        dpsi_BHdr.fill_ghosts();
        dump(&A2_BH);
        dump(&psi_BH);
    }
    void Initialize_q(double q_init) {
        const double x_center = (r_i + r_o) / 2.0;
        const double r_disk = (r_o - r_i) / 2.0;
        for (int i = 0; i < n_r; i++) {
            const double rl = psi_D.r(i);
            for (int j = 0; j < n_theta; j++) {
                const double sinthetal = psi_D.sintheta(j);
                const double costhetal = psi_D.costheta(j);
                for (int k = 0; k < n_phi; k++) {
                    const double xl = rl * sinthetal;
                    const double zl = rl * costhetal;
                    const double dist = sqrt((xl - x_center) * (xl - x_center) + zl * zl);
                    psi_D[i][j][k] = 0.0;
                    lapse_D[i][j][k] = 0.0;
                    if (dist < r_disk) {
                        q[i][j][k] = q_init * (1.0 - dist * dist / (r_disk * r_disk));
                    } else {
                        q[i][j][k] = 0.0;
                    }
                }
            }
        }
        dump(&q);
    }
    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_rt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_rp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_tt_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_tp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double h_pp_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double phi_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = r;
        return 0.0;
        // return log(sbn->conFactor(abs(r_C),M));
    }
    //================================================
    // Analytical solution for connection coefficients
    //================================================
    double lam_r_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = true; return 0;
    }
    double lam_t_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = true; return 0;
    }
    double lam_p_analytical(double r, double theta, double phi, double t,
        bool& done) {
        done = true; return 0;
    }
    //================================================
    // Analytical solution for extrinsic curvature
    //================================================
    double a_rr_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = r;
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        const double psi = 1.0;   // sbn->conFactor(abs(r_C),M);
        const double R_C = psi * psi * r_C;
        const double factor = -C / (R_C * R_C * R_C);
        // compute extrinsic curvature in cartesian coordinates
        tensor a_cart(factor * (3.0 * n_x * n_x - 1.0),
            factor * (3.0 * n_x * n_y),
            factor * (3.0 * n_x * n_z),
            factor * (3.0 * n_y * n_y - 1.0),
            factor * (3.0 * n_y * n_z),
            factor * (3.0 * n_z * n_z - 1.0));
        // transform to spherical polar coordinates
        tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
        // and return rr component
        return a[0][0];
    };
    double a_rt_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        const double psi = 1.0;  //sbn->conFactor(abs(r_C),M);
        const double R_C = psi * psi * r_C;
        const double factor = -C / (R_C * R_C * R_C);
        // compute extrinsic curvature in cartesian coordinates
        tensor a_cart(factor * (3.0 * n_x * n_x - 1.0),
            factor * (3.0 * n_x * n_y),
            factor * (3.0 * n_x * n_z),
            factor * (3.0 * n_y * n_y - 1.0),
            factor * (3.0 * n_y * n_z),
            factor * (3.0 * n_z * n_z - 1.0));
        // transform to spherical polar coordinates
        tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
        // and return rt component
        return a[0][1] / r;
    };
    double a_rp_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        const double psi = 1.0; // sbn->conFactor(abs(r_C),M);
        const double R_C = psi * psi * r_C;
        const double factor = -C / (R_C * R_C * R_C);
        // compute extrinsic curvature in cartesian coordinates
        tensor a_cart(factor * (3.0 * n_x * n_x - 1.0),
            factor * (3.0 * n_x * n_y),
            factor * (3.0 * n_x * n_z),
            factor * (3.0 * n_y * n_y - 1.0),
            factor * (3.0 * n_y * n_z),
            factor * (3.0 * n_z * n_z - 1.0));
        // transform to spherical polar coordinates
        tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
        // and return rp component
        return a[0][2] / (r * sin(theta));
    };
    double a_tt_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        const double psi = 1.0; // sbn->conFactor(abs(r_C),M);
        const double R_C = psi * psi * r_C;
        const double factor = -C / (R_C * R_C * R_C);
        // compute extrinsic curvature in cartesian coordinates
        tensor a_cart(factor * (3.0 * n_x * n_x - 1.0),
            factor * (3.0 * n_x * n_y),
            factor * (3.0 * n_x * n_z),
            factor * (3.0 * n_y * n_y - 1.0),
            factor * (3.0 * n_y * n_z),
            factor * (3.0 * n_z * n_z - 1.0));
        // transform to spherical polar coordinates
        tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
        // and return tt component
        return a[1][1] / (r * r);
    };
    double a_tp_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;

        const double psi = 1.0; // sbn->conFactor(abs(r_C),M);
        const double R_C = psi * psi * r_C;
        const double factor = -C / (R_C * R_C * R_C);
        // compute extrinsic curvature in cartesian coordinates
        tensor a_cart(factor * (3.0 * n_x * n_x - 1.0),
            factor * (3.0 * n_x * n_y),
            factor * (3.0 * n_x * n_z),
            factor * (3.0 * n_y * n_y - 1.0),
            factor * (3.0 * n_y * n_z),
            factor * (3.0 * n_z * n_z - 1.0));
        // transform to spherical polar coordinates
        tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
        // and return tp component
        return a[1][2] / (r * r * sin(theta));
    };
    double a_pp_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        const double psi = 1.0; //  sbn->conFactor(abs(r_C),M);
        const double R_C = psi * psi * r_C;
        const double factor = -C / (R_C * R_C * R_C);
        // compute extrinsic curvature in cartesian coordinates
        tensor a_cart(factor * (3.0 * n_x * n_x - 1.0),
            factor * (3.0 * n_x * n_y),
            factor * (3.0 * n_x * n_z),
            factor * (3.0 * n_y * n_y - 1.0),
            factor * (3.0 * n_y * n_z),
            factor * (3.0 * n_z * n_z - 1.0));
        // transform to spherical polar coordinates
        tensor a = Cartesian_to_Spherical(a_cart, r, theta, phi);
        // and return rr component
        return a[2][2] / (r * r * sin(theta) * sin(theta));
    };
    double K_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for gauge
    //================================================
    double lapse_analytical(double r, double theta, double phi, double t) {
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        return 1.0; // sbn->lapse(abs(r_C),M);
    };
    double shift_r_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        // compute shift vector in cartesian coordinates
        const double shift_norm = 0.0; //  sbn->shift(abs(r_C),M);
        vect beta_cart(shift_norm * n_x,
            shift_norm * n_y,
            shift_norm * n_z);
        // transform to spherical polar coordinates
        vect beta = Cartesian_to_Spherical_upper(beta_cart, r, theta, phi);
        // return r component
        return beta[0];
    };
    double shift_t_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        // compute shift vector in cartesian coordinates
        const double shift_norm = 1.0; // sbn->shift(abs(r_C),M);
        vect beta_cart(shift_norm * n_x,
            shift_norm * n_y,
            shift_norm * n_z);
        // transform to spherical polar coordinates
        vect beta = Cartesian_to_Spherical_upper(beta_cart, r, theta, phi);
        // return r component
        return beta[1] * r;
    };
    double shift_p_analytical(double r, double theta, double phi, double t) {
        // find cartesian coordinates
        const double x = r * sin(theta) * cos(phi);
        const double y = r * sin(theta) * sin(phi);
        const double z = r * cos(theta);
        const double r_C = sqrt((x) * (x)+(y) * (y)+(z) * (z));
        // set up normal vector
        const double n_x = (x) / r_C;
        const double n_y = (y) / r_C;
        const double n_z = (z) / r_C;
        // compute shift vector in cartesian coordinates
        const double shift_norm = 1.0; // sbn->shift(abs(r_C),M);
        vect beta_cart(shift_norm * n_x,
            shift_norm * n_y,
            shift_norm * n_z);
        // transform to spherical polar coordinates
        vect beta = Cartesian_to_Spherical_upper(beta_cart, r, theta, phi);
        // return r component
        return beta[2] * r * sin(theta);
    };
    //================================================
    // Analytical solution for auxiliary functions
    //================================================
    double Theta_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double B_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double B_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double B_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    //================================================
    // Analytical solution for hydro
    //================================================
    double rho_0_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double P_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double v_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double v_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double v_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for scalar field
    //================================================
    double sf_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double pi_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    //================================================
    // Analytical solution for Maxwell OR for Dual Maxwell:
    // For Maxwell need e_p and a_p only, for Dual Maxwell all a_i and as_i
    //================================================
    double e_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double a_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double as_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double as_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double as_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };  //================================================
    // Analytical solution for radiation
    //================================================
    double E_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double F_0_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    }
    double F_r_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double F_t_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };
    double F_p_analytical(double r, double theta, double phi, double t) {
        return 0.0;
    };


    //================================================
    // Dump function
    //================================================
    void dump(gf3d* fct, int step = 0) {
        ofstream outfile_slice;
        ostringstream filename_slice;
        // slice files 
        filename_slice << "output/" << fct->Name() << "_indata_slice_"
            << step << "_" << n_r - 2 * N_g
            << "_" << n_theta - 2 * N_g << ends;
        outfile_slice.open(filename_slice.str().c_str());
        for (int i = N_g; i < n_r; i++) {
            for (int j = N_g; j < n_theta - N_g; j++)
                outfile_slice << setprecision(16) << setw(24) << fct->r(i)
                << setprecision(16) << setw(24) << fct->theta(j)
                << setprecision(16) << setw(24) << (*fct)(i, j, N_g) << endl;
            outfile_slice << endl;
        }
        outfile_slice.close();
        // ray files 
        ofstream outfile_rays;
        ostringstream filename_rays;
        filename_rays << "output/" << fct->Name() << "_indata_rays_"
            << step << "_" << n_r - 2 * N_g
            << "_" << n_theta - 2 * N_g << ends;
        outfile_rays.open(filename_rays.str().c_str());
        for (int i = N_g; i < n_r; i++) {
            outfile_rays << setprecision(16) << setw(24) << fct->r(i)
                << setprecision(16) << setw(24) << (*fct)(i, 0.0, N_g)
                << setprecision(16) << setw(24) << (*fct)(i, PI / 2., N_g) << endl;
        }
        outfile_rays.close();
    }


    //================================================
    // Test shift solver
    //================================================
    double TestShiftSolver() {
        // For testing purposes, set up RHS_S:
        for (int i = N_g; i < n_r - N_g; i++) {
            const double rl = psi_D.r(i);
            for (int j = N_g; j < n_theta - N_g; j++) {
                const double sintheta = psi_D.sintheta(j);
                const double costheta = psi_D.costheta(j);
                const double sin2theta = sintheta * sintheta;
                const double cottheta = costheta / sintheta;
                for (int k = N_g; k < n_phi - N_g; k++) {
                    // RHS_r[i][j][k] = - 8.*M*(4.*M + rl) / (3*(M+rl)*(M+rl)*(M+rl)*(M+rl));
                    RHS_r[i][j][k] = 0.0; // - rl * exp(-rl*rl);
                    RHS_t[i][j][k] = rl * exp(-rl * rl) * sintheta * costheta;
                    RHS_p[i][j][k] = rl * exp(-rl * rl) * sintheta * costheta;

                    u_r[i][j][k] = 0.0;
                    u_t[i][j][k] = 0.0; //  rl / ((1.0 + rl)*(1.0 + rl)); // 0.0; // 0.0; 
                    u_p[i][j][k] = rl / ((1.0 + rl) * (1.0 + rl));
                }
            }
        }
#ifndef NoEllSolver
        cout << " DISK: calling Setu..." << endl;
        veclaplacian->Setu(u_r, u_t, u_p);
        cout << " DISK: setting up solver..." << endl;
        veclaplacian->SetupSolver();
        cout << " DISK: calling SetRHS... " << endl;
        veclaplacian->SetRHS(RHS_r, RHS_t, RHS_p);
        int max_it = 500;
        int num_it;
        double tol = 1.0e-8;
        cout << " DISK: solving..." << endl;
        veclaplacian->Solve(max_it, num_it, tol);
        cout << " DISK: extracting solution..." << endl;
        veclaplacian->GetSolution(beta_r, beta_t, beta_p);
#endif

        double res_norm = Momentum_Residual(u_r, u_t, u_p);
        cout << res_norm << endl;
        dump(&beta_r);
        dump(&beta_t);
        dump(&beta_p);
        dump(&res_r);
        dump(&res_t);
        dump(&res_p);
        dump(&RHS_r);
        dump(&RHS_t);
        dump(&RHS_p);
        dump(&u_r);
        dump(&u_t);
        dump(&u_p);
        exit(0);
        return res_norm;
    }
};

