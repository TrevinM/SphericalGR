#include "ConstraintSolver.h"
#include "Container.h"

//
// Constructor
//
ConstraintSolver::ConstraintSolver(Grid* grid_i, state* s_i, curvature* curve_i,
    Matter* matter_i) :
    grid(grid_i), s(s_i), curve(curve_i) {
    N_g = grid->N_ghosts();
    N_r = grid->N_r_tot();
    N_theta = grid->N_theta_tot();
    N_phi = grid->N_phi_tot();
    PI = acos(1.0);
    //
    // setup gridfunctions
    //
    int gf_counter = 2000;
    u.setup(grid, 1, "u", gf_counter++, +1, +1, +1);
    v.setup(grid, 1, "v", gf_counter++, +1, +1, +1);
    psi.setup(grid, 1, "psi", gf_counter++, +1, +1, +1);
    A2.setup(grid, 1, "A2", gf_counter++, +1, +1, +1);
    D2_psi.setup(grid, 1, "D2_psi", gf_counter++, +1, +1, +1);
    HamRes.setup(grid, 1, "Ham", gf_counter++, +1, +1, +1);
    delta_psi.setup(grid, 1, "delta_psi", gf_counter++, +1, +1, +1);
}

ConstraintSolver::~ConstraintSolver() {
    cout << " destructing CONSTRAINTSOLVER ... " << endl;
}

double ConstraintSolver::SolveHamiltonian(double tol, int itmax, bool precollapsed_lapse) {
    //
    // first, do some prep work... 
    // 
    for (int i = 0; i < N_r; i++)
        for (int j = 0; j < N_theta; j++)
            for (int k = 0; k < N_phi; k++) {
                psi[i][j][k] = exp(s->phi(i, j, k));
                //
                // CHECK: fix this for non-zero extrinsic curvature!!!
                // 
                if (s->a_rr(i, j, k) != 0.0) {
                    cout << " CONSTRAINTSOLVER: allow for non-zero K_ij! " << endl;
                    return 0.0;
                }
                A2[i][j][k] = 0.0;
            }
    //
    // Find initial residual
    //
    Compute_uv();
    double res = HamiltonianResidual();
    cout << " CONSTRAINTSOLVER: initial maximum residual = " << res << endl;
#ifdef NoEllSolver
    cout << " CONSTRAINTSOLVER: Cannot solve constraints without elliptic solvers..." << endl;
#else
    cout << " CONSTRAINTSOLVER: solving Hamiltonian... " << endl;
    //
    // Now iterate...
    //
    int it = 0;
    while ((it < itmax) && (res > tol)) {
        it++;
        //
        // setup elliptic solver
        //    
        ellsolver = new FlatEllSolver3D(grid);
        ellsolver->SetupSolver(1.0, u);
        ellsolver->SetRHS(1.0, HamRes);
        //
        // solve for correction
        //
        int num_it = 0;
        ellsolver->Solve(1000, num_it, 1.0e-4 * res);
        ellsolver->GetSolution(delta_psi);
        // 
        // add correction to current psi
        //
        const double factor = 0.8;
        for (int i = 0; i < N_r; i++)
            for (int j = 0; j < N_theta; j++)
                for (int k = 0; k < N_phi; k++)
                    psi[i][j][k] += factor * delta_psi(i, j, k);
        //
        // compute new u and v and re-evaluate residual
        //
        Compute_uv();
        res = HamiltonianResidual();
        cout << " CONSTRAINTSOLVER: residual = " << res << " after "
            << it << " iterations." << endl;
    }
#endif
    //
    // finally: compute phi...
    // 
    for (int i = 0; i < N_r; i++)
        for (int j = 0; j < N_theta; j++)
            for (int k = 0; k < N_phi; k++) {
                const double psi_l = psi(i, j, k);
                s->phi[i][j][k] = log(psi_l);
                if (precollapsed_lapse) s->lapse[i][j][k] = 1.0 / (psi_l * psi_l);
            }
    s->phi.fill_ghosts();
    s->lapse.fill_ghosts();
    return res;
}

double ConstraintSolver::HamiltonianResidual() {
    double res_max = 0.0;
    int i_max = 0;
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = grid->r(i);
        for (int j = N_g; j < N_theta - N_g; j++) {
            const double thetal = grid->theta(j);
            const double stl = grid->sintheta(j);
            const double st2 = stl * stl;
            const double ctl = grid->costheta(j);
            const double cottheta = ctl / stl;
            for (int k = N_g; k < N_phi - N_g; k++) {
                //
                // compute physical connection functions
                //
                const double Gam_r = s->lam_r(i, j, k)
                    - curve->gup_tt(i, j, k) * rl - curve->gup_pp(i, j, k) * rl * st2;
                const double Gam_t = s->lam_t(i, j, k) / rl
                    + 2.0 * curve->gup_rt(i, j, k) / rl
                    - curve->gup_pp(i, j, k) * stl * ctl;
                const double Gam_p = s->lam_p(i, j, k) / (rl * stl)
                    + 2.0 * curve->gup_rp(i, j, k) / rl
                    + 2.0 * curve->gup_tp(i, j, k) * cottheta;
                //
                // now compute covariant Laplace operator of psi
                //
                D2_psi[i][j][k] =
                    curve->gup_rr(i, j, k) * psi.ddr(i, j, k) +
                    curve->gup_tt(i, j, k) * psi.ddtheta(i, j, k) +
                    curve->gup_pp(i, j, k) * psi.ddphi(i, j, k) +
                    2.0 * (curve->gup_rt(i, j, k) * psi.drdtheta(i, j, k) +
                        curve->gup_rp(i, j, k) * psi.drdphi(i, j, k) +
                        curve->gup_tp(i, j, k) * psi.dthetadphi(i, j, k))
                    - Gam_r * psi.dr(i, j, k)
                    - Gam_t * psi.dtheta(i, j, k)
                    - Gam_p * psi.dphi(i, j, k);
                HamRes[i][j][k] = v(i, j, k) - D2_psi(i, j, k);
                const double absres = fabs(HamRes(i, j, k));
                if (absres > res_max) {
                    res_max = absres;
                    i_max = i;
                }
            }
        }
    }
    cout << " CONSTRAINTSOLVER: Found maximum residual at i = "
        << i_max << endl;
    return res_max;
};

int ConstraintSolver::Compute_uv() {
    for (int i = 0; i < N_r; i++)
        for (int j = 0; j < N_theta; j++)
            for (int k = 0; k < N_phi; k++) {
                const double psi_l = psi(i, j, k);
                const double psi4 = psi_l * psi_l * psi_l * psi_l;
                const double psi5 = psi4 * psi_l;
                // u:
                u[i][j][k] = -curve->trace_R(i, j, k) / 8.0
                    - 5.0 * psi4 * s->K(i, j, k) * s->K(i, j, k) / 12.0
                    + 5.0 * psi4 * A2(i, j, k) / 8.0
                    + 10.0 * PI * psi4 * Container::matter->adm_sources->rho_ADM(i, j, k);
                // v:
                v[i][j][k] = psi_l * curve->trace_R(i, j, k) / 8.0
                    + psi5 * s->K(i, j, k) * s->K(i, j, k) / 12.0
                    - psi5 * A2(i, j, k) / 8.0
                    - 2.0 * PI * psi5 * Container::matter->adm_sources->rho_ADM(i, j, k);
            }
    return 0;
};