//================================================
//
// File that contains routines for horizon finder
//
//================================================
// #include "Einstein.h"
#include "HorizonFinder.h"
#include "tensors.h"
#include <ctime>


//================================================
// Constructor
//================================================
HorizonFinder::HorizonFinder(Grid* grid_i, state* s, curvature* c, auxiliary* aux,
    const char* file_stem) :
    gup_rr_3d(&(c->gup_rr)), gup_rt_3d(&(c->gup_rt)), gup_rp_3d(&(c->gup_rp)),
    gup_tt_3d(&(c->gup_tt)), gup_tp_3d(&(c->gup_tp)), gup_pp_3d(&(c->gup_pp)),
    DG_r_rr_3d(&(c->DG_r_rr)), DG_r_rt_3d(&(c->DG_r_rt)), DG_r_rp_3d(&(c->DG_r_rp)),
    DG_r_tt_3d(&(c->DG_r_tt)), DG_r_tp_3d(&(c->DG_r_tp)), DG_r_pp_3d(&(c->DG_r_pp)),
    DG_t_rr_3d(&(c->DG_t_rr)), DG_t_rt_3d(&(c->DG_t_rt)), DG_t_rp_3d(&(c->DG_t_rp)),
    DG_t_tt_3d(&(c->DG_t_tt)), DG_t_tp_3d(&(c->DG_t_tp)), DG_t_pp_3d(&(c->DG_t_pp)),
    DG_p_rr_3d(&(c->DG_p_rr)), DG_p_rt_3d(&(c->DG_p_rt)), DG_p_rp_3d(&(c->DG_p_rp)),
    DG_p_tt_3d(&(c->DG_p_tt)), DG_p_tp_3d(&(c->DG_p_tp)), DG_p_pp_3d(&(c->DG_p_pp)),
    phi_3d(&(s->phi)), phi_r_3d(&(aux->dphi_dr)), phi_t_3d(&(aux->dphi_dt)), phi_p_3d(&(aux->dphi_dp)),
    a_rr_3d(&(s->a_rr)), a_rt_3d(&(s->a_rt)), a_rp_3d(&(s->a_rp)),
    a_tt_3d(&(s->a_tt)), a_tp_3d(&(s->a_tp)), a_pp_3d(&(s->a_pp)),
    K_3d(&(s->K)),
    foundhorizon(false), last_step(-1), a_friedmann(0.0), grid(grid_i) {
    cout << " HORIZONFINDER: Constructing HorizonFinder... " << endl;
    PI = acos(-1.0);
    N_g = grid->N_ghosts();
    N_theta = grid->N_theta_tot();
    N_phi = grid->N_phi_tot();
    //
    // Set up r_min and r_max
    //
    SetRMaxMin();
    //
    // allocate surface functions
    //
    // set up level surface function
    //
    h.setup(grid);
    h_old.setup(grid);
    //
    // set up 2D metric functions
    //
    gup_rr.setup(grid);
    gup_rt.setup(grid);
    gup_rp.setup(grid);
    gup_tt.setup(grid);
    gup_tp.setup(grid);
    gup_pp.setup(grid);
    m_rr.setup(grid);
    m_rt.setup(grid);
    m_rp.setup(grid);
    m_tt.setup(grid);
    m_tp.setup(grid);
    m_pp.setup(grid);
    //
    // set up 2D connection coefficients
    //
    Gamma_r_rr.setup(grid);
    Gamma_r_rt.setup(grid);
    Gamma_r_rp.setup(grid);
    Gamma_r_tt.setup(grid);
    Gamma_r_tp.setup(grid);
    Gamma_r_pp.setup(grid);

    Gamma_t_rr.setup(grid);
    Gamma_t_rt.setup(grid);
    Gamma_t_rp.setup(grid);
    Gamma_t_tt.setup(grid);
    Gamma_t_tp.setup(grid);
    Gamma_t_pp.setup(grid);

    Gamma_p_rr.setup(grid);
    Gamma_p_rt.setup(grid);
    Gamma_p_rp.setup(grid);
    Gamma_p_tt.setup(grid);
    Gamma_p_tp.setup(grid);
    Gamma_p_pp.setup(grid);
    //
    // set up 2D extrinsic curvature
    //
    A_rr.setup(grid);
    A_rt.setup(grid);
    A_rp.setup(grid);
    A_tt.setup(grid);
    A_tp.setup(grid);
    A_pp.setup(grid);
    K.setup(grid);
    //
    // set up conformal exponent and derivatives
    //
    phi_c.setup(grid);
    phi_r.setup(grid);
    phi_t.setup(grid);
    phi_p.setup(grid);
    //
    // horizon functions
    //
    lambda.setup(grid);
    expansion.setup(grid);
    res.setup(grid);
    rhs_grid.setup(grid);
    factor.setup(grid);
    integrand.setup(grid);
    //
    // tendicity and vorticity stuff
    //
    E_rr.setup(grid);
    E_rt.setup(grid);
    E_rp.setup(grid);
    E_tt.setup(grid);
    E_tp.setup(grid);
    E_pp.setup(grid);
    B_rr.setup(grid);
    B_rt.setup(grid);
    B_rp.setup(grid);
    B_tt.setup(grid);
    B_tp.setup(grid);
    B_pp.setup(grid);
    tendicity.setup(grid);
    vorticity.setup(grid);
    //
    // read parameters from input file
    //
    int read_error = ReadInput();
    //
    // if error occured, replace input parameters with default values
    //
    if (read_error != 0) {
        eta = -1.0;
        tol_ell = 1.e-9;
        tol_exp = 1.e-4;
        max_iter_ell = 10000;
        max_iter_exp = 10;
        find_steps = 0;
        find_times = 1.0;
        mass_guess = 1.0;
        //      strncpy(file_stem,"Horizon",sizeof(file_stem) - 1);
        //      file_stem[sizeof(file_stem)-1] = 0;
        //      string name = "Horizon";
        //      file_stem = name.c_str();
    }
    cout << " HORIZONFINDER: using parameters: " << endl;
    cout << "   find_steps    = " << find_steps << endl;
    cout << "   find_times    = " << find_times << endl;
    cout << "   eta          = " << eta << endl;
    cout << "   tol_ell      = " << tol_ell << endl;
    cout << "   max_iter_ell = " << max_iter_ell << endl;
    cout << "   tol_exp      = " << tol_exp << endl;
    cout << "   max_iter_exp = " << max_iter_exp << endl;
    cout << "   mass_guess   = " << mass_guess << endl;
    cout << "   using file names starting with '" << file_stem << "'" << endl;
    use_time_step_criterion = (find_steps > 0);
    next_time = 0.0;
    next_step = 0;
    //
    // set up 2D elliptic solver
    //
#ifndef NoEllSolver    
    ellsolver = new EllSolver2D(grid);
    ellsolver->SetupSolver(eta);
#endif
    //
    // finally create a monitor file...
    //
    const int N_r = K_3d->dim1();
    //    const int c = K_3d->log_factor();
    ostringstream monfilename;
    monfilename << "output/" << file_stem << "_" << N_r - 2 * N_g << "_"
        << N_theta - 2 * N_g << ".hor_mon" << ends;
    monitorfile.open(monfilename.str().c_str());
    monitorfile.setf(ios::left);
    time_t clocktime;
    struct tm* currenttime;
    time(&clocktime);
    currenttime = localtime(&clocktime);
    monitorfile << "# File created on " << asctime(currenttime);
    monitorfile << "# " << setw(14) << "time" <<
        setw(24) << "Irr mass" <<
        setw(16) << "Spin J" <<
        setw(16) << "Kerr mass M" <<
        setw(16) << "J / M^2" <<
        setw(18) << "Lin Mom P_z" <<
        setw(16) << "Eq. Circumf." <<
        setw(16) << "Pol. Circumf." <<
        setw(16) << "Coord. Eq." <<
        setw(16) << "Coord. Pole" <<
        setw(16) << "Tend. M^2 Eq." <<
        setw(16) << "Tend. M^2 Pole" <<
        setw(16) << "Accretion rate" << endl;
    monitorfile << "#===============================================================================================================================================" << endl;
    //
    // ... and a surface file
    //
    // NOTE: will assume axisymmetry in current implementation
    //
    ostringstream surfacefilename;
    surfacefilename << "output/" << file_stem << "_" << N_r - 2 * N_g << "_"
        << N_theta - 2 * N_g << ".hor_surface" << ends;
    surfacefile.open(surfacefilename.str().c_str());
    surfacefile.setf(ios::right);
    surfacefile << "# File created on " << asctime(currenttime);
    surfacefile << "# " << setw(14) << "time" <<
        setw(16) << "theta" <<
        setw(16) << "h(theta)" <<
        setw(16) << "expansion" <<
        setw(16) << "tendicity*M^2" <<
        setw(16) << "vorticity*M^2" << endl;
    surfacefile << "#===============================================================================================" << endl;
};


//================================================
//
// Read input from file "AH_finder_Input"
//
//================================================
int HorizonFinder::ReadInput() {
    int error = 0;
    ifstream infile;
    infile.open("AH_finder_Input");
    if (!infile) {
        cerr << " HORIZONFINDER: Can't open input file AH_finder_Input " << endl;
        cerr << " HORIZONFINDER: Will use default values " << endl;
        error = 1;
        return error;
    }
    char buf[500], c;
    infile.get(buf, 500, '='); infile.get(c); infile >> find_steps;
    infile.get(buf, 500, '='); infile.get(c); infile >> find_times;
    infile.get(buf, 500, '='); infile.get(c); infile >> eta;
    infile.get(buf, 500, '='); infile.get(c); infile >> tol_ell;
    infile.get(buf, 500, '='); infile.get(c); infile >> max_iter_ell;
    infile.get(buf, 500, '='); infile.get(c); infile >> tol_exp;
    infile.get(buf, 500, '='); infile.get(c); infile >> max_iter_exp;
    infile.get(buf, 500, '='); infile.get(c); infile >> mass_guess;
    //  infile.get(buf,500,'='); infile.get(c); infile.get(c); 
    //  infile.get(file_stem,64);
    if (infile.eof()) {
        cerr << " HORIZONFINDER: Error reading input file AH_finder_Input " << endl;
        cerr << " HORIZONFINDER: Will use default values " << endl;
        error = 2;
    }
    return error;
}
//================================================
//
// Time to find horizon?
//
//================================================
bool HorizonFinder::TimeToFindHorizon(int currenttimestep, Doub currenttime) {
    timestep = currenttimestep;
    phystime = currenttime;
    bool search_for_horizon = false;
    if (find_times < 0.0) {
        return false;
    }
    if (use_time_step_criterion) {
        // evaluate time step criterion
        if (timestep >= next_step)
            search_for_horizon = true;
    } else {
        // otherwise evaluate time criterion
        if (phystime >= next_time)
            search_for_horizon = true;
    };
    // if (timestep == last_step)
    //   search_for_horizon = true;
    return search_for_horizon;
}
//================================================
//
// Find horizon
// (Note: state s and fluxes needed for accretion rate only...)
//
//================================================
bool HorizonFinder::FindHorizon(int currenttimestep, Doub currenttime,
    state* s, Fluxes* fluxes,
    Doub adm_mass, Doub mom_guess, Doub a) {
    a_friedmann = a;
    timestep = currenttimestep;
    phystime = currenttime;
    if (timestep > last_step) {
        //    next_step += find_steps;
        //    next_time += find_times;
        next_step = currenttimestep + find_steps;
        next_time = currenttime + find_times;
    }
    last_step = timestep;       // remember this in case we call this routine again at same timestep...
    //  cout << " HORIZONFINDER: Looking for horizon at time t = " << phystime << ", Friedmann expansion factor a = " << a_friedmann << endl;
    //
    // if we haven't found a horizon yet, initialize h
    //
    // CHECK...
    if (!foundhorizon) {
        cout << " HORIZONFINDER: Initial guess M = " << adm_mass << ", P = " << mom_guess << endl;
        InitialGuess(adm_mass, mom_guess);
    }
    for (int j = 0; j < N_theta; j++)
        for (int k = 0; k < N_phi; k++)
            h_old[j][k] = h(j, k);
    //
    // compute initial expansion
    //   
    Project(h);
    double exp = Expansion();
    int its = 0;
    //  cout << " HORIZONFINDER: Initial L2 Norm of expansion = " << exp << endl;
    const double constant = -1.0;
    const double xi = 0.5;
#ifndef NoEllSolver
    while ((exp > tol_exp) && (its < max_iter_exp)) {
        its++;
        //
        // compute rhs side for elliptic solver
        //
        for (int j = N_g; j < N_theta - N_g; j++)
            for (int k = N_g; k < N_phi - N_g; k++) {
                rhs_grid[j][k] = constant * expansion(j, k) / factor(j, k) + h.Laplace(j, k) + eta * h(j, k);
            }
        ellsolver->SetRHS(rhs_grid);
        //
        // solve
        //
        int num_it = 0;
        ellsolver->Solve(max_iter_ell, num_it, tol_ell);
        //
        // get new h
        //    
        ellsolver->GetSolution(h);
        for (int j = 0; j < N_theta; j++)
            for (int k = 0; k < N_phi; k++) {
                h[j][k] = xi * h(j, k) + (1.0 - xi) * h_old(j, k);
                h_old[j][k] = h(j, k);
            }
        //
        // project grid functions on new surface h...
        //
        Project(h);
        //
        // ... and re-compute expansion
        //
        exp = Expansion();
        //    cout << " HORIZONFINDER: expansion = " << exp << " after " << its << " steps " 
        //	 << " (used " << num_it << " trilinos iterations)" << endl;
    }
#else
    cout << " HORIZONFINDER: Cannot search for horizons without elliptic solver. " << endl;
#endif
    // CHECK!!! kludge to avoid tiny black holes sometimes found early in evolution  
    if (exp < tol_exp && h(PI / 2.0, 2) > r_min && h(PI / 2.0, 2) < r_max) {
        // if ( exp < tol_exp ) {
        //
        // found horizon!
        // 
        foundhorizon = true;
        const double Mass = HorizonMass();
        const double AngMom = Spin();
        const double LinMom = LinearMomentum();
        double accretion = 0;
        if (fluxes != NULL)
            accretion = AccretionRate(s->lapse, fluxes);
        const double EqCir = EquatorialCircumference();
        const double PolCir = PolarCircumference();
        cout << " HORIZONFINDER: Found horizon -- mass = " << Mass
            << ", spin = " << AngMom
            << ", linear momentum " << LinMom << endl;
        cout << " HORIZONFINDER: Equatorial and polar circumferences = " << EqCir
            << " and  " << PolCir << endl;
        DumpHorizon();
        Note(phystime, Mass, AngMom, LinMom, EqCir, PolCir,
            h(PI / 2.0, 2), h(0.0, 2), tendicity(PI / 2.0, 2), tendicity(0.0, 2), accretion);
        return true;
    } else {
        //
        // did not find horizon...
        //
        return false;
    }
};
//================================================
//
// Compute Expansion
//
// exp = \bar \lambda \bar m^{ij} h_,{ij} - \bar m^{ij} \bar s_k \bar \Gamma^k_{ij} + 4 \bar s^k \partial_k \phi - (2/3) psi^2 K
//           + psi^2 \bar A_{ij} \bar s^i \bar s^k
//
//================================================
double HorizonFinder::Expansion() {

    int j_test = -4;
    int k_test = 11;
    for (int j = N_g; j < N_theta - N_g; j++)
        for (int k = N_g; k < N_phi - N_g; k++) {
            bool test = ((j == j_test) && (k == k_test));
            //
            // assign inverse conformal metric \bar \gamma^{ij}
            //
            tensor gup(gup_rr(j, k), gup_rt(j, k), gup_rp(j, k),
                gup_tt(j, k), gup_tp(j, k), gup_pp(j, k));
            //
            // first compute normal vector \bar s^i and \bar s_i
            //
            vect part_tau(1, -h.dtheta(j, k), -h.dphi(j, k));
            Doub mag = 0.0;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    mag += gup[a][b] * part_tau[a] * part_tau[b];
            lambda[j][k] = 1.0 / (sqrt(mag));
            vect s_low(lambda(j, k) * part_tau[0], lambda(j, k) * part_tau[1], lambda(j, k) * part_tau[2]);
            vect s;
            for (int a = 0; a < 3; a++) {
                s[a] = 0.0;
                for (int b = 0; b < 3; b++)
                    s[a] += gup[a][b] * s_low[b];
            }
            if (test) {
                cout << " h = " << h(j, k) << endl;
                s.print();
                s_low.print();
                gup.print();
            }
            //
            // compute inverse projected metric \bar m^{ij}
            //
            tensor m(gup[0][0] - s[0] * s[0],
                gup[0][1] - s[0] * s[1],
                gup[0][2] - s[0] * s[2],
                gup[1][1] - s[1] * s[1],
                gup[1][2] - s[1] * s[2],
                gup[2][2] - s[2] * s[2]);
            //      if (test) m.print();
            //
            // remember \bar m^{\theta\theta} for FindHorizon()
            // 
            factor[j][k] = -lambda(j, k) * m[1][1];
            //
            // also: invert m and remember m_tt and m_pp for circumferences...
            //
            tensor g = gup.inverse();
            m_pp[j][k] = g[2][2] - s_low[2] * s_low[2];
            m_tt[j][k] = g[1][1] - s_low[1] * s_low[1];
            //
            // NOTE: will rescale sin^2 theta from m_pp - make interpolation more accurate
            //
            const double theta = m_pp.theta(j);
            const double sintheta = sin(theta);
            m_pp[j][k] = m_pp(j, k) / (sintheta * sintheta);
            //
            // compute 2D Laplace operator on h:
            //
            expansion[j][k] = -lambda(j, k) * (m[1][1] * h.ddtheta(j, k) +
                2.0 * m[1][2] * h.dthetadphi(j, k) +
                m[2][2] * h.ddphi(j, k));
            //      if (test) cout << " Laplace of h = " << expansion(j,k) << endl;
            //
            // store gradient of phi in vector
            // 
            vect D_phi(phi_r[j][k], phi_t[j][k], phi_p[j][k]);
            //      if (test) D_phi.print();
            //
            // store extrinsic curvature in tensor
            // 
            tensor A(A_rr[j][k], A_rt[j][k], A_rp[j][k], A_tt[j][k], A_tp[j][k], A_pp[j][k]);
            //
            // store connection coefficients \bar \Gamma^i_{jk}
            // 
            rank3tens Gamma(Gamma_r_rr[j][k], Gamma_r_rt[j][k], Gamma_r_rp[j][k], Gamma_r_tt[j][k], Gamma_r_tp[j][k], Gamma_r_pp[j][k],
                Gamma_t_rr[j][k], Gamma_t_rt[j][k], Gamma_t_rp[j][k], Gamma_t_tt[j][k], Gamma_t_tp[j][k], Gamma_t_pp[j][k],
                Gamma_p_rr[j][k], Gamma_p_rt[j][k], Gamma_p_rp[j][k], Gamma_p_tt[j][k], Gamma_p_tp[j][k], Gamma_p_pp[j][k]);
            //
            // now compute remaining expansion terms...
            // Note cosmological expansion term
            //
            Doub psi2 = exp(2.0 * phi_c[j][k]) * a_friedmann;
            expansion[j][k] -= 2.0 * psi2 * K[j][k] / 3.0;
            for (int a = 0; a < 3; a++) {
                expansion[j][k] += 4.0 * s[a] * D_phi[a];
                for (int b = 0; b < 3; b++) {
                    expansion[j][k] += psi2 * A[a][b] * s[a] * s[b];
                    for (int c = 0; c < 3; c++) {
                        expansion[j][k] -= m[a][b] * s_low[c] * Gamma[c][a][b];
                    }
                }
            }
        }
    //
    // return integral over expansion^2
    // 
    return sqrt(SurfaceIntegral2(expansion, unitsphere));
};
//================================================
// Compute expansion for surface h, return local exansion and normal
//================================================
void HorizonFinder::Expansion(gf2d& h_S, gf2d& exp_S, gf2d& nr, gf2d& nt, gf2d& np) {
    Project(h_S);
    int j_test = -2;
    int k_test = 2;

    for (int j = N_g; j < N_theta - N_g; j++)
        for (int k = N_g; k < N_phi - N_g; k++) {
            bool test = ((j == j_test) && (k == k_test));
            //
            // assign inverse conformal metric \bar \gamma^{ij}
            //
            tensor gup(gup_rr(j, k), gup_rt(j, k), gup_rp(j, k),
                gup_tt(j, k), gup_tp(j, k), gup_pp(j, k));
            //
            // first compute normal vector \bar s^i and \bar s_i
            //
            vect part_tau(1, -h_S.dtheta(j, k), -h_S.dphi(j, k));
            Doub mag = 0.0;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    mag += gup[a][b] * part_tau[a] * part_tau[b];
            lambda[j][k] = 1.0 / (sqrt(mag));
            vect s_low(lambda(j, k) * part_tau[0], lambda(j, k) * part_tau[1], lambda(j, k) * part_tau[2]);
            vect s;
            for (int a = 0; a < 3; a++) {
                s[a] = 0.0;
                for (int b = 0; b < 3; b++)
                    s[a] += gup[a][b] * s_low[b];
            }
            if (test) {
                cout << " h = " << h_S(j, k) << endl;
                s.print();
                s_low.print();
                gup.print();
            }
            //
            // compute inverse projected metric \bar m^{ij}
            //
            tensor m(gup[0][0] - s[0] * s[0],
                gup[0][1] - s[0] * s[1],
                gup[0][2] - s[0] * s[2],
                gup[1][1] - s[1] * s[1],
                gup[1][2] - s[1] * s[2],
                gup[2][2] - s[2] * s[2]);
            //      if (test) m.print();
            //
            // remember \bar m^{\theta\theta} for FindHorizon()
            // 
            factor[j][k] = -lambda(j, k) * m[1][1];
            //
            // also: invert m and remember m_tt and m_pp for circumferences...
            //
            tensor g = gup.inverse();
            m_pp[j][k] = g[2][2] - s_low[2] * s_low[2];
            m_tt[j][k] = g[1][1] - s_low[1] * s_low[1];
            //
            // NOTE: will rescale sin^2 theta from m_pp - make interpolation more accurate
            //
            const double theta = m_pp.theta(j);
            const double sintheta = sin(theta);
            m_pp[j][k] = m_pp(j, k) / (sintheta * sintheta);
            //
            // compute 2D Laplace operator on h:
            //
            exp_S[j][k] = -lambda(j, k) * (m[1][1] * h_S.ddtheta(j, k) +
                2.0 * m[1][2] * h_S.dthetadphi(j, k) +
                m[2][2] * h_S.ddphi(j, k));
            //      if (test) cout << " Laplace of h = " << exp_S(j,k) << endl;
            //
            // store gradient of phi in vector
            // 
            vect D_phi(phi_r[j][k], phi_t[j][k], phi_p[j][k]);
            //      if (test) D_phi.print();
            //
            // store extrinsic curvature in tensor
            // 
            tensor A(A_rr[j][k], A_rt[j][k], A_rp[j][k], A_tt[j][k], A_tp[j][k], A_pp[j][k]);
            //
            // store connection coefficients \bar \Gamma^i_{jk}
            // 
            rank3tens Gamma(Gamma_r_rr[j][k], Gamma_r_rt[j][k], Gamma_r_rp[j][k], Gamma_r_tt[j][k], Gamma_r_tp[j][k], Gamma_r_pp[j][k],
                Gamma_t_rr[j][k], Gamma_t_rt[j][k], Gamma_t_rp[j][k], Gamma_t_tt[j][k], Gamma_t_tp[j][k], Gamma_t_pp[j][k],
                Gamma_p_rr[j][k], Gamma_p_rt[j][k], Gamma_p_rp[j][k], Gamma_p_tt[j][k], Gamma_p_tp[j][k], Gamma_p_pp[j][k]);
            //
            // now compute remaining exp_S terms...
            //
            Doub psi2 = exp(2.0 * phi_c[j][k]) * a_friedmann;
            exp_S[j][k] -= 2.0 * psi2 * K[j][k] / 3.0;
            for (int a = 0; a < 3; a++) {
                exp_S[j][k] += 4.0 * s[a] * D_phi[a];
                for (int b = 0; b < 3; b++) {
                    exp_S[j][k] += psi2 * A[a][b] * s[a] * s[b];
                    for (int c = 0; c < 3; c++) {
                        exp_S[j][k] -= m[a][b] * s_low[c] * Gamma[c][a][b];
                    }
                }
            }
            nr[j][k] = s[0] / psi2;
            nt[j][k] = s[1] / psi2;
            np[j][k] = s[2] / psi2;
        }
};

//================================================
// Project 3D functions onto 2D surfaces
//================================================
void HorizonFinder::Project(gf2d& h_S) {
    for (int j = N_g; j < N_theta - N_g; j++)
        for (int k = N_g; k < N_phi - N_g; k++) {
            double r = h_S(j, k);
            double r2 = r * r;
            double theta = gup_rr.theta(j);
            double sintheta = sin(theta);
            double costheta = cos(theta);
            double sin2theta = sintheta * sintheta;
            gup_rr[j][k] = (*gup_rr_3d)(r, j, k);
            gup_rt[j][k] = (*gup_rt_3d)(r, j, k);
            gup_rp[j][k] = (*gup_rp_3d)(r, j, k);
            gup_tt[j][k] = (*gup_tt_3d)(r, j, k);
            gup_tp[j][k] = (*gup_tp_3d)(r, j, k);
            gup_pp[j][k] = (*gup_pp_3d)(r, j, k);
            //
            // NOTE: 3D data are Delta Gammas, will add spherical 
            //       connection coefficients here so that we end up
            //       with connection associated with conformal metric
            //
            Gamma_r_rr[j][k] = (*DG_r_rr_3d)(r, j, k);
            Gamma_r_rt[j][k] = (*DG_r_rt_3d)(r, j, k);
            Gamma_r_rp[j][k] = (*DG_r_rp_3d)(r, j, k);
            Gamma_r_tt[j][k] = (*DG_r_tt_3d)(r, j, k) - r;
            Gamma_r_tp[j][k] = (*DG_r_tp_3d)(r, j, k);
            Gamma_r_pp[j][k] = (*DG_r_pp_3d)(r, j, k) - r * sin2theta;

            Gamma_t_rr[j][k] = (*DG_t_rr_3d)(r, j, k);
            Gamma_t_rt[j][k] = (*DG_t_rt_3d)(r, j, k) + 1.0 / r;
            Gamma_t_rp[j][k] = (*DG_t_rp_3d)(r, j, k);
            Gamma_t_tt[j][k] = (*DG_t_tt_3d)(r, j, k);
            Gamma_t_tp[j][k] = (*DG_t_tp_3d)(r, j, k);
            Gamma_t_pp[j][k] = (*DG_t_pp_3d)(r, j, k) - sintheta * costheta;

            Gamma_p_rr[j][k] = (*DG_p_rr_3d)(r, j, k);
            Gamma_p_rt[j][k] = (*DG_p_rt_3d)(r, j, k);
            Gamma_p_rp[j][k] = (*DG_p_rp_3d)(r, j, k) + 1.0 / r;
            Gamma_p_tt[j][k] = (*DG_p_tt_3d)(r, j, k);
            Gamma_p_tp[j][k] = (*DG_p_tp_3d)(r, j, k) + costheta / sintheta;
            Gamma_p_pp[j][k] = (*DG_p_pp_3d)(r, j, k);
            //
            // NOTE: will only deal with unscaled quantities in 
            //       horizon finder, so rescale extrinsic curvature here
            //
            A_rr[j][k] = (*a_rr_3d)(r, j, k);
            A_rt[j][k] = (*a_rt_3d)(r, j, k) * r;
            A_rp[j][k] = (*a_rp_3d)(r, j, k) * r * sintheta;
            A_tt[j][k] = (*a_tt_3d)(r, j, k) * r2;
            A_tp[j][k] = (*a_tp_3d)(r, j, k) * r2 * sintheta;
            A_pp[j][k] = (*a_pp_3d)(r, j, k) * r2 * sin2theta;
            K[j][k] = (*K_3d)(r, j, k);

            phi_c[j][k] = (*phi_3d)(r, j, k);
            phi_r[j][k] = (*phi_r_3d)(r, j, k);
            phi_t[j][k] = (*phi_t_3d)(r, j, k);
            phi_p[j][k] = (*phi_p_3d)(r, j, k);
            //
            // finally project electric and magnetic part of Weyl
            //
            if (EB_assigned) {
                E_rr[j][k] = (*E_rr_3d)(r, j, k);
                E_rt[j][k] = (*E_rt_3d)(r, j, k);
                E_rp[j][k] = (*E_rp_3d)(r, j, k);
                E_tt[j][k] = (*E_tt_3d)(r, j, k);
                E_tp[j][k] = (*E_tp_3d)(r, j, k);
                E_pp[j][k] = (*E_pp_3d)(r, j, k);
                B_rr[j][k] = (*B_rr_3d)(r, j, k);
                B_rt[j][k] = (*B_rt_3d)(r, j, k);
                B_rp[j][k] = (*B_rp_3d)(r, j, k);
                B_tt[j][k] = (*B_tt_3d)(r, j, k);
                B_tp[j][k] = (*B_tp_3d)(r, j, k);
                B_pp[j][k] = (*B_pp_3d)(r, j, k);
            }
        }
};
//================================================
// Dump horizon
//================================================
void HorizonFinder::DumpHorizon() {
    //
    // compute tendicity and vorticity
    //
    bool have_tendvort = TendicityVorticity();
    //
    // if we have computed tendicity and vorticity...
    //
    int k = 2;
    const double Mass = HorizonMass();
    surfacefile.setf(ios::right);
    if (have_tendvort) {
        for (int j = N_g; j < N_theta - N_g; j++) {
            double tl = h.theta(j);
            //      double pl = h.phi(k);
            surfacefile << setprecision(8) << setw(16) << phystime
                << setw(16) << tl
                << setw(16) << h[j][k]
                << setw(16) << expansion[j][k]
                << setw(16) << tendicity[j][k] * Mass * Mass
                << setw(16) << vorticity[j][k] * Mass * Mass
                << endl;
        }
    } else {
        for (int j = N_g; j < N_theta - N_g; j++) {
            double tl = h.theta(j);
            surfacefile << setprecision(8) << setw(16) << phystime
                << setw(16) << tl
                << setw(16) << h[j][k]
                << setw(16) << expansion[j][k]
                << endl;
        }
    }
    surfacefile << endl;
};
//================================================
//
// Horizon Diagnostics
//
//================================================
double HorizonFinder::HorizonMass() {
    return sqrt(SurfaceIntegral() / (16.0 * PI));
};
//
// Angular momentum
// 
double HorizonFinder::Spin() {
    // first compute spin integrand
    for (int j = N_g; j < N_theta - N_g; j++)
        for (int k = N_g; k < N_phi - N_g; k++) {
            const double psil = exp(phi_c(j, k));
            const double psi4 = psil * psil * psil * psil * a_friedmann * a_friedmann;
            const double psim4 = 1.0 / psi4;
            //
            // assign inverse of *physical* spatial metric
            //
            tensor gup(psim4 * gup_rr(j, k), psim4 * gup_rt(j, k), psim4 * gup_rp(j, k),
                psim4 * gup_tt(j, k), psim4 * gup_tp(j, k), psim4 * gup_pp(j, k));
            //
            // first compute normal vector s^i and s_i
            //
            vect part_tau(1, -h.dtheta(j, k), -h.dphi(j, k));
            Doub mag = 0.0;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    mag += gup[a][b] * part_tau[a] * part_tau[b];
            lambda[j][k] = 1.0 / (sqrt(mag));
            vect s_low(lambda(j, k) * part_tau[0], lambda(j, k) * part_tau[1], lambda(j, k) * part_tau[2]);
            vect s;
            for (int a = 0; a < 3; a++) {
                s[a] = 0.0;
                for (int b = 0; b < 3; b++)
                    s[a] += gup[a][b] * s_low[b];
            }
            //
            // compute *physical* spatial metric
            //
            tensor g = gup.inverse();
            //
            // now compute integrand for spin integral (see, e.g. eq. (25) in Dreyer etal., PRD 67, 024018 (2003)
            //
            // CHECK: is there a 1/3 missing???
            //
            const double Krp = psi4 * A_rp(j, k) + g[0][2] * K(j, k);
            const double Ktp = psi4 * A_tp(j, k) + g[1][2] * K(j, k);
            const double Kpp = psi4 * A_pp(j, k) + g[2][2] * K(j, k);
            integrand[j][k] = s[0] * Krp + s[1] * Ktp + s[2] * Kpp;
        }
    return SurfaceIntegral(integrand) / (8.0 * PI);
}
//
// Linear momentum
//
double HorizonFinder::LinearMomentum() {
    // first compute linear momentum integrand
    for (int j = N_g; j < N_theta - N_g; j++) {
        const double sintheta = grid->sintheta(j);
        const double costheta = grid->costheta(j);
        for (int k = N_g; k < N_phi - N_g; k++) {
            const double psil = exp(phi_c(j, k));
            const double psi4 = psil * psil * psil * psil * a_friedmann * a_friedmann;
            const double psim4 = 1.0 / psi4;
            //
            // assign inverse of *physical* spatial metric
            //
            tensor gup(psim4 * gup_rr(j, k), psim4 * gup_rt(j, k), psim4 * gup_rp(j, k),
                psim4 * gup_tt(j, k), psim4 * gup_tp(j, k), psim4 * gup_pp(j, k));
            //
            // first compute normal vector s^i and s_i
            //
            vect part_tau(1.0, -h.dtheta(j, k), -h.dphi(j, k));
            Doub mag = 0.0;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    mag += gup[a][b] * part_tau[a] * part_tau[b];
            lambda[j][k] = 1.0 / (sqrt(mag));
            vect s_low(lambda(j, k) * part_tau[0],
                lambda(j, k) * part_tau[1],
                lambda(j, k) * part_tau[2]);
            vect s;
            for (int a = 0; a < 3; a++) {
                s[a] = 0.0;
                for (int b = 0; b < 3; b++)
                    s[a] += gup[a][b] * s_low[b];
            }
            //
            // compute *physical* spatial metric
            //
            tensor g = gup.inverse();
            //
            // now compute integrand for linear momentum integral
            // (see, e.g., eq. (3) in Krishnan, Lousto & Zlochower,
            // PRD 76, 081501(R) (2007)
            //
            // First define auxiliary variable Kij as (K_ij - gamma_ij K)
            const double Kl = K(j, k);
            const double Krr = psi4 * A_rr(j, k) - 2.0 * g[0][0] * Kl / 3.0;
            const double Krt = psi4 * A_rt(j, k) - 2.0 * g[0][1] * Kl / 3.0;
            const double Krp = psi4 * A_rp(j, k) - 2.0 * g[0][2] * Kl / 3.0;
            const double Ktt = psi4 * A_tt(j, k) - 2.0 * g[1][1] * Kl / 3.0;
            const double Ktp = psi4 * A_tp(j, k) - 2.0 * g[1][2] * Kl / 3.0;
            // const double Kpp = psi4 * A_pp(j,k) - 2.0 * g[2][2] * Kl / 3.0;
            //
            const double pr = Krr * s[0] + Krt * s[1] + Krp * s[2];
            const double pt = Krt * s[0] + Ktt * s[1] + Ktp * s[2];
            //
            integrand[j][k] = pr * costheta - pt * sintheta / h[j][k];
        }
    }
    return SurfaceIntegral(integrand) / (8.0 * PI);
}
//
// compute accretion rate (for hydro only...)
//
double HorizonFinder::AccretionRate(gf3d& lapse, Fluxes* fluxes) {
    // first compute accretion rate integrand
    for (int j = N_g; j < N_theta - N_g; j++) {
        //    const double sintheta = grid->sintheta(j);
        //    const double costheta = grid->costheta(j);
        for (int k = N_g; k < N_phi - N_g; k++) {
            const double psil = exp(phi_c(j, k));
            const double psi4 = psil * psil * psil * psil * a_friedmann * a_friedmann;
            const double psim4 = 1.0 / psi4;
            //
            // assign inverse of *physical* spatial metric
            //
            tensor gup(psim4 * gup_rr(j, k), psim4 * gup_rt(j, k), psim4 * gup_rp(j, k),
                psim4 * gup_tt(j, k), psim4 * gup_tp(j, k), psim4 * gup_pp(j, k));
            //
            // first compute normal vector s_i
            //
            vect part_tau(1, -h.dtheta(j, k), -h.dphi(j, k));
            Doub mag = 0.0;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    mag += gup[a][b] * part_tau[a] * part_tau[b];
            lambda[j][k] = 1.0 / (sqrt(mag));
            vect s_low(lambda(j, k) * part_tau[0],
                lambda(j, k) * part_tau[1],
                lambda(j, k) * part_tau[2]);
            // interpolate fluxes to location of horizon
            // [fluxes are computed together with ADM Sources; *not* rescaled]
            vect j_up(fluxes->j_r(h(j, k), j, k),
                fluxes->j_t(h(j, k), j, k),
                fluxes->j_p(h(j, k), j, k));
            //
            // now compute \alpha j^i s_i
            //
            integrand[j][k] = 0.0;
            for (int a = 0; a < 3; a++) {
                integrand[j][k] -= lapse(h(j, k), j, k) * j_up[a] * s_low[a];
            }
        }
    }
    return SurfaceIntegral(integrand);
}
//
// Tendicity and Vorticiy...
//
bool HorizonFinder::TendicityVorticity() {
    if (EB_assigned) {
        for (int j = N_g; j < N_theta - N_g; j++)
            for (int k = N_g; k < N_phi - N_g; k++) {
                const double psil = exp(phi_c(j, k));
                const double psi4 = psil * psil * psil * psil * a_friedmann * a_friedmann;
                const double psim4 = 1.0 / psi4;
                //
                // assign inverse of *physical* spatial metric
                //
                tensor gup(psim4 * gup_rr(j, k), psim4 * gup_rt(j, k), psim4 * gup_rp(j, k),
                    psim4 * gup_tt(j, k), psim4 * gup_tp(j, k), psim4 * gup_pp(j, k));
                //
                // first compute normal vector s^i and s_i
                //
                vect part_tau(1, -h.dtheta(j, k), -h.dphi(j, k));
                Doub mag = 0.0;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        mag += gup[a][b] * part_tau[a] * part_tau[b];
                lambda[j][k] = 1.0 / (sqrt(mag));
                vect s_low(lambda(j, k) * part_tau[0], lambda(j, k) * part_tau[1], lambda(j, k) * part_tau[2]);
                vect s;
                for (int a = 0; a < 3; a++) {
                    s[a] = 0.0;
                    for (int b = 0; b < 3; b++)
                        s[a] += gup[a][b] * s_low[b];
                }
                tensor E(E_rr(j, k), E_rt(j, k), E_rp(j, k), E_tt(j, k), E_tp(j, k), E_pp(j, k));
                tensor B(B_rr(j, k), B_rt(j, k), B_rp(j, k), B_tt(j, k), B_tp(j, k), B_pp(j, k));
                //
                // now compute contractions of (physical) normal vector s^i with E_{ij} and B_{ij}
                //
                double tend = 0.0;
                double vort = 0.0;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        tend += E[a][b] * s[a] * s[b];
                        vort += B[a][b] * s[a] * s[b];
                    }
                tendicity[j][k] = tend;
                vorticity[j][k] = vort;
            }
        tendicity.fill_ghosts();
        vorticity.fill_ghosts();
        return true;
    } else {
        return false;
    }
}
//================================================
// Compute equatorial circumference
//================================================
double HorizonFinder::EquatorialCircumference() {
    // for interpolation:
    m_pp.fill_ghosts();
    phi_c.fill_ghosts();
    double circumference = 0.0;
    const double pi_half = PI / 2.0;
    double dphi = h.dphi();           // assuming uniform grid...
    for (int k = N_g; k < N_phi - N_g; k++) {  // loop around equator
        // NOTE: on equation sin(theta) = 1 - no need to multiply with (sin(theta))^2
        const double m_phiphi = m_pp(pi_half, k);  // interpolate surface functions to equator
        const double psi4 = exp(4.0 * phi_c(pi_half, k)) * a_friedmann * a_friedmann;
        circumference += sqrt(psi4 * m_phiphi) * dphi;
    };
    return circumference;
};
//================================================
// Compute polar circumference
//================================================
//
// REMEMBER: Can't do this, since we don't have coordinates that are adapted to surface, i.e.
// surface does not coincide with surface of constant r
//
// double HorizonFinder::PolarCircumference(int k) {
//   double circumference = 0.0;
//   double dtheta = h.dtheta();                  // assuming uniform grid...
//   const int K = 2 + (k + N_phi/2 - 4) % (N_phi - 4);  // index of location phi + pi
//   for (int j = N_g; j < N_theta - N_g; j++) {  // loop from north to south pole
//     circumference += ( sqrt(exp(4.0*phi_c(j,k))*m_tt(j,k)) + 
// 		       sqrt(exp(4.0*phi_c(j,K))*m_tt(j,K)) ) * dtheta; 
//   };
//   return circumference;
// };    
double HorizonFinder::PolarCircumference(int k) {
    double circumference = 0.0;
    // double d_y = grid->delta_y();
  // index of location phi + pi:
    const int K = (k - N_g + (N_phi - 2 * N_g) / 2) % (N_phi - 2 * N_g) + N_g;
    for (int j = N_g; j < N_theta - N_g; j++) {
        {
            //
            // compute vector Theta^i for "front" side k (see eq. (C.10) in Baumgarte & Shapiro 2010)
            //
            double r = h(j, k);
            vect Theta(h.dtheta(j, k) / r, 1.0 / r, 0.0);
            //
            // assign inverse conformal metric \bar \gamma^{ij}
            //
            tensor gup(gup_rr(j, k), gup_rt(j, k), gup_rp(j, k),
                gup_tt(j, k), gup_tp(j, k), gup_pp(j, k));
            //
            // invert inverse conformal metric to find conformal metric
            //
            tensor g = gup.inverse();
            //
            // compute dot products between Theta and Phi
            //
            double TT = 0.0;
            //
            double psi4 = exp(4.0 * phi_c(j, k)) * a_friedmann * a_friedmann;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) {
                    TT += psi4 * g[a][b] * Theta[a] * Theta[b];
                }
            //
            // now integrate circumference
            // 
            circumference += r * sqrt(TT) * grid->delta_theta(j);
        }
        {
            //
            // now repeat whole thing for "back" side K
            //
            double r = h(j, K);
            vect Theta(h.dtheta(j, K) / r, 1.0 / r, 0.0);
            //
            // assign inverse conformal metric \bar \gamma^{ij}
            //
            tensor gup(gup_rr(j, K), gup_rt(j, K), gup_rp(j, K),
                gup_tt(j, K), gup_tp(j, K), gup_pp(j, K));
            //
            // invert inverse conformal metric to find conformal metric
            //
            tensor g = gup.inverse();
            //
            // compute dot products between Theta and Phi
            //
            double TT = 0.0;
            //
            double psi4 = exp(4.0 * phi_c(j, K)) * a_friedmann * a_friedmann;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) {
                    TT += psi4 * g[a][b] * Theta[a] * Theta[b];
                }
            //
            // now integrate circumference
            // 
            circumference += r * sqrt(TT) * grid->delta_theta(j);
        }
    }
    //    
#ifdef EQSYMMETRY
    circumference *= 2.0;
#endif
    return circumference;
};
//================================================
//
// Surface integrals
//
// follow Appendix C in Baumgarte & Shapiro, 2010 
//
//================================================
double HorizonFinder::SurfaceIntegral() {
    double integral = 0.0;
    if (N_theta == 2 * N_g + 2) {
        // spherical symmetry...
        double r = h(N_g, N_g);
        double r2 = r * r;
        double sintheta = sin(gup_rr.theta(N_g));
        double sin2theta = sintheta * sintheta;
        double psi4 = exp(4.0 * phi_c(N_g, N_g)) * a_friedmann * a_friedmann;
        double h_tt = 1.0 / (gup_tt(N_g, N_g) * r2);
        double h_pp = 1.0 / (gup_pp(N_g, N_g) * r2 * sin2theta);
        //    integral = 4.0*PI*psi4*r2;
        //    cout << "SURFACE_INTEGRAL assuming spherical symmetry " << endl;
        integral = 4.0 * PI * psi4 * r2 * sqrt(h_tt * h_pp);  // CHECK: is this right???
    } else {
        // no symmetry...
        double dphi = grid->delta_phi();
        for (int j = N_g; j < N_theta - N_g; j++) {
            double dtheta = grid->delta_theta(j);
            for (int k = N_g; k < N_phi - N_g; k++)
                integral += SurfaceElement(j, k) * dtheta * dphi;
        }
#ifdef EQSYMMETRY
        integral *= 2.0;
#endif
    }
    return integral;
}
double HorizonFinder::SurfaceIntegral(gf2d& function) {
    double integral = 0.0;
    if (N_theta == 2 * N_g + 2) {
        // spherical symmetry...
        double r = h(N_g, N_g);
        double r2 = r * r;
        double sintheta = sin(gup_rr.theta(N_g));
        double sin2theta = sintheta * sintheta;
        double psi4 = exp(4.0 * phi_c(N_g, N_g)) * a_friedmann * a_friedmann;
        double h_tt = 1.0 / (gup_tt(N_g, N_g) * r2);
        double h_pp = 1.0 / (gup_pp(N_g, N_g) * r2 * sin2theta);
        integral = 4.0 * PI * psi4 * r2 * function(N_g, N_g) * sqrt(h_tt * h_pp);   // CHECK: is this right???
    } else {
        // no symmetry...
        double dphi = grid->delta_phi();
        for (int j = N_g; j < N_theta - N_g; j++) {
            double dtheta = grid->delta_theta(j);
            for (int k = N_g; k < N_phi - N_g; k++)
                integral += SurfaceElement(j, k) * function(j, k) * dtheta * dphi;
        }
#ifdef EQSYMMETRY
        integral *= 2.0;
#endif
    }
    return integral;
}
double HorizonFinder::SurfaceIntegral2(gf2d& function, int type) {
    // if (type == proper) 
    //   cout << " proper integral " << endl;
    // else 
    //   cout << " integral over unit sphere " << endl;
    double integral = 0.0;
    if (N_theta == 2 * N_g + 2) {
        // spherical symmetry...
        double r = h(N_g, N_g);
        double r2 = r * r;
        double sintheta = sin(gup_rr.theta(N_g));
        double sin2theta = sintheta * sintheta;
        double psi4 = exp(4.0 * phi_c(N_g, N_g)) * a_friedmann * a_friedmann;
        double h_tt = 1.0 / (gup_tt(N_g, N_g) * r2);
        double h_pp = 1.0 / (gup_pp(N_g, N_g) * r2 * sin2theta);
        if (type == proper)
            integral = 4.0 * PI * psi4 * r2 * function(N_g, N_g) * function(N_g, N_g) * sqrt(h_tt * h_pp);   // CHECK: is this right???
        else if (type == unitsphere)
            integral = 4.0 * PI * psi4 * function(N_g, N_g) * function(N_g, N_g) * sqrt(h_tt * h_pp);   // CHECK: is this right???
        else
            cerr << " Unknown type " << type <<
            " in HORIZONFINDER:SurfaceIntegral2()! " << endl;
    } else {
        // no symmetry...
        double dphi = grid->delta_phi();
        for (int j = N_g; j < N_theta - N_g; j++) {
            double dtheta = grid->delta_theta(j);
            for (int k = N_g; k < N_phi - N_g; k++)
                integral += SurfaceElement(j, k, type) *
                function(j, k) * function(j, k) * dtheta * dphi;
        }
#ifdef EQSYMMETRY
        integral *= 2.0;
#endif
    }
    return integral;
}
//================================================
// Surface element (for use in surface integrals)
//================================================
inline double HorizonFinder::SurfaceElement(int j, int k, int type) {
    //
    // define vectors Theta and Phi
    //
    double r = h(j, k);
    double theta = gup_rr.theta(j);
    double sintheta = sin(theta);
    double rst = r * sintheta;
    vect Theta(h.dtheta(j, k) / r, 1.0 / r, 0.0);
    vect Phi(h.dphi(j, k) / rst, 0.0, 1.0 / rst);
    //
    // assign inverse conformal metric \bar \gamma^{ij}
    //
    tensor gup(gup_rr(j, k), gup_rt(j, k), gup_rp(j, k),
        gup_tt(j, k), gup_tp(j, k), gup_pp(j, k));
    //
    // invert inverse conformal metric to find conformal metric
    //
    tensor g = gup.inverse();
    //
    // compute dot products between Theta and Phi
    //
    double TT = 0.0;
    double TP = 0.0;
    double PP = 0.0;
    //
    double psi4 = exp(4.0 * phi_c(j, k)) * a_friedmann * a_friedmann;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            TT += psi4 * g[a][b] * Theta[a] * Theta[b];
            TP += psi4 * g[a][b] * Theta[a] * Phi[b];
            PP += psi4 * g[a][b] * Phi[a] * Phi[b];
        }
    //
    // now return surface element (see eq. (C.9) in Baumgarte & Shapiro, 2010)
    // 
    if (type == proper)
        return sqrt(TT * PP - TP * TP) * r * rst;
    else if (type == unitsphere)
        return sqrt(TT * PP - TP * TP) * sintheta;
    else {
        cerr << " Unknown type " << type << " in HORIZONFINDER:SurfaceElement()! " << endl;
        return sqrt(TT * PP - TP * TP) * r * rst;
    }
};


void HorizonFinder::Note(double time, double mass, double AngMom, double LinMom,
    double eqcir, double polcir,
    double h_eq, double h_pole,
    double tend_eq, double tend_pole, double accretion) {
    monitorfile.setf(ios::left);
    const double j_irr = AngMom / (mass * mass);
    const double m_kerr = mass * sqrt(1.0 + j_irr * j_irr / 4.0);
    monitorfile << setw(16) << time
        << setw(24) << setprecision(10) << mass
        << setw(16) << AngMom
        << setw(16) << m_kerr
        << setw(16) << AngMom / (m_kerr * m_kerr)
        << setw(18) << LinMom
        << setw(16) << eqcir
        << setw(16) << polcir
        << setw(16) << h_eq
        << setw(16) << h_pole
        << setw(16) << tend_eq * mass * mass
        << setw(16) << tend_pole * mass * mass
        << setw(16) << accretion
        << endl;
};

bool HorizonFinder::AssignEB(gf3d* E_rr_p, gf3d* E_rt_p, gf3d* E_rp_p, gf3d* E_tt_p, gf3d* E_tp_p, gf3d* E_pp_p,
    gf3d* B_rr_p, gf3d* B_rt_p, gf3d* B_rp_p, gf3d* B_tt_p, gf3d* B_tp_p, gf3d* B_pp_p) {
    EB_assigned = true;
    E_rr_3d = E_rr_p;
    E_rt_3d = E_rt_p;
    E_rp_3d = E_rp_p;
    E_tt_3d = E_tt_p;
    E_tp_3d = E_tp_p;
    E_pp_3d = E_pp_p;
    B_rr_3d = B_rr_p;
    B_rt_3d = B_rt_p;
    B_rp_3d = B_rp_p;
    B_tt_3d = B_tt_p;
    B_tp_3d = B_tp_p;
    B_pp_3d = B_pp_p;
    return EB_assigned;
};

double HorizonFinder::HorizonLocation(int j, int k) {
    if (foundhorizon)
        return h[j][k];
    else
        return -1.0;
};

bool HorizonFinder::InsideHorizon(double r, int j, int k) {
    if (foundhorizon)
        return (r < h(j, k));
    else
        return false;
};

void HorizonFinder::InitialGuess(double M_init, double P_z) {
    //    cout << " HORIZONFINDER: Initializing horizon surface..." << endl;
    const double r_init = 0.5 * M_init;  // makes sense in isotropic coordinates
    for (int j = 0; j < N_theta; j++) {
        double ctl = h.costheta(j);
        for (int k = 0; k < N_phi; k++)
            //	h[j][k] = r_init * (1.0 + 0.1 * sin(h.theta(j)) * cos(h.phi(k)));
            // h[j][k] = r_init * (1.0 - 0.1*sin(h.theta(j) ) );
            h[j][k] = r_init * (1. - P_z * ctl / (8.0 * M_init));
    }
};

