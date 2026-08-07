// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Horizon finder
//
//================================================

#ifndef HORIZONFINDER_H
#define HORIZONFINDER_H

#include "EllSolver2D.h"

#include "nr3.h"
#include "surfacefunction.h"
#include "gridfunction.h"
// #include "Monitor.h"
#include "Grid.h"
#include "State.h"
#include "Curvature.h"
#include "Auxiliary.h"
#include "Fluxes.h"
#include <ctime>

enum { proper, unitsphere };

//
//================================================
// 
//================================================
//
class HorizonFinder {
private:
    // function h in level surface function tau = r - h(theta,phi)
    gf2d h, h_old;
    // projected metric
    gf2d gup_rr, gup_rt, gup_rp, gup_tt, gup_tp, gup_pp;
    // induced metric
    gf2d m_rr, m_rt, m_rp, m_tt, m_tp, m_pp;
    // conformal Christoffel symbols
    gf2d Gamma_r_rr, Gamma_r_rt, Gamma_r_rp, Gamma_r_tt, Gamma_r_tp, Gamma_r_pp;
    gf2d Gamma_t_rr, Gamma_t_rt, Gamma_t_rp, Gamma_t_tt, Gamma_t_tp, Gamma_t_pp;
    gf2d Gamma_p_rr, Gamma_p_rt, Gamma_p_rp, Gamma_p_tt, Gamma_p_tp, Gamma_p_pp;
    // extrinsic curvature
    gf2d A_rr, A_rt, A_rp, A_tt, A_tp, A_pp, K;
    // conformal exponent and derivatives
    gf2d phi_c, phi_r, phi_t, phi_p;
    // quantity needed for angular and linear momentum integrals
    gf2d integrand;
    // 
    gf2d lambda, expansion, res, rhs_grid, factor;
    // 
    // pointers to 3D grid functions
    //
    gf3d* gup_rr_3d, * gup_rt_3d, * gup_rp_3d;
    gf3d* gup_tt_3d, * gup_tp_3d, * gup_pp_3d;
    gf3d* DG_r_rr_3d, * DG_r_rt_3d, * DG_r_rp_3d, * DG_r_tt_3d, * DG_r_tp_3d, * DG_r_pp_3d;
    gf3d* DG_t_rr_3d, * DG_t_rt_3d, * DG_t_rp_3d, * DG_t_tt_3d, * DG_t_tp_3d, * DG_t_pp_3d;
    gf3d* DG_p_rr_3d, * DG_p_rt_3d, * DG_p_rp_3d, * DG_p_tt_3d, * DG_p_tp_3d, * DG_p_pp_3d;
    gf3d* phi_3d, * phi_r_3d, * phi_t_3d, * phi_p_3d;
    gf3d* a_rr_3d, * a_rt_3d, * a_rp_3d;
    gf3d* a_tt_3d, * a_tp_3d, * a_pp_3d;
    gf3d* K_3d;
    //
    // functions needed for vorticity and tencidity
    //
    bool EB_assigned;
    gf3d* E_rr_3d, * E_rt_3d, * E_rp_3d, * E_tt_3d, * E_tp_3d, * E_pp_3d;
    gf3d* B_rr_3d, * B_rt_3d, * B_rp_3d, * B_tt_3d, * B_tp_3d, * B_pp_3d;
    gf2d E_rr, E_rt, E_rp, E_tt, E_tp, E_pp;
    gf2d B_rr, B_rt, B_rp, B_tt, B_tp, B_pp;
    gf2d tendicity, vorticity;
    //
    int N_theta, N_phi;   // *include* ghost zones
    int N_g;
    bool foundhorizon;
    int timestep;
    Doub phystime;
#ifndef NoEllSolver
    EllSolver2D* ellsolver;
#endif
    //
    // parameters for solver (read in from file "AH_finder_Input")
    // 
    double tol_exp;    // expansion tolerance
    double tol_ell;    // tolerance for elliptic solver
    int max_iter_ell;  // maximum number of iterations in elliptic solver
    int max_iter_exp;  // maximum number of iterations for expansion
    double eta;      // parameter for expansion equation: \eta h 
    // added to both sides
    int find_steps;   // timesteps or times after which horizon should be found
    double find_times;
    int next_step, last_step;
    double next_time;
    double mass_guess;
    double a_friedmann;      // current value of cosmological expansion coefficient
    bool use_time_step_criterion;  // use step counter to decide whether to 
    // search for horizon; otherwise use time
//  char file_stem[64];
    ofstream monitorfile;
    ofstream surfacefile;
    double PI;
    double r_min, r_max;
    Grid* grid;
public:
    //================================================
    // Constructor
    //================================================
    HorizonFinder(Grid* grid_i, state* s, curvature* c, auxiliary* aux,
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
    // Destructor
    //================================================
    ~HorizonFinder() {
#ifndef NoEllSolver
        delete ellsolver;
#endif
        monitorfile.close();
        surfacefile.close();
    };
    //================================================
    // Read input
    //================================================
    int ReadInput();
    //================================================
    // time to find horizon?
    //================================================
    bool TimeToFindHorizon(int timestep, Doub time);
    //================================================
    // search for horizon (return true if horizon found)
    //================================================
    bool FindHorizon(int timestep, Doub time,
        state* s, Fluxes* fluxes,
        Doub adm_mass, Doub mom_guess, Doub a);
    //================================================
    // Horizon diagnostics
    //================================================
    double HorizonMass();
    double Spin();
    double LinearMomentum();
    double AccretionRate(gf3d& lapse, Fluxes* fluxes);
    //================================================
    // Adjust parameters
    //================================================
    double SetFindTimes(double time) { return find_times = time; };
    double SetFindSteps(int step) {
        find_steps = step;
        use_time_step_criterion = (find_steps > 0);
        return find_steps;
    };
    double SetMassGuess(double mass) { return mass_guess = mass; }
    void SetRMaxMin() {
        r_min = (*phi_t_3d).r(N_g);
        const int N_r = phi_t_3d->dim1();
        r_max = (*phi_t_3d).r(N_r - N_g);
    }
private:
    //================================================
    // Compute expansion
    //================================================
    double Expansion();
    //================================================
    // Compute expansion for surface h, return local exansion and normal
    //================================================
public:
    void Expansion(gf2d& h, gf2d& exp_at_h, gf2d& nr, gf2d& nt, gf2d& np);
    //================================================
    // Initial guess for horizon location
    // (use approximate result from Dennison et.al., PRD 74, 064016 (2006), eq. (24))
    //================================================
private:
    void InitialGuess() { InitialGuess(mass_guess); };
    void InitialGuess(double M_init, double P_z = 0) {
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
    //================================================
    // Project 3D functions onto 2D surfaces
    //================================================
    void Project(gf2d& h);
    //================================================
    // Dump horizon
    //================================================
    void DumpHorizon();
    void Note(double time, double mass, double AngMom, double LinMom,
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
    //================================================
    // Surface integrals
    //================================================
    double SurfaceIntegral();
    double SurfaceIntegral(gf2d& function);
    double SurfaceIntegral2(gf2d& function, int type = proper);
    inline double SurfaceElement(int j, int k, int type = proper);
    double EquatorialCircumference();
    double PolarCircumference() { return PolarCircumference(N_g); };
    double PolarCircumference(int k);
    //================================================
    // Routines for tendicity and vorticity
    //================================================
public:
    bool AssignEB(gf3d* E_rr_p, gf3d* E_rt_p, gf3d* E_rp_p, gf3d* E_tt_p, gf3d* E_tp_p, gf3d* E_pp_p,
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
private:
    bool TendicityVorticity();
    //================================================
    // Utility to decide whether point is inside horizon
    //================================================
public:
    double HorizonLocation(int j, int k) {
        if (foundhorizon)
            return h[j][k];
        else
            return -1.0;
    };
    bool InsideHorizon(double r, int j, int k) {
        if (foundhorizon)
            return (r < h(j, k));
        else
            return false;
    };
};

#endif   /* HORIZONFINDER_H */

