// Tell emacs that this is -*-c++-*- mode
//
//================================================
// 
// Initial for electromagnetic waves in electrovacuum, expressed in terms of two
// vector potentials A and *A with
//
//       B = curl A   and E = curl *A
//
// Will generate divergence-free A and A* by writing
//       A = A_tor   + curl V_A
//       *A = *A_tor + curl *V_A
// where only non-zero components of A_tor, A_tor, V_A and *V_A
// are the phi components, which
// depend on r and theta only (i.e. all purely toroidal).  Then:
//       B_pol = curl A_tor         B_tor = curl curl V_A
//       E_pol = curl *A_tor        E_tor = curl curl *V_A
//================================================
//
#ifndef DUALEMWAVE_H
#define DUALEMWAVE_H

#include "InData.h"



class DualEMWave : public InData {
private:
#ifndef NoEllSolver
    FlatEllSolver3D* laplace;
    VecLaplacian* veclaplacian;
#endif
    gf3d psi, rho, res, u, delta_psi;
    // all vector and tensor components rescaled...
    gf3d a_r, a_t, a_p;   // components of vector potential A (indices down)  
    gf3d as_r, as_t, as_p;  // components of vector potential *A (indices down)
    gf3d s_r, s_t, s_p;     // components of momentum densities (upstairs)
    gf3d W_r, W_t, W_p;     // vector potential for ext.  curv. (indices up)
    gf3d del_W_r, del_W_t, del_W_p;  // corrections...
    gf3d res_r, res_t, res_p;  // (upstairs...)
    // now extrinsic curvature: use *initial data* rescaling here,
    // and *upper* (geometrically) rescaled indices, i.e.
    // \bar A^{ij} = \psi^10 A^{ij}
    // need to (a) lower indices, and (b) convert to BSSN rescaling in
    // functions a_ij_analytical...
    gf3d A_rr, A_rt, A_rp, A_tt, A_tp, A_pp, A2;
    // gf3d v_a_p, v_as_p;     // phi component of tor. vectors that generate pol. parts of A and *A
    int n_r, n_theta, n_phi;
    int N_g;
    int n_psi;  // power of psi in vectors A and V_A, ie. A \propto psi^n etc
    int max_it;   // parameters for elliptic solver
    double tol;
    double a1_amp, a2_amp, a3_amp, as1_amp, as2_amp, as3_amp;   // parameters for initial data
    int n_a, n_v_a, n_as, n_v_as;       // powers of r in respective seed functions
    double r0;
    int flat;
    int dual_sign;
    double PI;
    bool all_clear;
    ostringstream indata_name;
public:
    //================================================
    // Constructor
    //================================================
    DualEMWave(char* indata_input, Grid* grid_i, Cosmology* cosmology);
    
    //================================================
    // Destructor
    //================================================
    ~DualEMWave() {};
    string Name() {
        return indata_name.str();
    };
    //================================================
    // Initializer
    //================================================
    bool Initialize(gf3d& fct);

    //================================================
    // Solve constraints
    //================================================
    double Solve_Constraints(bool verbose = false);

    //================================================
    // Total residual
    //================================================
    double Residual();

    //================================================
    // Analytical solution for h_{ij}
    //================================================
    double h_rr_analytical(double rl, double thetal, double phil, double t) override;
    double h_rt_analytical(double r, double theta, double phi, double t) override;
    double h_rp_analytical(double r, double theta, double phi, double t) override;
    double h_tt_analytical(double rl, double thetal, double phil, double t) override;
    double h_tp_analytical(double r, double theta, double phi, double t) override;
    double h_pp_analytical(double rl, double thetal, double phil, double t) override;
    double phi_analytical(double rl, double thetal, double phil, double tl) override;
    //================================================
    // Analytical solution for connection coefficients
    //================================================
    double lam_r_analytical(double r, double theta, double phi, double t,
        bool& done) override;
    double lam_t_analytical(double r, double theta, double phi, double t,
        bool& done) override;
    double lam_p_analytical(double r, double theta, double phi, double t,
        bool& done) override;
    //================================================
    // Analytical solution for extrinsic curvature
    // NOTE: in initialize we compute \bar A_{ij} = psi^6 \tilde A_{ij},
    // now need BSSN rescaling of extrinsic curvature...
    //================================================
    double a_rr_analytical(double r, double theta, double phi, double t) override;
    double a_rt_analytical(double r, double theta, double phi, double t) override;
    double a_rp_analytical(double r, double theta, double phi, double t) override;
    double a_tt_analytical(double r, double theta, double phi, double t) override;
    double a_tp_analytical(double r, double theta, double phi, double t) override;
    double a_pp_analytical(double r, double theta, double phi, double t) override;
    double K_analytical(double r, double theta, double phi, double t) override;
    //================================================
    // Analytical solution for gauge
    //================================================
    double lapse_analytical(double rl, double thetal, double phil, double t) override;
    double shift_r_analytical(double r, double theta, double phi, double t) override;
    double shift_t_analytical(double r, double theta, double phi, double t) override;
    double shift_p_analytical(double r, double theta, double phi, double t) override;
    //================================================
    // Analytical solution for auxiliary functions
    //================================================
    double Theta_analytical(double r, double theta, double phi, double t) override;
    double B_r_analytical(double r, double theta, double phi, double t) override;
    double B_t_analytical(double r, double theta, double phi, double t) override;
    double B_p_analytical(double r, double theta, double phi, double t) override;
    //================================================
    // Analytical solution for hydro
    //================================================
    double rho_0_analytical(double rl, double thetal, double phil, double tl) override;
    double P_analytical(double r, double theta, double phi, double t) override;
    double v_r_analytical(double r, double theta, double phi, double t) override;
    double v_t_analytical(double r, double theta, double phi, double t) override;
    double v_p_analytical(double r, double theta, double phi, double t) override;
    //================================================
    // Analytical solution for scalar field
    //================================================
    double sf_analytical(double r, double theta, double phi, double t) override;
    double pi_analytical(double r, double theta, double phi, double t) override;
    //================================================
    // Analytical solution for Maxwell OR for Dual Maxwell:
    // For Maxwell need e_p and a_p only, for Dual Maxwell all a_i and as_i
    //================================================
    // returns correct function at t=0 only
    double e_p_analytical(double r, double theta, double phi, double t) override;
    //
    // provides flat-space analytical result at any time,
    // but only for r_0 = 0; see
    // eq. (16) in Knapp, Walker & Baumgarte, 2002
    //
    double a_r_analytical(double r, double theta, double phi, double t) override;
    double a_t_analytical(double r, double theta, double phi, double t) override;
    double a_p_analytical(double r, double theta, double phi, double t) override;
    double as_r_analytical(double r, double theta, double phi, double t) override;
    double as_t_analytical(double r, double theta, double phi, double t) override;
    double as_p_analytical(double r, double theta, double phi, double t) override;
    //================================================
    // Analytical solution for radiation
    //================================================
    double E_analytical(double r, double theta, double phi, double t) override;
    double F_0_analytical(double r, double theta, double phi, double t) override;
    double F_r_analytical(double r, double theta, double phi, double t) override;
    double F_t_analytical(double r, double theta, double phi, double t) override;
    double F_p_analytical(double r, double theta, double phi, double t) override;
    
private:
    //===============================================================
    // Compute fields (including rho_ADM and S^i)
    //===============================================================
#include "DUALMAXWELL/Compute_Fields.h"
    //===============================================================
    // Solve Hamiltonian constraint
    //===============================================================
#include "DUALMAXWELL/Solve_Hamiltonian.h"
    //===============================================================
    // Solve momentum constraints
    //===============================================================
#include "DUALMAXWELL/Solve_Momentum.h"
    //=============================================================== 
    // Compute \bar A^{ij} from vector potential W^i
    //=============================================================== 
    void Compute_Aij();

    //===============================================================
    // Compute physical curl: returns scaled eps^ijk ( D_j A_k - D_k A_j )
    // where D_i is covariant derivative with respect to reference metric,
    // but epsilon is physical (rescaled) 3D epsilon: psi^{-6}[ijk]
    //===============================================================
    void curl(gf3d& a_r, gf3d& a_t, gf3d& a_p,
        double& curl_a_r, double& curl_a_t, double& curl_a_p,
        int i, int j, int k);

    //================================================
    // Dump function
    //================================================
    void dump(gf3d* fct, int step = 0);
};

#endif