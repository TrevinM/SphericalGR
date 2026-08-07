// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary variables for GR
//
//================================================
#ifndef AUXILIARY_H
#define AUXILIARY_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class auxiliary {
public:
    // storage for Lie derivatives of tensor
    gf3d Lt_rr, Lt_rt, Lt_rp, Lt_tt, Lt_tp, Lt_pp;
    // integrand
    gf3d integrand;
    // divergence of shift
    gf3d div_shift;
    // X = e^{- 2 phi}
    gf3d X;
    // derivatives of phi -- needed for horizon finder
    gf3d dphi_dr, dphi_dt, dphi_dp;
    // more derivatives needed for photon tracker
    gf3d dlapse_dr, dlapse_dt;
    gf3d dshift_r_dr, dshift_t_dr, dshift_r_dt, dshift_t_dt;
    //
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in auxiliary
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
    int N_g, N_r, N_t, N_p;  // ghost zones, *total* number grid points

public:
    //===============================================
    // Constructor
    //===============================================
    auxiliary(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 18;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        int gf_counter = 0;
        //
        // Lie derivative storage
        //
        fct_list[gf_counter] = Lt_rr.setup(grid, 1, "Lt_rr", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = Lt_rt.setup(grid, 1, "Lt_rt", gf_counter, -1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = Lt_rp.setup(grid, 1, "Lt_rp", gf_counter, +1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = Lt_tt.setup(grid, 1, "Lt_tt", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = Lt_tp.setup(grid, 1, "Lt_tp", gf_counter, -1, +1, -1);
        gf_counter++;
        fct_list[gf_counter] = Lt_pp.setup(grid, 1, "Lt_pp", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // divergence of shift
        //
        fct_list[gf_counter] = div_shift.setup(grid, 1, "div_shift", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // integrand
        //
        fct_list[gf_counter] = integrand.setup(grid, 1, "integrand", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // inverse conformal factor
        //
        fct_list[gf_counter] = X.setup(grid, 1, "X", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // derivatives of phi (needed for horizon finder)
        //
        fct_list[gf_counter] = dphi_dr.setup(grid, 2, "dphi_dr", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = dphi_dt.setup(grid, 2, "dphi_dt", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = dphi_dp.setup(grid, 2, "dphi_dp", gf_counter, -1, -1, +1);
        gf_counter++;
        //
        // more derivatives for photon tracker
        // 
        fct_list[gf_counter] = dlapse_dr.setup(grid, 2, "dlapse_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = dlapse_dt.setup(grid, 2, "dlapse_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = dshift_r_dr.setup(grid, 2, "dshift_r_dr", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = dshift_t_dr.setup(grid, 2, "dshift_t_dr", gf_counter);
        gf_counter++;
        fct_list[gf_counter] = dshift_r_dt.setup(grid, 2, "dshift_r_dt", gf_counter);
        gf_counter++;
        fct_list[gf_counter] = dshift_t_dt.setup(grid, 2, "dshift_t_dt", gf_counter);
        gf_counter++;
        //    
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN AUXILIARY!!! " << endl;

        N_g = grid->N_ghosts();
        N_r = grid->N_r_tot();
        N_t = grid->N_theta_tot();
        N_p = grid->N_phi_tot();

    };
    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts() {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].fill_ghosts();
    };
    //===============================================
    // addition
    //===============================================
    void add(double factor, auxiliary* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void add(auxiliary* auxiliary1, double factor, auxiliary* auxiliary2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(auxiliary1->fct_list[i], factor, auxiliary2->fct_list[i]);
    };
    //===============================================
    // equals
    //===============================================
    void equals(auxiliary* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // equals
    //===============================================
    void equals(double number) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(number);
        }
    };
    //===============================================
    // check whether functions are finite
    //===============================================
    bool FINITE() {
        bool fine = true;
        for (int i = 0; i < N_fcts; i++) {
            if (!fct_list[i]->FINITE()) {
                cout << " AUXILIARY: function " << fct_list[i]->Name() << " in "
                    << name << " is not finite!" << endl;
                fine = false;
            }
        }
        return fine;
    };
    //===============================================
    // print list of all functions in auxiliary
    //===============================================
    void function_names() {
        cout << " AUXILIARY: List of all auxiliary functions in " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of auxiliary
    //===============================================
    const char* Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
        const char* suffix = "") {
        if (dump->time_to_dump(timestep) || strcmp(suffix, "")) {
            cout << " AUXILIARY: Dumping auxiliary functions in " << name << " at time t = " << time << endl;
            for (int i = 0; i < N_dump; i++) {
                dump->dump(time, prop_time, timestep, dump_list[i], suffix);
                dump->slice(time, prop_time, timestep, dump_list[i], suffix);
            }
        }
    };
    //===============================================
    // Find all functions to be dumped
    //===============================================
    int assemble_dump_list(const char* dump_list_file) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << " AUXILIARY: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " AUXILIARY: reading dump list from file " << dump_list_file << endl;
        N_dump = 0;
        char fct_name[64];
        infile >> fct_name;
        while (!infile.eof()) {
            for (int i = 0; i < N_fcts; i++) {
                if (!strcmp(fct_name, fct_list[i]->Name())) {
                    cout << " ... found function name " << fct_name << endl;
                    dump_list[N_dump] = fct_list[i]->Address();
                    //	  cout << " function " << fct_name << " : " <<  dump_list[N_dump] << "  " << fct_list[i] << "  " << fct_list[i]->Address() << endl;
                    N_dump++;
                }
            }
            infile >> fct_name;
        }
        if (N_dump > N_fcts) {
            cerr << " AUXILIARY: found too many grid functions in assemble_dump_list() " << endl;
            N_dump = N_fcts;
        }
        return N_dump;
    };
    //===============================================
    // Destructor
    //===============================================
    ~auxiliary() {
        delete fct_list;
        delete dump_list;
        cout << " AUXILIARY: ... closing auxiliary functions " << name << "... " << endl;
    };
    //================================================
    //
    // Compute divergence of shift \bar D_i \beta^i
    // 
    //================================================
    void DivShift(state* s, curvature* c) {
        //
        // NOTE: routine assumes that determinant has already been computed at
        // correct time level!
        //
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            for (int j = N_g; j < N_t - N_g; j++) {
                const double rst = rl * grid->sintheta(j);
                const double ctl = grid->costheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    div_shift[i][j][k] = s->shift_r.dr(i, j, k) + s->shift_t.dtheta(i, j, k) / rl + s->shift_p.dphi(i, j, k) / rst +
                        2.0 * s->shift_r(i, j, k) / rl + s->shift_t(i, j, k) * ctl / rst +
                        (s->shift_r(i, j, k) * c->det.dr(i, j, k) + s->shift_t(i, j, k) * c->det.dtheta(i, j, k) / rl +
                            s->shift_p(i, j, k) * c->det.dphi(i, j, k) / rst) / (2.0 * c->det(i, j, k));
                }
            }
        }
        //  div_shift.fill_ghosts();    // not needed...
    };
    //================================================
    //
    // Lie derivative of scalar - two different version:
    //    one returns Lie derivative at one gridpoint, the other stores
    //    the Lie derivative in a grid function
    //
    // NOTE: uses upwind differencing for advective derivatives  
    //
    //================================================  
    double Lie(gf3d& b_r, gf3d& b_t, gf3d& b_p, gf3d& scalar,
        int i, int j, int k) {
        const double rl = grid->r(i);
        const double stl = grid->sintheta(j);
        return b_r(i, j, k) * scalar.dr(i, j, k, b_r(i, j, k)) +
            b_t(i, j, k) / rl * scalar.dtheta(i, j, k, b_t(i, j, k)) +
            b_p(i, j, k) / (rl * stl) * scalar.dphi(i, j, k, b_p(i, j, k));
    };
    void Lie(gf3d& b_r, gf3d& b_t, gf3d& b_p, gf3d& scalar, gf3d& Lscalar) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    Lscalar[i][j][k] = b_r(i, j, k) * scalar.dr(i, j, k, b_r(i, j, k)) +
                        b_t(i, j, k) / rl * scalar.dtheta(i, j, k, b_t(i, j, k)) +
                        b_p(i, j, k) / (rl * stl) * scalar.dphi(i, j, k, b_p(i, j, k));
                }
            }
        }
    };
    //================================================
    //
    // Lie derivative of contravariant rank-1 tensor (commutator)
    // NOTE: uses upwind differencing for advective derivatives of lambda (but not of shift)
    //
    //================================================  
    void Lie(gf3d& b_r, gf3d& b_t, gf3d& b_p,
        gf3d& l_r, gf3d& l_t, gf3d& l_p,
        gf3d& Ll_r, gf3d& Ll_t, gf3d& Ll_p) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                const double st2 = stl * stl;
                const double ctl = grid->costheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    // 
                    // construct physical variables from vector variables
                    // 
                    const double beta_r = b_r(i, j, k);
                    const double beta_t = b_t(i, j, k) / rl;
                    const double beta_p = b_p(i, j, k) / (rl * stl);
                    const double lamb_r = l_r(i, j, k);
                    const double lamb_t = l_t(i, j, k) / rl;
                    const double lamb_p = l_p(i, j, k) / (rl * stl);
                    //
                    // compute partial derivatives of beta and lambda in terms of code
                    // variables b and l.  Notation:  b^i_{,j} = b_i_j
                    //
                    const double beta_r_r = b_r.dr(i, j, k);
                    const double beta_r_t = b_r.dtheta(i, j, k);
                    const double beta_r_p = b_r.dphi(i, j, k);
                    const double beta_t_r = b_t.dr(i, j, k) / rl - b_t(i, j, k) / r2;
                    const double beta_t_t = b_t.dtheta(i, j, k) / rl;
                    const double beta_t_p = b_t.dphi(i, j, k) / rl;
                    const double beta_p_r = b_p.dr(i, j, k) / (rl * stl) - b_p(i, j, k) / (r2 * stl);
                    const double beta_p_t = b_p.dtheta(i, j, k) / (rl * stl) - b_p(i, j, k) * ctl / (rl * st2);
                    const double beta_p_p = b_p.dphi(i, j, k) / (rl * stl);

                    const double lamb_r_r = l_r.dr(i, j, k, beta_r);
                    const double lamb_r_t = l_r.dtheta(i, j, k, beta_t);
                    const double lamb_r_p = l_r.dphi(i, j, k, beta_p);
                    const double lamb_t_r = l_t.dr(i, j, k, beta_r) / rl - l_t(i, j, k) / r2;
                    const double lamb_t_t = l_t.dtheta(i, j, k, beta_t) / rl;
                    const double lamb_t_p = l_t.dphi(i, j, k, beta_p) / rl;
                    const double lamb_p_r = l_p.dr(i, j, k, beta_r) / (rl * stl) - l_p(i, j, k) / (r2 * stl);
                    const double lamb_p_t = l_p.dtheta(i, j, k, beta_t) / (rl * stl) - l_p(i, j, k) * ctl / (rl * st2);
                    const double lamb_p_p = l_p.dphi(i, j, k, beta_p) / (rl * stl);
                    //
                    // now compute Lie derivative
                    //
                    Ll_r[i][j][k] =
                        beta_r * lamb_r_r + beta_t * lamb_r_t + beta_p * lamb_r_p -
                        lamb_r * beta_r_r - lamb_t * beta_r_t - lamb_p * beta_r_p;

                    Ll_t[i][j][k] =
                        beta_r * lamb_t_r + beta_t * lamb_t_t + beta_p * lamb_t_p -
                        lamb_r * beta_t_r - lamb_t * beta_t_t - lamb_p * beta_t_p;

                    Ll_p[i][j][k] =
                        beta_r * lamb_p_r + beta_t * lamb_p_t + beta_p * lamb_p_p -
                        lamb_r * beta_p_r - lamb_t * beta_p_t - lamb_p * beta_p_p;
                    //
                    // finally rescale coefficients
                    //
                    Ll_t[i][j][k] *= rl;
                    Ll_p[i][j][k] *= rl * stl;
                }
            }
        }
    };
    //================================================
    //
    // Lie derivative of covariant rank-1 tensor 
    // NOTE: uses upwind differencing for advective derivatives of a (but not of shift)
    //
    //================================================  
    void Lie_covariant(gf3d& b_r, gf3d& b_t, gf3d& b_p,
        gf3d& a_r, gf3d& a_t, gf3d& a_p,
        gf3d& Ll_r, gf3d& Ll_t, gf3d& Ll_p) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                const double st2 = stl * stl;
                const double ctl = grid->costheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    // 
                    // construct physical variables from vector variables
                    // 
                    const double beta_r = b_r(i, j, k);
                    const double beta_t = b_t(i, j, k) / rl;
                    const double beta_p = b_p(i, j, k) / (rl * stl);
                    const double A_r = a_r(i, j, k);
                    const double A_t = a_t(i, j, k) * rl;
                    const double A_p = a_p(i, j, k) * (rl * stl);
                    //
                    // compute partial derivatives of beta and A in terms of code
                    // variables b and a.  Notation:  b^i_{,j} = b_i_j
                    //
                    const double beta_r_r = b_r.dr(i, j, k);
                    const double beta_r_t = b_r.dtheta(i, j, k);
                    const double beta_r_p = b_r.dphi(i, j, k);
                    const double beta_t_r = b_t.dr(i, j, k) / rl - b_t(i, j, k) / r2;
                    const double beta_t_t = b_t.dtheta(i, j, k) / rl;
                    const double beta_t_p = b_t.dphi(i, j, k) / rl;
                    const double beta_p_r = b_p.dr(i, j, k) / (rl * stl) - b_p(i, j, k) / (r2 * stl);
                    const double beta_p_t = b_p.dtheta(i, j, k) / (rl * stl) - b_p(i, j, k) * ctl / (rl * st2);
                    const double beta_p_p = b_p.dphi(i, j, k) / (rl * stl);

                    const double A_r_r = a_r.dr(i, j, k, beta_r);
                    const double A_r_t = a_r.dtheta(i, j, k, beta_t);
                    const double A_r_p = a_r.dphi(i, j, k, beta_p);
                    const double A_t_r = a_t.dr(i, j, k, beta_r) * rl + a_t(i, j, k);
                    const double A_t_t = a_t.dtheta(i, j, k, beta_t) * rl;
                    const double A_t_p = a_t.dphi(i, j, k, beta_p) * rl;
                    const double A_p_r = a_p.dr(i, j, k, beta_r) * (rl * stl) + a_p(i, j, k) * stl;
                    const double A_p_t = a_p.dtheta(i, j, k, beta_t) * (rl * stl) + a_p(i, j, k) * ctl * rl;
                    const double A_p_p = a_p.dphi(i, j, k, beta_p) * (rl * stl);
                    //
                    // now compute Lie derivative
                    //
                    Ll_r[i][j][k] =
                        beta_r * A_r_r + beta_t * A_r_t + beta_p * A_r_p +
                        A_r * beta_r_r + A_t * beta_t_r + A_p * beta_p_r;

                    Ll_t[i][j][k] =
                        beta_r * A_t_r + beta_t * A_t_t + beta_p * A_t_p +
                        A_r * beta_t_t + A_t * beta_t_t + A_p * beta_p_t;

                    Ll_p[i][j][k] =
                        beta_r * A_p_r + beta_t * A_p_t + beta_p * A_p_p +
                        A_r * beta_r_p + A_t * beta_t_p + A_p * beta_p_p;
                    //
                    // finally rescale coefficients
                    //
                    Ll_t[i][j][k] /= rl;
                    Ll_p[i][j][k] /= rl * stl;
                }
            }
        }
    };

    //================================================
    //
    // Lie derivative of covariant rank-2 tensors
    //
    // NOTE: we have to treat metric and extrinsic curvature separately, because
    //         for metric need to add 1.0 to some metric components... 
    // NOTE: uses upwind differencing for advective derivatives  
    // NOTE: returns *rescaled* tensor; i.e. automatically divides rt by r etc.
    //
    //================================================  
    void Lie_Metric(gf3d& b_r, gf3d& b_t, gf3d& b_p,
        gf3d& t_rr, gf3d& t_rt, gf3d& t_rp,
        gf3d& t_tt, gf3d& t_tp, gf3d& t_pp) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                const double st2 = stl * stl;
                const double ctl = grid->costheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    // 
                    // construct physical variables from vector variable
                    // 
                    const double beta_r = b_r(i, j, k);
                    const double beta_t = b_t(i, j, k) / rl;
                    const double beta_p = b_p(i, j, k) / (rl * stl);
                    //
                    // compute partial derivatives of beta in terms of code
                    // variable b.  Notation:  b^i_{,j} = b_i_j
                    //
                    const double beta_r_r = b_r.dr(i, j, k);
                    const double beta_r_t = b_r.dtheta(i, j, k);
                    const double beta_r_p = b_r.dphi(i, j, k);

                    const double beta_t_r = b_t.dr(i, j, k) / rl - b_t(i, j, k) / r2;
                    const double beta_t_t = b_t.dtheta(i, j, k) / rl;
                    const double beta_t_p = b_t.dphi(i, j, k) / rl;

                    const double beta_p_r = b_p.dr(i, j, k) / (rl * stl) - b_p(i, j, k) / (r2 * stl);
                    const double beta_p_t = b_p.dtheta(i, j, k) / (rl * stl) - b_p(i, j, k) * ctl / (rl * st2);
                    const double beta_p_p = b_p.dphi(i, j, k) / (rl * stl);
                    //
                    // compute physical variables from tensor variable
                    //
                    const double tens_rr = 1.0 + t_rr(i, j, k);
                    const double tens_rt = t_rt(i, j, k) * rl;
                    const double tens_tr = t_rt(i, j, k) * rl;
                    const double tens_rp = t_rp(i, j, k) * rl * stl;
                    const double tens_pr = t_rp(i, j, k) * rl * stl;
                    const double tens_tt = (1.0 + t_tt(i, j, k)) * r2;
                    const double tens_tp = t_tp(i, j, k) * r2 * stl;
                    const double tens_pt = t_tp(i, j, k) * r2 * stl;
                    const double tens_pp = (1.0 + t_pp(i, j, k)) * r2 * st2;
                    //
                    // Now compute components of Lie derivative
                    //
                    Lt_rr[i][j][k] =
                        beta_r * t_rr.dr(i, j, k, beta_r) +
                        beta_t * t_rr.dtheta(i, j, k, beta_t) +
                        beta_p * t_rr.dphi(i, j, k, beta_p) +
                        2.0 * (tens_rr * beta_r_r + tens_rt * beta_t_r + tens_rp * beta_p_r);
                    Lt_rt[i][j][k] =
                        beta_r * (t_rt.dr(i, j, k, beta_r) * rl + t_rt(i, j, k)) +
                        beta_t * t_rt.dtheta(i, j, k, beta_t) * rl +
                        beta_p * t_rt.dphi(i, j, k, beta_p) * rl +
                        tens_rt * beta_r_r + tens_tt * beta_t_r + tens_pt * beta_p_r +
                        tens_rr * beta_r_t + tens_tr * beta_t_t + tens_pr * beta_p_t;
                    Lt_rp[i][j][k] =
                        beta_r * (t_rp.dr(i, j, k, beta_r) * rl * stl + t_rp(i, j, k) * stl) +
                        beta_t * (t_rp.dtheta(i, j, k, beta_t) * rl * stl + t_rp(i, j, k) * rl * ctl) +
                        beta_p * t_rp.dphi(i, j, k, beta_p) * rl * stl +
                        tens_rp * beta_r_r + tens_tp * beta_t_r + tens_pp * beta_p_r +
                        tens_rr * beta_r_p + tens_tr * beta_t_p + tens_pr * beta_p_p;
                    Lt_tt[i][j][k] =
                        beta_r * (t_tt.dr(i, j, k, beta_r) * r2 + 2.0 * (1.0 + t_tt(i, j, k)) * rl) +
                        beta_t * t_tt.dtheta(i, j, k, beta_t) * r2 +
                        beta_p * t_tt.dphi(i, j, k, beta_p) * r2 +
                        2.0 * (tens_rt * beta_r_t + tens_tt * beta_t_t + tens_pt * beta_p_t);
                    Lt_tp[i][j][k] =
                        beta_r * (t_tp.dr(i, j, k, beta_r) * r2 * stl + 2.0 * t_tp(i, j, k) * rl * stl) +
                        beta_t * (t_tp.dtheta(i, j, k, beta_t) * r2 * stl + t_tp(i, j, k) * r2 * ctl) +
                        beta_p * t_tp.dphi(i, j, k, beta_p) * r2 * stl +
                        tens_rp * beta_r_t + tens_tp * beta_t_t + tens_pp * beta_p_t +
                        tens_rt * beta_r_p + tens_tt * beta_t_p + tens_pt * beta_p_p;
                    Lt_pp[i][j][k] =
                        beta_r * (t_pp.dr(i, j, k, beta_r) * r2 * st2 + 2.0 * (1.0 + t_pp(i, j, k)) * rl * st2) +
                        beta_t * (t_pp.dtheta(i, j, k, beta_t) * r2 * st2 + 2.0 * (1.0 + t_pp(i, j, k)) * r2 * stl * ctl) +
                        beta_p * t_pp.dphi(i, j, k, beta_p) * r2 * st2 +
                        2.0 * (tens_rp * beta_r_p + tens_tp * beta_t_p + tens_pp * beta_p_p);
                    //
                    // finally rescale components
                    //
                    Lt_rt[i][j][k] /= rl;
                    Lt_rp[i][j][k] /= rl * stl;
                    Lt_tt[i][j][k] /= r2;
                    Lt_tp[i][j][k] /= r2 * stl;
                    Lt_pp[i][j][k] /= r2 * st2;
                }
            }
        }
    };

    void Lie_ExCurv(gf3d& b_r, gf3d& b_t, gf3d& b_p,
        gf3d& t_rr, gf3d& t_rt, gf3d& t_rp,
        gf3d& t_tt, gf3d& t_tp, gf3d& t_pp) {
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            for (int j = N_g; j < N_t - N_g; j++) {
                const double stl = grid->sintheta(j);
                const double st2 = stl * stl;
                const double ctl = grid->costheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    // 
                    // construct physical variables from vector variable
                    // 
                    const double beta_r = b_r(i, j, k);
                    const double beta_t = b_t(i, j, k) / rl;
                    const double beta_p = b_p(i, j, k) / (rl * stl);
                    //
                    // compute partial derivatives of beta in terms of code
                    // variable b.  Notation:  b^i_{,j} = b_i_j
                    //
                    const double beta_r_r = b_r.dr(i, j, k);
                    const double beta_r_t = b_r.dtheta(i, j, k);
                    const double beta_r_p = b_r.dphi(i, j, k);

                    const double beta_t_r = b_t.dr(i, j, k) / rl - b_t(i, j, k) / r2;
                    const double beta_t_t = b_t.dtheta(i, j, k) / rl;
                    const double beta_t_p = b_t.dphi(i, j, k) / rl;

                    const double beta_p_r = b_p.dr(i, j, k) / (rl * stl) - b_p(i, j, k) / (r2 * stl);
                    const double beta_p_t = b_p.dtheta(i, j, k) / (rl * stl) - b_p(i, j, k) * ctl / (rl * st2);
                    const double beta_p_p = b_p.dphi(i, j, k) / (rl * stl);
                    //
                    // compute physical variables from tensor variable
                    //
                    const double tens_rr = t_rr(i, j, k);
                    const double tens_rt = t_rt(i, j, k) * rl;
                    const double tens_tr = t_rt(i, j, k) * rl;
                    const double tens_rp = t_rp(i, j, k) * rl * stl;
                    const double tens_pr = t_rp(i, j, k) * rl * stl;
                    const double tens_tt = t_tt(i, j, k) * r2;
                    const double tens_tp = t_tp(i, j, k) * r2 * stl;
                    const double tens_pt = t_tp(i, j, k) * r2 * stl;
                    const double tens_pp = t_pp(i, j, k) * r2 * st2;
                    //
                    // Now compute components of Lie derivative
                    //
                    Lt_rr[i][j][k] =
                        beta_r * t_rr.dr(i, j, k, beta_r) +
                        beta_t * t_rr.dtheta(i, j, k, beta_t) +
                        beta_p * t_rr.dphi(i, j, k, beta_p) +
                        2.0 * (tens_rr * beta_r_r + tens_rt * beta_t_r + tens_rp * beta_p_r);
                    Lt_rt[i][j][k] =
                        beta_r * (t_rt.dr(i, j, k, beta_r) * rl + t_rt(i, j, k)) +
                        beta_t * t_rt.dtheta(i, j, k, beta_t) * rl +
                        beta_p * t_rt.dphi(i, j, k, beta_p) * rl +
                        tens_rt * beta_r_r + tens_tt * beta_t_r + tens_pt * beta_p_r +
                        tens_rr * beta_r_t + tens_tr * beta_t_t + tens_pr * beta_p_t;
                    Lt_rp[i][j][k] =
                        beta_r * (t_rp.dr(i, j, k, beta_r) * rl * stl + t_rp(i, j, k) * stl) +
                        beta_t * (t_rp.dtheta(i, j, k, beta_t) * rl * stl + t_rp(i, j, k) * rl * ctl) +
                        beta_p * t_rp.dphi(i, j, k, beta_p) * rl * stl +
                        tens_rp * beta_r_r + tens_tp * beta_t_r + tens_pp * beta_p_r +
                        tens_rr * beta_r_p + tens_tr * beta_t_p + tens_pr * beta_p_p;
                    Lt_tt[i][j][k] =
                        beta_r * (t_tt.dr(i, j, k, beta_r) * r2 + 2.0 * t_tt(i, j, k) * rl) +
                        beta_t * t_tt.dtheta(i, j, k, beta_t) * r2 +
                        beta_p * t_tt.dphi(i, j, k, beta_p) * r2 +
                        2.0 * (tens_rt * beta_r_t + tens_tt * beta_t_t + tens_pt * beta_p_t);
                    Lt_tp[i][j][k] =
                        beta_r * (t_tp.dr(i, j, k, beta_r) * r2 * stl + 2.0 * t_tp(i, j, k) * rl * stl) +
                        beta_t * (t_tp.dtheta(i, j, k, beta_t) * r2 * stl + t_tp(i, j, k) * r2 * ctl) +
                        beta_p * t_tp.dphi(i, j, k, beta_p) * r2 * stl +
                        tens_rp * beta_r_t + tens_tp * beta_t_t + tens_pp * beta_p_t +
                        tens_rt * beta_r_p + tens_tt * beta_t_p + tens_pt * beta_p_p;
                    Lt_pp[i][j][k] =
                        beta_r * (t_pp.dr(i, j, k, beta_r) * r2 * st2 + 2.0 * t_pp(i, j, k) * rl * st2) +
                        beta_t * (t_pp.dtheta(i, j, k, beta_t) * r2 * st2 + 2.0 * t_pp(i, j, k) * r2 * stl * ctl) +
                        beta_p * t_pp.dphi(i, j, k, beta_p) * r2 * st2 +
                        2.0 * (tens_rp * beta_r_p + tens_tp * beta_t_p + tens_pp * beta_p_p);
                    //
                    // finally rescale components
                    //
                    Lt_rt[i][j][k] /= rl;
                    Lt_rp[i][j][k] /= rl * stl;
                    Lt_tt[i][j][k] /= r2;
                    Lt_tp[i][j][k] /= r2 * stl;
                    Lt_pp[i][j][k] /= r2 * st2;
                }
            }
        }
    };
    //
    //================================================
    // Compute auxiliary function X
    //================================================
    // 
    void Compute_X(state* s, double a_friedmann) {
        for (int i = 0; i < N_r; i++)
            for (int j = 0; j < N_t; j++)
                for (int k = 0; k < N_p; k++) {
                    X[i][j][k] = exp(-2.0 * s->phi(i, j, k)) / a_friedmann;
                }
    }
    //
    //================================================
    // Compute derivatives of phi
    //================================================
    // 
    void Compute_phi_derivs(state* s) {
        for (int i = N_g; i < N_r - N_g; i++)
            for (int j = N_g; j < N_t - N_g; j++)
                for (int k = N_g; k < N_p - N_g; k++) {
                    dphi_dr[i][j][k] = s->phi.dr(i, j, k);
                    dphi_dt[i][j][k] = s->phi.dtheta(i, j, k);
                    dphi_dp[i][j][k] = s->phi.dphi(i, j, k);
                }
        dphi_dr.fill_ghosts();
        dphi_dt.fill_ghosts();
        dphi_dp.fill_ghosts();
    }
    //
    //================================================
    // Compute derivatives for photon tracker
    //================================================
    // 
    void Compute_r_derivs(state* s, curvature* c) {
        s->phi.dr(&dphi_dr);
        s->phi.dtheta(&dphi_dt);
        s->lapse.dr(&dlapse_dr);
        s->lapse.dtheta(&dlapse_dt);
        s->shift_r.dr(&dshift_r_dr);
        s->shift_t.dr(&dshift_t_dr);
        s->shift_r.dtheta(&dshift_r_dt);
        s->shift_t.dtheta(&dshift_t_dt);
    }
};


#endif  /* AUXILIARY_H */
