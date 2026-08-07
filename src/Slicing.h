// Tell emacs that this is -*-c++-*- mode
//================================================
// Classes for slicing conditions
//================================================
#ifndef SLICING_H
#define SLICING_H

#include "nr3.h"
#include "tensors.h"
#include "Cosmology.h"
#include "Grid.h"
#include "State.h"
#include "Bona_Masso_function.h"
//
//================================================
// Base class - doesn't do much...
//================================================
//
class Slicing {
protected:
    //  int slicing_type;
    // pointers to all functions that will ever be used...
    Grid* grid;
    Cosmology* cosmology;
    int N_g, N_r, N_t, N_p;
    double eta_KO;
public:
    Slicing(Grid* grid_i, Cosmology* cosmology_i, double eta_KO_i) :
        grid(grid_i), cosmology(cosmology_i), eta_KO(eta_KO_i) {
        N_g = grid->N_ghosts();
        N_r = grid->N_r_tot();
        N_t = grid->N_theta_tot();
        N_p = grid->N_phi_tot();
    };
    virtual ~Slicing() {};
    //================================================
    // Virtual functions that will be overwritten by derived
    // classes below
    //================================================
    virtual void dot_lapse(state* c, state* derivs, double t = 0.0) = 0;
    virtual string Name() = 0;
};
//
//================================================
// Geodesic slicing
//================================================
//
class Geodesic : public Slicing {
public:
    Geodesic(Grid* grid_i, Cosmology* cosmology_i,
        double eta_KO_i) : Slicing(grid_i, cosmology_i, eta_KO_i) {
        cout << " SLICING: setting up geodesic slicing..." << endl;
        //    slicing_type = geodesic;
    };
    ~Geodesic() {};
    string Name() { return "geodesic slicing"; }
    //================================================
    // provide time derivative
    //================================================
    void dot_lapse(state* c, state* derivs, double t = 0.0) {
        for (int i = 0; i < N_r; i++)
            for (int j = 0; j < N_t; j++)
                for (int k = 0; k < N_p; k++) {
                    derivs->lapse[i][j][k] = 0.0;
                }
    };
};
//
//================================================
// (non-advective) 1+log slicing
//================================================
//
class OnePlusLog : public Slicing {
public:
    OnePlusLog(Grid* grid_i, Cosmology* cosmology_i,
        double eta_KO_i) : Slicing(grid_i, cosmology_i, eta_KO_i) {
        cout << " SLICING: setting up (non-advective) 1+log slicing..." << endl;
        //   slicing_type = onepluslog;
    };
    ~OnePlusLog() {};
    string Name() { return "(non-advective) 1+log slicing"; }
    //================================================
    // provide time derivative
    //================================================
    void dot_lapse(state* c, state* derivs, double t = 0.0) {
        const double K_0 = (*cosmology).K0(t);
        for (int i = N_g; i < N_r - N_g; i++)
            for (int j = N_g; j < N_t - N_g; j++)
                for (int k = N_g; k < N_p - N_g; k++) {
                    derivs->lapse[i][j][k] = -2.0 * c->lapse(i, j, k) *
                        (c->K(i, j, k) - K_0) + eta_KO * c->lapse.KO(i, j, k);
                }
    };
};
//
//================================================
// Advective 1+log slicing
//================================================
//
class Advective_OnePlusLog : public Slicing {
public:
    Advective_OnePlusLog(Grid* grid_i, Cosmology* cosmology_i,
        double eta_KO_i) :
        Slicing(grid_i, cosmology_i, eta_KO_i) {
        cout << " SLICING: setting up advective 1+log slicing..." << endl;
        //   slicing_type = advective_onepluslog; 
    };
    ~Advective_OnePlusLog() {};
    string Name() { return "advective 1+log slicing"; }
    //================================================
  // provide time derivative
  //================================================
    void dot_lapse(state* c, state* derivs, double t = 0.0) {
        const double K_0 = (*cosmology).K0(t);
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            for (int j = N_g; j < N_t - N_g; j++) {
                const double sintheta = grid->sintheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    const double shift_r = c->shift_r(i, j, k);
                    const double shift_t = c->shift_t(i, j, k) / rl;
                    const double shift_p = c->shift_p(i, j, k) / (rl * sintheta);
                    derivs->lapse[i][j][k] = -2.0 * c->lapse(i, j, k) *
                        (c->K(i, j, k) - K_0) + eta_KO * c->lapse.KO(i, j, k) +
                        +shift_r * c->lapse.dr(i, j, k, shift_r)
                        + shift_t * c->lapse.dtheta(i, j, k, shift_t)
                        + shift_p * c->lapse.dphi(i, j, k, shift_p);
                }
            }
        }
    };
};
//
//================================================
// harmonic slicing
//================================================
//
class Harmonic : public Slicing {
public:
    Harmonic(Grid* grid_i, Cosmology* cosmology_i,
        double eta_KO_i) : Slicing(grid_i, cosmology_i, eta_KO_i) {
        cout << " SLICING: setting up harmonic slicing..." << endl;
        //   slicing_type = harmonic;
    };
    ~Harmonic() {};
    string Name() { return "harmonic slicing"; }
    //================================================
    // provide time derivative
    //================================================
    void dot_lapse(state* c, state* derivs, double t = 0.0) {
        const double K_0 = (*cosmology).K0(t);
        for (int i = N_g; i < N_r - N_g; i++) {
            for (int j = N_g; j < N_t - N_g; j++) {
                for (int k = N_g; k < N_p - N_g; k++) {
                    derivs->lapse[i][j][k] = -c->lapse(i, j, k) * c->lapse(i, j, k) *
                        (c->K(i, j, k) - K_0) + eta_KO * c->lapse.KO(i, j, k);
                }
            }
        }
    }
};
//
//================================================
// Generalized Advective 1+log slicing:
// allow f(alpha) to be different from 2/alpha
//================================================
//
class BonaMasso : public Slicing {
private:
    BonaMasso_f* bonamasso_f;
public:
    BonaMasso(Grid* grid_i, Cosmology* cosmology_i,
        double eta_KO_i) :
        Slicing(grid_i, cosmology_i, eta_KO_i) {
        cout << " SLICING: setting up (generalized) Bona-Masso slicing..." << endl;
        int bona_masso_type = 1;
        double bona_masso_parameter = 2.0;
        ifstream infile;
        infile.open("Bona_Masso_Input");
        if (!infile) {
            cerr << " SLICING: Can't open input file Bona_Masso_Input " << endl;
            cerr << " SLICING: Will use 1+log instead..." << endl;
        } else {
            int n = 1000;   // size of buffer
            char buf[n], c;
            infile.get(buf, n, '='); infile.get(c); infile >> bona_masso_type;
            infile.get(buf, n, '='); infile.get(c); infile >> bona_masso_parameter;
        }
        if (bona_masso_type == 1) {
            bonamasso_f = new GenOnePlusLog(bona_masso_parameter);
        } else if (bona_masso_type == 2) {
            bonamasso_f = new Ken();
        } else if (bona_masso_type == 3) {
            bonamasso_f = new GaugeShockAvoid(bona_masso_parameter);
        } else if (bona_masso_type == 4) {
            bonamasso_f = new GaugeShockAvoid_lin(bona_masso_parameter);
        } else {
            cerr << " SLICING: Unknown Bona-Masso type " << bona_masso_type << endl;
            cerr << " SLICING: Will use standard 1+log... " << endl;
            bona_masso_parameter = 2.0;
            bonamasso_f = new GenOnePlusLog(bona_masso_parameter);
        }
    };
    ~BonaMasso() {};
    //
    string Name() { return bonamasso_f->Name(); }
    //
    void dot_lapse(state* c, state* derivs, double t = 0.0) {
        const double K_0 = (*cosmology).K0(t);
        for (int i = N_g; i < N_r - N_g; i++) {
            const double rl = grid->r(i);
            for (int j = N_g; j < N_t - N_g; j++) {
                const double sintheta = grid->sintheta(j);
                for (int k = N_g; k < N_p - N_g; k++) {
                    const double lapse_l = c->lapse(i, j, k);
                    const double f = (*bonamasso_f)(lapse_l);
                    const double shift_r = c->shift_r(i, j, k);
                    const double shift_t = c->shift_t(i, j, k) / rl;
                    const double shift_p = c->shift_p(i, j, k) / (rl * sintheta);
                    derivs->lapse[i][j][k] = -lapse_l * lapse_l * f *
                        (c->K(i, j, k) - K_0) + eta_KO * c->lapse.KO(i, j, k) +
                        +shift_r * c->lapse.dr(i, j, k, shift_r)
                        + shift_t * c->lapse.dtheta(i, j, k, shift_t)
                        + shift_p * c->lapse.dphi(i, j, k, shift_p);
                }
            }
        }
    };
};

#endif  /* SLICING_H */
