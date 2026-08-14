// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class that contains matter...
//
// class Matter is a base class, that provides a common interface for
// matter model.  
//
//================================================

#ifndef MATTER_H
#define MATTER_H

#include "InData.h"
#include "nr3.h"
#include "gridfunction.h"
#include "Cosmology.h"
#include "Grid.h"
#include "dumper.h"
#include "InData.h"
#include "Curvature.h"
#include "Diagnostics.h"
#include "ADM_Source_Terms.h"
#include "Fluxes.h"
#include <ctime>

class Matter {
public:
    Grid* grid;
    dumper* dump;
    InData* indata;
    Cosmology* cosmology;
    //  matter_state * m_state;
    int cowling;    // Cowling approximation 
    // (0: no, 1: fix gravity, 2: fix matter)  
    ADM_Source_Terms* adm_sources;
    Fluxes* fluxes;
    int N_g, N_r, N_t, N_p;   // number of ghost and grid points
    double PI;
public:
    //================================================
    // Constructor
    //================================================
    Matter(Grid* grid_i, dumper* dump_i, InData* indata_i,
        int cowling_i, Cosmology* cosmology_i) :
        grid(grid_i), dump(dump_i), indata(indata_i), cowling(cowling_i),
        cosmology(cosmology_i) {
        //
        N_g = grid->N_ghosts();
        N_r = grid->N_r_tot();
        N_t = grid->N_theta_tot();
        N_p = grid->N_phi_tot();
        PI = acos(-1.0);
        fluxes = NULL;
    };
    //================================================
    // Destructor
    //================================================
    virtual ~Matter() {
        cout << " MATTER: destructing base class Matter " << endl;
    }
    //================================================
    // Name
    //================================================
    virtual const char* Name() = 0;
    //================================================
    // Initialize
    //================================================
    virtual void Initialize(state* s, curvature* c, diagnostics* d) = 0;
    //================================================
    // Compute ADM Sources
    //================================================
    virtual void ADM_Sources(state* s, curvature* c) = 0;
    //===============================================
    // Compute RHS sides for matter equations
    //===============================================
    virtual void Compute_RHS(state* s, curvature* c, double time = 0.0) = 0;
    //===============================================
    // Start and finish RK steps
    //===============================================
    virtual void Start_RK() = 0;
    virtual void Finish_RK(state* s, curvature* c, double dt) = 0;
    //===============================================
    // Update matter state (adds dt * derivs to update)
    //===============================================
    virtual void Update(double dt) = 0;
    //===============================================
    // Compute intermediate state (computes inter = last + dt * derivs )
    //===============================================
    virtual void Compute_inter(double dt) = 0;
    //===============================================
    // Compute Diagnostics 
    //===============================================
    virtual double Compute_Diagnostics(state* s, curvature* c) = 0;
    //===============================================
    // Regridding etc
    //===============================================
    virtual int Regrid(VecDoub r_new) = 0;
    virtual double RegridCriterion() = 0;
    virtual gf3d* MatterField() = 0;
    //================================================
    // Note
    //================================================
    virtual void Note(int time_step, double time, double tau_c) = 0;
    //================================================
    // dump grid functions
    //================================================
    virtual void dump_fcts(double time = 0.0, double prop_time = 0.0,
        int timestep = 0, const char* suffix = "") = 0;
};

#endif  /* MATTER_H */
