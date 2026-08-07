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

//================================================
//
// Very simple derived class: vacuum
//
//================================================
class Vacuum : public Matter {
public:
    Vacuum(Grid* grid_i, dumper* dump_i, InData* indata_i, int cowling_i, Cosmology* cosmology_i) :
        Matter(grid_i, dump_i, indata_i, cowling_i, cosmology_i) {
        cout << " MATTER: setting up vacuum... " << endl;
        //
        adm_sources = new ADM_Source_Terms(grid, dump, "vacuum_souces");
        adm_sources->assemble_dump_list("Dump_List");
    };
    ~Vacuum() { cout << " MATTER: destructing derived class Vacuum " << endl; };
    const char* Name() { return "vacuum"; }
    //================================================
    // initialize: nuttin' to do...
    //================================================
    void Initialize(state* s, curvature* c, diagnostics* d) {};
    //===============================================
    // Compute RHS sides for matter equations
    //===============================================
    void Compute_RHS(state* s, curvature* c, double time) {};
    //===============================================
    // Start and finish RK steps
    //===============================================
    void Start_RK() {};
    void Finish_RK(state* s, curvature* c, double dt) {};
    //===============================================
    // Update matter state (adds dt * derivs to update)
    //===============================================
    void Update(double dt) {};
    //===============================================
    // Compute intermediate state (computes inter = last + dt * derivs )
    //===============================================
    void Compute_inter(double dt) {};
    //===============================================
    // Compute Diagnostics 
    //===============================================
    double Compute_Diagnostics(state* s, curvature* c) { return 0.0; };
    //===============================================
    // Regridding etc
    //===============================================
    int Regrid(VecDoub r_new) { return 0; };
    double RegridCriterion() { return 0.0; };
    gf3d* MatterField() { return NULL; };
    //================================================
    // Note
    //================================================
    void Note(int time_step, double time, double tau_c) {};
    //===============================================
    // dump grid functions
    //===============================================
    void dump_fcts(double time, double prop_time, int timestep,
        const char* suffix = "") {
        adm_sources->dump_fcts(time, prop_time, timestep);
    };
    //================================================
    // ADM sources
    //================================================
    void ADM_Sources(state* s, curvature* c) {
        adm_sources->rho_ADM.equals(0.0);
        adm_sources->S_r.equals(0.0);
        adm_sources->S_t.equals(0.0);
        adm_sources->S_p.equals(0.0);
        adm_sources->S_rr.equals(0.0);
        adm_sources->S_rt.equals(0.0);
        adm_sources->S_rp.equals(0.0);
        adm_sources->S_tt.equals(0.0);
        adm_sources->S_tp.equals(0.0);
        adm_sources->S_pp.equals(0.0);
        adm_sources->trace_S.equals(0.0);
    };
};


#include "Maxwell.h"
#include "DualMaxwell.h"
#include "ScalarField.h"
#include "Hydro.h"
#include "RadHydro.h"

#endif  /* MATTER_H */
