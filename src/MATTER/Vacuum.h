#ifndef VACUUM_H
#define VACUUM_H

#include "Matter.h"
#include "Matter_State.h"

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


class vacuum_state : public matter_state {
public:
    gf3d e_p, a_p;    // lower phi-components of fields, rescaled
    //===============================================
    // Constructor
    //===============================================
public:
    vacuum_state(Grid* grid_i, dumper* dump_i, const char* name_i) :
        matter_state(grid_i, dump_i, name_i) {
        cout << " MATTER: creating matter state for vacuum... " << endl;
        N_fcts = 0;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        //
        int gf_counter = 0;
        //
    }
    ~vacuum_state() {
    };
};










#endif