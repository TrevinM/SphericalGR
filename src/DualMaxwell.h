// Tell emacs that this is -*-c++-*- mode
#ifndef DUALMAXWELL_H
#define DUALMAXWELL_H


#include "DualMaxwell_State.h"
#include "DualMaxwell_Aux.h"
#include "Monitor.h"
#include "CheckPoint.h"

//================================================
//
// DualMaxwell
//
//================================================
class DualMaxwell : public Matter {
public:
    double rho_center, rho_c_max, rho_max, rho_max_MAX, drhoddr;    // diagnostics 
    int lapse_i, lapse_j, lapse_k, rho_i, rho_j, rho_k;
    dualmaxwell_state* last, * derivs, * inter, * updates;
    dualmaxwell_aux* aux;
    Monitor* monitor;
    CheckPoint* checkpoint;
    ofstream monitorfile;
    double eta_KO;   // Kreiss-Oliger coefficient
    int char_OB;    // decides how Sommerfeld BCs are implemented
    // 0: derivatives; 1: characteristic interpolation
public:
    DualMaxwell(Grid* grid_i, dumper* dump_i, InData* indata_i,
        int cowling_i, int char_OB_i, Cosmology* cosmology_i,
        Monitor* monitor_i, double eta_KO_i, CheckPoint* checkpoint_i);
    ~DualMaxwell();
    const char* Name() { return "DualMaxwell's equations"; }

    //================================================
    // initialize
    //================================================
    void Initialize(state* s, curvature* c, diagnostics* d);

    //===============================================
    // Compute RHS sides for matter equations
    //===============================================
    void Compute_RHS(state* s, curvature* c, double time = 0.0);
    void dot_a_as(dualmaxwell_state* m, state* s, double time = 0.0);
    void dot_as(dualmaxwell_state *m, state *s, curvature *c, double time = 0.0);

    //===============================================
    // Start and finish RK steps
    //===============================================
    void Start_RK();
    void Finish_RK(state* s, curvature* c, double dt);

    //===============================================
    // Update matter state (adds dt * derivs to update)
    //===============================================
    void Update(double dt);

    //===============================================
    // Compute intermediate state (computes inter = last + dt * derivs )
    //===============================================
    void Compute_inter(double dt);

    //===============================================
    // Regridding etc
    //===============================================
    int Regrid(VecDoub r_new);
    double RegridCriterion();
    gf3d* MatterField();

    //================================================
    // Note
    //================================================
    void Note(int time_step, double time, double tau_c);

    //================================================
    // dump grid functions
    //================================================
    void dump_fcts(double time, double prop_time, int timestep,
        const char* suffix);
        
    //================================================
    // Diagnostics
    //================================================
    double Compute_Diagnostics(state* s, curvature* c);

    //================================================
    // ADM sources
    //================================================
    void ADM_Sources(state* s, curvature* c);

    //================================================
    // Compute curl
    //================================================
    void curl(gf3d& a_r, gf3d& a_t, gf3d& a_p,
        double& curl_a_r, double& curl_a_t, double& curl_a_p,
        int i, int j, int k, state* s);
        
    //================================================
    // Compute Lie derivative of covariant vector
    //================================================  
    void Lie_covariant(gf3d& b_r, gf3d& b_t, gf3d& b_p,
        gf3d& a_r, gf3d& a_t, gf3d& a_p,
        gf3d& Ll_r, gf3d& Ll_t, gf3d& Ll_p);
};




#endif