// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of matter, i.e. collection
// of grid functions for dynamical matter variables
//
//================================================
#ifndef RADHYDRO_STATE_H
#define RADHYDRO_STATE_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class radhydro_state {
public:
    //
    // conserved fluid variables
    //  det_3D = \sqrt( \gamma / \hat \gamma )
    //
    //   D = det_3D W \rho_0
    //   S_r = det_3D W^2 \rho_0 h v_r
    //   S_t = det_3D W^2 \rho_0 h v_t / r 
    //   S_p = det_3D W^2 \rho_0 h v_p / (r sin \theta) 
    //   tau = det_3D (W^2 \rho_0 h - p) - W \rho_0)
    //
    //   tau_rad = det_3D * ( (4/3) W^2 E + 2 alpha W F^0 - E/3)
    //   S_rad_r = det_3D * lapse * ( (4/3) E u^t u_r + F^0 u_r + F_r u^t) 
    //   S_rad_t = det_3D * lapse * ( (4/3) E u^t u_t + F^0 u_t + F_t u^t) / r
    //   S_rad_p = det_3D * lapse * ( (4/3) E u^t u_p + F^0 u_p + F_p u^t) / (r sin \theta)
    //  NOTE: above definition of S_rad_i is consistent with
    //  S_rad_i = det_3D * ( (4/3) E W^2 v_a + W f_a + f W v_a )
    // 
    gf3d D, tau, S_r, S_t, S_p;
    gf3d tau_rad, S_rad_r, S_rad_t, S_rad_p;
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in radhydro_state
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    radhydro_state(Grid* grid_i, dumper* dump_i, const char* name_i);
    
    //===============================================
    // check properties
    //===============================================
    void check_props();

    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts();

    //===============================================
    // addition
    //===============================================
    void add(double factor, radhydro_state* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(radhydro_state* radhydro_state1, double factor, radhydro_state* radhydro_state2);

    //===============================================
    // equals
    //===============================================
    void equals(radhydro_state* rhs);
    void equals(double number);
    //===============================================
    // check whether states are equal
    //===============================================
    bool IsEqualTo(radhydro_state* s);

    //===============================================
    // print list of all functions in radhydro_state
    //===============================================
    void function_names();

    //===============================================
    // return name of radhydro_state
    //===============================================
    const char* Name();

    //===============================================
    // Dump grid functions
    //===============================================
    void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
        const char* suffix = "");

    //===============================================
    // Find all functions to be dumped
    //===============================================
    int assemble_dump_list(const char* dump_list_file);

    //===============================================
    // Regrid
    //===============================================
    void Regrid(VecDoub r_new);

    //===============================================
    // Destructor
    //===============================================
    ~radhydro_state();

};


#endif  /* RADHYDRO_STATE_H */
