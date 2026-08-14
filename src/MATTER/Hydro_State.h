// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of matter, i.e. collection
// of grid functions for dynamical matter variables
//
//================================================
#ifndef HYDRO_STATE_H
#define HYDRO_STATE_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class hydro_state {
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
    gf3d D, tau, S_r, S_t, S_p;
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in hydro_state
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    hydro_state(Grid* grid_i, dumper* dump_i, const char* name_i);

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
    void add(double factor, hydro_state* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(hydro_state* hydro_state1, double factor, hydro_state* hydro_state2);

    //===============================================
    // equals
    //===============================================
    void equals(hydro_state* rhs);
    void equals(double number);

    //===============================================
    // check whether states are equal
    //===============================================
    bool IsEqualTo(hydro_state* s);
    
    //===============================================
    // print list of all functions in hydro_state
    //===============================================
    void function_names();

    //===============================================
    // return name of hydro_state
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
    ~hydro_state();

};


#endif  /* HYDRO_STATE_H */
