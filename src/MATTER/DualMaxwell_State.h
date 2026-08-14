// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of matter, i.e. collection
// of grid functions for dynamical matter variables
//
//================================================
#ifndef DUALMAXWELL_STATE_H
#define DUALMAXWELL_STATE_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class dualmaxwell_state {
public:
    //
    // Note: will add functions to dualmaxwell_state through derived classes
    //
    //  rescaled vector potentials for B and E, indices downstairs
    //
    gf3d a_r, a_t, a_p, as_r, as_t, as_p;
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in dualmaxwell_state
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    dualmaxwell_state(Grid* grid_i, dumper* dump_i, const char* name_i);

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
    void add(double factor, dualmaxwell_state* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(dualmaxwell_state* dualmaxwell_state1, double factor, dualmaxwell_state* dualmaxwell_state2);

    //===============================================
    // addition
    //===============================================
    void equals(dualmaxwell_state* rhs);

    //===============================================
    // apply characteristic outer boundaries
    //===============================================
    void char_OB(dualmaxwell_state* last, double dt);

    //===============================================
    // print list of all functions in dualmaxwell_state
    //===============================================
    void function_names();

    //===============================================
    // return name of dualmaxwell_state
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
    ~dualmaxwell_state();
    
};


#endif  /* DUALMAXWELL_STATE_H */
