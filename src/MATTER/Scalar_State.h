// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of matter, i.e. collection
// of grid functions for dynamical matter variables
//
//================================================
#ifndef SCALAR_STATE_H
#define SCALAR_STATE_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class scalar_state {
public:
    //
    // Note: will add functions to scalar_state through derived classes
    //
    gf3d sf, pi;
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in scalar_state
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    scalar_state(Grid* grid_i, dumper* dump_i, const char* name_i);

    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts();

    //===============================================
    // addition
    //===============================================
    void add(double factor, scalar_state* rhs);
    void add(scalar_state* scalar_state1, double factor, scalar_state* scalar_state2);

    //===============================================
    // addition
    //===============================================
    void equals(scalar_state* rhs);

    //===============================================
    // print list of all functions in scalar_state
    //===============================================
    void function_names();

    //===============================================
    // return name of scalar_state
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
    int assemble_dump_list(const char* dump_list_file) ;

    //===============================================
    // Regrid
    //===============================================
    void Regrid(VecDoub r_new);

    //===============================================
    // Destructor
    //===============================================
    ~scalar_state();

};


#endif  /* SCALAR_STATE_H */
