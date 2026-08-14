// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of matter, i.e. collection
// of grid functions for dynamical matter variables
//
//================================================
#ifndef MAXWELL_STATE_H
#define MAXWELL_STATE_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class maxwell_state {
public:
    //
    // Note: will add functions to maxwell_state through derived classes
    //
    gf3d e_p, a_p;
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in maxwell_state
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    maxwell_state(Grid* grid_i, dumper* dump_i, const char* name_i);

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
    void add(double factor, maxwell_state* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(maxwell_state* maxwell_state1, double factor, maxwell_state* maxwell_state2);

    //===============================================
    // addition
    //===============================================
    void equals(maxwell_state* rhs);

    //===============================================
    // apply characteristic outer boundaries
    //===============================================
    void char_OB(maxwell_state* last, double dt);

    //===============================================
    // print list of all functions in maxwell_state
    //===============================================
    void function_names();

    //===============================================
    // return name of maxwell_state
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
    ~maxwell_state();
    
};


#endif  /* MAXWELL_STATE_H */
