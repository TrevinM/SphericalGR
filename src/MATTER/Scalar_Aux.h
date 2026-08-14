// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary functions for scalar
//
//================================================
#ifndef SCALAR_AUX_H
#define SCALAR_AUX_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class scalar_aux {
public:
    //
    // Note: will add functions to scalar_aux through derived classes
    //
    gf3d Omega;
    //
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in scalar_aux
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    scalar_aux(Grid* grid_i, dumper* dump_i, const char* name_i);

    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts();

    //===============================================
  // addition
  //===============================================
    void add(double factor, scalar_aux* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(scalar_aux* scalar_aux1, double factor, scalar_aux* scalar_aux2);

    //===============================================
    // addition
    //===============================================
    void equals(scalar_aux* rhs);

    //===============================================
    // print list of all functions in scalar_aux
    //===============================================
    void function_names();

    //===============================================
    // return name of scalar_aux
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
    // Destructor
    //===============================================
    ~scalar_aux();
    
};


#endif  /* SCALAR_AUX_H */
