// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary functions for maxwell
//
//================================================
#ifndef MAXWELL_AUX_H
#define MAXWELL_AUX_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class maxwell_aux {
public:
    //
    // Note: will add functions to maxwell_aux through derived classes
    //
    gf3d A_xi, A2;
    //
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in maxwell_aux
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    maxwell_aux(Grid* grid_i, dumper* dump_i, const char* name_i);

    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts();

    //===============================================
  // addition
  //===============================================
    void add(double factor, maxwell_aux* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(maxwell_aux* maxwell_aux1, double factor, maxwell_aux* maxwell_aux2);

    //===============================================
    // addition
    //===============================================
    void equals(maxwell_aux* rhs);

    //===============================================
    // print list of all functions in maxwell_aux
    //===============================================
    void function_names();

    //===============================================
    // return name of maxwell_aux
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
    ~maxwell_aux();
    
};


#endif  /* MAXWELL_AUX_H */
