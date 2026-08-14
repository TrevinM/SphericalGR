// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary functions for dualmaxwell
// Definitions in DualMaxwell_State.C
//
//================================================
#ifndef DUALMAXWELL_AUX_H
#define DUALMAXWELL_AUX_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class dualmaxwell_aux {
public:
    //
    // Note: will add functions to dualmaxwell_aux through derived classes
    //
    gf3d A_xi, A2;
    gf3d E_r, E_t, E_p, B_r, B_t, B_p;  // all rescaled, indices up 
    //
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in dualmaxwell_aux
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    dualmaxwell_aux(Grid* grid_i, dumper* dump_i, const char* name_i);
    
    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts();

    //===============================================
  // addition
  //===============================================
    void add(double factor, dualmaxwell_aux* rhs);
    
    //===============================================
    // addition 
    //===============================================
    void add(dualmaxwell_aux* dualmaxwell_aux1, double factor, dualmaxwell_aux* dualmaxwell_aux2);

    //===============================================
    // addition
    //===============================================
    void equals(dualmaxwell_aux* rhs);

    //===============================================
    // print list of all functions in dualmaxwell_aux
    //===============================================
    void function_names();

    //===============================================
    // return name of dualmaxwell_aux
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
    ~dualmaxwell_aux();
    
};


#endif  /* DUALMAXWELL_AUX_H */
