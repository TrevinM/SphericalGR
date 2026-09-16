// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of system, i.e. collection
// of grid functions for dynamical variables
//
// Reference: Baumgarte, Montero, Cordero-Carrion & Mueller, PRD 87, 044026 (2013)  (BMCM)
//
//================================================
#ifndef STATE_H
#define STATE_H

#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class StateBase {
public:
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in state
protected:
    gf3d** dump_list;  // list of all grid functions to be dumped
    const char* name;

    //===============================================
    // returns the prefix to be used in cout 
    //===============================================
    virtual const char* CoutName() = 0;

public:
    //===============================================
    // Constructor
    //===============================================
    StateBase(const char* name_i);

    //===============================================
    // return name of state
    //===============================================
    const char* Name();

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
    void add(double factor, StateBase* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(StateBase* state1, double factor, StateBase* state2);

    //===============================================
    // equals
    //===============================================
    void equals(StateBase* rhs);

    //===============================================
    // equals
    //===============================================
    void equals(double number);

    //===============================================
    // check whether states are equal
    //===============================================
    bool IsEqualTo(StateBase* s);

    //===============================================
    // apply characteristic outer boundaries
    //===============================================
    void char_OB(StateBase* last, double dt);

    //===============================================
    // check whether functions are finite
    //===============================================
    bool FINITE();

    //===============================================
    // print list of all functions in state
    //===============================================
    void function_names();

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
    virtual ~StateBase();
};


class state : public StateBase {
public:
    // metric components (rescaled as in (BMCM.20))
    gf3d phi, h_rr, h_rt, h_rp, h_tt, h_tp, h_pp;
    // extrinsic curvature components (see (BMCM.21))
    gf3d K, a_rr, a_rt, a_rp, a_tt, a_tp, a_pp;
    // connection functions (see (BMCM.22))
    gf3d lam_r, lam_t, lam_p;
    // Theta (for Z4)
    gf3d Theta;
    // lapse and shift (rescaled!), indices upstairs
    gf3d lapse, shift_r, shift_t, shift_p;
    gf3d B_r, B_t, B_p;
    // 
    Grid* grid;
public:
    //===============================================
    // Constructor
    //===============================================
    state(Grid* grid_i, dumper* dump_i, const char* name_i);
    ~state();

private:
    const char* CoutName() override;

};

#endif  /* STATE_H */
