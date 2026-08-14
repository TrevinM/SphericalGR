// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary functions for hydro
//
//================================================
#ifndef HYDRO_AUX_H
#define HYDRO_AUX_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class hydro_aux {
public:
    //
    // primitive hydro variables
    // 
    gf3d rho_0, p, eps;
    gf3d v_r, v_t, v_p;    // fluid three-velocity, indices *up*-stairs
    //
   // gamma factor
   //
    gf3d W;
    //
    // function values at interfaces
    //
    gf3d rho_L, rho_R;
    gf3d eps_L, eps_R;
    gf3d v_r_L, v_r_R, v_t_L, v_t_R, v_p_L, v_p_R;
    //
    // fluxes
    // 
    gf3d f_D_r, f_D_t, f_D_p;
    gf3d f_S_r_r, f_S_r_t, f_S_r_p;
    gf3d f_S_t_r, f_S_t_t, f_S_t_p;
    gf3d f_S_p_r, f_S_p_t, f_S_p_p;
    gf3d f_tau_r, f_tau_t, f_tau_p;
    //
    // flux through spheres of constant coordinate radius, sound speed,
    // and temperature ( = k_B T / (m_B c^2) )
    //
    gf3d flux, sound_speed, temperature, ut_low;
    gf3d hydro_errors;
    //
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in hydro_aux
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    hydro_aux(Grid* grid_i, dumper* dump_i, const char* name_i);

    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts();

    //===============================================
  // addition
  //===============================================
    void add(double factor, hydro_aux* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(hydro_aux* hydro_aux1, double factor, hydro_aux* hydro_aux2);

    //===============================================
    // addition
    //===============================================
    void equals(hydro_aux* rhs);

    //===============================================
    // print list of all functions in hydro_aux
    //===============================================
    void function_names();

    //===============================================
    // return name of hydro_aux
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
    ~hydro_aux();

};


#endif  /* HYDRO_AUX_H */
