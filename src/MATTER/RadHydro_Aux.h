// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary functions for radhydro
//
//================================================
#ifndef RADHYDRO_AUX_H
#define RADHYDRO_AUX_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class radhydro_aux {
public:
    //
    // primitive hydro variables
    // 
    gf3d rho_0, p, eps;
    gf3d v_r, v_t, v_p;    // fluid three-velocity, indices *up*-stairs
    // NOTE: uses Valencia convention, v^i = \gamma^i_a u^a / W
    //
    // primitive radiation variables
    // 
    gf3d E, f, f_r, f_t, f_p;   // indices up-stairs; f = - n_a F^a = \alpha F^0
    gf3d FF;   // magnitude squared of F^\alpha
    gf3d udotF;   // u_a F^a = W ( v_i f^i - f ): should vanish!  
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
    gf3d E_L, E_R, f_L, f_R;
    gf3d f_r_L, f_r_R, f_t_L, f_t_R, f_p_L, f_p_R;
    //
    // fluxes for fluid
    // 
    gf3d f_D_r, f_D_t, f_D_p;
    gf3d f_S_r_r, f_S_r_t, f_S_r_p;
    gf3d f_S_t_r, f_S_t_t, f_S_t_p;
    gf3d f_S_p_r, f_S_p_t, f_S_p_p;
    gf3d f_tau_r, f_tau_t, f_tau_p;
    //
    // fluxes for radiation
    // 
    gf3d f_S_rad_r_r, f_S_rad_r_t, f_S_rad_r_p;
    gf3d f_S_rad_t_r, f_S_rad_t_t, f_S_rad_t_p;
    gf3d f_S_rad_p_r, f_S_rad_p_t, f_S_rad_p_p;
    gf3d f_tau_rad_r, f_tau_rad_t, f_tau_rad_p;
    //
    // flux through spheres of constant coordinate radius, sound speed,
    // and temperature ( = k_B T / (m_B c^2) )
    //
    gf3d flux, sound_speed, temperature;
    //
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in radhydro_aux
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    radhydro_aux(Grid* grid_i, dumper* dump_i, const char* name_i);
    
    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts();

    //===============================================
  // addition
  //===============================================
    void add(double factor, radhydro_aux* rhs);

    //===============================================
    // addition 
    //===============================================
    void add(radhydro_aux* radhydro_aux1, double factor, radhydro_aux* radhydro_aux2);

    //===============================================
    // addition
    //===============================================
    void equals(radhydro_aux* rhs);

    //===============================================
    // print list of all functions in radhydro_aux
    //===============================================
    void function_names();

    //===============================================
    // return name of radhydro_aux
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
    ~radhydro_aux();

};


#endif  /* RADHYDRO_AUX_H */
