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
    radhydro_aux(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 64;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        //
        int gf_counter = 0;
        //
        fct_list[gf_counter] = rho_0.setup(grid, 1, "rho_0", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = p.setup(grid, 1, "p", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = eps.setup(grid, 1, "eps", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = v_r.setup(grid, 1, "v_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = v_t.setup(grid, 1, "v_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = v_p.setup(grid, 1, "v_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = E.setup(grid, 1, "E", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f.setup(grid, 1, "f", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_r.setup(grid, 1, "f_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_t.setup(grid, 1, "f_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_p.setup(grid, 1, "f_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = FF.setup(grid, 1, "FF", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = udotF.setup(grid, 1, "udotF", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // set up gamma factor
        // 
        fct_list[gf_counter] = W.setup(grid, 1, "W", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // set up interface functions
        // NOTE: store functions at interfaces, so we should never use fill
        // ghosts of these functions with fill_ghosts.  Will set up parities
        // anyway...
        // 
        fct_list[gf_counter] = rho_L.setup(grid, 1, "rho_L", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = rho_R.setup(grid, 1, "rho_R", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = eps_L.setup(grid, 1, "eps_L", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = eps_R.setup(grid, 1, "eps_R", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = v_r_L.setup(grid, 1, "v_r_L", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = v_r_R.setup(grid, 1, "v_r_R", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = v_t_L.setup(grid, 1, "v_t_L", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = v_t_R.setup(grid, 1, "v_t_R", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = v_p_L.setup(grid, 1, "v_p_L", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = v_p_R.setup(grid, 1, "v_p_R", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = E_L.setup(grid, 1, "E_L", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = E_R.setup(grid, 1, "E_R", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_L.setup(grid, 1, "f_L", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_R.setup(grid, 1, "f_R", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_r_L.setup(grid, 1, "f_r_L", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_r_R.setup(grid, 1, "f_r_R", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_t_L.setup(grid, 1, "f_t_L", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_t_R.setup(grid, 1, "f_t_R", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_p_L.setup(grid, 1, "f_p_L", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_p_R.setup(grid, 1, "f_p_R", gf_counter, -1, -1, +1);
        gf_counter++;
        //
        // set up fluxes
        //
        // IMPORTANT NOTE:
        //
        // fluxes are defined at cell interfaces, not cell centers.
        // We will store fluxes at interface i - 1/2 into index i
        //
        // E.g. for r-direction, assuming N_g = 3 ghost zones:
        //
        //                      r=0
        //     |  x  |  x  |  x  |  x | 
        // i=     0     1     2     3      for normal grid functions
        // i=        1     2     3    4    for fluxes
        //
        //
        fct_list[gf_counter] = f_D_r.setup(grid, 1, "f_D_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_D_t.setup(grid, 1, "f_D_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_D_p.setup(grid, 1, "f_D_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_r_r.setup(grid, 1, "f_S_r_r", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_r_t.setup(grid, 1, "f_S_r_t", gf_counter, -1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_r_p.setup(grid, 1, "f_S_r_p", gf_counter, +1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_t_r.setup(grid, 1, "f_S_t_r", gf_counter, -1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_t_t.setup(grid, 1, "f_S_t_t", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_t_p.setup(grid, 1, "f_S_t_p", gf_counter, -1, +1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_p_r.setup(grid, 1, "f_S_p_r", gf_counter, +1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_p_t.setup(grid, 1, "f_S_p_t", gf_counter, -1, +1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_p_p.setup(grid, 1, "f_S_p_p", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_tau_r.setup(grid, 1, "f_tau_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_tau_t.setup(grid, 1, "f_tau_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_tau_p.setup(grid, 1, "f_tau_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_r_r.setup(grid, 1, "f_S_rad_r_r", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_r_t.setup(grid, 1, "f_S_rad_r_t", gf_counter, -1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_r_p.setup(grid, 1, "f_S_rad_r_p", gf_counter, +1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_t_r.setup(grid, 1, "f_S_rad_t_r", gf_counter, -1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_t_t.setup(grid, 1, "f_S_rad_t_t", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_t_p.setup(grid, 1, "f_S_rad_t_p", gf_counter, -1, +1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_p_r.setup(grid, 1, "f_S_rad_p_r", gf_counter, +1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_p_t.setup(grid, 1, "f_S_rad_p_t", gf_counter, -1, +1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_S_rad_p_p.setup(grid, 1, "f_S_rad_p_p", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_tau_rad_r.setup(grid, 1, "f_tau_rad_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = f_tau_rad_t.setup(grid, 1, "f_tau_rad_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = f_tau_rad_p.setup(grid, 1, "f_tau_rad_p", gf_counter, -1, -1, +1);
        gf_counter++;
        //
        // accretion rate and sound speed
        //
        fct_list[gf_counter] = flux.setup(grid, 1, "flux", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = sound_speed.setup(grid, 1, "sound_speed", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = temperature.setup(grid, 1, "temperature", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN RADHYDRO_AUX!!! " << endl;
    };
    //===============================================
    // fill ghosts
    //===============================================
    void fill_ghosts() {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].fill_ghosts();
    };
    //===============================================
  // addition
  //===============================================
    void add(double factor, radhydro_aux* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void add(radhydro_aux* radhydro_aux1, double factor, radhydro_aux* radhydro_aux2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(radhydro_aux1->fct_list[i], factor, radhydro_aux2->fct_list[i]);
    };
    //===============================================
    // addition
    //===============================================
    void equals(radhydro_aux* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // print list of all functions in radhydro_aux
    //===============================================
    void function_names() {
        cout << " RADHYDRO_AUX: List of all functions in radhydro_aux " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of radhydro_aux
    //===============================================
    const char* Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
        const char* suffix = "") {
        if (dump->time_to_dump(timestep)) {
            cout << " RADHYDRO_AUX: Dumping functions in radhydro_aux "
                << name << " at time t = " << time << endl;
            for (int i = 0; i < N_dump; i++) {
                dump->dump(time, prop_time, timestep, dump_list[i], suffix);
                dump->slice(time, prop_time, timestep, dump_list[i], suffix);
            }
        }
    };
    //===============================================
    // Find all functions to be dumped
    //===============================================
    int assemble_dump_list(const char* dump_list_file) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << " RADHYDRO_AUX: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " RADHYDRO_AUX: reading dump list from file " << dump_list_file << endl;
        N_dump = 0;
        char fct_name[64];
        infile >> fct_name;
        while (!infile.eof()) {
            for (int i = 0; i < N_fcts; i++) {
                if (!strcmp(fct_name, fct_list[i]->Name())) {
                    cout << " ... found function name " << fct_name << endl;
                    dump_list[N_dump] = fct_list[i]->Address();
                    N_dump++;
                }
            }
            infile >> fct_name;
        }
        return N_dump;
    };
    //===============================================
    // Destructor
    //===============================================
    ~radhydro_aux() {
        delete fct_list;
        delete dump_list;
        cout << " RADHYDRO_AUX: ... closing radhydro_aux " << name << "... " << endl;
    };
};


#endif  /* RADHYDRO_AUX_H */
