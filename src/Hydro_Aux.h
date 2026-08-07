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
    hydro_aux(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 37;
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
        //
        // accretion rate and sound speed
        //
        fct_list[gf_counter] = flux.setup(grid, 1, "flux", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = sound_speed.setup(grid, 1, "sound_speed", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = temperature.setup(grid, 1, "temperature", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = ut_low.setup(grid, 1, "ut_low", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = hydro_errors.setup(grid, 1, "hydro_errors", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN HYDRO_AUX!!! " << endl;
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
    void add(double factor, hydro_aux* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void add(hydro_aux* hydro_aux1, double factor, hydro_aux* hydro_aux2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(hydro_aux1->fct_list[i], factor, hydro_aux2->fct_list[i]);
    };
    //===============================================
    // addition
    //===============================================
    void equals(hydro_aux* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // print list of all functions in hydro_aux
    //===============================================
    void function_names() {
        cout << " HYDRO_AUX: List of all functions in hydro_aux " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of hydro_aux
    //===============================================
    const char* Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
        const char* suffix = "") {
        if (dump->time_to_dump(timestep)) {
            cout << " HYDRO_AUX: Dumping functions in hydro_aux "
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
            cerr << " HYDRO_AUX: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " HYDRO_AUX: reading dump list from file " << dump_list_file << endl;
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
    ~hydro_aux() {
        delete fct_list;
        delete dump_list;
        cout << " HYDRO_AUX: ... closing hydro_aux " << name << "... " << endl;
    };
};


#endif  /* HYDRO_AUX_H */
