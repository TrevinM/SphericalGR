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
    maxwell_state(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 2;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        //
        int gf_counter = 0;
        //
        fct_list[gf_counter] = e_p.setup(grid, 1, "e_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = a_p.setup(grid, 1, "a_p", gf_counter, -1, -1, +1);
        gf_counter++;
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN MAXWELL_STATE!!! " << endl;
    };
    //===============================================
    // check properties
    //===============================================
    void check_props() {
        for (int i = 0; i < N_fcts; i++) {
            (*fct_list)[i].constants();
        }
        //      (*fct_list)[i].fill_ghosts();
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
    void add(double factor, maxwell_state* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void add(maxwell_state* maxwell_state1, double factor, maxwell_state* maxwell_state2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(maxwell_state1->fct_list[i], factor, maxwell_state2->fct_list[i]);
    };
    //===============================================
    // addition
    //===============================================
    void equals(maxwell_state* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // apply characteristic outer boundaries
    //===============================================
    void char_OB(maxwell_state* last, double dt) {
        for (int i = 0; i < N_fcts; i++)
            fct_list[i]->fill_outerboundary(*last->fct_list[i], dt);
    }
    //===============================================
    // print list of all functions in maxwell_state
    //===============================================
    void function_names() {
        cout << " MAXWELL_STATE: List of all functions in maxwell_state " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of maxwell_state
    //===============================================
    const char* Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
        const char* suffix = "") {
        if (dump->time_to_dump(timestep) || strcmp(suffix, "")) {
            cout << " MAXWELL_STATE: Dumping functions in maxwell_state "
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
            cerr << " MAXWELL_STATE: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " MAXWELL_STATE: reading dump list from file " << dump_list_file << endl;
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
    // Regrid
    //===============================================
    void Regrid(VecDoub r_new) {
        for (int i = 0; i < N_fcts; i++)
            fct_list[i]->Regrid(r_new);
    }
    //===============================================
    // Destructor
    //===============================================
    ~maxwell_state() {
        delete fct_list;
        delete dump_list;
        cout << " MAXWELL_STATE: ... closing maxwell_state " << name << "... " << endl;
    };
};


#endif  /* MAXWELL_STATE_H */
