// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing entire state of matter, i.e. collection
// of grid functions for dynamical matter variables
//
//================================================
#ifndef MATTER_STATE_H
#define MATTER_STATE_H


#include "gridfunction.h"
#include "Grid.h"
#include "dumper.h"

class matter_state {
public:
    //
    // Note: will add functions to matter_state through derived classes
    //
    int N_fcts, N_dump;
    gf3d** fct_list;   // list of all grid functions in matter_state
    gf3d** dump_list;  // list of all grid functions to be dumped
    Grid* grid;
    dumper* dump;
    const char* name;
public:
    //===============================================
    // Constructor
    //===============================================
    matter_state(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        // doesn't do much...
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
    void add(double factor, matter_state* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void add(matter_state* matter_state1, double factor, matter_state* matter_state2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(matter_state1->fct_list[i], factor, matter_state2->fct_list[i]);
    };
    //===============================================
    // addition
    //===============================================
    void equals(matter_state* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // print list of all functions in matter_state
    //===============================================
    void function_names() {
        cout << " MATTER_STATE: List of all functions in matter_state " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of matter_state
    //===============================================
    const char* Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0) {
        if (dump->time_to_dump(timestep)) {
            cout << " MATTER_STATE: Dumping functions in matter_state " << name << " at time t = " << time << endl;
            for (int i = 0; i < N_dump; i++)
                dump->dump(time, prop_time, timestep, dump_list[i]);
        }
    };
    //===============================================
    // Find all functions to be dumped
    //===============================================
    int assemble_dump_list(const char* dump_list_file) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << " MATTER_STATE: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " MATTER_STATE: reading dump list from file " << dump_list_file << endl;
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
    ~matter_state() {
        delete fct_list;
        delete dump_list;
        cout << " MATTER_STATE: ... closing matter_state " << name << "... " << endl;
    };
};

// class maxwell_state : public matter_state {
// public:
//     gf3d e_p, a_p;    // lower phi-components of fields, rescaled
//     //===============================================
//     // Constructor
//     //===============================================
//     maxwell_state(Grid* grid_i, dumper* dump_i, const char* name_i) :
//         matter_state(grid_i, dump_i, name_i) {
//         cout << " MATTER: creating matter state for Maxwell... " << endl;
//         N_fcts = 2;
//         fct_list = new gf3d * [N_fcts];
//         dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
//         N_dump = 0;                      // set in assemble_dump_list
//         //
//         int gf_counter = 2;
//         //
//         fct_list[gf_counter] = e_p.setup(grid, 1, "e_p", gf_counter, -1, -1, 1);
//         gf_counter++;
//         fct_list[gf_counter] = a_p.setup(grid, 1, "a_p", gf_counter, -1, -1, 1);
//         gf_counter++;
//     }
//     ~maxwell_state() {
//     };
// };


#endif  /* MATTER_STATE_H */
