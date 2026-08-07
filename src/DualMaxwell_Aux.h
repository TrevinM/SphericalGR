// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing auxiliary functions for dualmaxwell
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
    dualmaxwell_aux(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 8;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        //
        int gf_counter = 0;
        //
        fct_list[gf_counter] = A_xi.setup(grid, 1, "A_xi", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = A2.setup(grid, 1, "A2", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = E_r.setup(grid, 1, "E_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = E_t.setup(grid, 1, "E_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = E_p.setup(grid, 1, "E_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = B_r.setup(grid, 1, "B_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = B_t.setup(grid, 1, "B_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = B_p.setup(grid, 1, "B_p", gf_counter, -1, -1, +1);
        gf_counter++;
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN DUALMAXWELL_AUX!!! " << endl;
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
    void add(double factor, dualmaxwell_aux* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void add(dualmaxwell_aux* dualmaxwell_aux1, double factor, dualmaxwell_aux* dualmaxwell_aux2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(dualmaxwell_aux1->fct_list[i], factor, dualmaxwell_aux2->fct_list[i]);
    };
    //===============================================
    // addition
    //===============================================
    void equals(dualmaxwell_aux* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // print list of all functions in dualmaxwell_aux
    //===============================================
    void function_names() {
        cout << " DUALMAXWELL_AUX: List of all functions in dualmaxwell_aux " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of dualmaxwell_aux
    //===============================================
    const char* Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void dump_fcts(double time = 0.0, double prop_time = 0.0, int timestep = 0,
        const char* suffix = "") {
        if (dump->time_to_dump(timestep)) {
            cout << " DUALMAXWELL_AUX: Dumping functions in dualmaxwell_aux "
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
            cerr << " DUALMAXWELL_AUX: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " DUALMAXWELL_AUX: reading dump list from file " << dump_list_file << endl;
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
    ~dualmaxwell_aux() {
        delete fct_list;
        delete dump_list;
        cout << " DUALMAXWELL_AUX: ... closing dualmaxwell_aux " << name << "... " << endl;
    };
};


#endif  /* DUALMAXWELL_AUX_H */
