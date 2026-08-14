#include "DualMaxwell_State.h"

//===============================================
// Constructor
//===============================================
dualmaxwell_state::dualmaxwell_state(Grid* grid_i, dumper* dump_i, const char* name_i) :
    grid(grid_i), dump(dump_i), name(name_i) {
    N_fcts = 6;
    fct_list = new gf3d * [N_fcts];
    dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
    N_dump = 0;                      // set in assemble_dump_list
    //
    int gf_counter = 0;
    //
    fct_list[gf_counter] = a_r.setup(grid, 1, "a_r", gf_counter, -1, +1, -1);
    gf_counter++;
    fct_list[gf_counter] = a_t.setup(grid, 1, "a_t", gf_counter, +1, -1, +1);
    gf_counter++;
    fct_list[gf_counter] = a_p.setup(grid, 1, "a_p", gf_counter, -1, -1, +1);
    gf_counter++;
    fct_list[gf_counter] = as_r.setup(grid, 1, "as_r", gf_counter, -1, +1, -1);
    gf_counter++;
    fct_list[gf_counter] = as_t.setup(grid, 1, "as_t", gf_counter, +1, -1, +1);
    gf_counter++;
    fct_list[gf_counter] = as_p.setup(grid, 1, "as_p", gf_counter, -1, -1, +1);
    gf_counter++;
    //
    // sanity check
    //
    if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN DUALMAXWELL_STATE!!! " << endl;
};

//===============================================
// check properties
//===============================================
void dualmaxwell_state::check_props() {
    for (int i = 0; i < N_fcts; i++) {
        (*fct_list)[i].constants();
    }
    //      (*fct_list)[i].fill_ghosts();
};
//===============================================
// fill ghosts
//===============================================
void dualmaxwell_state::fill_ghosts() {
    for (int i = 0; i < N_fcts; i++)
        (*fct_list)[i].fill_ghosts();
};
//===============================================
// addition
//===============================================
void dualmaxwell_state::add(double factor, dualmaxwell_state* rhs) {
    for (int i = 0; i < N_fcts; i++)
        (*fct_list)[i].add(factor, rhs->fct_list[i]);
};
//===============================================
// addition 
//===============================================
void dualmaxwell_state::add(dualmaxwell_state* dualmaxwell_state1, double factor, dualmaxwell_state* dualmaxwell_state2) {
    for (int i = 0; i < N_fcts; i++)
        (*fct_list)[i].add(dualmaxwell_state1->fct_list[i], factor, dualmaxwell_state2->fct_list[i]);
};
//===============================================
// addition
//===============================================
void dualmaxwell_state::equals(dualmaxwell_state* rhs) {
    for (int i = 0; i < N_fcts; i++) {
        fct_list[i]->equals(rhs->fct_list[i]);
    }
};
//===============================================
// apply characteristic outer boundaries
//===============================================
void dualmaxwell_state::char_OB(dualmaxwell_state* last, double dt) {
    for (int i = 0; i < N_fcts; i++)
        fct_list[i]->fill_outerboundary(*last->fct_list[i], dt);
}
//===============================================
// print list of all functions in dualmaxwell_state
//===============================================
void dualmaxwell_state::function_names() {
    cout << " DUALMAXWELL_STATE: List of all functions in dualdualmaxwell_state " << name << ": " << endl;
    for (int i = 0; i < N_fcts; i++)
        cout << "      " << fct_list[i]->Name() << endl;
};
//===============================================
// return name of dualmaxwell_state
//===============================================
const char* dualmaxwell_state::Name() { return name; };
//===============================================
// Dump grid functions
//===============================================
void dualmaxwell_state::dump_fcts(double time, double prop_time, int timestep,
    const char* suffix) {
    if (dump->time_to_dump(timestep) || strcmp(suffix, "")) {
        cout << " DUALMAXWELL_STATE: Dumping functions in dualmaxwell_state "
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
int dualmaxwell_state::assemble_dump_list(const char* dump_list_file) {
    ifstream infile;
    infile.open(dump_list_file);
    if (!infile) {
        cerr << " DUALMAXWELL_STATE: can't oppen " << dump_list_file << " for input." << endl;
        return 0;
    }
    cout << " DUALMAXWELL_STATE: reading dump list from file " << dump_list_file << endl;
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
void dualmaxwell_state::Regrid(VecDoub r_new) {
    for (int i = 0; i < N_fcts; i++)
        fct_list[i]->Regrid(r_new);
}
//===============================================
// Destructor
//===============================================
dualmaxwell_state::~dualmaxwell_state() {
    delete fct_list;
    delete dump_list;
    cout << " DUALMAXWELL_STATE: ... closing dualmaxwell_state " << name << "... " << endl;
};
