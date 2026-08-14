#include "RadHydro_State.h"

    //===============================================
    // Constructor
    //===============================================
    radhydro_state::radhydro_state(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 9;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        //
        int gf_counter = 0;
        //
        fct_list[gf_counter] = D.setup(grid, 1, "D", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = tau.setup(grid, 1, "tau", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = S_r.setup(grid, 1, "S_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = S_t.setup(grid, 1, "S_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = S_p.setup(grid, 1, "S_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = tau_rad.setup(grid, 1, "tau_rad", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = S_rad_r.setup(grid, 1, "S_rad_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = S_rad_t.setup(grid, 1, "S_rad_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = S_rad_p.setup(grid, 1, "S_rad_p", gf_counter, -1, -1, +1);
        gf_counter++;
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN RADHYDRO_STATE!!! " << endl;
    };
    //===============================================
    // check properties
    //===============================================
    void radhydro_state::check_props() {
        for (int i = 0; i < N_fcts; i++) {
            (*fct_list)[i].constants();
        }
    };
    //===============================================
    // fill ghosts
    //===============================================
    void radhydro_state::fill_ghosts() {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].fill_ghosts();
    };
    //===============================================
    // addition
    //===============================================
    void radhydro_state::add(double factor, radhydro_state* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void radhydro_state::add(radhydro_state* radhydro_state1, double factor, radhydro_state* radhydro_state2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(radhydro_state1->fct_list[i], factor, radhydro_state2->fct_list[i]);
    };
    //===============================================
    // equals
    //===============================================
    void radhydro_state::equals(radhydro_state* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    void radhydro_state::equals(double number) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(number);
        }
    };
    //===============================================
    // check whether states are equal
    //===============================================
    bool radhydro_state::IsEqualTo(radhydro_state* s) {
        bool fct_equal = true;
        bool state_equal = true;
        for (int i = 0; i < N_fcts; i++) {
            fct_equal = fct_list[i]->IsEqualTo(s->fct_list[i]);
            if (!fct_equal) {
                cout << " RADHYDRO: Functions " << fct_list[i]->Name()
                    << " in radhydro_states "
                    << Name() << " and " << s->Name() << " are different!" << endl;
                state_equal = false;
            }
        }
        if (state_equal == true)
            cout << " RADHYDRO: radhydro_states " << Name() << " and " << s->Name()
            << " are equal! " << endl;
        return state_equal;
    };
    //===============================================
    // print list of all functions in radhydro_state
    //===============================================
    void radhydro_state::function_names() {
        cout << " RADHYDRO_STATE: List of all functions in radhydro_state " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of radhydro_state
    //===============================================
    const char* radhydro_state::Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void radhydro_state::dump_fcts(double time, double prop_time, int timestep,
        const char* suffix) {
        if (dump->time_to_dump(timestep) || strcmp(suffix, "")) {
            cout << " RADHYDRO_STATE: Dumping functions in radhydro_state "
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
    int radhydro_state::assemble_dump_list(const char* dump_list_file) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << " RADHYDRO_STATE: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " RADHYDRO_STATE: reading dump list from file " << dump_list_file << endl;
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
        if (N_dump > N_fcts) {
            cerr << " RADHYDRO: found too many grid functions in assemble_dump_list() " << endl;
            N_dump = N_fcts;
        }
        return N_dump;
    };
    //===============================================
    // Regrid
    //===============================================
    void radhydro_state::Regrid(VecDoub r_new) {
        for (int i = 0; i < N_fcts; i++)
            fct_list[i]->Regrid(r_new);
    }
    //===============================================
    // Destructor
    //===============================================
    radhydro_state::~radhydro_state() {
        delete fct_list;
        delete dump_list;
        cout << " RADHYDRO_STATE: ... closing radhydro_state " << name << "... " << endl;
    };
