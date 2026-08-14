#include "Maxwell_Aux.h"
    
    //===============================================
    // Constructor
    //===============================================
    maxwell_aux::maxwell_aux(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 2;
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
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN MAXWELL_AUX!!! " << endl;
    };
    //===============================================
    // fill ghosts
    //===============================================
    void maxwell_aux::fill_ghosts() {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].fill_ghosts();
    };
    //===============================================
  // addition
  //===============================================
    void maxwell_aux::add(double factor, maxwell_aux* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void maxwell_aux::add(maxwell_aux* maxwell_aux1, double factor, maxwell_aux* maxwell_aux2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(maxwell_aux1->fct_list[i], factor, maxwell_aux2->fct_list[i]);
    };
    //===============================================
    // addition
    //===============================================
    void maxwell_aux::equals(maxwell_aux* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // print list of all functions in maxwell_aux
    //===============================================
    void maxwell_aux::function_names() {
        cout << " MAXWELL_AUX: List of all functions in maxwell_aux " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of maxwell_aux
    //===============================================
    const char* maxwell_aux::Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void maxwell_aux::dump_fcts(double time, double prop_time, int timestep,
        const char* suffix) {
        if (dump->time_to_dump(timestep)) {
            cout << " MAXWELL_AUX: Dumping functions in maxwell_aux "
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
    int maxwell_aux::assemble_dump_list(const char* dump_list_file) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << " MAXWELL_AUX: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " MAXWELL_AUX: reading dump list from file " << dump_list_file << endl;
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
    maxwell_aux::~maxwell_aux() {
        delete fct_list;
        delete dump_list;
        cout << " MAXWELL_AUX: ... closing maxwell_aux " << name << "... " << endl;
    };