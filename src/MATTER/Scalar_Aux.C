#include "Scalar_Aux.h"
    
    //===============================================
    // Constructor
    //===============================================
    scalar_aux::scalar_aux(Grid* grid_i, dumper* dump_i, const char* name_i) :
        grid(grid_i), dump(dump_i), name(name_i) {
        N_fcts = 1;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        //
        int gf_counter = 0;
        //
        fct_list[gf_counter] = Omega.setup(grid, 1, "Omega", gf_counter, -1, -1, 1);
        gf_counter++;
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN SCALAR_AUX!!! " << endl;
    };
    //===============================================
    // fill ghosts
    //===============================================
    void scalar_aux::fill_ghosts() {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].fill_ghosts();
    };
    //===============================================
  // addition
  //===============================================
    void scalar_aux::add(double factor, scalar_aux* rhs) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(factor, rhs->fct_list[i]);
    };
    //===============================================
    // addition 
    //===============================================
    void scalar_aux::add(scalar_aux* scalar_aux1, double factor, scalar_aux* scalar_aux2) {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].add(scalar_aux1->fct_list[i], factor, scalar_aux2->fct_list[i]);
    };
    //===============================================
    // addition
    //===============================================
    void scalar_aux::equals(scalar_aux* rhs) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(rhs->fct_list[i]);
        }
    };
    //===============================================
    // print list of all functions in scalar_aux
    //===============================================
    void scalar_aux::function_names() {
        cout << " SCALAR_AUX: List of all functions in scalar_aux " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };
    //===============================================
    // return name of scalar_aux
    //===============================================
    const char* scalar_aux::Name() { return name; };
    //===============================================
    // Dump grid functions
    //===============================================
    void scalar_aux::dump_fcts(double time, double prop_time, int timestep,
        const char* suffix) {
        if (dump->time_to_dump(timestep)) {
            cout << " SCALAR_AUX: Dumping functions in scalar_aux "
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
    int scalar_aux::assemble_dump_list(const char* dump_list_file) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << " SCALAR_AUX: can't oppen " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << " SCALAR_AUX: reading dump list from file " << dump_list_file << endl;
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
    scalar_aux::~scalar_aux() {
        delete fct_list;
        delete dump_list;
        cout << " SCALAR_AUX: ... closing scalar_aux " << name << "... " << endl;
    };
