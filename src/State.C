    #include "State.h"
    #include "Container.h"
    
    //===============================================
    // Constructor
    //===============================================
    state::state(Grid* grid_i, dumper* dump_i, const char* name_i) : StateBase(name_i), grid(grid_i) {
        N_fcts = 25;
        fct_list = new gf3d * [N_fcts];
        dump_list = new gf3d * [N_fcts];   // allow for N_fcts, but restrict loops to N_dump...
        N_dump = 0;                      // set in assemble_dump_list
        int gf_counter = 0;
        //
        // metric components, g_{ij} = exp(4 \phi) (1 + h_{ij}), RESCALED
        //
        fct_list[gf_counter] = phi.setup(grid, 1, "phi", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = h_rr.setup(grid, 1, "h_rr", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = h_rt.setup(grid, 1, "h_rt", gf_counter, -1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = h_rp.setup(grid, 1, "h_rp", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = h_tt.setup(grid, 1, "h_tt", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = h_tp.setup(grid, 1, "h_tp", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = h_pp.setup(grid, 1, "h_pp", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // extrinsic curvature components, A_{ij} = exp(4 \phi) a_{ij}, RESCALED
        //
        double K_alpha_speed = 1.0;
        fct_list[gf_counter] = K.setup(grid, 2, "K", gf_counter, +1, +1, +1, K_alpha_speed);
        gf_counter++;
        fct_list[gf_counter] = a_rr.setup(grid, 1, "a_rr", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = a_rt.setup(grid, 1, "a_rt", gf_counter, -1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = a_rp.setup(grid, 1, "a_rp", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = a_tt.setup(grid, 1, "a_tt", gf_counter, +1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = a_tp.setup(grid, 1, "a_tp", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = a_pp.setup(grid, 1, "a_pp", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // connection functions, \Lambda^i, RESCALED
        //
        fct_list[gf_counter] = lam_r.setup(grid, 2, "lam_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = lam_t.setup(grid, 2, "lam_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = lam_p.setup(grid, 2, "lam_p", gf_counter, -1, -1, -1);
        gf_counter++;
        //
        // Theta function for Z4
        //
        fct_list[gf_counter] = Theta.setup(grid, 2, "Theta", gf_counter, +1, +1, +1);
        gf_counter++;
        //
        // lapse and shift
        // 
        double background = 1.0;
        fct_list[gf_counter] = lapse.setup(grid, 1, "lapse", gf_counter, +1, +1, +1,
            K_alpha_speed, background);
        gf_counter++;
        fct_list[gf_counter] = shift_r.setup(grid, 2, "shift_r", gf_counter, -1, +1, +1);
        gf_counter++;
        fct_list[gf_counter] = shift_t.setup(grid, 2, "shift_t", gf_counter, +1, -1, -1);
        gf_counter++;
        fct_list[gf_counter] = shift_p.setup(grid, 2, "shift_p", gf_counter, -1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = B_r.setup(grid, 2, "B_r", gf_counter, -1, +1, -1);
        gf_counter++;
        fct_list[gf_counter] = B_t.setup(grid, 2, "B_t", gf_counter, +1, -1, +1);
        gf_counter++;
        fct_list[gf_counter] = B_p.setup(grid, 2, "B_p", gf_counter, -1, -1, +1);
        gf_counter++;
        //
        // sanity check
        //
        if (gf_counter != N_fcts) cerr << " WRONG FUNCTION COUNT IN STATE!!! " << endl;
    };

    const char* state::CoutName() {
        return " STATE";
    }

    state::~state() {
        std::cout << CoutName() << "... closing state " << name << " ..." << endl; 
    }


    //===============================================
    // Constructor
    //===============================================
    StateBase::StateBase(const char* name_i) : name(name_i) {};

    //===============================================
    // return name of state
    //===============================================
    const char* StateBase::Name() { return name; };

    //===============================================
    // check properties
    //===============================================
    void StateBase::check_props() {
        for (int i = 0; i < N_fcts; i++) {
            (*fct_list)[i].constants();
        }
    }
    //===============================================
    // fill ghosts
    //===============================================
    void StateBase::fill_ghosts() {
        for (int i = 0; i < N_fcts; i++)
            (*fct_list)[i].fill_ghosts();
    };
    //===============================================
    // addition
    //===============================================
    void StateBase::add(double factor, StateBase* rhs) {
        if (typeid(*this) == typeid(*rhs))
            for (int i = 0; i < N_fcts; i++)
                (*fct_list)[i].add(factor, rhs->fct_list[i]);
        else {
            cout << CoutName() << ": Couldn't add " << typeid(*rhs).name() << " to " 
                << typeid(*this).name() << " " << name << endl;
            exit(1);
        }
    };
    //===============================================
    // addition 
    //===============================================
    void StateBase::add(StateBase* state1, double factor, StateBase* state2) {
        if (typeid(*this) == typeid(*state1) && typeid(*this) == typeid(*state2))
            for (int i = 0; i < N_fcts; i++)
                (*fct_list)[i].add(state1->fct_list[i], factor, state2->fct_list[i]);
        else {
            cout << CoutName() << ": Couldn't add " << typeid(*state1).name() << " and " << typeid(*state2).name()
                << " to " << typeid(*this).name() << " " << name << endl;
            exit(1);
        }
    };

    //===============================================
    // equals
    //===============================================
    void StateBase::equals(StateBase* rhs) {
        if (typeid(*this) == typeid(*rhs))
            for (int i = 0; i < N_fcts; i++) {
                fct_list[i]->equals(rhs->fct_list[i]);
            }
        else {
            cout << CoutName() << ": Couldn't equate " << typeid(*rhs).name() << " to " 
                << typeid(*this).name() << " " << name << endl;
            exit(1);
        }

    };

    //===============================================
    // equals
    //===============================================
    void StateBase::equals(double number) {
        for (int i = 0; i < N_fcts; i++) {
            fct_list[i]->equals(number);
        }
    };

    //===============================================
    // check whether states are equal
    //===============================================
    bool StateBase::IsEqualTo(StateBase* s) {
        if (typeid(*this) == typeid(s)) {
            bool fct_equal = true;
            bool state_equal = true;
            for (int i = 0; i < N_fcts; i++) {
                fct_equal = fct_list[i]->IsEqualTo(s->fct_list[i]);
                if (!fct_equal) {
                    cout << CoutName() << ": Functions " << fct_list[i]->Name() << " in states "
                        << Name() << " and " << s->Name() << " are different!" << endl;
                    state_equal = false;
                }
            }
            if (state_equal == true)
                cout << CoutName() << ": States " << Name() << " and " << s->Name()
                << " are equal! " << endl;
            return state_equal;
        } else {
            cout << CoutName() << ": Couldn't comapre " << typeid(*s).name() << " to " 
                << typeid(*this).name() << " " << name << endl;
            exit(1);
        }
    };
    //===============================================
    // apply characteristic outer boundaries
    //===============================================
    void StateBase::char_OB(StateBase* last, double dt) {
        for (int i = 0; i < N_fcts; i++)
            fct_list[i]->fill_outerboundary(*last->fct_list[i], dt);
    }

    //===============================================
    // check whether functions are finite
    //===============================================
    bool StateBase::FINITE() {
        bool fine = true;
        for (int i = 0; i < N_fcts; i++) {
            if (!fct_list[i]->FINITE()) {
                cout << CoutName() << ": Function " << fct_list[i]->Name() << " in state "
                    << name << " is not finite!" << endl;
                fine = false;
            }
        }
        return fine;
    };

    //===============================================
    // print list of all functions in state
    //===============================================
    void StateBase::function_names() {
        cout << CoutName() << ": List of all functions in state " << name << ": " << endl;
        for (int i = 0; i < N_fcts; i++)
            cout << "      " << fct_list[i]->Name() << endl;
    };

    //===============================================
    // Dump grid functions
    //===============================================
    void StateBase::dump_fcts(double time, double prop_time, int timestep,
        const char* suffix) {
        if (Container::dump->time_to_dump(timestep) || strcmp(suffix, "")) {
            cout << CoutName() << ": Dumping functions in state " << name
                << " at time t = " << time << endl;
            for (int i = 0; i < N_dump; i++) {
                Container::dump->dump(time, prop_time, timestep, dump_list[i], suffix);
                Container::dump->slice(time, prop_time, timestep, dump_list[i], suffix);
            }
        }
    };
    
    //===============================================
    // Find all functions to be dumped
    //===============================================
    int StateBase::assemble_dump_list(const char* dump_list_file) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << CoutName() << ": can't open " << dump_list_file << " for input." << endl;
            return 0;
        }
        cout << CoutName() << ": reading dump list from file " << dump_list_file << endl;
        N_dump = 0;
        char fct_name[64];
        infile >> fct_name;
        while (!infile.eof()) {
            for (int i = 0; i < N_fcts; i++) {
                if (!strcmp(fct_name, fct_list[i]->Name())) {
                    cout << " ... found function name " << fct_name << endl;
                    dump_list[N_dump] = fct_list[i]->Address();
                    //	  cout << " function " << fct_name << " : " <<  dump_list[N_dump] << "  " << fct_list[i] << "  " << fct_list[i]->Address() << endl;

                    N_dump++;
                }
            }
            infile >> fct_name;
        }
        if (N_dump > N_fcts) {
            cerr << CoutName() << ": found too many grid functions in assemble_dump_list() " << endl;
            N_dump = N_fcts;
        }
        return N_dump;
    };
    //===============================================
    // Regrid
    //===============================================
    void StateBase::Regrid(VecDoub r_new) {
        for (int i = 0; i < N_fcts; i++)
            fct_list[i]->Regrid(r_new);
    }
    //===============================================
    // Destructor
    //===============================================
    StateBase::~StateBase() {
        delete fct_list;
        delete dump_list;
    };
