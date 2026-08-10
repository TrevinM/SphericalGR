//================================================================
//
// Code for scalar field in Spherical Coordinates
//
//===============================================================
#include "Container.h"
#include "Read_Input.h"
//================================================================
//
// Main code
//
//================================================================

int main(int argc, char* argv[]) {
    cout << "==================================================" << endl;
    cout << "=======       Running SphericalGR         ========" << endl;
    cout << "==================================================" << endl;
    int error = 0;
    //================================================================
    // Read Input
    //================================================================
    double t_max = 0.0;
    double eta_KO = 0.0;  // coefficient for Kreiss-Oliger
    int dump_step = 0;
    int space_type = 0; // parameter for asymptotic space type
    int indata_type = 0;
    char indata_input[64];
    int matter_type = 0;
    int eos_type = 0;
    int slicing_type = 0;
    int gauge_type = 0;
    int get_rid_of_output = 0;
    char monitor_filename[64];
    int note_step = 0;
    int chkpt_step = 0;
    int read_from_chkpt = 0;
    int sigma = 0;      // switch between Lagrangian and Eulerian formulation 
    int cowling = 0;    // 0: no cowling, 1: fix gravity, 2: fix matter
    bool solve_constraints, rescale_metric, track_photons;
    bool write_profiles, extract_waves;
    int z4 = 0;       // switch between BSSN and Z4
    double kappa_11, kappa_12, kappa_2, kappa_ric; // parameters for Z4
    int RK_order = 0;  // order of Runge-Kutta (3 or 4)
    int char_OB = 0;   // switch for implementation of Sommerfeld BCs
    error = Read_Input(argc, argv, get_rid_of_output,
        t_max, eta_KO, sigma, cowling,
        z4, kappa_11, kappa_12, kappa_2,
        kappa_ric, RK_order, char_OB,
        dump_step, monitor_filename, note_step,
        read_from_chkpt, chkpt_step,
        space_type, indata_type, indata_input,
        matter_type, eos_type, slicing_type, gauge_type,
        solve_constraints, rescale_metric, track_photons,
        write_profiles, extract_waves);
    if (error) {
        cerr << " Read_Input detected error!" << endl;
        return error;
    }
    //================================================================
    // get rid of old output
    //================================================================
    if (get_rid_of_output == 1) {
        system("rm -f *_rays_*\n");
        system("rm -f *_slice_*\n");
        system("rm -f *Photon*.mon\n");
        system("rm -f Particles_*_*_*\n");
    }
#ifdef AXISYMMETRY
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
    cout << "!! Running code assuming AXISYMMETRY !! " << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
#endif /* AXISYMMETRY */
#ifdef EQSYMMETRY
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
    cout << "!! Running code assuming EQSYMMETRY  !! " << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
#endif /* EQSYMMETRY */
    //================================================================
    // Allocate Grid Class
    //================================================================
    cout << " GR: allocating Container::grid..." << endl;
    Container::grid = new Grid();
    //================================================================
    // Allocate checkpoint class
    //================================================================
    Container::checkpoint = new CheckPoint(read_from_chkpt, chkpt_step, Container::grid);
    //================================================================
    // Allocate Cosmology Class
    //================================================================
    cout << " GR: allocating Container::cosmology..." << endl;
    if (space_type == 1) {
        Container::cosmology = new Minkowski();
    } else if (space_type == 2) {
        Container::cosmology = new DeSitter();
    } else if (space_type == 3) {
        Container::cosmology = new Radiation();
    } else {
        cerr << " GR: no Container::cosmology of type " << space_type << "!!" << endl;
        return 1;
    }
    //================================================================
    // Allocate Equation of State
    //================================================================
    cout << " GR: allocating EOS..." << endl;
    if (eos_type == 0) {
        Container::eos = NULL;
    } else if (eos_type == 1) {
        Container::eos = new polytrope((char*)"Polytrope_Input");
    } else if (eos_type == 2) {
        Container::eos = new gamma_law((char*)"Polytrope_Input");
    } else if (eos_type == 3) {
        Container::eos = new piece_polytrope((char*)"Piece_Polytrope_Input");
    } else if (eos_type == 4) {
        Container::eos = new ideal_gas();
    } else if (eos_type == 5) {
        Container::eos = new gas_radiation((char*)"GasRadiation_Input");
    } else {
        cout << " GR: Unknown EOS for eos_type = " << eos_type << endl;
    }
    //================================================================
    // Allocate InData Class
    //================================================================
    cout << " GR: allocating indata..." << endl;
    if (Container::checkpoint->ReadFromChkpt()) {
        Container::indata = new ReadFromCheckPoint(read_from_chkpt, Container::grid, Container::cosmology);
    } else if (indata_type == 1) {
        Container::indata = new LinWave(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 2) {
        Container::indata = new Schwarzschild(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 3) {
        Container::indata = new Flat(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 4) {
        Container::indata = new TOV(indata_input, Container::eos, Container::grid, Container::cosmology);
    } else if (indata_type == 5) {
        Container::indata = new Trumpet(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 6) {
        Container::indata = new RNS(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 7) {
        Container::indata = new Brill(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 8) {
        Container::indata = new Kerr(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 9) {
        Container::indata = new BowenYork(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 10) {
        Container::indata = new KerrSchild(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 11) {
        Container::indata = new KenTrumpet(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 12) {
        Container::indata = new Brill_Lindquist(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 13) {
        Container::indata = new Shock(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 14) {
        Container::indata = new EvansColeman(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 15) {
        Container::indata = new Rad_TOV(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 16) {
        Container::indata = new OS(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 17) {
        Container::indata = new RotPerfectFluid(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 18) {
        Container::indata = new Bondi(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 19) {
        Container::indata = new Choptuik(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 20) {
        Container::indata = new EMWave(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 21) {
        Container::indata = new RadHydroShockTest(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 22) {
        Container::indata = new TOV_BH(indata_input, Container::eos, Container::grid, Container::cosmology);
    } else if (indata_type == 23) {
        Container::indata = new SMS_TOV(indata_input, Container::eos, Container::grid, Container::cosmology);
    } else if (indata_type == 24) {
        Container::indata = new GaugeWave(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 25) {
        Container::indata = new WindTunnel(indata_input, Container::eos, Container::grid, Container::cosmology);
    } else if (indata_type == 26) {
        Container::indata = new ShibataWave(indata_input, Container::grid, Container::cosmology);
    } else if (indata_type == 27) {
        Container::indata = new Disk(indata_input, Container::eos, Container::grid, Container::cosmology);
    } else if (indata_type == 28) {
        Container::indata = new DualEMWave(indata_input, Container::grid, Container::cosmology);
    } else {
        cerr << " No such Indata type! " << endl;
        return 1;
    }
    // tell initial data whether initial metric should be rescaled
    Container::indata->RescaleMetric(rescale_metric);
    //================================================================
    // Allocate Slicing Class
    //================================================================
    if (slicing_type == 1) {
        Container::slicing = new Geodesic(Container::grid, Container::cosmology, eta_KO);
    } else if (slicing_type == 2) {
        Container::slicing = new OnePlusLog(Container::grid, Container::cosmology, eta_KO);
    } else if (slicing_type == 3) {
        Container::slicing = new Advective_OnePlusLog(Container::grid, Container::cosmology, eta_KO);
    } else if (slicing_type == 4) {
        Container::slicing = new Harmonic(Container::grid, Container::cosmology, eta_KO);
        // } else if (slicing_type == 5) {
        //   slicing = new Maximal(Container::grid);  
        // } else if (slicing_type == 6) {
        //   slicing = new Maximal_so(Container::grid);
            //  } else if (slicing_type == 7) {
            //    slicing = new KenLog(Container::grid, Container::cosmology);
    } else if (slicing_type == 8) {
        Container::slicing = new BonaMasso(Container::grid, Container::cosmology, eta_KO);
    } else {
        cerr << " No such Slicing type! " << endl;
        return 1;
    }
    //================================================================
    // Allocate Gauge Class
    //================================================================
    if (gauge_type == 1) {
        Container::gauge = new Constant_Shift(Container::grid);
    } else if (gauge_type == 2) {
        Container::gauge = new Gamma_Driver(Container::grid);
    } else if (gauge_type == 3) {
        Container::gauge = new Advective_Gamma_Driver(Container::grid);
    } else if (gauge_type == 4) {
        Container::gauge = new Jena_Gamma_Driver(Container::grid);
    } else if (gauge_type == 5) {
        Container::gauge = new Advective_Jena_Gamma_Driver(Container::grid);
        // } else if (gauge_type == 6) {
        //   gauge = new Mod_Gamma_Driver(Container::grid,eta);
    } else if (gauge_type == 7) {
        Container::gauge = new Covariant_Advective_Jena_Gamma_Driver(Container::grid);
    } else if (gauge_type == 8) {
        Container::gauge = new Self_Sim_Shift(Container::grid);
    } else {
        cerr << " No such gauge type! " << endl;
        return 1;
    }
    //================================================================
    // Allocate dumper Class
    //================================================================
    cout << " GR: allocating dumper..." << endl;
    Container::dump = new dumper(dump_step, Container::indata, Container::slicing, Container::gauge);
    //================================================================
    // Allocate monitor Class
    //================================================================
    cout << " GR: allocating monitor..." << endl;
    Container::monitor = new Monitor(note_step, monitor_filename, Container::grid, Container::indata,
        Container::slicing, Container::gauge, Container::cosmology,
        eta_KO, sigma,
        z4, kappa_11, kappa_12, kappa_2, kappa_ric,
        RK_order, char_OB);
    //================================================================
    // Allocate profile Class
    //================================================================
    if (write_profiles) {
        cout << " GR: allocating profiler..." << endl;
        Container::profiles = new Profiles(monitor_filename, Container::grid, note_step, Container::indata,
            Container::slicing, Container::gauge, Container::cosmology);
    } else {
        Container::profiles = NULL;
    }
    //================================================================
    // Allocate wave extraction class
    //================================================================
    if (extract_waves) {
        cout << " GR: allocating wave extraction..." << endl;
        Container::waves = new WaveExtraction(monitor_filename, Container::grid, note_step, Container::indata,
            Container::slicing, Container::gauge, Container::cosmology);
    } else {
        Container::waves = NULL;
    }
    //================================================================
    // Allocate photon Class
    //================================================================
    if (track_photons) {
        cout << " GR: allocating photons ..." << endl;
        Container::photons = new Photons(monitor_filename, Container::grid);
    } else {
        Container::photons = NULL;
    }
    //================================================================
    // Allocate manager Class
    //================================================================
    cout << " GR: allocating manager..." << endl;
    Container::manager = new Manager(matter_type,
        sigma, cowling, eta_KO,
        z4, kappa_11, kappa_12, kappa_2, kappa_ric,
        RK_order, char_OB, solve_constraints);
    //================================================================
    // Initialize
    //================================================================
    if (t_max > Container::grid->r_max())
        t_max = Container::grid->r_max();

    Container::manager->Set_t_max(t_max);
    bool success = Container::manager->Initialize();
    //  manager.Test_Indices();
    //  einstein.Test_Gridfunction();
    //  manager.Test_Ricci_for_Schwarzschild();
    //  einstein.Test_Lambda_InData();
    // manager.Test_Derivatives();
    //  einstein.Test_Transformation_Vector();
    //  manager.Test_Flat_Metric();
    //  einstein.Test_DivShift();
    //  einstein.Test_Interpolation();
    //================================================================
    // Integrate to time t_max
    //================================================================
    if (success) Container::manager->Integrate(t_max);
    //
    if (Container::eos != NULL) delete Container::eos;
    if (Container::profiles != NULL) delete Container::profiles;
    if (Container::waves != NULL) delete Container::waves;
    delete Container::gauge;
    delete Container::slicing;
    delete Container::cosmology;
    delete Container::grid;
    return error;
}
