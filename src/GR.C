//================================================================
//
// Code for scalar field in Spherical Coordinates
//
//===============================================================
#include "GR.h"
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
    cout << " GR: allocating grid..." << endl;
    grid = new Grid();
    //================================================================
    // Allocate checkpoint class
    //================================================================
    CheckPoint checkpoint(read_from_chkpt, chkpt_step, grid);
    //================================================================
    // Allocate Cosmology Class
    //================================================================
    cout << " GR: allocating cosmology..." << endl;
    if (space_type == 1) {
        cosmology = new Minkowski();
    } else if (space_type == 2) {
        cosmology = new DeSitter();
    } else if (space_type == 3) {
        cosmology = new Radiation();
    } else {
        cerr << " GR: no cosmology of type " << space_type << "!!" << endl;
        return 1;
    }
    //================================================================
    // Allocate Equation of State
    //================================================================
    cout << " GR: allocating EOS..." << endl;
    if (eos_type == 0) {
        eos = NULL;
    } else if (eos_type == 1) {
        eos = new polytrope((char*)"Polytrope_Input");
    } else if (eos_type == 2) {
        eos = new gamma_law((char*)"Polytrope_Input");
    } else if (eos_type == 3) {
        eos = new piece_polytrope((char*)"Piece_Polytrope_Input");
    } else if (eos_type == 4) {
        eos = new ideal_gas();
    } else if (eos_type == 5) {
        eos = new gas_radiation((char*)"GasRadiation_Input");
    } else {
        cout << " GR: Unknown EOS for eos_type = " << eos_type << endl;
    }
    //================================================================
    // Allocate InData Class
    //================================================================
    cout << " GR: allocating indata..." << endl;
    if (checkpoint.ReadFromChkpt()) {
        indata = new ReadFromCheckPoint(read_from_chkpt, grid, cosmology);
    } else if (indata_type == 1) {
        indata = new LinWave(indata_input, grid, cosmology);
    } else if (indata_type == 2) {
        indata = new Schwarzschild(indata_input, grid, cosmology);
    } else if (indata_type == 3) {
        indata = new Flat(indata_input, grid, cosmology);
    } else if (indata_type == 4) {
        indata = new TOV(indata_input, eos, grid, cosmology);
    } else if (indata_type == 5) {
        indata = new Trumpet(indata_input, grid, cosmology);
    } else if (indata_type == 6) {
        indata = new RNS(indata_input, grid, cosmology);
    } else if (indata_type == 7) {
        indata = new Brill(indata_input, grid, cosmology);
    } else if (indata_type == 8) {
        indata = new Kerr(indata_input, grid, cosmology);
    } else if (indata_type == 9) {
        indata = new BowenYork(indata_input, grid, cosmology);
    } else if (indata_type == 10) {
        indata = new KerrSchild(indata_input, grid, cosmology);
    } else if (indata_type == 11) {
        indata = new KenTrumpet(indata_input, grid, cosmology);
    } else if (indata_type == 12) {
        indata = new Brill_Lindquist(indata_input, grid, cosmology);
    } else if (indata_type == 13) {
        indata = new Shock(indata_input, grid, cosmology);
    } else if (indata_type == 14) {
        indata = new EvansColeman(indata_input, grid, cosmology);
    } else if (indata_type == 15) {
        indata = new Rad_TOV(indata_input, grid, cosmology);
    } else if (indata_type == 16) {
        indata = new OS(indata_input, grid, cosmology);
    } else if (indata_type == 17) {
        indata = new RotPerfectFluid(indata_input, grid, cosmology);
    } else if (indata_type == 18) {
        indata = new Bondi(indata_input, grid, cosmology);
    } else if (indata_type == 19) {
        indata = new Choptuik(indata_input, grid, cosmology);
    } else if (indata_type == 20) {
        indata = new EMWave(indata_input, grid, cosmology);
    } else if (indata_type == 21) {
        indata = new RadHydroShockTest(indata_input, grid, cosmology);
    } else if (indata_type == 22) {
        indata = new TOV_BH(indata_input, eos, grid, cosmology);
    } else if (indata_type == 23) {
        indata = new SMS_TOV(indata_input, eos, grid, cosmology);
    } else if (indata_type == 24) {
        indata = new GaugeWave(indata_input, grid, cosmology);
    } else if (indata_type == 25) {
        indata = new WindTunnel(indata_input, eos, grid, cosmology);
    } else if (indata_type == 26) {
        indata = new ShibataWave(indata_input, grid, cosmology);
    } else if (indata_type == 27) {
        indata = new Disk(indata_input, eos, grid, cosmology);
    } else if (indata_type == 28) {
        indata = new DualEMWave(indata_input, grid, cosmology);
    } else {
        cerr << " No such Indata type! " << endl;
        return 1;
    }
    // tell initial data whether initial metric should be rescaled
    indata->RescaleMetric(rescale_metric);
    //================================================================
    // Allocate Slicing Class
    //================================================================
    if (slicing_type == 1) {
        slicing = new Geodesic(grid, cosmology, eta_KO);
    } else if (slicing_type == 2) {
        slicing = new OnePlusLog(grid, cosmology, eta_KO);
    } else if (slicing_type == 3) {
        slicing = new Advective_OnePlusLog(grid, cosmology, eta_KO);
    } else if (slicing_type == 4) {
        slicing = new Harmonic(grid, cosmology, eta_KO);
        // } else if (slicing_type == 5) {
        //   slicing = new Maximal(grid);  
        // } else if (slicing_type == 6) {
        //   slicing = new Maximal_so(grid);
            //  } else if (slicing_type == 7) {
            //    slicing = new KenLog(grid, cosmology);
    } else if (slicing_type == 8) {
        slicing = new BonaMasso(grid, cosmology, eta_KO);
    } else {
        cerr << " No such Slicing type! " << endl;
        return 1;
    }
    //================================================================
    // Allocate Gauge Class
    //================================================================
    if (gauge_type == 1) {
        gauge = new Constant_Shift(grid);
    } else if (gauge_type == 2) {
        gauge = new Gamma_Driver(grid);
    } else if (gauge_type == 3) {
        gauge = new Advective_Gamma_Driver(grid);
    } else if (gauge_type == 4) {
        gauge = new Jena_Gamma_Driver(grid);
    } else if (gauge_type == 5) {
        gauge = new Advective_Jena_Gamma_Driver(grid);
        // } else if (gauge_type == 6) {
        //   gauge = new Mod_Gamma_Driver(grid,eta);
    } else if (gauge_type == 7) {
        gauge = new Covariant_Advective_Jena_Gamma_Driver(grid);
    } else if (gauge_type == 8) {
        gauge = new Self_Sim_Shift(grid);
    } else {
        cerr << " No such gauge type! " << endl;
        return 1;
    }
    //================================================================
    // Allocate dumper Class
    //================================================================
    cout << " GR: allocating dumper..." << endl;
    dump = new dump(dump_step, indata, slicing, gauge);
    //================================================================
    // Allocate monitor Class
    //================================================================
    cout << " GR: allocating monitor..." << endl;
    monitor = new monitor(note_step, monitor_filename, grid, indata,
        slicing, gauge, cosmology,
        eta_KO, sigma,
        z4, kappa_11, kappa_12, kappa_2, kappa_ric,
        RK_order, char_OB);
    //================================================================
    // Allocate profile Class
    //================================================================
    if (write_profiles) {
        cout << " GR: allocating profiler..." << endl;
        profiles = new Profiles(monitor_filename, grid, note_step, indata,
            slicing, gauge, cosmology);
    } else {
        profiles = NULL;
    }
    //================================================================
    // Allocate wave extraction class
    //================================================================
    if (extract_waves) {
        cout << " GR: allocating wave extraction..." << endl;
        waves = new WaveExtraction(monitor_filename, grid, note_step, indata,
            slicing, gauge, cosmology);
    } else {
        waves = NULL;
    }
    //================================================================
    // Allocate photon Class
    //================================================================
    if (track_photons) {
        cout << " GR: allocating photons ..." << endl;
        photons = new Photons(monitor_filename, grid);
    } else {
        photons = NULL;
    }
    //================================================================
    // Allocate manager Class
    //================================================================
    cout << " GR: allocating manager..." << endl;
    Manager manager(grid, &dump, &monitor, &checkpoint, profiles,
        waves, indata, eos,
        slicing, gauge, cosmology, photons, matter_type,
        sigma, cowling, eta_KO,
        z4, kappa_11, kappa_12, kappa_2, kappa_ric,
        RK_order, char_OB, solve_constraints);
    //================================================================
    // Initialize
    //================================================================
    if (t_max > grid->r_max())
        t_max = grid->r_max();

    manager.Set_t_max(t_max);
    bool success = manager.Initialize();
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
    if (success) manager.Integrate(t_max);
    //
    if (eos != NULL) delete eos;
    if (profiles != NULL) delete profiles;
    if (waves != NULL) delete waves;
    delete gauge;
    delete slicing;
    delete cosmology;
    delete grid;
    return error;
}
