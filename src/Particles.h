// Tell emacs that this is -*-c++-*- mode
// 
//================================================
//
// Contains all methods for fluid particles
//
//================================================
//
#ifndef PARTICLES_H
#define PARTICLES_H
#include "gridfunction.h"
//
//================================================
//================================================
//
// first: an individual particle
//
//================================================
//================================================
//
class Particle {
private:
    double last[4];
    double inter[4];
    double updates[4];
    double derivs[4];
    int number;
    int N_g;
public:
    //============================================
    // constructor
    //============================================ 
    Particle(double r_i, double t_i, double p_i,
        //	   gf3d * lapse_i, gf3d * W_i,
        //	   gf3d * shift_r_i, gf3d * shift_t_i, gf3d * shift_p_i, 
        //	   gf3d * v_r_i, gf3d * v_t_i, gf3d * v_p_i,
        Grid* grid, int number_i) :
        //    lapse(lapse_i), W(W_i),
        //    shift_r(shift_r_i), shift_t(shift_t_i), shift_p(shift_p_i), 
        //    v_r(v_r_i), v_t(v_t_i), v_p(v_p_i), 
        number(number_i) {
        N_g = grid->N_ghosts();
        last[0] = inter[0] = 0.0;  // tau
        last[1] = inter[1] = r_i;  // r
        last[2] = inter[2] = t_i;  // theta
        last[3] = inter[3] = p_i;  // phi
        for (int i = 0; i < 4; i++) {
            inter[i] = last[i];
            updates[i] = last[i];
            derivs[i] = 0.0;
        }
    };
    //============================================
    // destructor
    //============================================
    ~Particle() {};
    //============================================
    // get location
    //============================================
    void Get_Location(double& tau, double& r, double& theta, double& phi) {
        tau = last[0];
        r = last[1];
        theta = last[2];
        phi = last[3];
    }
    //============================================
    // return address and number
    //============================================
    Particle* Address() { return this; };
    int Number() { return number; };
    //============================================
    // compute time derivatives for inter_pos
    //============================================
    void Compute_RHS(gf3d* lapse, gf3d* W,
        gf3d* shift_r, gf3d* shift_t, gf3d* shift_p,
        gf3d* v_r, gf3d* v_t, gf3d* v_p) {
#ifndef AXISYMMETRY
        cout << " PARTICLE: OOOUUUCCCHHH : need 3D interpolater!!! " << endl;
#endif
        const double rl = inter[1];
        const double tl = inter[2];
        const double W_l = (*W)(rl, tl, N_g);
        const double lapse_l = (*lapse)(rl, tl, N_g);
        const double v_r_l = (*v_r)(rl, tl, N_g);
        const double v_t_l = (*v_t)(rl, tl, N_g) / rl;
        const double v_p_l = (*v_p)(rl, tl, N_g) / (rl * sin(tl));
        const double shift_r_l = (*shift_r)(rl, tl, N_g);
        const double shift_t_l = (*shift_t)(rl, tl, N_g) / rl;
        const double shift_p_l = (*shift_p)(rl, tl, N_g) / (rl * sin(tl));
        derivs[0] = lapse_l / W_l;
        if (!isfinite(derivs[0])) cout << " OOOOPPPS " << lapse_l
            << "  " << W_l << endl;
        derivs[1] = lapse_l * v_r_l - shift_r_l;
        derivs[2] = lapse_l * v_t_l - shift_t_l;
        derivs[3] = lapse_l * v_p_l - shift_p_l;
    };
    //============================================
    // Runge Kutta stuff...
    //============================================ 
    void Start_RK() { for (int i = 0; i < 4; i++) inter[i] = last[i]; }
    void Finish_RK() { for (int i = 0; i < 4; i++) inter[i] = last[i] = updates[i]; }
    void Update(double dt) { for (int i = 0; i < 4; i++) updates[i] += dt * derivs[i]; }
    void Compute_inter(double dt) { for (int i = 0; i < 4; i++) inter[i] = last[i] + dt * derivs[i]; }
};
//
//================================================
//================================================
//
// Then: a whole bunch of particles...
//
//================================================
//================================================
//
class Particles {
private:
    Particle** tracers;
    int N_particles;
    gf3d** fct_list;
    int N_fcts;           // actual number of functions to be tracked
    int N_fcts_max;       // maximum possible number - set below
    int N_r, N_theta, N_g;
public:
    //============================================
    // constructor
    //============================================ 
    Particles(int N_particles_i, double R_max, Grid* grid) :
        N_particles(N_particles_i) {
        cout << " PARTICLES: initializing " << N_particles << " fluid tracers... " << endl;
        tracers = new Particle * [N_particles];
        double r_init = 0.0;
        double theta_init = 0.0;
        double phi_init = 0.0;
        for (int i = 0; i < N_particles; i++) {
            Init_Position(i, R_max, r_init, theta_init, phi_init);
            tracers[i] = new Particle(r_init, theta_init, phi_init, grid, i);
        }
        //
        // deal with list of functions - then call "add_to_dump_list..." to 
        // populate
        // 
        N_fcts = 0;
        N_fcts_max = 20;
        fct_list = new gf3d * [N_fcts_max];
        //
        // create file stem for output
        //
        N_r = grid->N_r_tot();
        N_theta = grid->N_theta_tot();
        N_g = grid->N_ghosts();
    }
    //============================================
    // destructor
    //============================================ 
    ~Particles() {
        for (int i = 0; i < N_particles; i++) {
            delete tracers[i];
        }
        delete tracers;
        delete fct_list;
        cout << " PARTICLES: ... closing fluid tracers... " << endl;
    }
    //============================================
    // compute time derivatives for inter_pos
    //============================================
    void Compute_RHS(gf3d* lapse, gf3d* W,
        gf3d* shift_r, gf3d* shift_t, gf3d* shift_p,
        gf3d* v_r, gf3d* v_t, gf3d* v_p) {
        for (int i = 0; i < N_particles; i++) {
            tracers[i]->Compute_RHS(lapse, W,
                shift_r, shift_t, shift_p,
                v_r, v_t, v_p);
        }
    }
    //============================================
    // Runge Kutta stuff
    //============================================ 
    void Start_RK() { for (int i = 0; i < N_particles; i++) tracers[i]->Start_RK(); };
    void Finish_RK() { for (int i = 0; i < N_particles; i++) tracers[i]->Finish_RK(); };
    void Update(double dt) { for (int i = 0; i < N_particles; i++) tracers[i]->Update(dt); };
    void Compute_inter(double dt) { for (int i = 0; i < N_particles; i++) tracers[i]->Compute_inter(dt); };
    //============================================
    // Assemble dump list
    //============================================ 
    template <class bundle>
    int add_to_dump_list(const char* dump_list_file, bundle* s) {
        ifstream infile;
        infile.open(dump_list_file);
        if (!infile) {
            cerr << " PARTICLES: can't open " << dump_list_file
                << " for input." << endl;
            return 0;
        }
        cout << " PARTICLES: looking for files listed in "
            << dump_list_file << " in " << s->Name() << endl;
        char fct_name[64];
        infile >> fct_name;
        while (!infile.eof() && N_fcts < N_fcts_max) {
            for (int i = 0; i < s->N_fcts; i++) {
                if (!strcmp(fct_name, s->fct_list[i]->Name())) {
                    cout << " ... found match for function name " << fct_name
                        << " in " << s->Name() << endl;
                    fct_list[N_fcts] = s->fct_list[i]->Address();
                    N_fcts++;
                }
            }
            infile >> fct_name;
        }
        if (N_fcts == N_fcts_max) {
            cout << " PARTICLES: ran out of room for functions to track - "
                << " increase N_fcts_max! " << endl;
        }
        return N_fcts;
    }
    //============================================
    // Gather info and dump...
    //============================================ 
    void Dump(double t = 0, double tau = 0, int timestep = 0,
        const char* suffix = "") {
        //
        // create output file
        //
        ofstream outfile;
        ostringstream filename;
        filename << "Particles_" << N_r - 2 * N_g
            << "_" << N_theta - 2 * N_g
            << "_" << setfill('0') << setw(8) << timestep
            << suffix << ends;
        outfile.open(filename.str().c_str());
        if (!outfile) cerr << " Could not open file " << filename.str().c_str()
            << " in Particles::Dump " << endl;
        outfile.setf(ios::right);
        int digs = 14;
        outfile << "# Particle data at t = " << t
            << " and tau_c = " << tau << endl;
        outfile << "# " << setw(digs - 2) << " tau "
            << setw(digs) << " r "
            << setw(digs) << " theta "
            << setw(digs) << " phi ";
        for (int fct_count = 0; fct_count < N_fcts; fct_count++)
            outfile << setw(digs) << fct_list[fct_count]->Name();
        outfile << endl;
        int characters = (4 + N_fcts) * digs;
        outfile << "#" << setfill('=') << setw(characters) << "=" << endl;
        outfile << setfill(' ');
        //
        // now loop over all particles
        //
        double tau_p, r, theta, phi;
        for (int i = 0; i < N_particles; i++) {
            tracers[i]->Get_Location(tau_p, r, theta, phi);
            outfile << setw(digs) << tau_p
                << setw(digs) << r
                << setw(digs) << theta
                << setw(digs) << phi;
            //
            // now find values of all functions in function list
            // 
            for (int fct_count = 0; fct_count < N_fcts; fct_count++) {
#ifndef AXISYMMETRY
                cout << " FIX ME!!! " << endl;
#endif
                double fct_value = (*fct_list[fct_count])(r, theta, N_g);
                outfile << setw(digs) << fct_value;
            }
            outfile << endl;
        }
        outfile.close();
    }
    //============================================
    // provide initial positions - for now uniform, one line in radius...
    //============================================ 
    void Init_Position(int i, double R_max,
        double& r_init, double& theta_init, double& phi_init) {
        const double PI = acos(-1.0);
        r_init = double(i + 0.5) / double(N_particles) * R_max;
        theta_init = PI / 2.0;
        phi_init = 0.0;
    };
    //============================================
    // Find radius of outermost particle
    //============================================ 
    double RadiusOuterParticle() {
        double tau_p, r, theta, phi;
        tracers[N_particles - 1]->Get_Location(tau_p, r, theta, phi);
        return r;
    };
};


#endif    /* PARTICLES_H */    

