// Tell emacs that this is -*-c++-*- mode

#include "DualMaxwell_State.h"
#include "DualMaxwell_Aux.h"
#include "Monitor.h"
#include "CheckPoint.h"

//================================================
//
// DualMaxwell
//
//================================================
class DualMaxwell : public Matter {
public:
    double rho_center, rho_c_max, rho_max, rho_max_MAX, drhoddr;    // diagnostics 
    int lapse_i, lapse_j, lapse_k, rho_i, rho_j, rho_k;
    dualmaxwell_state* last, * derivs, * inter, * updates;
    dualmaxwell_aux* aux;
    Monitor* monitor;
    CheckPoint* checkpoint;
    ofstream monitorfile;
    double eta_KO;   // Kreiss-Oliger coefficient
    int char_OB;    // decides how Sommerfeld BCs are implemented
    // 0: derivatives; 1: characteristic interpolation
public:
    DualMaxwell(Grid* grid_i, dumper* dump_i, InData* indata_i,
        int cowling_i, int char_OB_i, Cosmology* cosmology_i,
        Monitor* monitor_i, double eta_KO_i, CheckPoint* checkpoint_i) :
        Matter(grid_i, dump_i, indata_i, cowling_i, cosmology_i),
        char_OB(char_OB_i), monitor(monitor_i), eta_KO(eta_KO_i),
        checkpoint(checkpoint_i) {
        cout << " DUALMAXWELL: setting up dualmaxwell... " << endl;
        // 
        // create states for dynamical variables
        // 
        cout << " DUALMAXWELL: setting up states... " << endl;
        last = new dualmaxwell_state(grid, dump, "dualmaxwell_last");
        derivs = new dualmaxwell_state(grid, dump, "dualmaxwell_derivs");
        inter = new dualmaxwell_state(grid, dump, "dualmaxwell_inter");
        updates = new dualmaxwell_state(grid, dump, "dualmaxwell_updates");
        // create dump list for state last:
        last->assemble_dump_list("Dump_List");
        derivs->assemble_dump_list("Dump_List");
        // let checkpointer know about last
        checkpoint->CollectDynVariables(last);
        //
        // create ADM sources
        // 
        adm_sources = new ADM_Source_Terms(grid, dump, "dualmaxwell_souces");
        adm_sources->assemble_dump_list("Dump_List");
        //
        // create auxiliary variables
        // 
        aux = new dualmaxwell_aux(grid, dump, "dualmaxwell_aux");
        aux->assemble_dump_list("Dump_List");
        //
       // finally create a monitor file...
       //

        ostringstream monfilename;
        monfilename << "output/" << monitor->Filestem() << "_" << N_r - 2 * N_g << "_"
            << N_t - 2 * N_g << ".dualmaxwell_mon" << ends;
        monitorfile.open(monfilename.str().c_str());
        monitorfile.setf(ios::right);
        time_t clocktime;
        struct tm* currenttime;
        time(&clocktime);
        currenttime = localtime(&clocktime);
        monitorfile << "# File created on " << asctime(currenttime);
        monitorfile << "# DualMaxwell evolution " << endl;
        monitorfile << "# " << setw(16) << "time"
            << setw(18) << "pr time (r=0)"
            << setw(18) << "rho_ADM_c"
            << setw(18) << "rho_ADM_c_max"
            << setw(18) << "rho_ADM_max"
            << setw(18) << "rho_ADM_max_MAX"
            << setw(18) << "r(rho_ADM_max)"
            << setw(18) << "th(rho_ADM_max)"
            << setw(18) << "r(lapse_min)"
            << setw(18) << "th(lapse_min)"
            << endl;
        monitorfile << "#=======================================================================================================================================================================" << endl;
        //
        rho_center = rho_c_max = drhoddr = 0.0;
        rho_max = rho_max_MAX = 0.0;
    };
    ~DualMaxwell() {
        delete last;
        delete derivs;
        delete inter;
        delete updates;
        delete adm_sources;
        delete aux;
        monitorfile.close();
        cout << " DUALMAXWELL: destructing derived class DualMaxwell " << endl;
    };
    const char* Name() { return "DualMaxwell's equations"; }
    //================================================
    // initialize
    //================================================
    void Initialize(state* s, curvature* c, diagnostics* d);
    //===============================================
    // Compute RHS sides for matter equations
    //===============================================
    void Compute_RHS(state* s, curvature* c, double time = 0.0);
    void dot_a_as(dualmaxwell_state* m, state* s, double time = 0.0);
    //  void dot_as(dualmaxwell_state *m, state *s, curvature *c, double time = 0.0);
    //===============================================
    // Start and finish RK steps
    //===============================================
    void Start_RK() {
        inter->equals(last);
        updates->equals(last);
    };
    void Finish_RK(state* s, curvature* c, double dt) {
        updates->char_OB(last, dt);
        last->equals(updates);
        last->fill_ghosts();
        inter->equals(last);
        ADM_Sources(s, c);
    }
    //===============================================
    // Update matter state (adds dt * derivs to update)
    //===============================================
    void Update(double dt) {
        updates->add(dt, derivs);
    };
    //===============================================
    // Compute intermediate state (computes inter = last + dt * derivs )
    //===============================================
    void Compute_inter(double dt) {
        inter->add(last, dt, derivs);
        inter->fill_ghosts();
        inter->char_OB(last, dt);
    };
    //===============================================
    // Regridding etc
    //===============================================
    int Regrid(VecDoub r_new) { last->Regrid(r_new); return 0; };
    double RegridCriterion() {
        const Doub tiny = 1.e-12;
        // const Doub drhoddr = abs(a_p_o.ddr(N_g,N_g,N_g));
        Doub rho_l = fabs(rho_center);
        if (rho_l < 1.0) rho_l = 1.0;
        const Doub scale = sqrt(rho_l / (fabs(drhoddr) + tiny));
        const Doub grid_scale = grid->delta_r(N_g);
        // if ( grid_scale / scale > 0.15 ) {
        //   cout << " rho_center = " << rho_center 
        // 	   << " drhoddr = " << drhoddr 
        // 	   << " scale = " << scale 
        // 	   << " grid_scale = " << grid->delta_r(N_g) 
        // 	   << endl;
        // }
        return grid_scale / scale;
    };
    gf3d* MatterField() { return &last->a_p; };
    //================================================
    // Note
    //================================================
    void Note(int time_step, double time, double tau_c) {
        if (rho_center > rho_c_max) rho_c_max = rho_center;
        if (time_step % monitor->Note_Step() == 0) {
            monitorfile.setf(ios::left);
            monitorfile << setw(18) << time
                << setw(18) << tau_c
                << setprecision(10) << setw(18) << rho_center
                << setprecision(10) << setw(18) << rho_c_max
                << setprecision(10) << setw(18) << rho_max
                << setprecision(10) << setw(18) << rho_max_MAX
                << setprecision(10) << setw(18) << grid->r(rho_i)
                << setprecision(10) << setw(18) << grid->theta(rho_j)
                << setprecision(10) << setw(18) << grid->r(lapse_i)
                << setprecision(10) << setw(18) << grid->theta(lapse_j)
                << endl;
        }
    };
    //================================================
    // dump grid functions
    //================================================
    void dump_fcts(double time, double prop_time, int timestep,
        const char* suffix) {
        if (dump->time_to_dump(timestep) || strcmp(suffix, "")) {
            cout << " DUALMAXWELL: Dumping matter functions at time t = "
                << time << endl;
            last->dump_fcts(time, prop_time, timestep, suffix);
            adm_sources->dump_fcts(time, prop_time, timestep, suffix);
            aux->dump_fcts(time, prop_time, timestep, suffix);
        }
    }
    //================================================
    // Diagnostics
    //================================================
    double Compute_Diagnostics(state* s, curvature* c) {
        for (int i = 0; i < N_r; i++) {
            const double rl = grid->r(i);
            for (int j = 0; j < N_t; j++) {
                const double sinthetal = grid->sintheta(j);
                for (int k = 0; k < N_p; k++) {
                    // recall that gup is *not* rescaled...
                    aux->A2[i][j][k] = exp(-4.0 * s->phi(i, j, k)) * c->gup_pp(i, j, k) *
                        last->a_p(i, j, k) * last->a_p(i, j, k) * rl * rl * sinthetal * sinthetal;
                    const double g_pp = exp(4.0 * s->phi(i, j, k))
                        * (1.0 + s->h_pp(i, j, k));
                    aux->A_xi[i][j][k] = last->a_p(i, j, k) / sqrt(g_pp);
                }
            }
        }
        return aux->A2(0.0, N_g, N_g);
    };
    //================================================
    // ADM sources
    //================================================
    void ADM_Sources(state* s, curvature* c);
    //================================================
    // Compute curl
    //================================================
    void curl(gf3d& a_r, gf3d& a_t, gf3d& a_p,
        double& curl_a_r, double& curl_a_t, double& curl_a_p,
        int i, int j, int k, state* s);
    //================================================
    // Compute Lie derivative of covariant vector
    //================================================  
    void Lie_covariant(gf3d& b_r, gf3d& b_t, gf3d& b_p,
        gf3d& a_r, gf3d& a_t, gf3d& a_p,
        gf3d& Ll_r, gf3d& Ll_t, gf3d& Ll_p);
};




