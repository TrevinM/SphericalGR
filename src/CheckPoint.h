// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing routines that handle check-pointing
//================================================
//
#ifndef CHKPT_H
#define CHKPT_H

#include "Grid.h"
#include "State.h"
// #include "Matter.h"
#include "gridfunction.h"



class CheckPoint {
private:
    Grid* grid;
    gf3d** fct_list;
    int N_fcts, N_max;
    int start_from_chkpt, chkpt_step;
    bool read_from_chkpt;
    // const char * chkpt_file;
    double t_start, tau_c_start;
    double r_out_start;
public:
    //=========================================
    // Constructor
    //=========================================
    CheckPoint(int start_from_chkpt_i, int chkpt_step_i, Grid* grid_i) :
        start_from_chkpt(start_from_chkpt_i),
        chkpt_step(chkpt_step_i), grid(grid_i) {
        cout << " CHECKPOINT: creating checkpoint class..." << endl;
        if (start_from_chkpt > 0) {
            cout << " CHECKPOINT: will start from checkpoint files at time step "
                << start_from_chkpt << endl;
            read_from_chkpt = true;
        } else {
            start_from_chkpt = 0;
            read_from_chkpt = false;
        };
        N_fcts = 0;
        N_max = 75;
        fct_list = new gf3d * [N_max];
        //
        // read summary file
        //
        t_start = 0.0;
        tau_c_start = 0.0;
        if (start_from_chkpt != 0) {
            ReadSummaryFile(start_from_chkpt);
        }
    }
    //=========================================
    // Destructor
    //=========================================
    ~CheckPoint() {};
    //=========================================
    // checkpoint index
    //=========================================
    int CheckPointStep() { return chkpt_step; };
    //=========================================
    // Read from check point?
    //=========================================
    bool ReadFromChkpt() { return read_from_chkpt; };
    //=========================================
    // Initial times and outer radius
    //=========================================
    double TStart() { return t_start; }
    double TauStart() { return tau_c_start; }
    double ROutStart() { return r_out_start; }
    int TimeStepStart() { return start_from_chkpt; }
    //=========================================
    // Collect all dynamical variables
    //=========================================
    template <class dyn_vars>
    void CollectDynVariables(dyn_vars* in) {
        cout << " CHECKPOINT: Adding " << in->Name() << " to checkpoint list "
            << endl;
        for (int i = 0; i < in->N_fcts; i++) {
            if (N_fcts > N_max - 1) {
                cout << " CHECKPOINT: N_fcts > N_max -- need to increase N_max "
                    << endl;
                exit(0);
            }
            fct_list[N_fcts] = in->fct_list[i]->Address();
            // cout << " CHECKPOINT: " << N_fcts << "  " << fct_list[N_fcts]->Name()
            //	   << endl;
            N_fcts++;
        }
    }
    //=========================================
    // Write out check-point file
    //=========================================
    bool WriteCheckPoint(int timestep, double t, double tau_c) {
        // write checkpoint files
        cout << " CHECKPOINT: writing checkpoint for timestep = "
            << timestep << endl;
        for (int i = 0; i < N_fcts; i++)
            fct_list[i]->checkpoint(timestep);
        // write summary file
        WriteSummaryFile(timestep, t, tau_c);
        stringstream command1, command2;
        command1 << "tar -czf CHKPT_" << setfill('0') << setw(8)
            << timestep << ".tgz *_" << setfill('0') << setw(8)
            << timestep << ".cpt *Input *List code.tgz" << ends;
        cout << command1.str().c_str() << endl;
        system(command1.str().c_str());
        command2 << "rm" << " *_" << setfill('0') << setw(8)
            << timestep << ".cpt" << ends;
        cout << command2.str().c_str() << endl;
        system(command2.str().c_str());
        return true;
    }

    //=========================================
    // Write summary file
    //=========================================
    bool WriteSummaryFile(int timestep, double t, double tau_c) {
        ofstream outfile;
        ostringstream outfilename;
        outfilename << "CHKPT_Summary_" << setfill('0') << setw(8)
            << timestep << ".cpt" << ends;
        outfile.open(outfilename.str().c_str());
        if (!outfile) {
            cerr << " CHECKPOINT: Could not open checkpoint summary file "
                << outfilename.str().c_str() << endl;
            return false;
        }
        outfile << " Checkpoint at timestep = " << setprecision(16)
            << timestep << endl;
        outfile << "                   time = " << t << endl;
        outfile << "  proper time at center = " << tau_c << endl;
        outfile << "         outer boundary = " << grid->r_max() << endl;
        outfile.close();
        return true;
    }
    //=========================================
    // Read summary file
    //=========================================
    bool ReadSummaryFile(int timestep) {
        ifstream infile;
        stringstream infilename;
        infilename << "CHKPT_Summary_" << setfill('0') << setw(8)
            << timestep << ".cpt" << ends;
        infile.open(infilename.str().c_str());
        if (!infile) {
            cerr << " CHECKPOINT: Could not open checkpoint summary file "
                << infilename.str().c_str() << " for reading" << endl;
            return false;
            exit(0);
        }
        char buf[500], c;
        int index;
        infile.get(buf, 1000, '='); infile.get(c); infile >> index;
        infile.get(buf, 1000, '='); infile.get(c); infile >> t_start;
        infile.get(buf, 1000, '='); infile.get(c); infile >> tau_c_start;
        infile.get(buf, 1000, '='); infile.get(c); infile >> r_out_start;
        infile.close();
        if (index != timestep) {
            cout << " CHECKPOINT: fatal error in ReadSummaryFile: index != timestep "
                << endl;
            exit(0);
        }
        cout << " CHECKPOINT: found t_start = " << t_start << endl;
        return true;
    }
};
#endif  /* CHKPT */
