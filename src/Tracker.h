#ifndef TRACKER_H
#define TRACKER_H

#include <vector>
#include "Monitor.h"

//================================================
//
// Class that follows extrema of variables
//
//================================================

class Tracker {
private:
    ofstream monitorfile;

    int timestep_last_track;
    double r_min;
    double lapse_min, lapse_r, lapse_th;
    double rho_max, rho_r, rho_th;
    double S_p_max, S_p_r, S_p_th;
    double matter_max, matter_r, matter_th;

    int lapse_i, lapse_j, lapse_k;
    int rho_i, rho_j, rho_k;
    int S_p_i, S_p_j, S_p_k;
    int matter_i, matter_j, matter_k;

    int num_focus;
    double r_focus;

    // std::vector<double> lapse_r_past, lapse_th_past;
    // std::vector<double> rho_r_past, rho_th_past;
    // std::vector<double> S_p_r_past, S_p_th_past;
    // std::vector<double> matter_r_past, matter_th_past;
    // std::vector<double> delta;

    double avg;
    int worst;
public:
    // Constructor
    Tracker();

    // Public getters
    double rFocus();

    // Deconstructor
    ~Tracker();

    // Note pretracked variables
    void Note();

private:
    // Update variables before using them
    void Update();

};





#endif