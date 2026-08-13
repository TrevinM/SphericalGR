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
public:
    double lapse_min, lapse_r, lapse_th;
    double rho_max, rho_r, rho_th;
    double S_p_max, S_p_r, S_p_th;
    double a_p_max, a_p_r, a_p_th;
    int lapse_i, lapse_j, lapse_k;
    int rho_i, rho_j, rho_k;
    int S_p_i, S_p_j, S_p_k;
    int a_p_i, a_p_j, a_p_k;
private:
    ofstream monitorfile;
    std::vector<double> lapse_r_past, lapse_th_past;
    std::vector<double> rho_r_past, rho_th_past;
    std::vector<double> S_p_r_past, S_p_th_past;
    std::vector<double> delta;
    double avg;
    int worst;

public:
    // Constructor
    Tracker();

    // Deconstructor
    ~Tracker();

    // Method to be called every timestep
    void Execute();

    // Note pretracked variables
    void Note();

    double rFocus();
private:
    // Returns the average ignoring 1 outlier
    double Smooth(std::vector<double> past, double val);

};





#endif