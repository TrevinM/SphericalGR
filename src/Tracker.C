#include "Tracker.h"
#include "Container.h"
#include "Manager.h"


//================================================
// Constructor
//================================================
Tracker::Tracker() {
    ostringstream monfilename;
    monfilename << "output/" << Container::monitor->Filestem() << "_" << Manager::matter->N_r - 2 * Manager::matter->N_g << "_"
        << Manager::matter->N_t - 2 * Manager::matter->N_g << ".tracker_mon" << ends;
    monitorfile.open(monfilename.str().c_str());
    monitorfile.setf(ios::left);
    time_t clocktime;
    struct tm* currenttime;
    time(&clocktime);
    currenttime = localtime(&clocktime);
    monitorfile << "# File created on " << asctime(currenttime);
    monitorfile << "# " << setw(16) << "time"
        << setw(18) << "rho_ADM max"
        << setw(18) << "r"
        << setw(18) << "th"
        << setw(18) << "lapse min"
        << setw(18) << "r"
        << setw(18) << "th"
        << setw(18) << "S_p abs_max"
        << setw(18) << "r"
        << setw(18) << "th"
        << endl;
    monitorfile << "#=======================================================================================================================================================================" << endl;

}

//================================================
// Deconstructor
//================================================
Tracker::~Tracker() {
    cout << " TRACKER: destructing tracker" << endl;
}

//================================================
// Execute
//================================================
void Tracker::Execute() {
    lapse_min = Manager::last->lapse.min(lapse_i, lapse_j, lapse_k);
    rho_max = Manager::matter->adm_sources->rho_ADM.max(lapse_i, lapse_j, lapse_k);
    S_p_max = Manager::matter->adm_sources->S_p.abs_max(lapse_i, lapse_j, lapse_k);

    lapse_r = Smooth(lapse_r_past, Container::grid->r(lapse_i));
    lapse_th = Smooth(lapse_th_past, Container::grid->theta(lapse_j));
    rho_r = Smooth(rho_r_past, Container::grid->r(rho_i));
    rho_th = Smooth(rho_th_past, Container::grid->theta(rho_j));
    S_p_r = Smooth(S_p_r_past, Container::grid->r(S_p_i));
    S_p_th = Smooth(S_p_th_past, Container::grid->theta(S_p_j));
}

void Tracker::Note() {
    monitorfile.setf(ios::left);
    monitorfile << setw(18) << time << Container::manager->t
        << setprecision(10) << setw(18) << rho_max
        << setprecision(10) << setw(18) << rho_r
        << setprecision(10) << setw(18) << rho_th
        << setprecision(10) << setw(18) << lapse_min
        << setprecision(10) << setw(18) << lapse_r
        << setprecision(10) << setw(18) << lapse_th
        << setprecision(10) << setw(18) << S_p_max
        << setprecision(10) << setw(18) << S_p_r
        << setprecision(10) << setw(18) << S_p_th
        << endl;
}


double Tracker::Smooth(std::vector<double> past, double val) {
    if (past.size() < 4) {
        past.push_back(val);
        return 0.0;
    } else {
        past[0] = past[1];
        past[1] = past[2];
        past[2] = past[3];
        past[3] = val;

        avg = (past[0] + past[1] + past[2] + past[3])/4.;
        delta = {abs(past[0] - avg), abs(past[1] - avg),  abs(past[2] - avg), abs(past[3] - avg)};
        worst = 0;
        for ( int i=1 ; i < 4 ; i++ ) {
            if (delta[i] > delta[worst]) {
                worst = i;
            }
        }
        avg = 0;
        for ( int i=0 ; i < 4 ; i++ ) {
            if ( i != worst ) {
                avg += past[i]/3.;
            }
        }

        return avg;
    }
}