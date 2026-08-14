#include "Tracker.h"
#include "Container.h"
#include "Manager.h"
#include <typeinfo>
#include "Matter.h"

//================================================
// Constructor
//================================================
Tracker::Tracker() {
    ostringstream monfilename;
    monfilename << "output/" << Container::monitor->Filestem() << "_" << Container::matter->N_r - 2 * Container::matter->N_g << "_"
        << Container::matter->N_t - 2 * Container::matter->N_g << ".tracker_mon" << ends;
    monitorfile.open(monfilename.str().c_str());
    monitorfile.setf(ios::left);
    time_t clocktime;
    struct tm* currenttime;
    time(&clocktime);
    currenttime = localtime(&clocktime);
    monitorfile << "# File created on " << asctime(currenttime);
    monitorfile << "# " << setw(16) << "time"
        << setw(18) << "r_focus"
        << setw(18) << "rho_ADM max"
        << setw(18) << "(r,"
        << setw(18) << "th)"
        << setw(18) << "lapse"
        << setw(18) << "(r,"
        << setw(18) << "th)"
        << setw(18) << "S_p abs_max"
        << setw(18) << "(r,"
        << setw(18) << "th)"
        << setw(18) << "Matter max"
        << setw(18) << "(r,"
        << setw(18) << "th)"
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
// Updates tracked vales
//================================================
void Tracker::Update() {
    int timestep = Container::manager->Step();
    if (timestep != timestep_last_track) {
        timestep_last_track = timestep;
        r_min = Container::grid->r(Container::grid->N_ghosts());
        lapse_min = Manager::last->lapse.min(lapse_i, lapse_j, lapse_k);
        rho_max = Container::matter->adm_sources->rho_ADM.max(rho_i, rho_j, rho_k);
        S_p_max = Container::matter->adm_sources->S_p.abs_max(S_p_i, S_p_j, S_p_k);
        matter_max = Container::matter->MatterField()->abs_max(matter_i, matter_j, matter_k);

        lapse_r = Container::grid->r(lapse_i);
        lapse_th = Container::grid->theta(lapse_j);
        rho_r = Container::grid->r(rho_i);
        rho_th = Container::grid->theta(rho_j);
        S_p_r = Container::grid->r(S_p_i);
        S_p_th = Container::grid->theta(S_p_j);
        matter_r = Container::grid->r(matter_i);
        matter_th = Container::grid->theta(matter_j);

        num_focus = 0;
        r_focus = 0.;
        if (lapse_r > r_min) {
            r_focus += lapse_r;
            num_focus++;
        }
        if (rho_r > r_min) {
            r_focus += rho_r;
            num_focus++;
        }
        if (S_p_r > r_min) {
            r_focus += S_p_r;
            num_focus++;
        }
        if (matter_r > r_min) {
            r_focus += matter_r;
            num_focus++;
        }
        if (num_focus > 0){
            r_focus /= num_focus;
        }
    }
}

void Tracker::Note() {
    Update();
    monitorfile.setf(ios::left);
    monitorfile << setw(18) << Container::manager->Time()
        << setprecision(10) << setw(18) << rFocus()
        << setprecision(10) << setw(18) << rho_max
        << setprecision(10) << setw(18) << rho_r
        << setprecision(10) << setw(18) << rho_th
        << setprecision(10) << setw(18) << lapse_min
        << setprecision(10) << setw(18) << lapse_r
        << setprecision(10) << setw(18) << lapse_th
        << setprecision(10) << setw(18) << S_p_max
        << setprecision(10) << setw(18) << S_p_r
        << setprecision(10) << setw(18) << S_p_th
        << setprecision(10) << setw(18) << matter_max
        << setprecision(10) << setw(18) << matter_r
        << setprecision(10) << setw(18) << matter_th
        << endl;
}


// double Tracker::Smooth(std::vector<double> &past, double val) {
//     if (past.size() < 4) {
//         past.push_back(val);
//         return 0.0;
//     } else {
//         past[0] = past[1];
//         past[1] = past[2];
//         past[2] = past[3];
//         past[3] = val;

//         avg = (past[0] + past[1] + past[2] + past[3])/4.;
//         delta = {abs(past[0] - avg), abs(past[1] - avg),  abs(past[2] - avg), abs(past[3] - avg)};
//         worst = 0;
//         for ( int i=1 ; i < 4 ; i++ ) {
//             if (delta[i] > delta[worst]) {
//                 worst = i;
//             }
//         }
//         avg = 0;
//         for ( int i=0 ; i < 4 ; i++ ) {
//             if ( i != worst ) {
//                 avg += past[i]/3.;
//             }
//         }

//         return avg;
//     }
// }

double Tracker::rFocus() {
    Update();
    return r_focus;
}
