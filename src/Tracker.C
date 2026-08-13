#include "Tracker.h"
#include "Container.h"
#include "Manager.h"
#include <typeinfo>
#include "DualMaxwell.h"

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
        << setw(18) << "r smooth"
        << setw(18) << "th smooth"
        << setw(18) << "lapse min"
        << setw(18) << "r smooth"
        << setw(18) << "th smooth"
        << setw(18) << "S_p abs_max"
        << setw(18) << "r smooth"
        << setw(18) << "th smooth"
        << setw(18) << "a_p max"
        << setw(18) << "r"
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
    rho_max = Container::matter->adm_sources->rho_ADM.max(lapse_i, lapse_j, lapse_k);
    S_p_max = Container::matter->adm_sources->S_p.abs_max(lapse_i, lapse_j, lapse_k);
    if (typeid(Container::matter) == typeid(DualMaxwell)) {
        auto dual_em = static_cast<DualMaxwell*>(Container::matter);
        a_p_max = dual_em->last->a_p.max(a_p_i, a_p_j, a_p_k);
    } else {
        a_p_max = 0.0;
        a_p_i = a_p_j = a_p_k = Container::grid->N_ghost;
    }

    lapse_r = Smooth(lapse_r_past, Container::grid->r(lapse_i));
    lapse_th = Smooth(lapse_th_past, Container::grid->theta(lapse_j));
    rho_r = Smooth(rho_r_past, Container::grid->r(rho_i));
    rho_th = Smooth(rho_th_past, Container::grid->theta(rho_j));
    S_p_r = Smooth(S_p_r_past, Container::grid->r(S_p_i));
    S_p_th = Smooth(S_p_th_past, Container::grid->theta(S_p_j));
    a_p_r = Container::grid->r(a_p_i);
}

void Tracker::Note() {
    monitorfile.setf(ios::left);
    monitorfile << setw(18) << time << Container::manager->t
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
        << setprecision(10) << setw(18) << a_p_max
        << setprecision(10) << setw(18) << a_p_r
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

double Tracker::rFocus() {
    int num_focus = 0;
    double r_focus = 0.;
    if (lapse_r != 0.) {
        r_focus += lapse_r;
        num_focus++;
    }
    if (rho_r != 0.) {
        r_focus += rho_r;
        num_focus++;
    }
    if (S_p_r != 0.) {
        r_focus += S_p_r;
        num_focus++;
    }
    if (num_focus > 0){
        return r_focus/num_focus;
    } else {
        return 0.;
    }
}