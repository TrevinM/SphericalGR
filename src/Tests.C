// 
//================================================
//
// Tests Routines for Einstein class
//
//================================================
//
#include "Manager.h"
#include "InData.h"
#include "Transformations.h"

//#include <iostream>
//#include <cmath>
//#include <complex>
//#include <iomanip>
//#include <vector>
//#include <limits>
//#include <stdlib.h>
//#include <stdio.h>
//#include <time.h>
//#include <fcntl.h>
//#include <string.h>
//#include <ctype.h>
//using namespace std;
//using std::cerr;
//using std::endl;
//#include <fstream>
//using std::ofstream;
//#include <cstdlib> // for exit function
// This program outputs values to a file called SLy_EOS_Tables.dat


void Manager::Test_EOS() {
    double rho_0 = 0.003499997129395651;
    double epsilon = GR::eos->cold_eps(rho_0);
    double P = GR::eos->P(rho_0, epsilon);
    double temp = GR::eos->temperature(rho_0, epsilon);
    double sound_speed = GR::eos->sound_speed(rho_0, epsilon);
    cout << "rho_0 = " << rho_0 << " epsilon = " << epsilon << " P = "
        << P << " temp = " << temp << " sound_speed = " << sound_speed << endl;
    cout << "rho_0(P) = " << GR::eos->rho_0(P) << endl;

    ofstream outdata;
    outdata.open("PWP_EOS_Tables.txt"); // opens the file
    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    outdata << std::left << setw(16) << "# rho_0 " << std::left
        << setw(16) << "P " << std::left << setw(16)
        << "epsilon " << setw(16) << "sound speed" << endl;
    for (double rho_0 = 1.0e-5; rho_0 <= 1.0e-1; rho_0 *= 1.05) {
        double epsilon = GR::eos->cold_eps(rho_0);
        double P = GR::eos->P(rho_0, epsilon);
        outdata << std::left << setw(16) << rho_0 << std::left
            << setw(16) << P << std::left << setw(16)
            << epsilon << setw(16) << GR::eos->sound_speed(rho_0, epsilon) << endl;
    }
    outdata.close();

    return;
}


void Manager::Test_Indices() {
    cout << " TESTS: testing indices for r... " << endl;
    for (int i = 0; i < N_r; i++) {
        double rl = GR::grid->r(i);
        int index = GR::grid->i_ind(rl);
        cout << " i = " << i << " r = " << rl << " found index "
            << index << " for r = " << r[index] << endl;
    }
    cout << " TESTS: testing indices for theta... " << endl;
    for (int j = 0; j < N_t; j++) {
        double thetal = GR::grid->theta(j);
        int index = GR::grid->j_ind(thetal);
        cout << " j = " << j << " theta = " << thetal << " found index "
            << index << " for theta = " << theta[index] << endl;
    }
    cout << " TESTS: testing indices for phi... " << endl;
    for (int k = 0; k < N_p; k++) {
        double phil = GR::grid->phi(k);
        int index = GR::grid->k_ind(phil);
        cout << " k = " << k << " phi = " << phil << " found index "
            << index << " for phi = " << phi[index] << endl;
    }
}

// int Einstein::Test_Gridfunction() {
//   for (int i = N_g; i < N_r; i++)    
//     for (int j = N_g; j < N_theta-N_g; j++)
//       for (int k = N_g; k < N_phi-N_g; k++) {    
// 	lapse_o[i][j][k] = exp(r[i] * r[i]) * costheta[j] * costheta[j] * cos(2.0*GR::grid->phi(k)) * cos(2.0 * GR::grid->phi(k));
//       }
//   lapse_o.fill_ghosts();

//   cout << "\n WILL WRITE FIRST DERIVS INTO shift_i and SECOND DERIVS INTO R_ii \n" << endl;

//   //  dump->dump(0.0, 0.0, 0, &lapse_o);
//   //  dump->slice(0.0, 0.0, 0, &lapse_o);

//   for (int i = N_g; i < N_r-N_g; i++) {    
//     const double rl = r[i];
//     for (int j = N_g; j < N_theta-N_g; j++) {
//       const double costhetal = costheta[j];
//       const double sinthetal = sintheta[j];
//       for (int k = N_g; k < N_phi-N_g; k++) {    
// 	const double phil = GR::grid->phi(k);
// 	//
// 	const double f_r = 2.0*rl * exp(rl*rl) * costhetal * costhetal *  cos(2.0*phil) * cos(2.0*phil);
// 	const double f_t = - 2.0 * exp(rl*rl) * costhetal * sinthetal * cos(2.0*phil) * cos(2.0*phil);
// 	const double f_p = - 4.0 * exp(rl*rl) * costhetal * costhetal * sin(2.0*phil) * cos(2.0*phil);
// 	//
// 	const double f_rr = 2.0* (1.0 + 2.0 * rl*rl) * exp(rl*rl) * costhetal *costhetal * cos(2.0*phil) * cos(2.0*phil);
// 	const double f_rt = - 4.0 * rl * exp(rl*rl) * sinthetal * costhetal * cos(2.0*phil) * cos(2.0*phil);
// 	const double f_rp = - 8.0 * rl * exp(rl*rl) * costhetal * costhetal * sin(2.0*phil) * cos(2.0*phil);
// 	const double f_tt = exp(rl*rl) * (2.0 - 4.0 * costhetal*costhetal ) * cos(2.0*phil) * cos(2.0*phil);;
// 	const double f_tp = 8.0* exp(rl*rl) * sinthetal * costhetal  * sin(2.0*phil)  * cos(2.0*phil);;
// 	const double f_pp = 4.0 * exp(rl*rl) * costhetal * costhetal * (2.0 - 4.0 * cos(2.0*phil) * cos(2.0*phil));
// 	//
// 	if (i == N_g && j == N_theta/2 && k == N_g) {
// 	  cout << " Function at r = " << r[i] << ", theta = " << GR::grid->theta(j)	
// 	       << ", phi = " << phil << " is : " << lapse_o[i][j][k] 
// 	       << " = " << lapse_o(i,j,k) << endl;
// 	  cout << " Derivatives... " << endl;
// 	  cout << " df / dr =         " << lapse_o.dr(i,j,k) << " = " << f_r << endl;
// 	  cout << " df / dr =         " << lapse_o.dr(i,j,k,+1.0) << " = " << f_r << endl;
// 	  cout << " df / dr =         " << lapse_o.dr(i,j,k,-1.0) << " = " << f_r << endl;
// 	  cout << " df / dr =         " << lapse_o.dr_so(i,j,k) << " = " << f_r << endl;
// 	  cout << " df / dtheta =     " << lapse_o.dtheta(i,j,k) << " = " << f_t << endl;
// 	  cout << " df / dphi =       " << lapse_o.dphi(i,j,k) << " = " << f_p << endl;
// 	  cout << " d^2f / dr^2 =     " << lapse_o.ddr(i,j,k) << " = " << f_rr << endl;
// 	  cout << " d^2f / dr^2 =     " << lapse_o.ddr_so(i,j,k) << " = " << f_rr << endl;
// 	  cout << " df / drdtheta =   " << lapse_o.drdtheta(i,j,k) << " = " << f_rt << endl;
// 	  cout << " df / drdphi =     " << lapse_o.drdphi(i,j,k) << " = " << f_rp << endl;
// 	  cout << " d^2f / dtheta^2 = " << lapse_o.ddtheta(i,j,k) << " = " << f_tt << endl;
// 	  cout << " d^2f / dtheta^2 = " << lapse_o.ddtheta_so(i,j,k) << " = " << f_tt << endl;
// 	  cout << " df / dthetadphi = " << lapse_o.dthetadphi(i,j,k) << " = " << f_tp << endl;
// 	  cout << " df / dthetadphi = " << lapse_o.dthetadphi_so(i,j,k) << " = " << f_tp << endl;
// 	  cout << " df / dphi^2 =     " << lapse_o.ddphi(i,j,k) << " = " << f_pp << endl;
// 	}
// 	shift_r_o[i][j][k] = lapse_o.dr(i,j,k) - f_r;
// 	shift_t_o[i][j][k] = lapse_o.dtheta(i,j,k) - f_t;
// 	shift_p_o[i][j][k] = lapse_o.dphi(i,j,k) - f_p;
// 	//
// 	R_rr[i][j][k] = lapse_o.ddr(i,j,k) - f_rr;
// 	R_rt[i][j][k] = lapse_o.drdtheta(i,j,k) - f_rt;
// 	R_rp[i][j][k] = lapse_o.drdphi(i,j,k) - f_rp;
// 	R_tt[i][j][k] = lapse_o.ddtheta(i,j,k) - f_tt;
// 	R_tp[i][j][k] = lapse_o.dthetadphi(i,j,k) - f_tp;
// 	R_pp[i][j][k] = lapse_o.ddphi(i,j,k) - f_pp;
//       }
//     }
//   }
//   dump->dump(0.0, 0.0, 0, &shift_r_o);
//   dump->dump(0.0, 0.0, 0, &shift_t_o);
//   dump->dump(0.0, 0.0, 0, &shift_p_o);
//   dump->dump(0.0, 0.0, 0, &R_rr);
//   dump->dump(0.0, 0.0, 0, &R_rt);
//   dump->dump(0.0, 0.0, 0, &R_rp);
//   dump->dump(0.0, 0.0, 0, &R_tt);
//   dump->dump(0.0, 0.0, 0, &R_tp);
//   dump->dump(0.0, 0.0, 0, &R_pp);
//   dump->slice(0.0, 0.0, 0, &shift_r_o);
//   dump->slice(0.0, 0.0, 0, &shift_t_o);
//   dump->slice(0.0, 0.0, 0, &shift_p_o);
//   dump->slice(0.0, 0.0, 0, &R_rr);
//   dump->slice(0.0, 0.0, 0, &R_rt);
//   dump->slice(0.0, 0.0, 0, &R_rp);
//   dump->slice(0.0, 0.0, 0, &R_tt);
//   dump->slice(0.0, 0.0, 0, &R_tp);
//   dump->slice(0.0, 0.0, 0, &R_pp);
// }

//
//================================================
// Test connection coefficients and Ricci tensor for 
// Schwarzschild solution in isotropic coordinates
//================================================
//
void Manager::Test_Ricci_for_Schwarzschild() {
    const double M = 1.0;
    cout << " Testing Ricci for Schwarzschild... " << endl;
    //
    // first compute h_ij for Schwarzschild solution
    //
    for (int i = N_g; i < N_r - N_g; i++)
        for (int j = N_g; j < N_t - N_g; j++)
            for (int k = N_g; k < N_p - N_g; k++) {
                const double rl = GR::grid->r(i);
                double psi = 1.0 + M / (2.0 * rl);
                double psi4 = psi * psi * psi * psi;
                last->h_rr[i][j][k] = psi4 - 1.0;
                last->h_rt[i][j][k] = 0.0;
                last->h_rp[i][j][k] = 0.0;
                last->h_tt[i][j][k] = psi4 - 1.0;
                last->h_tp[i][j][k] = 0.0;
                last->h_pp[i][j][k] = psi4 - 1.0;
            }
    //
    // Now compute Connection
    // 
    int i = N_r / 2;
    int j = N_t / 3;
    int k = N_g;
    const double rl = GR::grid->r(i);
    double r2l = rl * rl;
    double psi = 1.0 + M / (2.0 * rl);
    double psi4 = psi * psi * psi * psi;
    double dpsidr = -M / (2.0 * r2l);

    curve->Compute_Metric_Derivatives(last);
    curve->Compute_Inverse_Metric(last);

    cout << "\n Results at r = " << rl << ", theta = " << GR::grid->theta(j) << ", phi = " << GR::grid->phi(k) << endl;

    cout << "\nINVERSE METRIC : " << endl;
    cout << setw(10) << "gup^{ij} " << setw(15) << "numerical" << setw(15) << "analytical" << endl;
    cout << setw(10) << "gup^{rr} " << setw(15) << curve->gup_rr(i, j, k) << setw(15) << 1.0 / psi4 << endl;
    cout << setw(10) << "gup^{rt} " << setw(15) << curve->gup_rt(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "gup^{rp} " << setw(15) << curve->gup_rp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "gup^{tt} " << setw(15) << curve->gup_tt(i, j, k) << setw(15) << 1.0 / psi4 / r2l << endl;
    cout << setw(10) << "gup^{tp} " << setw(15) << curve->gup_tp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "gup^{pp} " << setw(15) << curve->gup_pp(i, j, k) << setw(15) << 1.0 / (psi4 * r2l * sintheta[j] * sintheta[j]) << endl;

    curve->Compute_Connection();

    cout << "\nCONNECTION : (without the flat part) " << endl;
    cout << setw(10) << "Gam^i_{jk} " << setw(15) << "numerical" << setw(15) << "analytical" << endl;
    cout << setw(10) << "Gam^r_{rr} " << setw(15) << curve->DG_r_rr(i, j, k) << setw(15) << 2.0 * dpsidr / psi << endl;
    cout << setw(10) << "Gam^r_{rt} " << setw(15) << curve->DG_r_rt(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^r_{rp} " << setw(15) << curve->DG_r_rp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^r_{tt} " << setw(15) << curve->DG_r_tt(i, j, k) << setw(15) << -2.0 * r2l * dpsidr / psi << endl;
    cout << setw(10) << "Gam^r_{tp} " << setw(15) << curve->DG_r_tp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^r_{pp} " << setw(15) << curve->DG_r_pp(i, j, k) << setw(15) << -2.0 * r2l * dpsidr / psi * sintheta[j] * sintheta[j] << endl;

    cout << setw(10) << "Gam^t_{rr} " << setw(15) << curve->DG_t_rr(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^t_{rt} " << setw(15) << curve->DG_t_rt(i, j, k) << setw(15) << 2.0 * dpsidr / psi << endl;
    cout << setw(10) << "Gam^t_{rp} " << setw(15) << curve->DG_t_rp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^t_{tt} " << setw(15) << curve->DG_t_tt(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^t_{tp} " << setw(15) << curve->DG_t_tp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^t_{pp} " << setw(15) << curve->DG_t_pp(i, j, k) << setw(15) << 0.0 << endl;

    cout << setw(10) << "Gam^p_{rr} " << setw(15) << curve->DG_p_rr(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^p_{rt} " << setw(15) << curve->DG_p_rt(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^p_{rp} " << setw(15) << curve->DG_p_rp(i, j, k) << setw(15) << 2.0 * dpsidr / psi << endl;
    cout << setw(10) << "Gam^p_{tt} " << setw(15) << curve->DG_p_tt(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^p_{tp} " << setw(15) << curve->DG_p_tp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "Gam^p_{pp} " << setw(15) << curve->DG_p_pp(i, j, k) << setw(15) << 0.0 << endl;
    //
    // Before computing Ricci we first compute Lambda's
    //
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = GR::grid->r(i);
        for (int j = N_g; j < N_t - N_g; j++)
            for (int k = N_g; k < N_p - N_g; k++) {
                double psi = 1.0 + M / (2.0 * rl);
                double psi4 = psi * psi * psi * psi;
                double dpsidr = -M / (2.0 * rl * rl);
                last->lam_r[i][j][k] = -2.0 * dpsidr / (psi4 * psi);
                last->lam_t[i][j][k] = 0.0;
                last->lam_p[i][j][k] = 0.0;
            }
    }
    curve->Compute_Ricci(last);

    i = N_r / 2;
    j = N_t / 3;
    k = N_g;
    r2l = rl * rl;
    psi = 1.0 + M / (2.0 * rl);
    psi4 = psi * psi * psi * psi;
    dpsidr = -M / (2.0 * r2l);

    cout << "\nRICCI : " << endl;
    cout << setw(10) << "R_{ij} " << setw(15) << "numerical" << setw(15) << "analytical" << endl;
    cout << setw(10) << "R_{rr} " << setw(15) << curve->R_rr(i, j, k) << setw(15) << -2.0 * M / (psi * psi * r2l * rl) << endl;
    cout << setw(10) << "R_{rt} " << setw(15) << curve->R_rt(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "R_{rp} " << setw(15) << curve->R_rp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "R_{tt} " << setw(15) << curve->R_tt(i, j, k) << setw(15) << M / (psi * psi * rl) << endl;
    cout << setw(10) << "R_{tp} " << setw(15) << curve->R_tp(i, j, k) << setw(15) << 0.0 << endl;
    cout << setw(10) << "R_{pp} " << setw(15) << curve->R_pp(i, j, k) << setw(15) << M / (psi * psi * rl) * sintheta[j] * sintheta[j] << endl;
}
// //
// //==========================================================================
// // Test set-up of conformal connection functions in InData for Schwarzschild 
// //==========================================================================
// //
// int Einstein::Test_Lambda_InData() { 
//   cout << " Testing Conformal Connection Function setup for Schwarzschild... " << endl;
//   int i = N_r/2;
//   int j = N_theta/3;
//   int k = N_phi/4;
//   double M = 1.0;

//   double psi = 1.0 + M/(2.0*r[i]);
//   double psi4 = psi*psi*psi*psi;
//   double dpsidr = - M/(2.0*r[i]*r[i]);
//   cout << "Lambda^r = " << setw(16) << lam_r_o(i,j,k) << " analytical = " << setw(16) << - 2.0 * dpsidr/(psi4*psi) << endl; 
//   cout << "Lambda^t = " << setw(16) << lam_t_o(i,j,k) << " analytical = " << setw(16) << 0.0 << endl; 
//   cout << "Lambda^p = " << setw(16) << lam_p_o(i,j,k) << " analytical = " << setw(16) << 0.0 << endl; 
// }
// //
// //==========================================================================
// // Test derivatives at interfaces
// //==========================================================================
// //
void Manager::Test_Derivatives() {
    //
    // set up function in interior GR::grid
    //
    for (int i = N_g; i < N_r; i++) {
        const double rl = last->lapse.r(i);
        for (int j = N_g; j < N_t - N_g; j++) {
            const double stl = last->lapse.sintheta(j);
            const double st2 = stl * stl;
            for (int k = N_g; k < N_p - N_g; k++)
                last->lapse[i][j][k] = exp(-rl * rl) * (2.0 - 3.0 * st2);
        }
    }
    //
    // fill ghosts...
    //
    last->lapse.fill_ghosts();
    //
    // now compute derivatives and compare with analytical solutions...
    //
    last->shift_r.equals(0.0);
    last->shift_t.equals(0.0);
    last->h_rr.equals(0.0);
    last->h_rt.equals(0.0);
    last->h_tt.equals(0.0);
    for (int i = N_g; i < N_r - N_g; i++) {
        const double rl = last->lapse.r(i);
        for (int j = N_g; j < N_t - N_g; j++) {
            const double stl = last->lapse.sintheta(j);
            const double ctl = last->lapse.costheta(j);
            for (int k = N_g; k < N_p - N_g; k++) {
                last->shift_r[i][j][k] = last->lapse.dr(i, j, k, -1.0)
                    - (-2.0 * rl * last->lapse(i, j, k));
                last->shift_t[i][j][k] = last->lapse.dtheta(i, j, k, -1.0)
                    - (-6.0 * exp(-rl * rl) * stl * ctl);
                last->h_rr[i][j][k] = last->lapse.ddr(i, j, k)
                    - ((4.0 * rl * rl - 2.0) * last->lapse(i, j, k));
                last->h_rt[i][j][k] = last->lapse.drdtheta(i, j, k)
                    - (12.0 * rl * exp(-rl * rl) * stl * ctl);
                last->h_tt[i][j][k] = last->lapse.ddtheta(i, j, k)
                    - (-6.0 * exp(-rl * rl) * (ctl * ctl - stl * stl));
            }
        }
    }
    dump->slice(0.0, 0.0, 0, last->shift_r.Address());
    dump->slice(0.0, 0.0, 0, last->shift_t.Address());
    dump->slice(0.0, 0.0, 0, last->h_rr.Address());
    dump->slice(0.0, 0.0, 0, last->h_rt.Address());
    dump->slice(0.0, 0.0, 0, last->h_tt.Address());
    dump->dump(0.0, 0.0, 0, last->h_tt.Address());
    exit(0);
}
//   h_rr_o.fill_ghosts();
//   h_rt_o.fill_ghosts();
//   h_rp_o.fill_ghosts();
//   h_tt_o.fill_ghosts();
//   h_tp_o.fill_ghosts();
//   h_pp_o.fill_ghosts();
//   phi_o.fill_ghosts();

//   double dr = 1.e-3;
//   double dt = 1.e-3;
//   double dp = 1.e-3;
//   cout << " Testing derivatives ... " << endl;

//   cout << "\n Testing at CENTER \n" << endl;

//   int i = N_r/2;
//   int j = N_theta/2;
//   int k = N_phi/2;

//   i = 2;
//   double rl = r[i];
//   double tl = GR::grid->theta(j);
//   double pl = GR::grid->phi(k);  

//   cout << " at r = " << rl << ", theta = " << tl << ", phi = " << pl << endl;
//   cout << " at x = " << rl*sin(tl)*cos(pl) << ", y = " << rl*sin(tl)*sin(pl) << ", z = " << rl*cos(tl) << endl << endl;

//   //==========================================================================
//   // h_rr
//   //==========================================================================
//   cout << "\n h_rr \n" << endl;

//   double fct_l = GR::indata->h_rr_analytical(rl,tl,pl,t);
//   double fctpr = GR::indata->h_rr_analytical(rl+dr,tl,pl,t);
//   double fctmr = GR::indata->h_rr_analytical(rl-dr,tl,pl,t);
//   double fctpt = GR::indata->h_rr_analytical(rl,tl+dt,pl,t);
//   double fctmt = GR::indata->h_rr_analytical(rl,tl-dt,pl,t);
//   double fctpp = GR::indata->h_rr_analytical(rl,tl,pl+dp,t);
//   double fctmp = GR::indata->h_rr_analytical(rl,tl,pl-dp,t);

//   double fctprpt = GR::indata->h_rr_analytical(rl+dr,tl+dt,pl,t);
//   double fctprmt = GR::indata->h_rr_analytical(rl+dr,tl-dt,pl,t);
//   double fctmrpt = GR::indata->h_rr_analytical(rl-dr,tl+dt,pl,t);
//   double fctmrmt = GR::indata->h_rr_analytical(rl-dr,tl-dt,pl,t);

//   double fctprpp = GR::indata->h_rr_analytical(rl+dr,tl,pl+dp,t);
//   double fctprmp = GR::indata->h_rr_analytical(rl+dr,tl,pl-dp,t);
//   double fctmrpp = GR::indata->h_rr_analytical(rl-dr,tl,pl+dp,t);
//   double fctmrmp = GR::indata->h_rr_analytical(rl-dr,tl,pl-dp,t);

//   double fctptpp = GR::indata->h_rr_analytical(rl,tl+dt,pl+dp,t);
//   double fctptmp = GR::indata->h_rr_analytical(rl,tl+dt,pl-dp,t);
//   double fctmtpp = GR::indata->h_rr_analytical(rl,tl-dt,pl+dp,t);
//   double fctmtmp = GR::indata->h_rr_analytical(rl,tl-dt,pl-dp,t);


//   cout << setprecision(8);

//   const double t = 0.0;
//   double f_r  = (fctpr - fctmr)/(2.0*dr);
//   double f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   double f_t  = (fctpt - fctmt)/(2.0*dt);
//   double f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   double f_p  = (fctpp - fctmp)/(2.0*dp);
//   double f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   double f_rt = (fctprpt - fctprmt - fctmrpt + fctmrmt)/(4.0*dr*dt);
//   double f_rp = (fctprpp - fctprmp - fctmrpp + fctmrmp)/(4.0*dr*dp);
//   double f_tp = (fctptpp - fctptmp - fctmtpp + fctmtmp)/(4.0*dt*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_rr_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_rr_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_rr_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_rr_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_rr_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_rr_o.ddphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/drdt" << setw(16) << f_rt << setw(16) << h_rr_o.drdtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/drdp" << setw(16) << f_rp << setw(16) << h_rr_o.drdphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dtdp" << setw(16) << f_tp << setw(16) << h_rr_o.dthetadphi(i,j,k) << endl;


//   //==========================================================================
//   // h_rt
//   //==========================================================================
//   cout << "\n h_rt \n" << endl;

//   fct_l = GR::indata->h_rt_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_rt_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_rt_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_rt_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_rt_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_rt_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_rt_analytical(rl,tl,pl-dp,t);

//   fctprpt = GR::indata->h_rt_analytical(rl+dr,tl+dt,pl,t);
//   fctprmt = GR::indata->h_rt_analytical(rl+dr,tl-dt,pl,t);
//   fctmrpt = GR::indata->h_rt_analytical(rl-dr,tl+dt,pl,t);
//   fctmrmt = GR::indata->h_rt_analytical(rl-dr,tl-dt,pl,t);

//   fctprpp = GR::indata->h_rt_analytical(rl+dr,tl,pl+dp,t);
//   fctprmp = GR::indata->h_rt_analytical(rl+dr,tl,pl-dp,t);
//   fctmrpp = GR::indata->h_rt_analytical(rl-dr,tl,pl+dp,t);
//   fctmrmp = GR::indata->h_rt_analytical(rl-dr,tl,pl-dp,t);

//   fctptpp = GR::indata->h_rt_analytical(rl,tl+dt,pl+dp,t);
//   fctptmp = GR::indata->h_rt_analytical(rl,tl+dt,pl-dp,t);
//   fctmtpp = GR::indata->h_rt_analytical(rl,tl-dt,pl+dp,t);
//   fctmtmp = GR::indata->h_rt_analytical(rl,tl-dt,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   f_rt = (fctprpt - fctprmt - fctmrpt + fctmrmt)/(4.0*dr*dt);
//   f_rp = (fctprpp - fctprmp - fctmrpp + fctmrmp)/(4.0*dr*dp);
//   f_tp = (fctptpp - fctptmp - fctmtpp + fctmtmp)/(4.0*dt*dp);


//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_rt_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_rt_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_rt_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_rt_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_rt_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_rt_o.ddphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/drdt" << setw(16) << f_rt << setw(16) << h_rt_o.drdtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/drdp" << setw(16) << f_rp << setw(16) << h_rt_o.drdphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dtdp" << setw(16) << f_tp << setw(16) << h_rt_o.dthetadphi(i,j,k) << endl;

//   //==========================================================================
//   // h_rp
//   //==========================================================================
//   cout << "\n h_rp \n" << endl;

//   fct_l = GR::indata->h_rp_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_rp_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_rp_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_rp_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_rp_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_rp_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_rp_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_rp_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_rp_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_rp_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_rp_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_rp_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_rp_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // h_tt
//   //==========================================================================
//   cout << "\n h_tt \n" << endl;

//   fct_l = GR::indata->h_tt_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_tt_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_tt_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_tt_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_tt_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_tt_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_tt_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_tt_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_tt_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_tt_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_tt_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_tt_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_tt_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // h_tp
//   //==========================================================================
//   cout << "\n h_tp \n" << endl;

//   fct_l = GR::indata->h_tp_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_tp_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_tp_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_tp_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_tp_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_tp_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_tp_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_tp_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_tp_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_tp_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_tp_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_tp_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_tp_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // h_pp
//   //==========================================================================
//   cout << "\n h_pp \n" << endl;

//   fct_l = GR::indata->h_pp_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_pp_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_pp_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_pp_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_pp_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_pp_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_pp_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_pp_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_pp_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_pp_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_pp_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_pp_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_pp_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // phi
//   //==========================================================================
//   cout << "\n phi \n" << endl;

//   fct_l = GR::indata->phi_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->phi_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->phi_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->phi_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->phi_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->phi_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->phi_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << phi_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << phi_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << phi_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << phi_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << phi_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << phi_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // Now test on axis
//   //==========================================================================
//   cout << "\n Testing on AXIS \n" << endl;

//   i = N_r/2;
//   j = N_theta/2;
//   k = N_phi/2;

//   j = 2;
//   rl = r[i];
//   tl = GR::grid->theta(j);
//   pl = GR::grid->phi(k);  


//   cout << " at r = " << rl << ", theta = " << tl << ", phi = " << pl << endl;
//   cout << " at x = " << rl*sin(tl)*cos(pl) << ", y = " << rl*sin(tl)*sin(pl) << ", z = " << rl*cos(tl) << endl << endl;

//   //==========================================================================
//   // h_rr
//   //==========================================================================
//   cout << "\n h_rr \n" << endl;

//   fct_l = GR::indata->h_rr_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_rr_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_rr_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_rr_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_rr_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_rr_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_rr_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_rr_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_rr_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_rr_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_rr_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_rr_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_rr_o.ddphi(i,j,k) << endl;


//   //==========================================================================
//   // h_rt
//   //==========================================================================
//   cout << "\n h_rt \n" << endl;

//   fct_l = GR::indata->h_rt_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_rt_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_rt_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_rt_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_rt_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_rt_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_rt_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_rt_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_rt_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_rt_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_rt_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_rt_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_rt_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // h_rp
//   //==========================================================================
//   cout << "\n h_rp \n" << endl;

//   fct_l = GR::indata->h_rp_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_rp_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_rp_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_rp_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_rp_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_rp_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_rp_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_rp_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_rp_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_rp_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_rp_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_rp_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_rp_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // h_tt
//   //==========================================================================
//   cout << "\n h_tt \n" << endl;

//   fct_l = GR::indata->h_tt_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_tt_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_tt_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_tt_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_tt_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_tt_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_tt_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_tt_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_tt_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_tt_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_tt_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_tt_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_tt_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // h_tp
//   //==========================================================================
//   cout << "\n h_tp \n" << endl;

//   fct_l = GR::indata->h_tp_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_tp_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_tp_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_tp_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_tp_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_tp_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_tp_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_tp_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_tp_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_tp_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_tp_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_tp_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_tp_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // h_pp
//   //==========================================================================
//   cout << "\n h_pp \n" << endl;

//   fct_l = GR::indata->h_pp_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->h_pp_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->h_pp_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->h_pp_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->h_pp_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->h_pp_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->h_pp_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << h_pp_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << h_pp_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << h_pp_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << h_pp_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << h_pp_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << h_pp_o.ddphi(i,j,k) << endl;

//   //==========================================================================
//   // phi
//   //==========================================================================
//   cout << "\n phi \n" << endl;

//   fct_l = GR::indata->phi_analytical(rl,tl,pl,t);
//   fctpr = GR::indata->phi_analytical(rl+dr,tl,pl,t);
//   fctmr = GR::indata->phi_analytical(rl-dr,tl,pl,t);
//   fctpt = GR::indata->phi_analytical(rl,tl+dt,pl,t);
//   fctmt = GR::indata->phi_analytical(rl,tl-dt,pl,t);
//   fctpp = GR::indata->phi_analytical(rl,tl,pl+dp,t);
//   fctmp = GR::indata->phi_analytical(rl,tl,pl-dp,t);


//   cout << setprecision(8);

//   f_r  = (fctpr - fctmr)/(2.0*dr);
//   f_rr = (fctpr - 2.0*fct_l + fctmr)/(dr*dr);
//   f_t  = (fctpt - fctmt)/(2.0*dt);
//   f_tt = (fctpt - 2.0*fct_l + fctmt)/(dt*dt);
//   f_p  = (fctpp - fctmp)/(2.0*dp);
//   f_pp = (fctpp - 2.0*fct_l + fctmp)/(dp*dp);

//   cout << setw(10) << "df/dr" << setw(16) << f_r << setw(16) << phi_o.dr(i,j,k) << endl;
//   cout << setw(10) << "d2f/dr2" << setw(16) << f_rr << setw(16) << phi_o.ddr(i,j,k) << endl;
//   cout << setw(10) << "df/dt" << setw(16) << f_t << setw(16) << phi_o.dtheta(i,j,k) << endl;
//   cout << setw(10) << "d2f/dt2" << setw(16) << f_tt << setw(16) << phi_o.ddtheta(i,j,k) << endl;
//   cout << setw(10) << "df/dp" << setw(16) << f_p << setw(16) << phi_o.dphi(i,j,k) << endl;
//   cout << setw(10) << "d2f/dp2" << setw(16) << f_pp << setw(16) << phi_o.ddphi(i,j,k) << endl;

// };

// Doub Vx(Doub x, Doub y, Doub z) { return x*x; }
// Doub Vy(Doub x, Doub y, Doub z) { return y; }
// Doub Vz(Doub x, Doub y, Doub z) { return 1.0; }

// Doub Vr(Doub r, Doub t, Doub p) { 
//   const Doub x = r * sin(t) * cos(p);
//   const Doub y = r * sin(t) * sin(p);
//   const Doub z = r * cos(t);
//   vect V_cart(Vx(x,y,z),Vy(x,y,z),Vz(x,y,z));
//   vect V_sc;
//   V_sc = Cartesian_to_Spherical_upper(V_cart,r,t,p);
//   return V_sc[0];
// }

// Doub Vt(Doub r, Doub t, Doub p) { 
//   const Doub x = r * sin(t) * cos(p);
//   const Doub y = r * sin(t) * sin(p);
//   const Doub z = r * cos(t);
//   vect V_cart(Vx(x,y,z),Vy(x,y,z),Vz(x,y,z));
//   vect V_sc;
//   V_sc = Cartesian_to_Spherical_upper(V_cart,r,t,p);
//   return V_sc[1] * r;
// }

// Doub Vp(Doub r, Doub t, Doub p) { 
//   const Doub x = r * sin(t) * cos(p);
//   const Doub y = r * sin(t) * sin(p);
//   const Doub z = r * cos(t);
//   vect V_cart(Vx(x,y,z),Vy(x,y,z),Vz(x,y,z));
//   vect V_sc;
//   V_sc = Cartesian_to_Spherical_upper(V_cart,r,t,p);
//   return V_sc[2] * r * sin(t);
// }





// int Einstein::Test_Transformation_Vector() {
//   //
//   // Initialize vector
//   //
//   // Doub r = 1.0;
//   // Doub theta = 1.3;
//   // Doub phi = 4.5;
//   // Doub st = sin(theta);
//   // Doub ct = cos(theta);
//   // Doub sp = sin(phi);
//   // Doub cp = cos(phi);
//   // Doub x = r * st * cp;
//   // Doub y = r * st * sp;
//   // Doub z = r * ct;
//   // cout << " x = " << x << " y = " << y << " z = " << z << endl; 
//   // vect V_cart(Vx(x,y,z),Vy(x,y,z),Vz(x,y,z));
//   // vect V_sc;
//   // V_sc = Cartesian_to_Spherical_upper(V_cart,r,theta,phi);
//   // V_cart.print();
//   // V_sc.print();

//   int i_check = 5;
//   int j_check = 14;
//   int k_check = 12;
//   for (int i = N_g; i < N_r-N_g; i++)    
//     for (int j = N_g; j < N_theta-N_g; j++)
//       for (int k = N_g; k < N_phi-N_g; k++) {
//   	bool print = (i == 1 && j == j_check && k == k_check);
//   	// create vector (upper indices) 
//   	lam_r_o[i][j][k] = Vr(r[i],GR::grid->theta(j),GR::grid->phi(k));
//   	lam_t_o[i][j][k] = Vt(r[i],GR::grid->theta(j),GR::grid->phi(k));
//   	lam_p_o[i][j][k] = Vp(r[i],GR::grid->theta(j),GR::grid->phi(k));
//       }
//   lam_r_o.fill_ghosts();
//   lam_t_o.fill_ghosts();
//   lam_p_o.fill_ghosts();

//   const Doub dr = 1.e-6;
//   const Doub dt = 1.e-6;
//   const Doub dp = 1.e-6;

//   //
//   // Center
//   // 

//   int i = 1;
//   int j = j_check;
//   int k = k_check;

//   cout << endl << " Testing at CENTER " << endl << endl;
//   cout << " r = " << r[i] << " theta = " << GR::grid->theta(j) << " phi = " << GR::grid->phi(k) << endl;
//   Doub rl = r[i];
//   Doub tl = GR::grid->theta(j);
//   Doub pl = GR::grid->phi(k);
//   Doub Vrr_ana = (Vr(rl+dr,tl,pl) - Vr(rl-dr,tl,pl))/(2.0*dr);
//   Doub Vtr_ana = (Vt(rl+dr,tl,pl) - Vt(rl-dr,tl,pl))/(2.0*dr);
//   Doub Vpr_ana = (Vp(rl+dr,tl,pl) - Vp(rl-dr,tl,pl))/(2.0*dr);
//   cout << " V^r_{,r} = " << setw(16) << lam_r_o.dr(i,j,k) << " analytically: " 
//        << setw(16) << Vrr_ana << endl;
//   cout << " V^t_{,r} = " << setw(16) << lam_t_o.dr(i,j,k) << " analytically: " 
//        << setw(16) << Vtr_ana << endl;
//   cout << " V^p_{,r} = " << setw(16) << lam_p_o.dr(i,j,k) << " analytically: " 
//        << setw(16) << Vpr_ana << endl;

//   //
//   // Axis
//   // 

//   i = i_check;
//   j = 1;
//   k = k_check;

//   cout << endl << " Testing at AXIS " << endl << endl;
//   cout << " r = " << r[i] << " theta = " << GR::grid->theta(j) << " phi = " << GR::grid->phi(k) << endl;
//   rl = r[i];
//   tl = GR::grid->theta(j);
//   pl = GR::grid->phi(k);
//   Doub Vrt_ana = (Vr(rl,tl+dt,pl) - Vr(rl,tl-dt,pl))/(2.0*dt);
//   Doub Vtt_ana = (Vt(rl,tl+dt,pl) - Vt(rl,tl-dt,pl))/(2.0*dt);
//   Doub Vpt_ana = (Vp(rl,tl+dt,pl) - Vp(rl,tl-dt,pl))/(2.0*dt);
//   cout << " V^r_{,t} = " << setw(16) << lam_r_o.dtheta(i,j,k) << " analytically: " 
//        << setw(16) << Vrt_ana << endl;
//   cout << " V^t_{,t} = " << setw(16) << lam_t_o.dtheta(i,j,k) << " analytically: " 
//        << setw(16) << Vtt_ana << endl;
//   cout << " V^p_{,t} = " << setw(16) << lam_p_o.dtheta(i,j,k) << " analytically: " 
//        << setw(16) << Vpt_ana << endl;


//   //
//   // periodicity
//   // 

//   i = i_check;
//   j = j_check;
//   k = 1;

//   cout << endl << " Testing PERIODICITY " << endl << endl;
//   cout << " r = " << r[i] << " theta = " << GR::grid->theta(j) << " phi = " << GR::grid->phi(k) << endl;
//   rl = r[i];
//   tl = GR::grid->theta(j);
//   pl = GR::grid->phi(k);
//   Doub Vrp_ana = (Vr(rl,tl,pl+dp) - Vr(rl,tl,pl-dp))/(2.0*dp);
//   Doub Vtp_ana = (Vt(rl,tl,pl+dp) - Vt(rl,tl,pl-dp))/(2.0*dp);
//   Doub Vpp_ana = (Vp(rl,tl,pl+dp) - Vp(rl,tl,pl-dp))/(2.0*dp);
//   cout << " V^r_{,p} = " << setw(16) << lam_r_o.dphi(i,j,k) << " analytically: " 
//        << setw(16) << Vrp_ana << endl;
//   cout << " V^t_{,p} = " << setw(16) << lam_t_o.dphi(i,j,k) << " analytically: " 
//        << setw(16) << Vtp_ana << endl;
//   cout << " V^p_{,p} = " << setw(16) << lam_p_o.dphi(i,j,k) << " analytically: " 
//        << setw(16) << Vpp_ana << endl;



//
//=================================================================
// Test flat metric in funky coordinate system
//=================================================================
//
void Manager::Test_Flat_Metric() {
    cout << " Testing flat metric... " << endl;
    int i_check = (N_r - 4) / 2 + 2;
    int j_check = (N_t - 4) / 4 + 2;
    int k_check = (N_p - 4) / 4 + 2;

    const int i = i_check;
    const int j = j_check;
    const int k = k_check;

    const double rl = GR::grid->r(i);
    cout << " Results at r = " << rl << ", theta = " << GR::grid->theta(j) << " and phi = " << GR::grid->phi(k) << endl;

    const double st = GR::grid->sintheta(j);
    const double ct = GR::grid->costheta(j);
    const double sp = sin(GR::grid->phi(k));
    const double cp = cos(GR::grid->phi(k));

    const double A = 0.1;  // make sure these agree with paramters defined in InData.h
    const double B = 0.3;
    const double C = 0.5;
    const double x = rl * st * cp;
    const double y = rl * st * sp;
    const double z = rl * ct;
    const double f = 1.0 + A * x * x + 0.2;
    const double g = 1.0 + B * y * y + 0.5;
    const double h = 1.0 + C * z * z + 0.1;

    // define transformations
    tensor part_x_part_r(st * cp, st * sp, ct,
        rl * ct * cp, rl * ct * sp, -rl * st,
        -rl * st * sp, rl * st * cp, 0.0);
    tensor part_r_part_x(st * cp, st * sp, ct,
        ct * cp / rl, ct * sp / rl, -st / rl,
        -sp / (rl * st), cp / (rl * st), 0.0);

    cout << " h_ij: " << endl;
    tensor h_ij(last->h_rr(i, j, k), last->h_rt(i, j, k), last->h_rp(i, j, k),
        last->h_tt(i, j, k), last->h_tp(i, j, k), last->h_pp(i, j, k));
    h_ij.print();

    cout << " lam^i : " << endl;
    vect lam(last->lam_r(i, j, k), last->lam_t(i, j, k), last->lam_p(i, j, k));
    lam.print();
    //
    // cartesian components of vector Lambda^i
    //
    double Gamx = A * x / (f * f);
    double Gamy = B * y / (g * g);
    double Gamz = C * z / (h * h);
    vect lam_ana_cart(Gamx, Gamy, Gamz);
    // now transform to spherical coordinates
    vect lam_ana;
    for (int a = 0; a < 3; a++) {
        lam_ana[a] = 0.0;
        for (int b = 0; b < 3; b++)
            lam_ana[a] += part_r_part_x[a][b] * lam_ana_cart[b];
    }
    // properly rescale lamba:
    lam_ana[1] *= rl;
    lam_ana[2] *= rl * st;
    cout << " lam^i analytically : " << endl;
    lam_ana.print();
    //
    // derivatives of lambda
    //
    tensor D_lam_cart(A / (f * f) - 4.0 * A * A * x * x / (f * f * f), 0.0, 0.0, B / (g * g) - 4.0 * B * B * y * y / (g * g * g),
        0.0, C / (h * h) - 4.0 * C * C * z * z / (h * h * h));
    tensor D_lam_sc;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            D_lam_sc[a][b] = 0.0;
            for (int d = 0; d < 3; d++)
                for (int e = 0; e < 3; e++)
                    D_lam_sc[a][b] += part_x_part_r[a][d] * part_r_part_x[b][e] * D_lam_cart[d][e];
        }
    cout << " D_i Lam^j analytically: " << endl;
    D_lam_sc.print();


    cout << "==============================================" << endl;
    cout << "Inverse metric" << endl;
    tensor g_cart_up(1.0 / f, 0.0, 0.0, 1.0 / g, 0.0, 1.0 / h);
    tensor g_up_ana;
    // Now transform...
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            g_up_ana[a][b] = 0.0;
            for (int d = 0; d < 3; d++)
                for (int e = 0; e < 3; e++)
                    g_up_ana[a][b] += part_r_part_x[a][d] * part_r_part_x[b][e] * g_cart_up[d][e];
        }
    tensor g_up(curve->gup_rr(i, j, k), curve->gup_rt(i, j, k), curve->gup_rp(i, j, k),
        curve->gup_tt(i, j, k), curve->gup_tp(i, j, k), curve->gup_pp(i, j, k));
    cout << " from code: " << endl;
    g_up.print();
    cout << " analytically: " << endl;
    g_up_ana.print();
    cout << "==============================================" << endl;
    cout << " Connection functions..." << endl;
    double gamxxx = A * x / f;
    double gamyyy = B * y / g;
    double gamzzz = C * z / h;
    rank3tens Gam_cart(gamxxx, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, gamyyy, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, gamzzz);
    rank3tens Gam_sc;
    // Now transform...
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int c = 0; c < 3; c++) {
                Gam_sc[a][b][c] = 0.0;
                for (int d = 0; d < 3; d++)
                    for (int e = 0; e < 3; e++)
                        for (int f = 0; f < 3; f++)
                            Gam_sc[a][b][c] += part_r_part_x[a][d] * part_x_part_r[b][e] * part_x_part_r[c][f] * Gam_cart[d][e][f];
            }

    // double test = 0.5 * ( curve->gup_rr(i,j,k) * (2.0 * curve->Dt_e_rt(i,j,k) - curve->Dr_e_tt(i,j,k)) + 
    // 			curve->gup_rt(i,j,k) * curve->Dt_e_tt(i,j,k) +
    // 			curve->gup_rp(i,j,k) * (2.0 * curve->Dt_e_tp(i,j,k) - curve->Dp_e_tt(i,j,k)) );
    cout << "Gam^r_rr : " << setw(16) << curve->DG_r_rr(i, j, k) << setw(16) << Gam_sc[0][0][0] << endl;
    cout << "Gam^r_rt : " << setw(16) << curve->DG_r_rt(i, j, k) << setw(16) << Gam_sc[0][0][1] << endl;
    cout << "Gam^r_rp : " << setw(16) << curve->DG_r_rp(i, j, k) << setw(16) << Gam_sc[0][0][2] << endl;
    cout << "Gam^r_tt : " << setw(16) << curve->DG_r_tt(i, j, k) << setw(16) << Gam_sc[0][1][1] << endl;
    cout << "Gam^r_tp : " << setw(16) << curve->DG_r_tp(i, j, k) << setw(16) << Gam_sc[0][1][2] << endl;
    cout << "Gam^r_pp : " << setw(16) << curve->DG_r_pp(i, j, k) << setw(16) << Gam_sc[0][2][2] << endl;

    cout << "Gam^t_rr : " << setw(16) << curve->DG_t_rr(i, j, k) << setw(16) << Gam_sc[1][0][0] << endl;
    cout << "Gam^t_rt : " << setw(16) << curve->DG_t_rt(i, j, k) << setw(16) << Gam_sc[1][0][1] << endl;
    cout << "Gam^t_rp : " << setw(16) << curve->DG_t_rp(i, j, k) << setw(16) << Gam_sc[1][0][2] << endl;
    cout << "Gam^t_tt : " << setw(16) << curve->DG_t_tt(i, j, k) << setw(16) << Gam_sc[1][1][1] << endl;
    cout << "Gam^t_tp : " << setw(16) << curve->DG_t_tp(i, j, k) << setw(16) << Gam_sc[1][1][2] << endl;
    cout << "Gam^t_pp : " << setw(16) << curve->DG_t_pp(i, j, k) << setw(16) << Gam_sc[1][2][2] << endl;

    cout << "Gam^p_rr : " << setw(16) << curve->DG_p_rr(i, j, k) << setw(16) << Gam_sc[2][0][0] << endl;
    cout << "Gam^p_rt : " << setw(16) << curve->DG_p_rt(i, j, k) << setw(16) << Gam_sc[2][0][1] << endl;
    cout << "Gam^p_rp : " << setw(16) << curve->DG_p_rp(i, j, k) << setw(16) << Gam_sc[2][0][2] << endl;
    cout << "Gam^p_tt : " << setw(16) << curve->DG_p_tt(i, j, k) << setw(16) << Gam_sc[2][1][1] << endl;
    cout << "Gam^p_tp : " << setw(16) << curve->DG_p_tp(i, j, k) << setw(16) << Gam_sc[2][1][2] << endl;
    cout << "Gam^p_pp : " << setw(16) << curve->DG_p_pp(i, j, k) << setw(16) << Gam_sc[2][2][2] << endl;

    cout << "==============================================" << endl;
    tensor R(curve->R_rr(i, j, k), curve->R_rt(i, j, k), curve->R_rp(i, j, k),
        curve->R_tt(i, j, k), curve->R_tp(i, j, k), curve->R_pp(i, j, k));
    cout << "\nRicci:\n" << endl;
    R.print();

    cout << " Norms of Ricci : "
        << setw(12) << curve->R_rr.L2_norm()
        << setw(12) << curve->R_rt.L2_norm()
        << setw(12) << curve->R_rp.L2_norm()
        << setw(12) << curve->R_tt.L2_norm()
        << setw(12) << curve->R_tp.L2_norm()
        << setw(12) << curve->R_pp.L2_norm()
        << endl;

};

// int Einstein::Test_DivShift() {
//   //
//   // sanity??
//   //
//   cout << " Did you initialize with flat conformal metric?? " << endl;
//   int i_check = (N_r-4)/2 + 2;
//   int j_check = (N_theta-4)/2 + 2;
//   int k_check = (N_phi-4)/2 + 2;

//   for (int i = N_g; i < N_r-N_g; i++)    
//     for (int j = N_g; j < N_theta-N_g; j++)
//       for (int k = N_g; k < N_phi-N_g; k++) {
// 	//	bool print = (i == 1 && j == j_check && k == k_check);
//   	// create vector (upper indices)  (already rescaled!)
//   	shift_r_o[i][j][k] = Vr(rl,GR::grid->theta(j),GR::grid->phi(k));
//   	shift_t_o[i][j][k] = Vt(rl,GR::grid->theta(j),GR::grid->phi(k));
//   	shift_p_o[i][j][k] = Vp(rl,GR::grid->theta(j),GR::grid->phi(k));
//       }
//   DivShift(shift_r_o,shift_t_o,shift_p_o);
//   int i = i_check;
//   int j = j_check;
//   int k = k_check;

//   cout << " Results at r = " << rl << ", theta = " << GR::grid->theta(j) << " and phi = " << GR::grid->phi(k) << endl;
//   cout << " Determinant : " << det(i_check,j_check,k_check) << endl;
//   cout << " Divergence of shift : " << div_shift(i_check,j_check,k_check) << endl;
//   const double x = r[i_check] * sintheta[j_check] * cos(GR::grid->phi(k_check));
//   cout << " Analytically        : " << 2.0 * x + 1.0 << endl;
// };

// Doub Polynomial(Doub r, Doub A, Doub B, Doub C) {
//   return 1.0 + A * r + B * r*r + C * r*r*r;
// }

// int Einstein::Test_Interpolation() {
//   // choose coefficients...
//   Doub A = 0.3;
//   Doub B = 4.0;
//   Doub C = 60.0;
//   for (int i = N_g; i < N_r; i++)    
//     for (int j = N_g; j < N_theta-N_g; j++)
//       for (int k = N_g; k < N_phi-N_g; k++) {
// 	Doub r_l = r[i];
// 	phi_o[i][j][k] = Polynomial(r_l,A,B,C);
//       }
//   int j = N_g;
//   int k = N_g;
//   for (Doub r_test = 0.1; r_test < 0.7; r_test += 0.05) {
//     cout << " r = " << setw(10) << r_test << " Analytic = " << setw(16) << Polynomial(r_test,A,B,C) 
// 	 << " Numerical = " << setw(16) << phi_o(r_test,j,k) << endl;
//   }
// };

// int Einstein::Test_Regrid() {
//   // choose coefficients...
//   Doub A = 0.3;
//   Doub B = 4.0;
//   Doub C = 60.0;
//   for (int i = 0; i < N_r; i++)    
//     for (int j = 0; j < N_theta; j++)
//       for (int k = 0; k < N_phi; k++) {
// 	Doub r_l = r[i];
// 	shift_t_o[i][j][k] = Polynomial(r_l,A,B,C);
//       }
//   Doub r_max = GR::grid->r_max();
//   cout << " TEST_REGRID: initial r_max = " << r_max << endl;
//   int j = N_g;
//   int k = N_g;
//   for (int i = 0; i < N_r; i++) {
//     Doub r_l = r[i];
//     cout << i << " r = " << r_l << " fct = " << shift_t_o(i,j,k) << "   " <<  Polynomial(r_l,A,B,C) << endl;
//   }
//   shift_t_o.fill_ghosts();
//   for (int i = 0; i < N_r; i++) {
//     Doub r_l = r[i];
//     cout << i << " r = " << r_l << " fct = " << shift_t_o(i,j,k) << "   " <<  Polynomial(r_l,A,B,C) << endl;
//   }
//   Doub t_max = 10.;
//   Regrid(t_max);
//   for (int i = 0; i < N_r; i++) {
//     Doub r_l = r[i];
//     cout << i << " r = " << r_l << " fct = " << shift_t_o(i,j,k) << "   " <<  Polynomial(r_l,A,B,C) << endl;
//   }
// }

