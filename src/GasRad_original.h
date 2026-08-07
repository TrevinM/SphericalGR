// Tell emacs that this is -*-c++-*- mode
//================================================
// Derived from EOS.h: Equation of state for mixture of
// gas and radiation pressure
//================================================
#ifndef GASRAD_H
#define GASRAD_H

#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

//================================================
// 
// gas & radiation class derived from EOS
//
// Notes:
//
// - measure T in units of m_B, i.e. tau \equiv k_B T / m_B
// - measure entropy s in units of k_B
// - rescale everything else with respect to
//        K = a/3 ( 3s/(4 m_B a) )^{4/3}
//   (has units of length^{2/3})
// - radiation constant:
//        alpha = K^3 a m_B^4 / k_B^4
// - then entropy:
//        s_r = 4 alpha tau^3 / (3 rho_0)
//        s_g = ln( gamma tau^3 / rho_0^2 ) + 5
//   where
//        gamma = m_e^{3/2) m_B^{13/2} / 2 pi^{3} hbar^6
//              = 32 pi^3 m_e^{3/2) m_B^{13/2} / h^6
//   has units of length^{-4}
//================================================
//
class gas_radiation : public EOS {
private:
    double beta = 0.0;
    double K = 0.0;
    const double m_e_cm = 6.764e-56;  // see MTW, in units of cm
    const double m_B_cm = 1.2419e-52;
    const double PI = acos(-1.0);
    const double h_cm = 2. * PI * 2.612e-66;  // h in units of cm^2 (NOT hbar...)
    double s_init = 0.0;  // initial entroy s_bar (in units of k_B)
    const double eta = 1.0;  // number of particles contributing to radiation
    double alpha = 0.0;  // radiation constant in code units
    const double tol = 1.e-12;
    double gamma = 0.0;  // constant for gas entropy - see Notes above
    double tau_guess = 1.e-7;
    double rho_0_guess = 1e-4;
public:
    //===============================================================
    // Constructor
    //===============================================================
    gas_radiation(char* input_file) : EOS(0.0) {
        eos_type = gasrad;
        ifstream infile;
        infile.open(input_file);
        if (!infile) {
            cerr << " GASRAD: Can't open " << input_file
                << " for input. This is bad. " << endl;
        } else
            cout << " GASRAD: Reading parameters for gas & radiation from file "
            << input_file << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> beta;
        cout << " GASRAD: Setting up gas & radiation with beta = " << beta << endl;
        s_init = 8. / beta;
        // define constant K:
        const double factor1 = pow(3. / 4., 4. / 3.) * pow(15., 1. / 3.) / (6. * pow(PI, 5. / 3.));
        const double factor2 = h_cm / pow(m_B_cm, 4. / 3.);
        K = factor1 * factor2 * pow(s_init, 4. / 3.);
        name << "for gas and radiation mixture with initial beta = "
            << beta << " and K = " << K << " cm^{2/3} ";
        cout << " GASRAD: K = " << K << " cm^{2/3}" << endl;
        // define alpha \equiv a m_B^4 / k_B^4 (see Notes at top of file)
        const double m_o_h = m_B_cm / h_cm;
        const double m_o_h3 = m_o_h * m_o_h * m_o_h;
        const double alpha_cm = 8. * PI * PI * PI * PI * PI * m_o_h3 * m_B_cm / 15.;
        cout << " GASRAD: alpha = " << alpha_cm << " cm^{-2}" << endl;
        // and convert to code units and include eta:
        alpha = eta * K * K * K * alpha_cm;
        cout << " GASRAD: alpha = " << alpha << " (code units)" << endl;
        // define gamma (see Notes at top of file)
        const double m_o_h6 = m_o_h3 * m_o_h3;
        const double gamma_cm = 32. * PI * PI * PI * pow(m_e_cm, 3. / 2.) * sqrt(m_B_cm) * m_o_h6;
        gamma = K * K * K * K * K * K * gamma_cm;
        cout << " GASRAD: gamma = " << gamma << " (code units)" << endl;
        // const double P_test = 0.01;
        // cout << " Now run test.." << endl;
        // Test(P_test);
    }
    //===============================================================
    // Destructor
    //===============================================================
    ~gas_radiation() {};
    //===============================================================
    //
    // find pressure two ways:
    // first way (with two arguments): analytical solution
    // second way (with three arguments): iterative solution
    //
    //===============================================================
    double P(double rho_0, double epsilon) {
        if (rho_0 <= 0.0) return 0.0;
        //
        // first find tau by solving epsilon equation analytically,
        // then compute P from rho_0 and tau
        //
        const double tau = temperature(rho_0, epsilon);
        // const double check = epsilon - eps(rho_0, tau);
        // cout << " GASRAD: check solution in GasRad::P(rho_0, epsilon) = "
        // 	 << check << " tau = " << tau << endl;
        //
        // in case accuracy is not high enough follow up with iteration (see
        // NOTE in routine temperature(rho_0, epsilon) below)
        //
        // CHECK: better to add this iteration in temperature, so that it is
        // always executed?
        //
        return P(rho_0, epsilon, tau);
        // return press(rho_0, tau);

    };
    //===============================================================
    double P(double rho_0, double epsilon, double tau) {
        // tau passed in is assumed to be reasonable first guess for temperature,
        // return value is actual temperature
        //
        // first find tau:
        int itmax = 50;
        int it = 0;
        double delta_eps = eps(rho_0, tau) - epsilon;
        while ((fabs(delta_eps) > tol * fabs(epsilon)) and (it < itmax)) {
            it++;
            const double depsdtau = 4 * alpha * tau * tau * tau / rho_0 + 3.;
            tau -= delta_eps / depsdtau;
            delta_eps = eps(rho_0, tau) - epsilon;
        }
        if (it >= itmax) {
            cerr << "GASRAD: failed to converge in GasRad::P for rho_0 = "
                << rho_0 << " epsilon = " << epsilon << endl;
        }
        // const double check = epsilon - eps(rho_0, tau);
        // cout << " GASRAD: check solution in GasRad::P(rho_0, epsilon,tau) = "
        // 	 << check <<  " tau = " << tau << endl;
        // then compute and return pressure for given value of tau:
        return press(rho_0, tau);
    };
    //===============================================================
    //
    // find rho_0 given pressure P and *initial* entropy s_init
    // (use for initial data only!)
    //
    //===============================================================    
    double rho_0(double P) {
        if (P <= 0.0) return 0.0;
        // need to find rho_0 and tau simultaneously for given values of
        // entropy and P
        int itmax = 50;
        int it = 0;

        // get initial guess from limit of Gamma = 4/3 polytrope
        double tau = pow(3. * P / alpha, 1. / 4.);
        double rho_r = pow(P, 3. / 4.);   // use r instead of 0 to avoid confusion...
        double delta_P = press(rho_r, tau) - P;
        double delta_s = ent(rho_r, tau) - s_init;
        double delta = sqrt(delta_P * delta_P + delta_s * delta_s);
        // cout << " rho_0 = " << rho_r << " tau = " << tau << " delta = " << delta << endl;
        while ((fabs(delta) > tol * s_init) and (it < itmax)) {
            it++;
            const double dPdtau = 4. / 3. * alpha * tau * tau * tau + 2. * rho_r;
            const double dPdrho = 2.0 * tau;
            const double dsdtau = 4. * alpha * tau * tau / rho_r + 3. / tau;
            const double dsdrho = -4. * alpha * tau * tau * tau / (3. * rho_r * rho_r) - 2. / rho_r;
            const double det = dPdrho * dsdtau - dPdtau * dsdrho;
            rho_r -= (dsdtau * delta_P - dPdtau * delta_s) / det;
            tau -= (-dsdrho * delta_P + dPdrho * delta_s) / det;
            delta_P = press(rho_r, tau) - P;
            delta_s = ent(rho_r, tau) - s_init;
            delta = sqrt(delta_P * delta_P + delta_s * delta_s);
            // cout << " rho_0 = " << rho_r << " tau = " << tau
            // 	   << " delta = " << delta << endl;
        }
        if (it >= itmax) {
            cout << "GASRAD: failed to converge in GasRad::rho_0 for P = "
                << P << " s = " << s_init << endl;
        }
        return rho_r;
    };
    //===============================================================
    //
    // overwrite function rho_0_of_rho
    // (use for initial data only!)
    //
    //===============================================================    
    double rho_0_of_rho(double rho) {
        if (rho <= 0.0) return 0.0;
        // need to find rho_0 and tau simultaneously for given values of
        // entropy and rho
        int itmax = 50;
        int it = 0;
        // get initial guess from limit of Gamma = 4/3 polytrope
        double rho_r = rho;  // use subscript 'r' instead of '0' to avoid confusion with function
        double P = pow(rho_r, 4. / 3.);
        double tau = pow(3. * P / alpha, 1. / 4.);
        double delta_rho = dens(rho_r, tau) - rho;
        double delta_s = ent(rho_r, tau) - s_init;
        double delta = sqrt(delta_rho * delta_rho + delta_s * delta_s);
        // cout << " rho_0 = " << rho_r << " tau = " << tau << " delta = " << delta << endl;
        while ((fabs(delta) > tol * s_init) and (it < itmax)) {
            it++;
            const double drhodtau = 4. * alpha * tau * tau * tau + 3. * rho_r;
            const double drhodrho_0 = 1.0 + 3.0 * tau;
            const double dsdtau = 4. * alpha * tau * tau / rho_r + 3. / tau;
            const double dsdrho_0 = -4. * alpha * tau * tau * tau / (3. * rho_r * rho_r) - 2. / rho_r;
            const double det = drhodrho_0 * dsdtau - drhodtau * dsdrho_0;
            rho_r -= (dsdtau * delta_rho - drhodtau * delta_s) / det;
            tau -= (-dsdrho_0 * delta_rho + drhodrho_0 * delta_s) / det;
            delta_rho = dens(rho_r, tau) - rho;
            delta_s = ent(rho_r, tau) - s_init;
            delta = sqrt(delta_rho * delta_rho + delta_s * delta_s);
            // cout << " rho_0 = " << rho_r << " tau = " << tau
            // 	   << " delta = " << delta << endl;
        }
        if (it >= itmax) {
            cout << "GASRAD: failed to converge in GasRad::rho_0_of_rho for rho = "
                << rho << " s = " << s_init << endl;
        }
        return rho_r;
    };
    //===============================================================
    //
    // "cold" epsilon, i.e. for initial entropy 
    //
    //===============================================================    
    double cold_eps(double rho_0) {
        if (rho_0 <= 0.0) return 0.0;
        int itmax = 100;
        int it = 0;
        // iterate over s to find temperature...
        double tau = tau_guess;
        double delta_s = ent(rho_0, tau) - s_init;
        while ((fabs(delta_s) > tol * s_init) and (it < itmax)) {
            it++;
            const double dsdtau = 4 * alpha * tau * tau / rho_0 + 3. / tau;
            tau -= delta_s / dsdtau;
            delta_s = ent(rho_0, tau) - s_init;
            if (it > 90) {
                cout << "GASRAD: in cold_eps: it = " << it << " tau = "
                    << tau << " delta_s = " << delta_s
                    << " dsdtau = " << dsdtau << endl;
            }
        }
        if (it >= itmax) {
            cout << "GASRAD: failed to converge in GasRad::cold_eps for rho_0 = "
                << rho_0 << " -- found tau = " << tau << endl;
        }
        // cout << " For rho_0 = " << rho_0 << " found tau = " << tau << endl;
        // then compute and return epsilon for given value of tau:
        return eps(rho_0, tau);
    };
    //===============================================================
    //
    // entropy
    //
    // uses same algorith as P(rho_0, epsilon) above!
    //
    //===============================================================    
    double h(double rho_0, double epsilon) {
        if (rho_0 <= 0.0)
            return 0.0;
        else {
            //
            // first find tau by solving epsilon equation analytically,
            // then compute P from rho_0 and tau
            //
            const double tau = temperature(rho_0, epsilon);
            const double P = press(rho_0, tau);
            return 1.0 + epsilon + P / rho_0;
        };
    };
    //===============================================================
    //
    // find temperature (analytically from rho_0 and epsilon)
    //
    //===============================================================    
    // double temperature() { return tau_last; };
    double temperature(double rho_0, double epsilon) {
        if (rho_0 <= 0.0)
            return 0.0;
        else {
            const double beta3 = 3. * rho_0 / alpha;
            const double gamma4 = epsilon * rho_0 / alpha;
            const double Q = -4. * gamma4 / 3.;
            const double R = -beta3 * beta3 / 2.;
            const double A_arg = -R + sqrt(R * R - Q * Q * Q);
            const double A = pow(A_arg, 1. / 3.);
            const double x = A + Q / A;
            //
            // NOTE: numerical value for x may not be very accurate,
            // since A and Q/A may have
            // similar magnitude but have opposite sign.
            // 
            // cout << " A = " << A << " Q/A = " << Q/A << endl;
            if (x < 0.0) {
                cerr << " GASRAD: x < 0 in GasRad::temperature for rho_0 = "
                    << rho_0 << " epsilon = " << epsilon << endl;
                exit(0);
            }
            const double a = sqrt(x);
            const double T_disc = beta3 / (2. * a) - a * a / 4.;
            if (T_disc < 0.0) {
                cerr << " GASRAD: T_disc < 0 in GasRad::temperature for rho_0 = "
                    << rho_0 << " epsilon = " << epsilon << endl;
                exit(0);
            }
            // const double b2 = a*a/2 - beta3/(2.*a);
            // cout << " b2 = " << b2 << " a2/4 = " << a*a / 4 << endl;
            // const double q = - 0.5*(a + sqrt( a*a - 4*b2 ));
            // const double tau = b2 / q;
            double tau = -a / 2. + sqrt(T_disc);
            //
            // because of truncation error resulting from the x issue
            // above, may want to follow up with a few iteration steps...
            //
            int itmax = 50;
            int it = 0;
            double delta_eps = eps(rho_0, tau) - epsilon;
            //      cout << "it = " << it << " delta_eps = " << delta_eps << endl;
            while ((fabs(delta_eps) > tol * fabs(epsilon)) and (it < itmax)) {
                it++;
                const double depsdtau = 4 * alpha * tau * tau * tau / rho_0 + 3.;
                tau -= delta_eps / depsdtau;
                delta_eps = eps(rho_0, tau) - epsilon;
                // cout << "it = " << it << " delta_eps = " << delta_eps << endl;
            }
            if (it >= itmax) {
                cerr << "GASRAD: failed to converge in GasRad::temperature for rho_0 = "
                    << rho_0 << " epsilon = " << epsilon << endl;
            }
            return tau;
        };
    };
    //===============================================================    
    //
    // sound speed
    //
    //===============================================================    
    double sound_speed(double rho_0, double epsilon) {
        const double tau = temperature(rho_0, epsilon);
        const double ent_rad = 4. * alpha * tau * tau * tau / (3. * rho_0);
        const double sigma = (ent_rad + 2.) / (3. * ent_rad + 3);
        //
        const double eps_rad = alpha * tau * tau * tau * tau / rho_0;
        const double eps_gas = 3. * tau;
        const double eps = eps_rad + eps_gas;
        //
        const double P_rad = alpha * tau * tau * tau * tau / 3.0;
        const double P_gas = 2. * rho_0 * tau;
        const double P = P_rad + P_gas;
        // 
        const double H = 1 + eps + P / rho_0;
        const double temp = 2. * sigma * eps_rad + (1. + sigma) * eps_gas;
        return sqrt(2. * temp / (3. * H));
    };
    //===============================================================    
    //
    // A bunch of relics that should never be used...
    //
    //===============================================================    
    // double dPdrho_0(double rho_0, double epsilon)  {   
    //   cout << " in dPdrho_0!" << endl; exit(0); return 0.0;
    // };
    // double dPdepsilon(double rho_0, double epsilon) {
    //   cout << " in dPdepsilon!" << endl; exit(0); return 0.0; };
    // // 
    double epsilon_root(double tau) {
        cout << " in epsilon_root! " << endl;
        exit(0);
        return 0.0;
    };
    //===============================================================    
    //
    // extra functions
    //
    //===============================================================    
    inline double press(double rho_0, double tau) {
        const double P_rad = alpha * tau * tau * tau * tau / 3.0;
        const double P_gas = 2. * rho_0 * tau;
        return P_rad + P_gas;
    };
    inline double dens(double rho_0, double tau) {
        return rho_0 * (1.0 + eps(rho_0, tau));
    }
    inline double eps(double rho_0, double tau) {
        const double eps_rad = alpha * tau * tau * tau * tau / rho_0;
        const double eps_gas = 3. * tau;
        return eps_rad + eps_gas;
    };
    inline double ent(double rho_0, double tau) {
        const double ent_rad = 4. * alpha * tau * tau * tau / (3. * rho_0);
        const double ent_gas = log(gamma * tau * tau * tau / (rho_0 * rho_0)) + 5.0;
        return ent_rad + ent_gas;
    };
    //===============================================================    
    //
    // Test...
    //
    //===============================================================    
    double Test(double P_test) {
        cout << " GASRAD_TEST: start with P = " << P_test << " and s_init = " << s_init << endl;
        double rho_r = rho_0(P_test);
        cout << " GASRAD_TEST: rho_0 = " << rho_r << endl;
        double epsilon = cold_eps(rho_r);
        cout << " GASRAD_TEST: epsilon = " << epsilon << endl;
        double P_fin_1 = P(rho_r, epsilon);
        double P_fin_2 = P(rho_r, epsilon, 0.0);
        cout << " GASRAD_TEST: P_fin = " << P_fin_1 << " or " << P_fin_2 << endl;
        double rho = rho_r * (1 + epsilon);
        double rho_0_final = rho_0_of_rho(rho);
        cout << " GASRAD_TEST: rho_0_init = " << rho_r << " rho_0_final " << rho_0_final << endl;
        return P_fin_1;
    };
};


#endif /* GASRAD */
