// Tell emacs that this is -*-c++-*- mode
//================================================
//
// General TOV integrator: 
// integrates TOV equations for central rest_mass density rho_0 = rho_c,
// then provides TOV solution as function of isotropic radius for general
// EOS, provided by class EOS.h
//
// In constructor, provide
// - central *rest mass* density rho_c (= rho_0 at center)
// - dr_init: initial radial step size for integration
// - N_array: size of arrays in which interation results are stored in 
//   stellar interior; used for interpolation to arbitrary radii
//
// 
//
//================================================
#ifndef SMS_TOV_H
#define SMS_TOV_H

#include "../EOS.h"
#include "ludcmp.h"

class SMS_TOV_Solution {
private:
    int N_array;
    EOS* eos;
    double M;      // stellar gravitational mass
    double R;      // stellar areal radius
    double R_iso;  // stellar isotropic radius
    double C;      // constant of integration relating areal and isotropic radius
    double phi_c;  // central value of phi, where lapse = e^phi
    double s;      // entropy per baryon (in units of k_B)
    double P_current, rho_0_current;
    double rho_last, P_last, tau_last;
    double alpha, s_fac;   // radiation constant and factor in entropy
    double PI;
    VecDoub* rho_0_array, * lapse_array, * psi_array;
    VecDoub* r_array, * r_iso_array, * E_array, * p_gas_array;
    typedef void (SMS_TOV_Solution::* USR_FUN_PTR)(VecDoub_I& x,
        VecDoub_O& fvec,
        MatDoub_O& fjac);
    typedef double (SMS_TOV_Solution::* FCT_PTR)(double);
public:
    //================================================
    // Constructor: carries out integration and sets up arrays for interpolation
    //================================================
    SMS_TOV_Solution(double rho_c, double s_i,
        double dr_init, double N_array_i) :
        s(s_i), N_array(N_array_i) {
        PI = acos(-1.0);
        //
        // define natural constants in cgs units
        // 
        const double G = 6.67408e-8;
        const double c = 2.99792e10;
        const double h_cgs = 6.62607e-27;
        const double m_B_cgs = 1.67493e-24;
        const double m_e_cgs = 9.109384e-28;
        const double M_sun_cgs = 1.9891e33;
        //
        // convert to solar masses
        // 
        const double M_sun = G * M_sun_cgs / (c * c);  // solar mass in cm 
        const double h = G * h_cgs / (c * c * c * M_sun * M_sun);
        const double m_B = G * m_B_cgs / (c * c * M_sun);
        const double m_e = G * m_e_cgs / (c * c * M_sun);
        //
        // assemble \bar \alpha
        //
        const double PI_5 = PI * PI * PI * PI * PI;
        const double m_B_4 = m_B * m_B * m_B * m_B;
        const double h_3 = h * h * h;
        alpha = (8. * PI_5 / 15.) * m_B_4 / h_3;
        s_fac = 4.0 * (8.0 * PI * PI * PI) * pow(m_e, 1.5) * pow(m_B / h, 6) * sqrt(m_B);
        //
        // cout << " SMS_TOV_Solution: Natural constants: " << endl;
        // cout << "       h = " << h << " m_B = " << m_B << " m_e = " << m_e 
        // 	 << "       bar alpha = " << alpha << " s_fac = " << s_fac << endl;
        //
        // integrate a first time to find radius and mass
        //
        rho_last = rho_c;
        tau_last = 1.e-8;
        P_last = P_of_rho_0(rho_c);
        double eps = eps_of_rho_0(rho_c, tau_last);
        cout << " SMS_TOV: Central values: eps = " << eps << " P = " << P_last
            << " rho_0c = " << rho_0_of_P(P_last) << " tau = " << tau_last << endl;

        Integrate(rho_c, dr_init);
        //
        // FIX THIS!!!
        eps = eps_of_rho_0(rho_c, tau_last);
        cout << " SMS_TOV: rho_0_c = " << rho_c << ", rho_c = " << rho_c * (1.0 + eps)
            << "          M = " << M << ", R = " << R << ", R_iso = " << R_iso << endl;
        //
        // allocate arrays for interpolation
        //
        rho_0_array = new VecDoub(N_array);
        lapse_array = new VecDoub(N_array);
        psi_array = new VecDoub(N_array);
        r_array = new VecDoub(N_array);
        r_iso_array = new VecDoub(N_array);
        E_array = new VecDoub(N_array);
        p_gas_array = new VecDoub(N_array);
        //
        // Now integrate again to fill arrays
        //
        rho_last = rho_c;
        Integrate(rho_c);
        //    
    }
    //================================================
    // Destructor
    //================================================
    ~SMS_TOV_Solution() {
        delete rho_0_array;
        delete lapse_array;
        delete psi_array;
        delete r_array;
        delete r_iso_array;
        delete E_array;
        delete p_gas_array;
    }
    //================================================
    // lapse as function of isotropic radius
    //================================================
    double rho_0(double r_bar) {
        if (r_bar > R_iso) {
            //
            // in exterior
            //    
            return 0.0;
        } else {
            //
            // in interior
            //    
            return Interpolator(r_bar, rho_0_array);
        }
    };
    //================================================
    // lapse as function of isotropic radius
    //================================================
    double lapse(double r_bar) {
        if (r_bar > R_iso) {
            //
            // in exterior
            //    
            double help = 1.0 + M / (2.0 * r_bar);
            double r_areal = r_bar * help * help;
            return sqrt(1.0 - 2.0 * M / r_areal);
        } else {
            //
            // in interior
            //    
            return Interpolator(r_bar, lapse_array);
        }
    };
    //================================================
    // conformal factor as function of isotropic radius
    //================================================
    double psi(double r_bar) {
        if (r_bar > R_iso) {
            //
            // in exterior
            //    
            return 1.0 + M / (2.0 * r_bar);
        } else {
            //
            // in interior
            //    
            return Interpolator(r_bar, psi_array);
        }
    };
    //================================================
    // radiation energy density as function of isotropic radius
    //================================================
    double E(double r_bar) {
        if (r_bar > R_iso) {
            //
            // in exterior
            //    
            return 0.0;
        } else {
            //
            // in interior
            //    
            return Interpolator(r_bar, E_array);
        }
    };
    //================================================
    // gas pressure
    //================================================
    double p_gas(double r_bar) {
        if (r_bar > R_iso) {
            //
            // in exterior
            //    
            return 0.0;
        } else {
            //
            // in interior
            //    
            return Interpolator(r_bar, p_gas_array);
        }
    };
    //================================================
    // Interpolator, for any array
    //================================================
    double Interpolator(double r_bar, VecDoub* data_array) {
        int order = 4;   // use order points in interpolation
        //
        // find lower index
        //    
        int ind_low = 0;
        hunt(*r_iso_array, r_bar, ind_low);    // find index left to r_bar
        //    cout << " r_bar = " << r_bar << " ind_low = " << ind_low << " bracket: "
        //	 << (*r_iso_array)[ind_low] << " : " << (*r_iso_array)[ind_low + 1] << endl;
        int ind_min = ind_low - (order / 2 - 1); // lower index for interpolation 
        if (ind_min < 0) ind_min = 0;          // adjust for lower end
        if (ind_min > N_array - order) ind_min = N_array - order;  // adjust at upper end
        //
        // fill short arrays for interpolation
        //    
        VecDoub r_iso_short(order);
        VecDoub data_short(order);
        for (int i = 0; i < order; i++) {
            r_iso_short[i] = (*r_iso_array)[ind_min + i];
            data_short[i] = (*data_array)[ind_min + i];
        }
        //
        // use short arrays for interpolation
        //    
        double data;
        double d_data;
        polint(r_iso_short, data_short, r_bar, data, d_data);
        return data;
    }

    //================================================
    // First: Integrate TOV equations until pressure vanishes
    //================================================
    void Integrate(double rho_c, double dr_init) {
        VecDoub vars(4);
        VecDoub d_vars_dr(4);
        // initialize variables
        vars[0] = P_of_rho_0(rho_c);  // pressure
        vars[1] = 0.0;                // mass
        vars[2] = 0.0;                // phi
        vars[3] = 0.0;                // integral for isotropic radius
        double r = 0.0;
        double dr = dr_init;
        while (vars[0] > 0.0) {
            Derivs(r, vars, d_vars_dr);
            r += dr;
            vars[0] += d_vars_dr[0] * dr;
            vars[1] += d_vars_dr[1] * dr;
            vars[2] += d_vars_dr[2] * dr;
            vars[3] += d_vars_dr[3] * dr;
            //      cout << " at r = " << r << ": P = " << vars[0] << " rho_0 = " << rho_0_of_P(vars[0]) << endl;
        }
        M = vars[1];
        R = r;
        R_iso = 0.5 * (sqrt(R * (R - 2.0 * M)) + R - M);
        C = R_iso / R * exp(-vars[3]);
        phi_c = log(1.0 - 2.0 * M / R) / 2.0 - vars[2];
    };
    //================================================
    // Then: Integrate TOV equations to fill arrays
    //================================================
    void Integrate(double rho_c) {
        //
        VecDoub vars(4);
        VecDoub d_vars_dr(4);
        //
        // initialize integration variables
        //
        vars[0] = P_of_rho_0(rho_c);  // pressure
        vars[1] = 0.0;                // mass
        vars[2] = phi_c;              // phi
        vars[3] = 0.0;                // integral for isotropic radius
        //
        // fill central values of arrays
        //    
        (*rho_0_array)[0] = rho_c;
        (*lapse_array)[0] = exp(phi_c);
        (*r_array)[0] = 0.0;
        (*r_iso_array)[0] = 0.0;
        (*psi_array)[0] = sqrt(1.0 / C);
        (*E_array)[0] = alpha * tau_last * tau_last * tau_last * tau_last / 3.0;
        (*p_gas_array)[0] = 2.0 * rho_c * tau_last;
        //
        double r = 0.0;              // current radius
        double dr = R / (N_array - 1); // distance between array points
        int N = 100;                 // carry out N integraion steps between array points:
        double dr_int = dr / N;        // integration steps 
        for (int i = 1; i < N_array; i++) {  // loop over array points
            for (int j = 1; j <= N; j++) {      // integration loop
                Derivs(r, vars, d_vars_dr);
                r += dr_int;
                vars[0] += d_vars_dr[0] * dr_int;
                vars[1] += d_vars_dr[1] * dr_int;
                vars[2] += d_vars_dr[2] * dr_int;
                vars[3] += d_vars_dr[3] * dr_int;
                if (vars[0] < 0.0) vars[0] = 0.0;
            }
            //      cout << " radius at grid point " << i << " = " << rl << endl;
            //
            // fill array...
            //

            (*rho_0_array)[i] = rho_0_of_P(vars[0]);
            (*lapse_array)[i] = exp(vars[2]);
            (*r_array)[i] = r;
            (*r_iso_array)[i] = C * r * exp(vars[3]);
            (*psi_array)[i] = sqrt(r / (*r_iso_array)[i]);
            (*E_array)[i] = alpha * tau_last * tau_last * tau_last * tau_last / 3.0;
            (*p_gas_array)[i] = 2.0 * (*rho_0_array)[i] * tau_last;
        }
    };
    //================================================
    // Compute derivatives (TOV equations, including equation (4) for transformation to
    // isotropic coordinates)
    //================================================
    void Derivs(double r, VecDoub f, VecDoub& dfdr) {
        if (r == 0.0) {
            dfdr[0] = 0.0;
            dfdr[1] = 0.0;
            dfdr[2] = 0.0;
            dfdr[3] = 0.0;
        } else {
            const double P = f[0];
            const double m = f[1];
            //      const double phi = f[2];
            const double rho_0 = rho_0_of_P(P);
            const double eps = eps_of_rho_0(rho_0, tau_last);  // specific internal energy
            const double rho = rho_0 * (1.0 + eps);            // total mass-energy density
            dfdr[0] = -(rho + P) * (m + 4.0 * PI * r * r * r * P) / (r * r - 2.0 * m * r);
            dfdr[1] = 4.0 * PI * rho * r * r;
            dfdr[2] = -dfdr[0] / (rho + P);
            dfdr[3] = (1.0 - sqrt(1.0 - 2.0 * m / r)) / (r * sqrt(1.0 - 2.0 * m / r));
        }
    };
    //================================================
    // EOS for mix of plasma and radiation
    //================================================
    double P_of_rho_0(double rho_0) {
        rho_0_current = rho_0;
        // first find temperature from entropy
        double tau1 = 0.9 * tau_last;
        double tau2 = tau_last + 1.e-8;
        // make sure roots are bracketed 
        zbrac(&SMS_TOV_Solution::T_of_rho0_fun, tau1, tau2);
        // now find temperature
        const double tol = 1.e-12;
        tau_last = zbrent(&SMS_TOV_Solution::T_of_rho0_fun, tau1, tau2, tol * tau2);
        const double tau4 = tau_last * tau_last * tau_last * tau_last;
        return alpha * tau4 / 3.0 + 2.0 * rho_0 * tau_last;
        // 
    }
    double rho_0_of_P(double P) {
        P_current = P;
        const int n_trail = 100;
        const double tol = 1.e-8;
        VecDoub x(2);
        x[0] = tau_last;
        x[1] = rho_last;
        VecDoub f(2);
        MatDoub M(2, 2);
        rho0_of_P_fun(x, f, M);
        const double tolf = tol * (fabs(f[0]) + fabs(f[1]));
        int n_steps = mnewt(&SMS_TOV_Solution::rho0_of_P_fun, n_trail, x,
            tol, tolf);
        if (n_steps < 0) {
            cerr << "OUCH!!  Could not find solution in SMS_TOV_Solution::rho_0_of_P" << endl;
            exit(0);
        }
        // else 
        //   cout << " SMS_TOV_Solution: took " << n_steps << " steps in mnewt " 
        // 	   << endl;
        tau_last = x[0];
        rho_last = x[1];
        return rho_last;
    }
    double eps_of_rho_0(double rho_0, double tau) {
        return alpha * tau * tau * tau * tau / rho_0 + 3.0 * tau;
    }
    //================================================
    // User functions for minimizations...
    //================================================
    double T_of_rho0_fun(double tau) {
        const double rho = rho_0_current;
        const double tau3 = tau * tau * tau;
        return 4. * alpha * tau3 / (3.0 * rho) + log(s_fac * tau3 / (rho * rho)) + 5 - s;
    };
    void rho0_of_P_fun(VecDoub_I& x, VecDoub_O& fvec, MatDoub_O& fjac) {
        // F1(tau, rho_0) = P - a tau^4 / 3 - 2 rho_0 tau
        // F2(tau, rho_0) = s - ...
        // J11 = \partial_T J_1 = - 4 a tau^4 /3 - 2 rho_0
        // J12 = \partial_rho J_1 = = - 2 tau
        // J21 = \partial_T J_2 = ...
        // J22 = \partial_rho J_2 = ...
        const double tau = x[0];
        const double rho = x[1];  // rest-mass density
        const double tau3 = tau * tau * tau;
        const double tau4 = tau3 * tau;
        const double rho2 = rho * rho;
        const double F1 = P_current - alpha * tau4 / 3.0 - 2.0 * rho * tau;
        const double F2 = s - 4.0 * alpha * tau3 / (3.0 * rho)
            - log(s_fac * tau3 / rho2) - 5.0;
        const double J11 = -4.0 * alpha * tau3 / 3.0 - 2.0 * rho;
        const double J12 = -2.0 * tau;
        const double J21 = -4.0 * alpha * tau * tau / rho - 3.0 / tau;
        const double J22 = 4.0 * alpha * tau3 / (3.0 * rho2) + 2.0 / rho;
        fvec[0] = F1;
        fvec[1] = F2;
        fjac[0][0] = J11;
        fjac[0][1] = J12;
        fjac[1][0] = J21;
        fjac[1][1] = J22;
    }
    //================================================
    // Numerical recipes routines for interpolation
    //================================================
    void polint(VecDoub& xa, VecDoub& ya, const Doub x, Doub& y, Doub& dy) {
        int i, m, ns = 0;
        Doub den, dif, dift, ho, hp, w;

        int n = xa.size();
        VecDoub c(n), d(n);
        dif = fabs(x - xa[0]);
        for (i = 0;i < n;i++) {
            if ((dift = fabs(x - xa[i])) < dif) {
                ns = i;
                dif = dift;
            }
            c[i] = ya[i];
            d[i] = ya[i];
        }
        y = ya[ns--];
        for (m = 1;m < n;m++) {
            for (i = 0;i < n - m;i++) {
                ho = xa[i] - x;
                hp = xa[i + m] - x;
                w = c[i + 1] - d[i];
                if ((den = ho - hp) == 0.0) cerr << "Error in routine polint" << endl;
                den = w / den;
                d[i] = hp * den;
                c[i] = ho * den;
            }
            y += (dy = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
        }
    };
    void hunt(VecDoub& xx, const Doub x, int& jlo) {
        int jm, jhi, inc;
        bool ascnd;

        int n = xx.size();
        ascnd = (xx[n - 1] >= xx[0]);
        if (jlo < 0 || jlo > n - 1) {
            jlo = -1;
            jhi = n;
        } else {
            inc = 1;
            if (x >= xx[jlo] == ascnd) {
                if (jlo == n - 1) return;
                jhi = jlo + 1;
                while (x >= xx[jhi] == ascnd) {
                    jlo = jhi;
                    inc += inc;
                    jhi = jlo + inc;
                    if (jhi > n - 1) {
                        jhi = n;
                        break;
                    }
                }
            } else {
                if (jlo == 0) {
                    jlo = -1;
                    return;
                }
                jhi = jlo--;
                while (x < xx[jlo] == ascnd) {
                    jhi = jlo;
                    inc <<= 1;
                    if (inc >= jhi) {
                        jlo = -1;
                        break;
                    } else jlo = jhi - inc;
                }
            }
        }
        while (jhi - jlo != 1) {
            jm = (jhi + jlo) >> 1;
            if (x >= xx[jm] == ascnd)
                jlo = jm;
            else
                jhi = jm;
        }
        if (x == xx[n - 1]) jlo = n - 2;
        if (x == xx[0]) jlo = 0;
    };


    int mnewt(USR_FUN_PTR usrfun, const Int ntrial, VecDoub_IO& x,
        const Doub tolx, const Doub tolf) {
        Int i, n = x.size();
        VecDoub p(n), fvec(n);
        MatDoub fjac(n, n);
        for (Int k = 0;k < ntrial;k++) {
            (this->*usrfun)(x, fvec, fjac);
            Doub errf = 0.0;
            for (i = 0;i < n;i++) errf += abs(fvec[i]);
            if (errf <= tolf) return k;
            for (i = 0;i < n;i++) p[i] = -fvec[i];
            LUdcmp alu(fjac);
            alu.solve(p, p);
            Doub errx = 0.0;
            for (i = 0;i < n;i++) {
                errx += abs(p[i]);
                x[i] += p[i];
            }
            if (errx <= tolx) return k;
        }
        return -1;
    };
    //================================================
    // Numerical Recipes routines for rootfinding
    //================================================
    // zbrac
    //================================================
    Bool zbrac(FCT_PTR func, Doub& x1, Doub& x2) {
        const Int NTRY = 50;
        const Doub FACTOR = 1.6;
        if (x1 == x2) cerr << "EOS: Bad initial range in EOS::zbrac" << endl;
        Doub f1 = (this->*func)(x1);
        Doub f2 = (this->*func)(x2);
        for (Int j = 0;j < NTRY;j++) {
            if (f1 * f2 < 0.0) return true;
            if (abs(f1) < abs(f2))
                f1 = (this->*func)(x1 += FACTOR * (x1 - x2));
            else
                f2 = (this->*func)(x2 += FACTOR * (x2 - x1));
        }
        return false;
    }
    //================================================
    // zbrent
    //================================================
    Doub zbrent(FCT_PTR func, const Doub x1, const Doub x2, const Doub tol) {
        const Int ITMAX = 100;
        const Doub EPS = numeric_limits<Doub>::epsilon();
        Doub a = x1, b = x2, c = x2, d, e, fa = (this->*func)(a), fb = (this->*func)(b), fc, p, q, r, s, tol1, xm;
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
            cerr << "Root must be bracketed in zbrent" << endl;
        fc = fb;
        for (Int iter = 0;iter < ITMAX;iter++) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                c = a;
                fc = fa;
                e = d = b - a;
            }
            if (abs(fc) < abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            tol1 = 2.0 * EPS * abs(b) + 0.5 * tol;
            xm = 0.5 * (c - b);
            if (abs(xm) <= tol1 || fb == 0.0) return b;
            if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
                s = fb / fa;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0) q = -q;
                p = abs(p);
                Doub min1 = 3.0 * xm * q - abs(tol1 * q);
                Doub min2 = abs(e * q);
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }
            a = b;
            fa = fb;
            if (abs(d) > tol1)
                b += d;
            else
                b += SIGN(tol1, xm);
            fb = (this->*func)(b);
        }
        cerr << "EOS: Maximum number of iterations exceeded in zbrent" << endl;
        return b;
    }
};

#endif  /* SMS_TOV_H */
