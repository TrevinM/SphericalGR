// Tell emacs that this is -*-c++-*- mode
//================================================
// Classes for initial data
//================================================
#ifndef EOS_H
#define EOS_H

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

using namespace std;
enum { poly, gammalaw };

//
//================================================
// 
// Base class for equation of state - doesn't do much...
//
//================================================
//
class EOS {
protected:
    int eos_type;
    double last_rho_0, rho_current, tol;   // used for iteration in rho_0(rho)
    ostringstream name;
public:
    EOS() : last_rho_0(0.0), rho_current(0.0), tol(1.e-18) {};
    ~EOS() {};
    virtual double P(double rho_0, double epsilon) = 0;
    virtual double rho_0(double P) = 0;     // should be used for initial data only.
    virtual double dPdrho_0(double rho_0, double epsilon) = 0;
    virtual double dPdepsilon(double rho_0, double epsilon) = 0;
    virtual double h(double rho_0, double epsilon) = 0;
    virtual double cold_eps(double rho_0) = 0;
    double sound_speed(double rho_0, double epsilon) {
        const double p = P(rho_0, epsilon);
        const double h = 1.0 + epsilon + p / rho_0;
        const double cs2 = (dPdrho_0(rho_0, epsilon) +
            p * dPdepsilon(rho_0, epsilon) / (rho_0 * rho_0)) / h;
        if (!isfinite(cs2))
            cout << " EOS: ouch... cs2 messed up - h = " << h << " rho_0 " << rho_0 << endl;
        if (cs2 > 1.0)
            cout << " EOS: ouch... cs2 > 1.0 for rho_0 = " << rho_0 << " and epsilon = " << epsilon << endl;
        if (cs2 < 0.0)
            cout << " EOS: Ouch... cs2 < 0.0 for rho_0 = " << rho_0 << " and epsilon = " << epsilon << endl;
        return sqrt(abs(cs2));
    };
    //================================================
    // Return name of EOS
    //================================================
    const char* Name() {
        return name.str().c_str();
    };
    //================================================
    // Routines needed for initial data
    //================================================
    double rho_0_of_rho(double rho) {       // needed in TOV_BH
        if (rho < 0.0) cerr << "EOS: found negative rho in rho_0_of_rho!" << endl;
        if (rho == 0.0) return 0.0;
        rho_current = rho;
        if (last_rho_0 == 0) last_rho_0 = rho;
        double rho_0_1 = last_rho_0;
        double rho_0_2 = 0.9 * rho_0_1;
        zbrac(rho_0_1, rho_0_2);   // make sure that roots are bracketed!
        return last_rho_0 = zbrent(rho_0_1, rho_0_2, tol * rho);
    }
    double rho_0_root(double rho_0) {       // used for rootfinding in rho_0_of_rho
        double eps = cold_eps(rho_0);
        double rho = rho_0 * (1.0 + eps);
        return rho - rho_current;
    }
    //================================================
    // Numerical Recipes routines for rootfinding
    //================================================
    Bool zbrac(Doub& x1, Doub& x2) {
        const Int NTRY = 50;
        const Doub FACTOR = 1.6;
        if (x1 == x2) cerr << "EOS: Bad initial range in EOS::zbrac" << endl;
        Doub f1 = rho_0_root(x1);
        Doub f2 = rho_0_root(x2);
        for (Int j = 0;j < NTRY;j++) {
            if (f1 * f2 < 0.0) return true;
            if (abs(f1) < abs(f2))
                f1 = rho_0_root(x1 += FACTOR * (x1 - x2));
            else
                f2 = rho_0_root(x2 += FACTOR * (x2 - x1));
        }
        return false;
    }
    Doub zbrent(const Doub x1, const Doub x2, const Doub tol) {
        const Int ITMAX = 100;
        const Doub EPS = numeric_limits<Doub>::epsilon();
        Doub a = x1, b = x2, c = x2, d, e, fa = rho_0_root(a), fb = rho_0_root(b), fc, p, q, r, s, tol1, xm;
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
            fb = rho_0_root(b);
        }
        cerr << "EOS: Maximum number of iterations exceeded in zbrent" << endl;
    }


};
//
//================================================
// derived class for Polytrope
//
//    P = K * rho_0^Gamma
//
//================================================
//
class polytrope : public EOS {
private:
    double K;
    double Gamma;
    double n;
public:
    // Constructor
    polytrope(char* input_file) : EOS() {
        eos_type = poly;
        ifstream infile;
        infile.open(input_file);
        if (!infile) {
            cerr << " EOS: Can't open " << input_file
                << " for input. This is bad. " << endl;
        } else
            cout << " EOS: Reading EOS polytropic parameters from file "
            << input_file << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> K;
        infile.get(buf, 100, '='); infile.get(c); infile >> Gamma;
        cout << " EOS: Setting up polytropic EOS with Kappa = " << K
            << " and Gamma = " << Gamma << endl;
        n = 1.0 / (Gamma - 1.0);
        name << "polytropic EOS with K = " << K
            << " and Gamma = " << Gamma << ends;
    };
    // Destructor
    ~polytrope() {};
    //
    // Routines
    //
    inline double P(double rho_0, double epsilon) {
        if (Gamma == 2.0)
            return K * rho_0 * rho_0;
        else
            return K * pow(rho_0, Gamma);
    };
    inline double rho_0(double P) {
        ;
        return pow(P / K, 1.0 / Gamma);
    }
    inline double dPdrho_0(double rho_0, double epsilon) {
        return epsilon * (Gamma - 1.0);
        // if (Gamma == 2.0)
        //   return 2.0 * K * rho_0;
        // else 
        //   return Gamma * K * pow(rho_0,Gamma-1.0);
    };
    inline double dPdepsilon(double rho_0, double epsilon) {
        return rho_0 * (Gamma - 1.0);
    };
    inline double h(double rho_0, double epsilon) {
        if (Gamma == 2.0)
            return 1.0 + 2.0 * K * rho_0;
        else
            return 1.0 + (n + 1.0) * K * pow(rho_0, Gamma - 1.0);
    };
    inline double cold_eps(double rho_0) {
        if (Gamma == 2.0)
            return K * rho_0;
        else
            return K * pow(rho_0, Gamma - 1.0) / (Gamma - 1.0);
    };
};


//
//================================================
// derived class for Gamma-law EOS
//
//    P = rho_0 * epsilon ( \Gamma - 1 )
//
//================================================
//
class gamma_law : public EOS {
private:
    double K;  // needed only for cold part
    double Gamma;
    double n;
public:
    // Constructor
    gamma_law(char* input_file) : EOS() {
        eos_type = gammalaw;
        ifstream infile;
        infile.open(input_file);
        if (!infile) {
            cerr << " EOS: Can't open " << input_file
                << " for input. This is bad. " << endl;
        } else
            cout << " EOS: Reading EOS polytropic parameters from file "
            << input_file << endl;
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> K;
        infile.get(buf, 100, '='); infile.get(c); infile >> Gamma;
        n = 1.0 / (Gamma - 1.0);
        cout << " EOS: Setting up gamma-law EOS with Gamma = " << Gamma << endl;
        name << "Gamma-law EOS with Gamma = " << Gamma << ends;
    };
    // Destructor
    ~gamma_law() {};
    //
    // Routines
    //
    inline double P(double rho_0, double epsilon) {
        return rho_0 * epsilon * (Gamma - 1.0);
    };
    inline double rho_0(double P) {
        ;
        return pow(P / K, 1.0 / Gamma);
    }
    inline double dPdrho_0(double rho_0, double epsilon) {
        return epsilon * (Gamma - 1.0);
    };
    inline double dPdepsilon(double rho_0, double epsilon) {
        return rho_0 * (Gamma - 1.0);
    };
    inline double h(double rho_0, double epsilon) {
        return 1.0 + Gamma * epsilon;
    };
    inline double cold_eps(double rho_0) {
        if (Gamma == 2.0)
            return K * rho_0;
        else
            return K * pow(rho_0, Gamma - 1.0) / (Gamma - 1.0);
    };
};

#endif   /* EOS_H */
