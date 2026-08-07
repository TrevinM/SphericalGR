// Tell emacs that this is -*-c++-*- mode
//================================================
//
// TOV integrator: 
// integrates TOV equations for central rest_mass density rho_0 = rho_c,
// then provides TOV solution as function of isotropic radius
//
// In constructor, provide
// - polytropic constant kappa and polytropic index Gamma
// - central *rest mass* density rho_c (= rho_0 at center)
// - dr_init: initial radial step size for integration
// - N_array: size of arrays in which interation results are stored in 
//   stellar interior; used for interpolation to arbitrary radii
//
// 
//
//================================================
#ifndef TOV_H
#define TOV_H

class TOV_Solution {
private:
    double Kappa;
    double Gamma;
    int N_array;
    double M;      // stellar gravitational mass
    double R;      // stellar areal radius
    double R_iso;  // stellar isotropic radius
    double C;      // constant of integration relating areal and isotropic radius
    double phi_c;  // central value of phi, where lapse = e^phi
    double PI;
    VecDoub* rho_0_array, * lapse_array, * psi_array, * r_array, * r_iso_array;
public:
    //================================================
    // Constructor: carries out integration and sets up arrays for interpolation
    //================================================
    TOV_Solution(double Kappa_i, double Gamma_i, double rho_c,
        double dr_init, double N_array_i) :
        Kappa(Kappa_i), Gamma(Gamma_i), N_array(N_array_i) {
        PI = acos(-1.0);
        //
        // integrate a first time to find radius and mass
        //
        Integrate(rho_c, dr_init);
        //
        double eps = eps_of_rho_0(rho_c);
        cout << " TOV: rho_0_c = " << rho_c << ", rho_c = " << rho_c * (1.0 + eps)
            << ", M = " << M << ", R = " << R << ", R_iso = " << R_iso << endl;
        //
        // allocate arrays for interpolation
        //
        rho_0_array = new VecDoub(N_array);
        lapse_array = new VecDoub(N_array);
        psi_array = new VecDoub(N_array);
        r_array = new VecDoub(N_array);
        r_iso_array = new VecDoub(N_array);
        //
        // Now integrate again to fill arrays
        //
        Integrate(rho_c);
        //    
    }
    //================================================
    // Destructor
    //================================================
    ~TOV_Solution() { delete rho_0_array; delete lapse_array; delete psi_array; delete r_array; delete r_iso_array; }
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
    // lapse as function of isotropic radius
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
            const double eps = P / rho_0 / (Gamma - 1.0);        // specific internal energy
            const double rho = rho_0 * (1.0 + eps);            // total mass-energy density
            dfdr[0] = -(rho + P) * (m + 4.0 * PI * r * r * r * P) / (r * r - 2.0 * m * r);
            dfdr[1] = 4.0 * PI * rho * r * r;
            dfdr[2] = -dfdr[0] / (rho + P);
            dfdr[3] = (1.0 - sqrt(1.0 - 2.0 * m / r)) / (r * sqrt(1.0 - 2.0 * m / r));
        }
    };
    //================================================
    // Different versions of EOS
    //================================================
    double P_of_rho_0(double rho_0) { return Kappa * pow(rho_0, Gamma); }
    double rho_0_of_P(double P) { return pow(P / Kappa, 1.0 / Gamma); }
    double eps_of_rho_0(double rho_0) { return Kappa * pow(rho_0, Gamma - 1.0) / (Gamma - 1.0); }
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
                if ((den = ho - hp) == 0.0) toss("Error in routine polint");
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
};

#endif  /* TOV_H */
