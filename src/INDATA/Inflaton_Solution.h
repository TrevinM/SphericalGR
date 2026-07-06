#include "Potential.H"
#include <fstream>

#define DEBUG_MODE=1
#ifdef DEBUG_MODE
    #define DEBUG(x,i) if (i%100 == 0) cout << (*x)[i] << ", ";
    #define DEBUG2(x,y,i) DEBUG(x) DEBUG(y) << " || ";
#else
    #define DEBUG(x,i)
    #define DEBUG2(x,y,i)
#endif


// Scalar Blob IC solver using CTTK method
class ScalarBlob_Solution
{
  public:

    // Variables given
    const double N_r, tol_ham, tol_mom, max_iter, Dr;
    const VecDoub *sf;
    const VecDoub *R;
    double PI;
    Potential* potential;
    bool choose_K;
    ofstream output_file;
    double extraction_radius;
    double MS_mass;
    double search_step;
    double A_sf;
    double K_prop;

    ScalarBlob_Solution(double a_N_r, double a_Dr, double a_A_sf, VecDoub *a_sf, VecDoub *a_K, double a_K_prop, VecDoub* a_psi, VecDoub* a_R, double a_tol_ham, double a_tol_mom,
                        Potential *a_potential, double a_extraction_radius, double a_search_step, double a_max_iter) : N_r(a_N_r), Dr(a_Dr), sf(a_sf),
                        K(a_K), K_prop(a_K_prop), psi(a_psi), R(a_R), tol_ham(a_tol_ham), tol_mom(a_tol_mom), potential(a_potential), search_step(a_search_step),
                        output_file("./Plotting/ICSolverData/ICSolverData.txt"), extraction_radius(a_extraction_radius), max_iter(a_max_iter), A_sf(a_A_sf)
    {
        psi_copy = new VecDoub(N_r);
        Wr = new VecDoub(N_r);
        delta_Wr = new VecDoub(N_r);
        PI = acos(-1.0);
        rho = new VecDoub(N_r);
        arr = new VecDoub(N_r);
        att = new VecDoub(N_r);
        app = new VecDoub(N_r);
        phi = new VecDoub(N_r);
        Lap_psi = new VecDoub(N_r);
        delta_psi = new VecDoub(N_r);
        Q = new VecDoub(N_r);
        delta_Q = new VecDoub(N_r);
        R_areal = new VecDoub(N_r);

        for (int i = 0; i < N_r; ++i)
        {
            const double rm1 = 1./(*R)[i];
            (*Lap_psi)[i] = laplacian(psi, i);
            (*phi)[i] = log((*psi)[i]);
            (*delta_psi)[i] = 0.;
            (*Q)[i] = 0.;
            (*delta_Q)[i] = 0.;
        }

        for (int i = 0; i < N_r; ++i)
        {
            double dsf_dr = dr(sf, i);
            const double psil = (*psi)[i];
            const double psim4 = 1./psil/psil/psil/psil;
            (*rho)[i] = 0.5 * psim4 * dsf_dr * dsf_dr + potential->V((*sf)[i]);
        }

        for (int i = 0; i < N_r; i+=1)
        {
            (*arr)[i] = 0.;
            (*att)[i] = 0.;
            (*app)[i] = 0.;
            (*Wr)[i] = 0.;
            (*delta_Wr)[i] = 0.;
        }

        cout << "Central density is " << (*rho)[0] << endl;
        
        cout << "Solver set up, now solving" << endl;
        SolveConstraints();
        cout << "Solver finished!" << endl;
    }

    void SolveConstraints()
    {
        int iter = 0;

        // Big loop
        while (iter < max_iter) // && (abs(err_mom) > tol_mom || abs(err_ham) > tol_ham))
        {
            //calculate_max_errors();
            SolveHamiltonian(iter == 0);
            //cout << "After solving Ham on iteration " << iter << " psi outer is " << (*psi)[N_r-1] << endl;

            //cout << "DrDr sf = ";
            for (int i = 0; i < N_r; ++i)
            {
                double dsf_dr = dr(sf, i);
                //cout << drdr(sf, i) << ", ";
                const double psil = (*psi)[i];
                const double psim4 = 1./psil/psil/psil/psil;
                (*rho)[i] = 0.5 * psim4 * dsf_dr * dsf_dr + potential->V((*sf)[i]);
            }
            //cout << endl;

            SolveMomentum(iter);

            if (iter % (int)(max_iter/100) == 0)
            {
                pair<double, double> max_errors = calculate_max_errors();
                cout << "On iteration " << iter << " Ham error max = " << max_errors.first
                     << " and mom error max = " << max_errors.second << endl;
            }
            // cout << "Psi on iter " << iter << ": ";
            // for (int i = 0; i < N_r; i += N_r/200)
            // {
            //     cout << (*psi)[i] << " ";
            // }
            // cout << endl;
            ++iter;
        }

        // Calculate phi
        for (int i = 0; i < N_r; ++i)
        {
            (*phi)[i] = log((*psi)[i]);
        }

        VecDoub* A2 = new VecDoub(N_r);
        // Rescale Aij
        for (int i = 0; i < N_r; i++)
        {
            const double psil = (*psi)[i];
            const double psim6 = 1/psil/psil/psil/psil/psil/psil;
            // Rescaled
            (*arr)[i] = psim6 * (*arr)[i];
            (*att)[i] = psim6 * (*att)[i];
            (*app)[i] = psim6 * (*app)[i];
            (*R_areal)[i] = (*psi)[i]*(*psi)[i]*(*R)[i];
        }
        write_var(Wr, "Wr");
        write_var(arr, "arr");
        write_var(arr, "rho");

        // Calculate MS mass
        int i_ext = (extraction_radius-Dr/2.)/Dr;

        MS_mass = pow((*R)[i_ext],3) * pow((*psi)[i_ext],6) * (*K)[i_ext] * (*K)[i_ext] / 18. 
                    - 2 * pow((*R)[i_ext],3) * dr(psi, i_ext) * dr(psi, i_ext)
                    - 2 * (*R)[i_ext]*(*R)[i_ext] * (*psi)[i_ext] * dr(psi, i_ext);

        cout << "K_prop and MS mass (at " << extraction_radius << ") are " << K_prop << " " << MS_mass << endl;

    }

    void SolveHamiltonian(bool first_iter = false)
    {

        if (first_iter)
        {
            // Hard-coded params for solving Ham - MAKE THESE PROPER PARAMETERS
            const int max_iter_ham = 1000;
            const double step_size = 0.1;

            for (int i = 0; i < N_r; ++i)
            {
                (*psi_copy)[i] = (*psi)[i];
            }

            // Find delta psi at r0 that gives delta psi -> 0 at outer boundary
            
            // Find low and high (for delta_psi at r0) either side of root
            double low = 0.;
            double high = 0.;

            double asymptote_high = 1;
            double asymptote_low = 1;

            int counter = 0;
            while (asymptote_low*asymptote_high > 0 && counter++ < 1000)
            {
                for (int i = 0; i < N_r; ++i)
                {
                    (*psi)[i] = (*psi_copy)[i];
                }

                Integrate(delta_psi, pair<double, double>{high, 0.}, true);
                // Update psi
                for (int i = 0; i < N_r; ++i)
                {
                    (*psi)[i] += (*delta_psi)[i];
                }

                // cout << "psi:";
                // for (int i = 0; i < 10; ++i)
                // {
                //     cout << (*psi)[i] << ", ";
                // }
                // cout << endl;
                // for (int i = N_r-10; i < N_r; ++i)
                // {
                //     cout << (*psi)[i] << ", ";
                // }
                // cout << endl;

                // Solve for delta psi
                for (int iter_ham = 0; iter_ham < max_iter_ham; ++iter_ham)
                {
                    Integrate(delta_psi, pair<double, double>{0., 0.}, true);
                    // Update psi
                    for (int i = 0; i < N_r; ++i)
                    {
                        (*psi)[i] += step_size*(*delta_psi)[i];
                    }
                }
                asymptote_high = (*psi)[N_r-1] - 1.;

                for (int i = 0; i < N_r; ++i)
                {
                    (*psi)[i] = (*psi_copy)[i];
                }

                if (low + (*psi)[0] < 1)
                    low = 1 - (*psi)[0];

                Integrate(delta_psi, pair<double, double>{low, 0.}, true);
                // Update psi
                for (int i = 0; i < N_r; ++i)
                {
                    (*psi)[i] += (*delta_psi)[i];
                }
                // Solve for delta psi
                for (int iter_ham = 0; iter_ham < max_iter_ham; ++iter_ham)
                {
                    Integrate(delta_psi, pair<double, double>{0., 0.}, true);
                    // Update psi
                    //cout << "psi:";
                    for (int i = 0; i < N_r; ++i)
                    {
                        (*psi)[i] += step_size*(*delta_psi)[i];
                        // DEBUG(psi,i)
                        //cout << (*psi)[i] << ", ";
                    }
                    //cout << endl;
                }
                asymptote_low = (*psi)[N_r-1] - 1.;
                // cout << "Outer value is " << (*psi)[N_r-1] << endl;
                if (first_iter)
                    cout << "Starting values of " << low << " and " << high << " gives asymptotes of " << asymptote_low << " and " << asymptote_high << endl;
                if (!first_iter)
                    low -= search_step * (*psi)[0];
                high += search_step * (*psi)[0];
            }
            
            double tol = 0.000001;
            double asymptote = 1;

            while (abs(asymptote) > tol && high-low >= 1e-8)
            {
                for (int i = 0; i < N_r; ++i)
                {
                    (*psi)[i] = (*psi_copy)[i];
                }
                double mid = (high+low)/2;
                // Solve for delta psi
                Integrate(delta_psi, pair<double, double>{mid, 0.}, true);
                // Update psi
                for (int i = 0; i < N_r; ++i)
                {
                    (*psi)[i] += (*delta_psi)[i];
                }
                // Solve for delta psi
                for (int iter_ham = 0; iter_ham < max_iter_ham; ++iter_ham)
                {
                    Integrate(delta_psi, pair<double, double>{0., 0.}, true);
                    // Update psi
                    for (int i = 0; i < N_r; ++i)
                    {
                        (*psi)[i] += step_size*(*delta_psi)[i];
                    }
                }
                asymptote = (*psi)[N_r-1] - 1.;
                //cout << "Starting value of " << (*psi)[0] << " gives asymptote of " << asymptote << endl;
                //cout << "This means psi0 has been corrected to " << (*psi)[0] << endl;
                // BRING BACK
                if (asymptote > 0)
                {
                    if (asymptote_high < 0)
                        low = mid;
                    else
                        high = mid;
                }
                else
                    if (asymptote_high > 0)
                        low = mid;
                    else
                        high = mid;
            }
            if (high-low < 1e-8)
            {
                cout << "ROOT FINDING FAILED" << endl;
                abort();
            }
            cout << "Solution found for psi" << endl;
        }

        //cout << "sqrt(K) = ";
        for (int i = 0; i < N_r; ++i)
        {
            const double psil = (*psi)[i];
            const double psim5 = 1./psil/psil/psil/psil/psil;
            const double psim12 = psim5*psim5/psil/psil;

            const double rm1 = 1./(*R)[i];
            double A2 = (*arr)[i]*(*arr)[i] + (*att)[i]*(*att)[i] + (*app)[i]*(*app)[i]; // Raising with h (rescaled metric)
            
            double dsf_dr = dr(sf, i);
            const double psim4 = 1./psil/psil/psil/psil;
            double rho_grad = 0.5 * psim4 * dsf_dr * dsf_dr;
            double rho_pot = K_prop * potential->V((*sf)[i]);

            (*K)[i] = -sqrt(24.*PI*(rho_grad+rho_pot) + 3./2.*psim12*A2);
            //cout << 24.*PI*rho_grad << " " << 3./2.*psim12*A2 << " ||| ";
        }
        //cout << endl;
    }

    void SolveMomentum(int iter)
    {
        /*
        int max_iter_mom = 100;
        double step_size = 0.1;

        for (int iter_mom = 0; iter_mom < max_iter_mom; ++iter_mom)
        {
            Integrate(delta_Wr, pair<double, double>{0., 0.}, false, delta_Q);
            // Update Wr
            for (int i = 0; i < N_r; ++i)
            {
                (*Wr)[i] += step_size*(*delta_Wr)[i];
                (*Q)[i] += step_size*(*delta_Q)[i];
            }
        }
        */

        Integrate(Wr, pair<double, double>{0.,0.}, false, Q);

        // Shift Q and recalculate Wr
        double Q_0 = -(*Q)[N_r-1];
        double W_0 = Q_0*(*R)[0]/3.;

        /*
        Integrate(delta_Wr, pair<double, double>{W_0,Q_0}, false, delta_Q);
        // Update Wr
        for (int i = 0; i < N_r; ++i)
        {
            (*Wr)[i] += (*delta_Wr)[i];
            (*Q)[i] += (*delta_Q)[i];
        }

        for (int iter_mom = 0; iter_mom < max_iter_mom; ++iter_mom)
        {
            Integrate(delta_Wr, pair<double, double>{0., 0.}, false, delta_Q);
            // Update Wr
            for (int i = 0; i < N_r; ++i)
            {
                (*Wr)[i] += step_size*(*delta_Wr)[i];
                (*Q)[i] += step_size*(*delta_Q)[i];
            }
        }
        */
       
        Integrate(Wr, pair<double, double>{W_0,Q_0}, false, Q);
        
        // Update Aij
        for (int i = 0; i < N_r; ++i)
        {
            const double rm1 = 1./(*R)[i];
            const double drWr = dr(Wr, i, false);
            // Rescaled
            (*arr)[i] = 4./3.*drWr-4./3.*rm1*(*Wr)[i];
            (*att)[i] = 2./3.*rm1*(*Wr)[i]-2./3.*drWr;
            (*app)[i] = 2./3.*rm1*(*Wr)[i]-2./3.*drWr;
        }
    }

    pair<double, double> rhs(pair<double, double>& vars, int i, bool use_ham, bool use_mid)
    {
        if (use_ham)
        {
            return ham_rhs(vars, i, use_mid);
        }
        else
        {
            return mom_rhs(vars, i, use_mid);
        }
    }

    pair<double, double> ham_rhs(pair<double, double>& vars, int i, bool use_mid = false)
    {
        // Calculate RHS half-way to the right of point i if use_mid = true
        if (i >= N_r)
            cout << "i = " << i << " > N_r = " << N_r << "!!!" << endl;

        double psil = (*psi)[i];
        if (use_mid)
            psil = (psil + (*psi)[i+1])/2.;
        const double psi4 = psil*psil*psil*psil;
        const double psi5 = psi4*psil;
        const double psim7 = 1/psi5/psil/psil;
        const double psim8 = psim7/psil;
        
        double A2 = (*arr)[i]*(*arr)[i] + (*att)[i]*(*att)[i] + (*app)[i]*(*app)[i]; // Raising with h (rescaled metric)
        double A2_r = (*arr)[i+1]*(*arr)[i+1] + (*att)[i+1]*(*att)[i+1] + (*app)[i+1]*(*app)[i+1];
        if (use_mid)
            A2 = (A2 + A2_r) / 2.;
        //double v = psi5 * (*K)[i] * (*K)[i] / 12. - psim7/8. * A2 - 2. * PI * psi5 * (*rho)[i];
        //double u = -5. * psi4 * (*K)[i] * (*K)[i] / 12. - psim8*7./8. * A2 + 10. * PI * psi4 * (*rho)[i];

        // New version, with K handling gradient and A2 terms:
        double Vl = potential->V((*sf)[i]);
        if (use_mid)
            Vl = (Vl + potential->V((*sf)[i+1]))/2.;
        double v = - 2. * PI * psi5 * (1-K_prop) * Vl;
        double u = 10. * PI * psi4 * (1-K_prop) * Vl;

        double lap_l = laplacian(psi, i);
        if (use_mid)
            lap_l = (lap_l + laplacian(psi, i+1)) / 2.;
        double Ham_res = v - lap_l;

        double rl = (*R)[i];
        if (use_mid)
            rl = (rl + (*R)[i+1]) / 2.;
        double dr_vars_0 = vars.second;
        double dr_vars_1 = - (2. / rl) * vars.second - u * vars.first + Ham_res;
        //cout << ((*R)[i]) << ", " << vars.second << endl;
        //cout << - (2. / (*R)[i]) * vars.second << ", " <<  - u * vars.first << ", " << Ham_res << " || " << endl;
        return pair<double, double>{dr_vars_0,dr_vars_1};
    }

    pair<double, double> mom_rhs(pair<double, double>& vars, int i, bool use_mid = false)
    {
        // Calculate momentum constraint RHS
        double psil = (*psi)[i];
        if (use_mid)
            psil = (psil + (*psi)[i+1])/2.;
        const double psi6 = psil*psil*psil*psil*psil*psil;

        double rl = (*R)[i];
        if (use_mid)
            rl = (rl + (*R)[i+1]) / 2.;

        double drKl = dr(K, i);
        if (use_mid)
            drKl = (drKl + dr(K, i+1))/2.;
        return pair<double, double> {-2./rl*vars.first + vars.second, 1./2.*psi6*drKl /*- vector_laplacian(Wr, i)*/};
    }

    void Integrate(VecDoub* variable, pair<double, double> initial_data, bool use_ham, VecDoub* second_variable = nullptr)
    {
        pair<double, double> l = initial_data; // last last point

        pair<double, double> m = {0.,0.}; // mid point

        pair<double, double> r = {0.,0.}; // next point

        (*variable)[0] = l.first;
        if (second_variable != nullptr)
        {
            (*second_variable)[0] = l.second;
        }

        // // Euler method to get first point
        // pair<double, double> dydt = rhs(l2, 0, use_ham);

        // l1.first = l2.first+dydt.first*Dr;
        // l1.second = l2.second+dydt.second*Dr;

        // (*variable)[1] = l1.first;
        // if (second_variable != nullptr)
        // {
        //     (*second_variable)[1] = l1.second;
        // }

        // Integrate using leapfrog method:
        // (i points at the first unfilled point)

        for (int i = 0; i < N_r-1; i+=1)
        {
            // Estimate midpoint:

            pair<double, double> dvdr_l = rhs(l, i, use_ham, false);

            m.first = l.first + Dr/2. * dvdr_l.first;
            m.second = l.second + Dr/2. * dvdr_l.second;

            pair<double, double> dvdr_m = rhs(m, i, use_ham, true);

            r.first = l.first + dvdr_m.first*Dr;
            r.second = l.second + dvdr_m.second*Dr;

            (*variable)[i+1] = r.first;
            if (second_variable != nullptr)
            {
                (*second_variable)[i+1] = r.second;
            }

            l = r;

            //cout << " |||||| at i = " << i << " ||||||| " << endl;
            //cout << "l2: " << l2.first << " " << l2.second << endl;
            //cout << "l1: " << l1.first << " " << l1.second << endl;
            // pair<double, double> dydt_l1 = rhs(l1, i-1, use_ham);

            // pair<double, double> r1;
            // r1.first = l2.first + dydt_l1.first*2*Dr;
            // r1.second = l2.second + dydt_l1.second*2*Dr;

            // pair<double, double> dydt_r1 = rhs(r1, i, use_ham);
            // pair<double, double> r2;
            // r2.first = l1.first + dydt_r1.first*2*Dr;
            // r2.second = l1.second + dydt_r1.second*2*Dr;

            // cout << "drl1: " << dydt_l1.first << " " << dydt_l1.second << endl;
            // cout << "drr1: " << dydt_r1.first << " " << dydt_r1.second << endl;
            // cout << "r1: " << r1.first << " " << r1.second << endl;
            // cout << "r2: " << r2.first << " " << r2.second << endl;

            // (*variable)[i] = r1.first;
            // (*variable)[i+1] = r2.first;
            // if (second_variable != nullptr)
            // {
            //     (*second_variable)[i] = r1.second;
            //     (*second_variable)[i+1] = r2.second;
            // }

            // l2 = r1;
            // l1 = r2;
        }
    }

    double dr(const VecDoub *values, int i, bool isEven = true)
    {
        if (i == 0)
        {
            double ghost_val = (*values)[i];
            if (!isEven)
            {
                ghost_val *= -1;
            }
            return (-ghost_val+(*values)[i+1]) / (2*Dr);
        }
        else if (i == N_r-1)
        {
            return ((*values)[i-2]-4*(*values)[i-1]+3*(*values)[i]) / (2*Dr);
        }
        else
        {
            return (-(*values)[i-1]+(*values)[i+1]) / (2*Dr);
        }
        /*
        // Forward:
        if (i==N_r-1)
            return (-(*values)[i-1]+(*values)[i]) / Dr;
        else
            return (-(*values)[i]+(*values)[i+1]) / Dr;
        */
    }

    double drdr(const VecDoub *values, int i, bool isEven = true)
    {
        if (i == 0)
        {
            double ghost_val = (*values)[i];
            if (!isEven)
            {
                ghost_val *= -1;
            }

            return (ghost_val-2*(*values)[i]+(*values)[i+1]) / (Dr*Dr);
        }
        else if (i == N_r-1)
        {
            return (-(*values)[i-3]+4*(*values)[i-2]-5*(*values)[i-1]+2*(*values)[i]) / (Dr*Dr);
        }
        else
        {
            return ((*values)[i-1]-2*(*values)[i]+(*values)[i+1]) / (Dr*Dr);
        }

        /*
        // Forward:
        if (i >= N_r-2)
            return ((*values)[i-2]-2*(*values)[i-1]+(*values)[i]) / (Dr*Dr);
        else
            return ((*values)[i]-2*(*values)[i+1]+(*values)[i+2]) / (Dr*Dr);
        */
    }

    double laplacian(const VecDoub *values, int i, bool isEven = true)
    {
        return drdr(values, i, isEven) + 2. / (*R)[i] * dr(values, i, isEven);
    }

    double vector_laplacian(const VecDoub *values, int i, bool isEven = true)
    {
        return drdr(values, i, isEven) + 2. / (*R)[i] * dr(values, i, isEven) - 2. / (*R)[i] / (*R)[i] * (*values)[i];
    }

    pair<double, double> calculate_max_errors()
    {
        // Find max hamiltonian error
        double err_ham_here = 0.;
        double err_mom_here = 0.;
        double err_ham = 0;
        double err_mom = 0;
        //cout << "Ham error is ";
        for (int i = 0; i < N_r-2; ++i)
        {
            const double rm1 = 1./(*R)[i];
            const double psil = (*psi)[i];
            const double psim5 = 1/psil/psil/psil/psil/psil;
            const double psim12 = psim5*psim5/psil/psil;
            const double psi6 = psil*psil*psil*psil*psil*psil;
            double A2 = (*arr)[i]*(*arr)[i] + (*att)[i]*(*att)[i] + (*app)[i]*(*app)[i]; // Raising with h (rescaled metric)
            const double d2rWr = drdr(Wr, i, false);

            err_ham_here = abs((*K)[i]*(*K)[i] - (24.*PI*(*rho)[i] + 3./2.*psim12*A2 + 12.*psim5*laplacian(psi, i)));
            err_ham = max(err_ham_here, err_ham);

            err_mom_here = 4./3.*(d2rWr+2.*rm1*dr(Wr,i)-2.*rm1*rm1*(*Wr)[i]) - 2./3. * psi6*dr(K, i);
            err_mom = max(err_mom_here, err_mom);
            //cout << err_ham_here << ", " << d2rWr << " || ";
            //DEBUG(psi, i);
        }
        //cout << endl;

        return pair<double, double>{err_ham, err_mom};
    }

    double K_0(double r_bar) {
        return Interpolator(r_bar, K);
    }

    double a_rr(double r_bar) {
        return Interpolator(r_bar, arr);
    }

    double a_tt(double r_bar) {
        return Interpolator(r_bar, att);
    }

    double a_pp(double r_bar) {
        return Interpolator(r_bar, app);
    }

    double phi_0(double r_bar) {
        return Interpolator(r_bar, phi);
    }

    double Interpolator(double r_bar, VecDoub *data_array) 
    {
        if (r_bar > (*R)[R->size()-1])
        {
            std::cout << "1) Looking out at " << r_bar << std::endl;
        }
        int index_low = std::round((r_bar)/Dr)-1;
        if (index_low > R->size()-1)
        {
            std::cout << "2) Looking out at " << r_bar << " and index " << index_low << std::endl;
            return (*data_array)[R->size()-1];
        }
        // std::cout << std::setprecision(6) << r_bar << " " << (*R)[index_low] << " " << (*R)[index_low+1] << " " << Dr << " " << std::endl;
        int index_high = index_low + 1;
        return (*data_array)[index_low];//+((*data_array)[index_high]-(*data_array)[index_low])*(r_bar-index_low*Dr)/Dr;
    }

    void write_var(VecDoub *data_array, string name)
    {
        output_file << name << std::setprecision(15);
        for (int i = 0; i < N_r; ++i)
        {
            output_file << " " << (*data_array)[i];
        }
        output_file << endl;
    }

    ~ScalarBlob_Solution() 
    {
        delete K; delete rho; delete Wr;
        delete arr; delete att; delete app;
        delete phi; delete psi; delete Lap_psi; delete delta_psi; delete rho;
        delete psi_copy;
        // DELETE MORE STUFF
    }
  //private:

    // Variables calculated here
    VecDoub *K, *Wr, *delta_Wr, *phi, *psi, *Lap_psi, *delta_psi, *rho, *Q, *delta_Q, *R_areal;
    VecDoub *arr, *att, *app, *psi_copy;
};
