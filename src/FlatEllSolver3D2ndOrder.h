// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing flat, 2nd order 3D elliptic solver
//
//================================================

#ifndef NoEllSolver

//================================================
//
// Uses Trilinos software to solve covariant 3D Laplace
// operator in spherical polar coordinates
//
//  \nabla^2 f + factor * u * f = RHS
//
// where u is a function and factor a constant
//
// Most boundary conditions are determined by
// periodicity; these are implemented in index utilities.  The 
// only exception are outer boundary conditions in r - we assume
// fall-off with power 'fall-off'
//
//================================================

//================================================
//
// Important note concerning grid:
//
// n_r, n_theta and n_phi include N_g ghost_zones, which 
// are not needed for elliptic solver.  In interior
// grid, have N_g <= i < n_r-N_g (and similar for j and
// k).  However, want "superindex" II go from 0 to its
// maximum value, (n_r - 2 N_g)*(n_theta - 2 N_g)*(n_phi - 2 N_g) - 1.
//
//================================================

#ifndef _FLATELLSOLVER3D_H_
#define _FLATELLSOLVER3D_H_


#include <cstdlib>
// #include <cassert>
#include <string>
#include <vector>

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
//#include "Epetra_Version.h"
#include "AztecOO.h"

#ifdef EPETRA_MPI
#  include "mpi.h"
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

#include "nr3.h"
#include "surfacefunction.h"
#include "gridfunction.h"


using namespace std;

class FlatEllSolver3D {
private:
    Grid* grid;
    int fall_off;   // default set to 1 - change with SetFallOff()
    int n_r;
    int n_theta;
    int n_phi;
    int n_global;
    int N_g;      // number of ghost zones - hardcoded to three
    int NumMyElements;
    VecDoub* r;  // , *theta, *phi;
#ifdef EPETRA_MPI
    Epetra_MpiComm* Comm;
#else
    Epetra_SerialComm* Comm;
#endif
    Epetra_Map* Map;
    Epetra_Vector* rhs, * sol;
    Epetra_CrsMatrix* A;
    bool verbose;
public:
    //================================================
    // Constructor
    //================================================
    FlatEllSolver3D(Grid* grid_i, bool verbose_i = false) :
        grid(grid_i), fall_off(1), verbose(verbose_i) {
        N_g = grid->N_ghosts();
        n_r = grid->N_r_tot();
        n_theta = grid->N_theta_tot();
        n_phi = grid->N_phi_tot();
        r = grid->r();
        if (verbose) cout << " Constructing FlatEllSolver3D2ndOrder";
        // If Trilinos was built with MPI, initialize MPI, otherwise
        // initialize the serial "communicator" that stands in for MPI.
#ifdef EPETRA_MPI
        Comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
        Comm = new Epetra_SerialComm();
#endif
        //    cout << *Comm << endl;
        //
        // Create Epetra Map...
        // 
        const int NumGlobalElements = (n_r - 2 * N_g) * (n_theta - 2 * N_g) * (n_phi - 2 * N_g);
        if (verbose) cout << " for " << NumGlobalElements << "x"
            << NumGlobalElements << " matrix..." << endl;
        Map = new Epetra_Map(NumGlobalElements, 0, *Comm);
        n_global = NumGlobalElements;
        //
        // ... then allocate vectors
        // 
        rhs = new Epetra_Vector(*Map);
        sol = new Epetra_Vector(*Map);
        //
        A = NULL;
    };
    //================================================
    // Destructor
    //================================================
    ~FlatEllSolver3D() {
        if (verbose) cout << " Destructing FlatEllSolver3D... " << endl;
        delete Comm;
        delete Map;
        delete rhs;
        delete A;
    };
    //================================================
    // Make verbose
    //================================================
    bool SetVerbose(bool verb) { return verbose = verb; }
    //================================================
    // Change fall-off (set to 1 by default)
    //================================================
    int SetFallOff(int fall_off_i) { return fall_off = fall_off_i; }
    //================================================
    // Set up Solver 
    //
    // (note that DG_i's are \bar \gamma^{lm} \Delta \Gamma^i_{lm};
    //  will add flat contributions in spherical polar coordinates in
    //  SetupCovLaplace below)
    //================================================
    int SetupSolver(double factor, gf3d& u) {
        int ierr = 0;
        // even though rest of code is serial, we will allow for
        // parallel implementation here
        NumMyElements = Map->NumMyElements();  // # of elements on local proc.
        // create vector that contains (global) indices of local elements
        //
        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        //
        // in our application, all elements have 6 off-diagonal terms
        //
        std::vector<int> NumNz(NumMyElements);
        for (int i = 0; i < NumMyElements; i++)
            NumNz[i] = 7;
        //
        // create matrix
        // 
        if (A != NULL) delete A;
        A = new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]);
        ierr = SetupCovLaplace(factor, u);
        ierr = A->FillComplete();
        //    cout << *A << endl;
        assert(ierr == 0);
        return ierr;
    }
    //================================================
    // Set up Solver (max_it on input is maximum number of iterations, num_it
    //                number of iterations used; returns residual)
    //================================================
    double Solve(int max_it, int& num_it, double tol) {
        Epetra_LinearProblem linprob(A, sol, rhs);
        AztecOO Solver(linprob);
        ofstream aztec_output;
        aztec_output.open("output/Aztec_Output_CovEllSolver");
        Solver.SetOutputStream(aztec_output);
        //    Solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
        //    Solver.SetAztecOption(AZ_precond, AZ_ilut);
        //    Solver.SetAztecOption(AZ_output, AZ_summary);
        Solver.SetAztecOption(AZ_conv, AZ_rhs);
        Solver.Iterate(max_it, tol);
        num_it = Solver.NumIters();
        aztec_output.close();
        return Solver.TrueResidual();
    }
    //================================================
    // Set right hand side
    //================================================
    int SetRHS(gf3d& rhs_grid) { return SetRHS(1.0, rhs_grid); };
    int SetRHS(double factor, gf3d& rhs_grid) {
        //    cout << " Setting up RHS..." << endl;
        for (int II = 0; II < n_global; II++) { // loop over superindex
            (*rhs)[II] = factor * rhs_grid(i_ind(II), j_ind(II), k_ind(II));
        }
        return 1;
    };
    //================================================
    // Zero initial guess
    //================================================
    int ZeroInitialGuess() {
        //    cout << " Setting up RHS..." << endl;
        for (int II = 0; II < n_global; II++) { // loop over superindex
            (*sol)[II] = 0.0;
        }
        return 1;
    };
    //================================================
    // Set initial guess
    //================================================
    int SetInitialGuess(gf3d& guess) {
        //    cout << " Setting up initial guess..." << endl;
        for (int II = 0; II < n_global; II++) { // loop over superindex
            (*sol)[II] = guess(i_ind(II), j_ind(II), k_ind(II));
        }
        return 1;
    };
    //================================================
    // Get solution
    //================================================
    void GetSolution(gf3d& sol_grid) {
        for (int i = N_g; i < n_r - N_g; i++) // loop over grid indices
            for (int j = N_g; j < n_theta - N_g; j++)
                for (int k = N_g; k < n_phi - N_g; k++)
                    sol_grid[i][j][k] = (*sol)[II_ind(i, j, k)];
        //
        // now fill outer ghostzones with correct falloff
        //
        //    int i_in = n_r - 4;
        int i_in = n_r - N_g - 1;
        for (int j = N_g; j < n_theta - N_g; j++)
            for (int k = N_g; k < n_phi - N_g; k++)
                for (int i = n_r - N_g; i < n_r; i++)
                    sol_grid[i][j][k] = fall_off_factor((*r)[i_in], (*r)[i]) * sol_grid(i_in, j, k);
        // int i = n_r - 3;
        // sol_grid[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid(i_in,j,k);
        // i = n_r - 2;
        // sol_grid[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid(i_in,j,k);
        // i = n_r - 1;
        // sol_grid[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid(i_in,j,k);

        sol_grid.fill_ghosts();
    };
private:
    //================================================
    // Set up covariant Laplace operator
    //================================================
    int SetupCovLaplace(double factor, gf3d& u) {
        int ierr = 0;
        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        // assume equidistant grid in phi direction
        const double d_x = grid->delta_x();
        const double d_y = grid->delta_y();
        const double d_phi = grid->delta_phi();
        //
        // set up vectors for indices and corresponding matrix entries 
        // of off-diagonal entries
        //
        std::vector<int> Indices(7);
        std::vector<double> Values(7);
        // fill matrix with derivative terms
        //
        for (int ii = 0; ii < NumMyElements; ii++) { // loop over superindex
            int II = MyGlobalElements[ii];               // global superindex
            //
            // get grid indices
            //
            int i = i_ind(II);
            int j = j_ind(II);
            int k = k_ind(II);
            //
            const double rl = grid->r(i);
            const double r2 = rl * rl;
            const double sintheta = grid->sintheta(j);
            const double sin2theta = sintheta * sintheta;
            const double costheta = grid->costheta(j);
            const double cottheta = costheta / sintheta;
            //
            // find auxiliary functions
            //
            const double dxdr = grid->dxdr(i);
            const double ddxdr = grid->ddxdr(i);
            const double dydtheta = grid->dydtheta(j);
            const double ddydtheta = grid->ddydtheta(j);
            //
            // compute Gamma^i's from background
            //
            const double Gam_r = -2.0 / rl;
            const double Gam_t = -cottheta / r2;
            //
            // define factors for derivatives in 2nd order differencing
            //
            const double ddx_term = 1.0 / (d_x * d_x) * dxdr * dxdr;
            const double ddy_term = 1.0 / (d_y * d_y) * dydtheta * dydtheta / r2;
            const double ddphi_term = 1.0 / (d_phi * d_phi) / (r2 * sin2theta);
            const double dx_term = 1.0 / (2.0 * d_x) * (ddxdr - Gam_r * dxdr);
            const double dy_term = 1.0 / (2.0 * d_y) * (ddydtheta / r2 - Gam_t * dydtheta);
            //
            // prepare diagonal element
            //
            double diag = -2.0 * (ddx_term + ddy_term + ddphi_term)
                + factor * u(i, j, k);
            //
            // now go through grid, but ignore i+1 terms for now
            //
            // i-1, j, k:
            Values[1] = ddx_term - dx_term;
            Indices[1] = II_ind(i - 1, j, k);
            // i, j-1, k:
            Values[2] = ddy_term - dy_term;
            Indices[2] = II_ind(i, j - 1, k);
            // i, j, k-1:
            Values[3] = ddphi_term;
            Indices[3] = II_ind(i, j, k - 1);
            // i, j+1, k:
            Values[4] = ddy_term + dy_term;
            Indices[4] = II_ind(i, j + 1, k);
            // i, j, k+1:
            Values[5] = ddphi_term;
            Indices[5] = II_ind(i, j, k + 1);
            //
            // rest depends on where i is...
            //
            if (i < n_r - N_g - 1) {   // i+1 is in interior of grid...
                // i+1, j, k:
                Values[6] = ddx_term + dx_term;
                Indices[6] = II_ind(i + 1, j, k);
                // diagonal element:
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 7, &Values[0], &Indices[0]);
                assert(ierr == 0);
            } else if (i == n_r - N_g - 1) { // i+1 is in ghost zone
                const double iterm = ddx_term + dx_term;
                const double factor = fall_off_factor((*r)[i], (*r)[i + 1]);
                // i+1, j, k corrects diagonal term:
                diag += factor * iterm;
                // diagonal element:
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 6, &Values[0], &Indices[0]);
                assert(ierr == 0);
            } else {
                cout << " FLATELLSOLVER3D: Ooops - should never have reached i = " << i << endl;
            }
        }
        return ierr;
    };
    //================================================
    // index utilies: translate between 
    // - grid indices i, j and k, running from N_g to (n_r, n_theta or n_phi) - N_g
    // - superindex, running from 0 to (n_r - 2 N_g)(n_theta - 2 N_g)(n_phi - 2 N_g)
    //
    // - for grid indices in outer boundary will use fall-off conditions (not done in this routine)
    // - for grid indices in all other ghost zones will use periodicity conditions
    //   to return corresponding superindex for interior grid point
    //
    //================================================
    inline int II_ind(int i, int j, int k) {
#ifdef AXISYMMETRY
        k = N_g;
#endif
        // first use periodicity for i indices in inner ghost zones...
        if (i >= n_r - N_g)
            cerr << " FLATELLSOLVER3D: II_ind should never be called in fall-off region! " << endl;
        //
        // deal with i in inner ghost zone
        // 
        if (i < N_g) {
            i = 2 * N_g - 1 - i;
#ifndef EQSYMMETRY
            j = n_theta - j - 1;                       // theta -> pi - theta
#endif
#ifndef AXISYMMETRY
            // CHECK INDEX: should this be n_phi?  (replaced np with n_phi to make compiler happy...)
            k = (k - N_g + (n_phi - 2 * N_g) / 2) % (n_phi - 2 * N_g) + N_g;  // phi -> phi + pi
#endif
        }
        //     if (i == 0) {
        //       cerr << " FLATELLSOLVER3D: why is II_ind called for i == 0??? " << endl;
        //       i = 5;
        // #ifndef EQSYMMETRY
        //       j = n_theta - j - 1;                       // theta -> pi - theta
        // #endif
        // #ifndef AXISYMMETRY
        //       // CHECK INDEX: should this be n_phi?  (replaced np with n_phi to make compiler happy...)
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
        //     if (i == 1) {
        //       i = 4;
        // #ifndef EQSYMMETRY
        //       j = n_theta - j - 1;                       // theta -> pi - theta
        // #endif
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
        //     if (i == 2) {
        //       i = 3;
        // #ifndef EQSYMMETRY
        //       j = n_theta - j - 1;                       // theta -> pi - theta
        // #endif
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
            //
            // then use periodicity for j indices in ghost zones...
            //
        if (j < N_g) {
            j = 2 * N_g - 1 - j;
#ifndef AXISYMMETRY
            k = (k - N_g + (n_phi - 2 * N_g) / 2) % (n_phi - 2 * N_g) + N_g;  // phi -> phi + pi
#endif
        }
        //     if (j == 0) {
        //       j = 5;
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
        //     if (j == 1) {
        //       j = 4;
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
        //     if (j == 2) {
        //       j = 3;
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
        if (j >= n_theta - N_g) {
            int delta_j = j - (n_theta - N_g);
            j = n_theta - N_g - 1 - delta_j;
#ifndef AXISYMMETRY
            k = (k - N_g + (n_phi - 2 * N_g) / 2) % (n_phi - 2 * N_g) + N_g;  // phi -> phi + pi
#endif
        }
        //     if (j == n_theta-3) {
        //       j = n_theta-4;
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
        //     if (j == n_theta-2) {
        //       j = n_theta-5;
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
        //     if (j == n_theta-1) {
        //       j = n_theta-6;
        // #ifndef AXISYMMETRY
        //       k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
        // #endif
        //     }
            // ... and finally use periodicity for k indices in ghost zones...
        if (k < N_g) k = n_phi - 2 * N_g + k;
        if (k >= n_phi - N_g) k = N_g + (k - (n_phi - N_g));
        // if (k == 0) k = n_phi-6;
        // if (k == 1) k = n_phi-5;
        // if (k == 2) k = n_phi-4;
        // if (k == n_phi-1) k = 5;
        // if (k == n_phi-2) k = 4;
        // if (k == n_phi-3) k = 3;
        // now return corresponding superindex
        return (i - N_g) + (j - N_g) * (n_r - 2 * N_g) + (k - N_g) * (n_r - 2 * N_g) * (n_theta - 2 * N_g);
    };
    inline int i_ind(int II) { return (II % ((n_r - 2 * N_g) * (n_theta - 2 * N_g))) % (n_r - 2 * N_g) + N_g; };
    inline int j_ind(int II) { return (II % ((n_r - 2 * N_g) * (n_theta - 2 * N_g))) / (n_r - 2 * N_g) + N_g; };
    inline int k_ind(int II) { return II / ((n_r - 2 * N_g) * (n_theta - 2 * N_g)) + N_g; };
    //===============================================================
    // provide fall_off factor
    //===============================================================
    double fall_off_factor(double r_in, double r_out) {
        const double ratio = r_in / r_out;
        if (fall_off == 1) return ratio;
        else if (fall_off == 2) return ratio * ratio;
        else return pow(ratio, fall_off);
    };
public:
    //===============================================================
    //
    // As a test, solve the problem
    //
    //  \nabla^2 f + \lambda f = (A + \lambda) Y^{11}
    //
    //===============================================================

    void Test(double AA, gf3d& rhs_grid, gf3d& sol_grid, gf3d& res) {
        //
        // First set up Laplace operator
      //   // 
      //   double lambda = 1.0;
      //   SetupSolver(lambda);
      //   for (int j = N_g; j < n_theta-N_g; j++) // loop over grid indices
      //     for (int k = N_g; k < n_phi-N_g; k++)
      // 	rhs_grid[j][k] = ( AA + lambda ) * sin(sol_grid.theta(j)) * cos(sol_grid.phi(k));
      //   SetRHS(rhs_grid);
      //   //
      //   // now solve
      //   //
      //   int max_it = 100;
      //   int num_it = 0;
      //   double tol = 1.e-8;

      //   double res_tri = Solve(max_it,num_it,tol);
      //   cout << " Residual after " << num_it << " iteration steps = " << res_tri << endl;	  
      //   //
      //   // get solution
      //   //    
      //   GetSolution(sol_grid);
      //   sol_grid.fill_ghosts();
      //   //
      //   // Let's compute residual ourselvers...
      //   //
      //   for (int j = N_g; j < n_theta-N_g; j++) // loop over grid indices
      //     for (int k = N_g; k < n_phi-N_g; k++) {
      // 	res[j][k] = sol_grid.Laplace(j,k) + lambda*sol_grid(j,k) - rhs_grid(j,k);
      //     }
      //   cout << " Our residual = " << res.L2_norm() << endl;
      // }  
    };
};

#endif  /* _FLATELLSOLVER3D_H_ */

#endif  /* NoEllSolver */
