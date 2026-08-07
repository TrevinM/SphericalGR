// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing 2D elliptic solver
//
//================================================
// #define NoEllSolver

#ifndef NoEllSolver

//================================================
//
// Uses Trilinos software to solve flat 2D Laplace
// equation on unit sphere
//
//  \nabla^2 f + \lambda f 
//     = (\partial^2_\theta + \cot \theta \partial_\theta 
//       + \sin^{-2} \theta \partial^2_\phi + \lambda) f 
//     = ( (dx/dtheta)^2 \partial^2_x 
//       + ( dx/d\theta \cot theta + d^2/d\theta^2) partial_x  
//       + \sin^{-2} \theta \partial^2_\phi + \lambda) f 
//     =  RHS
//
// where \lambda is a constant.
//
// Boundary conditions are completely determined by
// periodicity; these are implemented in index utilities.
//
//================================================

//================================================
//
// Important note concerning grid:
//
// n_theta and n_phi include N_g ghost_zones, which 
// are not needed for elliptic solver.  In interior
// grid, have N_g <= j < n_theta-N_g (and similar for
// k).  However, want "superindex" II go from 0 to its
// maximum value, (n_theta - 2 N_g)*(n_phi - 2 N_g) - 1.
//
//================================================

#ifndef _ELLSOLVER2D_H_
#define _ELLSOLVER2D_H_


#include <cstdlib>
#include <cassert>
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
#include "Grid.h"

using namespace std;

class EllSolver2D {
private:
    int n_theta;
    int n_phi;
    int n_global;
    int N_g;
    int NumMyElements;
    //  VecDoub *theta, *phi;
#ifdef EPETRA_MPI
    Epetra_MpiComm* Comm;
#else
    Epetra_SerialComm* Comm;
#endif
    Epetra_Map* Map;
    Epetra_Vector* rhs, * sol;
    Epetra_CrsMatrix* A;
    Grid* grid;
public:
    //================================================
    // Constructor
    //================================================
    EllSolver2D(Grid* grid_i) : grid(grid_i) {
        N_g = grid->N_ghosts();
        n_theta = grid->N_theta_tot();
        n_phi = grid->N_phi_tot();
        // theta = grid->theta();
        // phi = grid->phi();
        cout << " Constructing EllSolver2D";
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
        const int NumGlobalElements = (n_theta - 2 * N_g) * (n_phi - 2 * N_g);
        cout << " for " << NumGlobalElements << "x"
            << NumGlobalElements << " matrix..." << endl;
        Map = new Epetra_Map(NumGlobalElements, 0, *Comm);
        n_global = NumGlobalElements;
        //
        // ... then allocate vectors
        // 
        rhs = new Epetra_Vector(*Map);
        sol = new Epetra_Vector(*Map);
        // for (int j = N_g; j < n_theta-N_g; j++) // loop over grid indices
        //   for (int k = N_g; k < n_phi-N_g; k++) 
        // 	cout << j << "  " << k << "  " << II_ind(j,k) << "  " 
        // 	     << j_ind(II_ind(j,k)) << "  " << k_ind(II_ind(j,k)) << endl;
        // even though rest of code is serial, we will allow for
        // parallel implementation here
        NumMyElements = Map->NumMyElements();  // # of elements on local proc.
        // create vector that contains (global) indices of local elements

        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        //
        // in our application, all elements have 8 off-diagonal terms
        //
        std::vector<int> NumNz(NumMyElements);
        for (int i = 0; i < NumMyElements; i++)
            NumNz[i] = 9;
        //
        // create matrix
        // 
        A = new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]);
    };
    //================================================
    // Destructor
    //================================================
    ~EllSolver2D() {
        cout << " Destructing EllSolver2D... " << endl;
        delete Comm;
        delete Map;
        delete rhs;
        delete sol;
        delete A;
    };
    //================================================
    // Set up Solver
    //================================================
    int SetupSolver(double lambda) {
        int ierr = 0;
        ierr = SetupFlatLaplace(lambda);
        ierr = A->FillComplete();
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
        aztec_output.open("Aztec_Output");
        //    Solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
        //    Solver.SetAztecOption(AZ_precond, AZ_ilut);
        //    Solver.SetAztecOption(AZ_output, AZ_summary);
        Solver.SetAztecOption(AZ_conv, AZ_rhs);
        Solver.SetOutputStream(aztec_output);
        Solver.Iterate(max_it, tol);
        num_it = Solver.NumIters();
        aztec_output.close();
        return Solver.TrueResidual();
    }

    //================================================
    // Set right hand side
    //================================================
    void SetRHS(gf2d& rhs_func) {
        //    cout << " Setting up RHS..." << endl;
        for (int II = 0; II < n_global; II++) { // loop over superindex
            (*rhs)[II] = rhs_func(j_ind(II), k_ind(II));
        }
    };
    //================================================
    // Get solution
    //================================================
    void GetSolution(gf2d& h) {
        for (int j = N_g; j < n_theta - N_g; j++) // loop over grid indices
            for (int k = N_g; k < n_phi - N_g; k++)
                h[j][k] = (*sol)[II_ind(j, k)];
        h.fill_ghosts();
    };
private:
    //================================================
    // Set up flat Laplace operator
    //================================================
    int SetupFlatLaplace(double lambda) {
        int ierr = 0;
        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        // assume equidistant grid in phi, but not in theta...
        const double d_phi = grid->delta_phi();
        const double d_y = grid->delta_y();
        // define factors for derivatives in 4th order differencing
        const double ddphi_fact = 1.0 / (12.0 * d_phi * d_phi);
        const double ddy_fact = 1.0 / (12.0 * d_y * d_y);
        const double dy_fact = 1.0 / (12.0 * d_y);
        //
        // set up vectors for indices and corresponding matrix entries 
        // of off-diagonal entries
        //
        std::vector<int> Indices(9);
        std::vector<double> Values(9);
        // fill matrix with derivative terms
        //
        for (int ii = 0; ii < NumMyElements; ii++) { // loop over superindex
            int II = MyGlobalElements[ii];               // global superindex
            //
            // get grid indices
            //
            int j = j_ind(II);
            int k = k_ind(II);
            const double sintheta = grid->sintheta(j);
            const double sin2theta = sintheta * sintheta;
            const double cottheta = grid->costheta(j) / sintheta;
            const double dydtheta = grid->dydtheta(j);
            const double ddydtheta = grid->ddydtheta(j);
            const double ddphi_term = ddphi_fact / sin2theta;
            const double ddy_term = ddy_fact * dydtheta * dydtheta;
            const double dy_term = dy_fact * (cottheta * dydtheta + ddydtheta);
            // k - 2:
            Values[0] = -ddphi_term;
            Indices[0] = II_ind(j, k - 2);
            // k - 1:
            Values[1] = 16.0 * ddphi_term;
            Indices[1] = II_ind(j, k - 1);
            // j - 2
            Values[2] = -ddy_term + dy_term;
            Indices[2] = II_ind(j - 2, k);
            // j - 1
            Values[3] = 16.0 * ddy_term - 8.0 * dy_term;
            Indices[3] = II_ind(j - 1, k);
            // diagonal element
            Values[4] = -30.0 * (ddy_term + ddphi_term) + lambda;
            Indices[4] = II;
            // j + 1
            Values[5] = 16.0 * ddy_term + 8.0 * dy_term;
            Indices[5] = II_ind(j + 1, k);
            // j + 2
            Values[6] = -ddy_term - dy_term;
            Indices[6] = II_ind(j + 2, k);
            // k + 1:
            Values[7] = 16.0 * ddphi_term;
            Indices[7] = II_ind(j, k + 1);
            // k + 2:
            Values[8] = -ddphi_term;
            Indices[8] = II_ind(j, k + 2);
            // 
            // now store matrix elements
            //
            ierr = A->InsertGlobalValues(MyGlobalElements[ii], 9, &Values[0], &Indices[0]);
            assert(ierr == 0);
        }
        return ierr;
    };
    //================================================
    // index utilies: translate between 
    // - grid indices j and k, running from N_g t0 (n_theta or n_phi) - N_g
    // - superindex, running from 0 to (n_theta - 2 N_g)(n_phi - 2 N_g)
    //
    // for grid indices in ghost zones will use periodicity conditions
    //   to return corresponding superindex for interior grid point
    //
    //================================================
    inline int II_ind(int j, int k) {
        // first use periodicity for j indices in ghost zones...
#ifdef AXISYMMETRY
        k = N_g;
#endif
        if (j < N_g) {
            //      cout << " j = " << j << " becomes j = " << 2*N_g -1 - j << endl;
            j = 2 * N_g - 1 - j;
#ifndef AXISYMMETRY
            k = (k - N_g + (n_phi - 2 * N_g) / 2) % (n_phi - 2 * N_g) + N_g;  // phi -> phi + pi
#endif
        }
        //     if (j == 0) {
        //       j = 5;
        // #ifndef AXISYMMETRY
        //       // CHECK INDEX: should this be n_phi?  (replaced np with n_phi to make compiler happy...)
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
                //
            // then use periodicity for j indices in ghost zones...
            //
        if (j >= n_theta - N_g) {
            int delta_j = j - (n_theta - N_g);
            //      cout << " j = " << j << " becomes j = " << n_theta - N_g - 1 - delta_j << endl;
            j = n_theta - N_g - 1 - delta_j;
#ifndef AXISYMMETRY
            k = (k - N_g + (n_phi - 2 * N_g) / 2) % (n_phi - 2 * N_g) + N_g;  // phi -> phi + pi
#endif
        }
        //    if (j == n_theta-3) {
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
           // then use periodicity for k indices in ghost zones...
        if (k < N_g) {
            k = n_phi - 2 * N_g + k;
        }
        if (k >= n_phi - N_g) {
            k = N_g + (k - (n_phi - N_g));
        }
        // if (k == 0) k = n_phi-6;
        // if (k == 1) k = n_phi-5;
        // if (k == 2) k = n_phi-4;
        // if (k == n_phi-1) k = 5;
        // if (k == n_phi-2) k = 4;
        // if (k == n_phi-3) k = 3;
        //
        // now return corresponding superindex
        return (j - N_g) + (k - N_g) * (n_theta - 2 * N_g);
    };
    inline int j_ind(int II) { return II % (n_theta - 2 * N_g) + N_g; };
    inline int k_ind(int II) { return II / (n_theta - 2 * N_g) + N_g; };
public:
    //===============================================================
    //
    // As a test, solve the problem
    //
    //  \nabla^2 f + \lambda f = (A + \lambda) Y^{11}
    //
    //===============================================================

    void Test(double AA, gf2d& rhs_grid, gf2d& sol_grid, gf2d& res) {
        //
        // First set up Laplace operator
        // 
        double lambda = 1.0;
        SetupSolver(lambda);
        for (int j = N_g; j < n_theta - N_g; j++) // loop over grid indices
            for (int k = N_g; k < n_phi - N_g; k++)
                rhs_grid[j][k] = (AA + lambda) * sin(sol_grid.theta(j)) * cos(sol_grid.phi(k));
        SetRHS(rhs_grid);
        //
        // now solve
        //
        int max_it = 100;
        int num_it = 0;
        double tol = 1.e-8;

        double res_tri = Solve(max_it, num_it, tol);
        cout << " Residual after " << num_it << " iteration steps = " << res_tri << endl;
        //
        // get solution
        //    
        GetSolution(sol_grid);
        sol_grid.fill_ghosts();
        //
        // Let's compute residual ourselvers...
        //
        for (int j = N_g; j < n_theta - N_g; j++) // loop over grid indices
            for (int k = N_g; k < n_phi - N_g; k++) {
                res[j][k] = sol_grid.Laplace(j, k) + lambda * sol_grid(j, k) - rhs_grid(j, k);
            }
        cout << " Our residual = " << res.L2_norm() << endl;
    }
};

#endif  /* _ELLSOLVER2D_H_ */

#endif /* NoEllSolver */
