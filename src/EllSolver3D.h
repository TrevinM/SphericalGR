// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing 3D elliptic solver
//
// OLD VERSION, with N_g hardcoded to 2!!
//
//================================================

#ifndef NoEllSolver

//================================================
//
// Uses Trilinos software to solve flat 3D Laplace
// operator in spherical polar coordinates
//
//  \nabla^2 f + factor * u * f 
//     = (\partial^2_r + 2 \partial_r / r + 
//        r^{-2} \partial^2_\theta + r^{-2} \cot \theta \partial_\theta + 
//        r^{-2} sin^{-2} \theta \partial^2_\phi + factor * u) f =  RHS
//
// or, alternatively, with extra first derivative term
//
//  \nabla^2 f + V^i D_i f + \lambda f = RHS 
//
// where u is a function, V^i a vector, and factor is a constant
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
//
//================================================

#ifndef _ELLSOLVER3D_H_
#define _ELLSOLVER3D_H_


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

class EllSolver3D {
private:
    int n_r;
    int n_theta;
    int n_phi;
    int n_global;
    int N_g;      // number of ghost zones - hardcoded to two
    int NumMyElements;
    VecDoub* r, * theta, * phi;
    int fall_off;   // default set to 1 - change with SetFallOff()
    bool verbose;
#ifdef EPETRA_MPI
    Epetra_MpiComm* Comm;
#else
    Epetra_SerialComm* Comm;
#endif
    Epetra_Map* Map;
    Epetra_Vector* rhs, * sol;
    Epetra_CrsMatrix* A;
public:
    //================================================
    // Constructor
    //================================================
    EllSolver3D(int n_r_i, int n_theta_i, int n_phi_i, VecDoub* r_i, VecDoub* theta_i, VecDoub* phi_i, bool verbose_i = false) :
        n_r(n_r_i), n_theta(n_theta_i), n_phi(n_phi_i), N_g(2),
        r(r_i), theta(theta_i), phi(phi_i),
        fall_off(1), verbose(verbose_i) {
        if (verbose) cout << " Constructing EllSolver3D";
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
        if (verbose)  cout << " for " << NumGlobalElements << "x"
            << NumGlobalElements << " matrix..." << endl;
        Map = new Epetra_Map(NumGlobalElements, 0, *Comm);
        n_global = NumGlobalElements;
        //
        // ... then allocate vectors
        // 
        rhs = new Epetra_Vector(*Map);
        sol = new Epetra_Vector(*Map);
        //
        // check for grid indices...
        //
        // for (int i = N_g; i < n_r-N_g; i++) // loop over grid indices
        //   for (int j = N_g; j < n_theta-N_g; j++) 
        // 	for (int k = N_g; k < n_phi-N_g; k++) 
        // 	  cout << i << "  " << j << "  " << k << "  " << II_ind(i,j,k) << "  " << i_ind(II_ind(i,j,k)) << "  "
        // 	       << j_ind(II_ind(i,j,k)) << "  " << k_ind(II_ind(i,j,k)) << endl;
        //
        A = NULL;
    };
    //================================================
    // Destructor
    //================================================
    ~EllSolver3D() {
        if (verbose) cout << " Destructing EllSolver3D... " << endl;
        delete Comm;
        delete Map;
        delete rhs;
        delete A;
    };
    //================================================
    // Change fall-off (set to 1 by default)
    //================================================
    int SetFallOff(int fall_off_i) { return fall_off = fall_off_i; }
    //================================================
    // Set up Solver (three versions, depending on whether first 
    // derivative term or metric coefficients are included)
    //================================================
    int SetupSolver(double factor, gf3d& u) {
        int ierr = 0;
        // even though rest of code is serial, we will allow for
        // parallel implementation here
        NumMyElements = Map->NumMyElements();  // # of elements on local proc.
        // create vector that contains (global) indices of local elements

        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        //
        // in our application, almost all elements have 12 off-diagonal terms
        //
        std::vector<int> NumNz(NumMyElements);
        for (int i = 0; i < NumMyElements; i++)
            NumNz[i] = 13;
        //
        // create matrix
        // 
        if (A != NULL) delete A;
        A = new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]);
        ierr = SetupFlatLaplace(factor, u);
        ierr = A->FillComplete();
        //    cout << *A << endl;
        assert(ierr == 0);
        return ierr;
    }
    int SetupSolver(gf3d& Vr, gf3d& Vt, gf3d& Vp, double factor, gf3d& u) {
        int ierr = 0;
        // even though rest of code is serial, we will allow for
        // parallel implementation here
        NumMyElements = Map->NumMyElements();  // # of elements on local proc.
        // create vector that contains (global) indices of local elements

        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        //
        // in our application, almost all elements have 12 off-diagonal terms
        //
        std::vector<int> NumNz(NumMyElements);
        for (int i = 0; i < NumMyElements; i++)
            NumNz[i] = 13;
        //
        // create matrix
        // 
        if (A != NULL) delete A;
        A = new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]);
        ierr = SetupFlatLaplace(Vr, Vt, Vp, factor, u);
        ierr = A->FillComplete();
        //    cout << *A << endl;
        assert(ierr == 0);
        return ierr;
    }
    int SetupSolver(gf3d& gup_rr, gf3d& gup_rt, gf3d& gup_rp,
        gf3d& gup_tt, gf3d& gup_tp, gf3d& gup_pp,
        gf3d& DG_r, gf3d& DG_t, gf3d& DG_p,
        gf3d& Vr, gf3d& Vt, gf3d& Vp, double factor, gf3d& u) {
        int ierr = 0;
        if (verbose) cout << " In SetupSolver for covariant Laplace operator " << endl;
        // even though rest of code is serial, we will allow for
        // parallel implementation here
        NumMyElements = Map->NumMyElements();  // # of elements on local proc.
        // create vector that contains (global) indices of local elements

        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        //
        // in our application, almost all elements have 24 off-diagonal terms
        //
        std::vector<int> NumNz(NumMyElements);
        for (int i = 0; i < NumMyElements; i++)
            NumNz[i] = 25;
        //
        // create matrix
        // 
        if (verbose) cout << " Setting up matrix... ";
        if (A != NULL) delete A;
        A = new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]);
        if (verbose) cout << "... done! " << endl;
        ierr = SetupFlatLaplace(gup_rr, gup_rt, gup_rp, gup_tt, gup_tp, gup_pp, DG_r, DG_t, DG_p, Vr, Vt, Vp, factor, u);
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
        aztec_output.open("Aztec_Output_EllSolver3D");
        // Solver.SetAztecOption(AZ_precond, AZ_ilut);
        //    Solver.SetAztecOption(AZ_solver, AZ_cgs);
        //    Solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
        //    Solver.SetAztecOption(AZ_conv,AZ_Anorm);
        Solver.SetOutputStream(aztec_output);
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
        int i_in = n_r - 3;
        for (int j = N_g; j < n_theta - N_g; j++)
            for (int k = N_g; k < n_phi - N_g; k++) {
                int i = n_r - 2;
                sol_grid[i][j][k] = fall_off_factor((*r)[i_in], (*r)[i]) * sol_grid(i_in, j, k);
                i = n_r - 1;
                sol_grid[i][j][k] = fall_off_factor((*r)[i_in], (*r)[i]) * sol_grid(i_in, j, k);
            }
        sol_grid.fill_ghosts();
    };
private:
    //================================================
    // Set up flat Laplace operator (again, three versions...)
    //
    // First version, without first-derivative terms
    //================================================
    int SetupFlatLaplace(double factor, gf3d& u) {
        int ierr = 0;
        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        // assume equidistant grid
        const double d_r = (*r)[1] - (*r)[0];
        const double d_phi = (*phi)[1] - (*phi)[0];
        const double d_theta = (*theta)[1] - (*theta)[0];
        // define factors for second derivatives in 4th order differencing
        const double ddr_fact = 1.0 / (12.0 * d_r * d_r);
        const double dr_fact = 1.0 / (12.0 * d_r);
        const double ddphi_fact = 1.0 / (12.0 * d_phi * d_phi);
        const double ddtheta_fact = 1.0 / (12.0 * d_theta * d_theta);
        const double dtheta_fact = 1.0 / (12.0 * d_theta);
        //
        // set up vectors for indices and corresponding matrix entries 
        // of off-diagonal entries
        //
        std::vector<int> Indices(13);
        std::vector<double> Values(13);
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
            const double rl = (*r)[i];
            const double r2 = rl * rl;
            const double sintheta = sin((*theta)[j]);
            const double sin2theta = sintheta * sintheta;
            const double cottheta = cos((*theta)[j]) / sintheta;
            //
            // prepare diagonal element
            //
            double diag = -30.0 * (ddr_fact + ddtheta_fact / r2 + ddphi_fact / (r2 * sin2theta)) + factor * u(i, j, k);
            // k - 2:
            Values[1] = -ddphi_fact / (r2 * sin2theta);
            Indices[1] = II_ind(i, j, k - 2);
            // k - 1:
            Values[2] = 16.0 * ddphi_fact / (r2 * sin2theta);
            Indices[2] = II_ind(i, j, k - 1);
            // j - 2
            Values[3] = -ddtheta_fact / r2 + cottheta / r2 * dtheta_fact;
            Indices[3] = II_ind(i, j - 2, k);
            // j - 1
            Values[4] = 16.0 * ddtheta_fact / r2 - 8.0 * cottheta / r2 * dtheta_fact;
            Indices[4] = II_ind(i, j - 1, k);
            // i - 2:
            Values[5] = -ddr_fact + 2.0 / rl * dr_fact;
            Indices[5] = II_ind(i - 2, j, k);
            // i - 1:
            Values[6] = 16.0 * ddr_fact - 8.0 * 2.0 / rl * dr_fact;
            Indices[6] = II_ind(i - 1, j, k);
            // j + 1
            Values[7] = 16.0 * ddtheta_fact / r2 + 8.0 * cottheta / r2 * dtheta_fact;
            Indices[7] = II_ind(i, j + 1, k);
            // j + 2
            Values[8] = -ddtheta_fact / r2 - cottheta / r2 * dtheta_fact;
            Indices[8] = II_ind(i, j + 2, k);
            // k + 1:
            Values[9] = 16.0 * ddphi_fact / (r2 * sin2theta);
            Indices[9] = II_ind(i, j, k + 1);
            // k + 2:
            Values[10] = -ddphi_fact / (r2 * sin2theta);
            Indices[10] = II_ind(i, j, k + 2);
            //
            // rest depends on where i is...
            //
            if (i < n_r - N_g - 2) {   // all neighbors are in interior of grid
                // i + 1:
                Values[11] = 16.0 * ddr_fact + 8.0 * 2.0 / rl * dr_fact;
                Indices[11] = II_ind(i + 1, j, k);
                // i + 2:
                Values[12] = -ddr_fact - 2.0 / rl * dr_fact;
                Indices[12] = II_ind(i + 2, j, k);
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 13, &Values[0], &Indices[0]);
                assert(ierr == 0);
            } else if (i == n_r - N_g - 2) {  // point i+2 is in ghost zone
                // i + 1:
                Values[11] = 16.0 * ddr_fact + 8.0 * 2.0 / rl * dr_fact;
                Indices[11] = II_ind(i + 1, j, k);
                // i + 2:
                double ip2_term = -ddr_fact - 2.0 / rl * dr_fact;
                Values[11] += fall_off_factor((*r)[i + 1], (*r)[i + 2]) * ip2_term;
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 12, &Values[0], &Indices[0]);
            } else if (i == n_r - N_g - 1) { // both upper neighbors in ghost zone
                double ip1_term = 16.0 * ddr_fact + 8.0 * 2.0 / rl * dr_fact;
                double ip2_term = -ddr_fact - 2.0 / rl * dr_fact;
                diag += fall_off_factor((*r)[i], (*r)[i + 2]) * ip2_term + fall_off_factor((*r)[i], (*r)[i + 1]) * ip1_term;
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 11, &Values[0], &Indices[0]);
            } else {
                cout << " ELLSOLVER3D: Ooops - should never have reached i = " << i << endl;
            }
        }
        return ierr;
    };
    //================================================
    // Set up flat Laplace operator (again, three versions...)
    //
    // Second version, with first-derivative terms
    //================================================
    int SetupFlatLaplace(gf3d& Vr, gf3d& Vt, gf3d& Vp, double factor, gf3d& u) {
        int ierr = 0;
        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        // assume equidistant grid
        const double d_r = (*r)[1] - (*r)[0];
        const double d_phi = (*phi)[1] - (*phi)[0];
        const double d_theta = (*theta)[1] - (*theta)[0];
        // define factors for second derivatives in 4th order differencing
        const double ddr_fact = 1.0 / (12.0 * d_r * d_r);
        const double dr_fact = 1.0 / (12.0 * d_r);
        const double ddphi_fact = 1.0 / (12.0 * d_phi * d_phi);
        const double dphi_fact = 1.0 / (12.0 * d_phi);
        const double ddtheta_fact = 1.0 / (12.0 * d_theta * d_theta);
        const double dtheta_fact = 1.0 / (12.0 * d_theta);
        //
        // set up vectors for indices and corresponding matrix entries 
        // of off-diagonal entries
        //
        std::vector<int> Indices(13);
        std::vector<double> Values(13);
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
            const double rl = (*r)[i];
            const double r2 = rl * rl;
            const double sintheta = sin((*theta)[j]);
            const double sin2theta = sintheta * sintheta;
            const double cottheta = cos((*theta)[j]) / sintheta;
            //
            // prepare diagonal element
            //
            double diag = -30.0 * (ddr_fact + ddtheta_fact / r2 + ddphi_fact / (r2 * sin2theta)) + factor * u(i, j, k);
            // k - 2:
            Values[1] = -ddphi_fact / (r2 * sin2theta) + Vp(i, j, k) * dphi_fact;
            Indices[1] = II_ind(i, j, k - 2);
            // k - 1:
            Values[2] = 16.0 * ddphi_fact / (r2 * sin2theta) - 8.0 * Vp(i, j, k) * dphi_fact;
            Indices[2] = II_ind(i, j, k - 1);
            // j - 2
            Values[3] = -ddtheta_fact / r2 + (cottheta / r2 + Vt(i, j, k)) * dtheta_fact;
            Indices[3] = II_ind(i, j - 2, k);
            // j - 1
            Values[4] = 16.0 * ddtheta_fact / r2 - 8.0 * (cottheta / r2 + Vt(i, j, k)) * dtheta_fact;
            Indices[4] = II_ind(i, j - 1, k);
            // i - 2:
            Values[5] = -ddr_fact + (2.0 / rl + Vr(i, j, k)) * dr_fact;
            Indices[5] = II_ind(i - 2, j, k);
            // i - 1:
            Values[6] = 16.0 * ddr_fact - 8.0 * (2.0 / rl + Vr(i, j, k)) * dr_fact;
            Indices[6] = II_ind(i - 1, j, k);
            // j + 1
            Values[7] = 16.0 * ddtheta_fact / r2 + 8.0 * (cottheta / r2 + Vt(i, j, k)) * dtheta_fact;
            Indices[7] = II_ind(i, j + 1, k);
            // j + 2
            Values[8] = -ddtheta_fact / r2 - (cottheta / r2 + Vt(i, j, k)) * dtheta_fact;
            Indices[8] = II_ind(i, j + 2, k);
            // k + 1:
            Values[9] = 16.0 * ddphi_fact / (r2 * sin2theta) + 8.0 * Vp(i, j, k) * dphi_fact;
            Indices[9] = II_ind(i, j, k + 1);
            // k + 2:
            Values[10] = -ddphi_fact / (r2 * sin2theta) - Vp(i, j, k) * dphi_fact;
            Indices[10] = II_ind(i, j, k + 2);
            //
            // rest depends on where i is...
            //
            if (i < n_r - N_g - 2) {   // all neighbors are in interior of grid
                // i + 1:
                Values[11] = 16.0 * ddr_fact + 8.0 * (2.0 / rl + Vr(i, j, k)) * dr_fact;
                Indices[11] = II_ind(i + 1, j, k);
                // i + 2:
                Values[12] = -ddr_fact - (2.0 / rl + Vr(i, j, k)) * dr_fact;
                Indices[12] = II_ind(i + 2, j, k);
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 13, &Values[0], &Indices[0]);
                assert(ierr == 0);
            } else if (i == n_r - N_g - 2) {  // point i+2 is in ghost zone
                // i + 1:
                Values[11] = 16.0 * ddr_fact + 8.0 * (2.0 / rl + Vr(i, j, k)) * dr_fact;
                Indices[11] = II_ind(i + 1, j, k);
                // i + 2:
                double ip2_term = -ddr_fact - (2.0 / rl + Vr(i, j, k)) * dr_fact;
                Values[11] += fall_off_factor((*r)[i + 1], (*r)[i + 2]) * ip2_term;
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 12, &Values[0], &Indices[0]);
            } else if (i == n_r - N_g - 1) { // both upper neighbors in ghost zone
                double ip1_term = 16.0 * ddr_fact + 8.0 * (2.0 / rl + Vr(i, j, k)) * dr_fact;
                double ip2_term = -ddr_fact - (2.0 / rl + Vr(i, j, k)) * dr_fact;
                diag += fall_off_factor((*r)[i], (*r)[i + 2]) * ip2_term + fall_off_factor((*r)[i], (*r)[i + 1]) * ip1_term;
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 11, &Values[0], &Indices[0]);
            } else {
                cout << " ELLSOLVER3D: Ooops - should never have reached i = " << i << endl;
            }
        }
        return ierr;
    };
    //================================================
    // Set up flat Laplace operator (again, three versions...)
    //
    // Third version: covariant Laplace operator, except that mixed
    // derivatives are treated only to second order
    //================================================
    int SetupFlatLaplace(gf3d& gup_rr, gf3d& gup_rt, gf3d& gup_rp,
        gf3d& gup_tt, gf3d& gup_tp, gf3d& gup_pp,
        gf3d& DG_r, gf3d& DG_t, gf3d& DG_p,
        gf3d& Vr, gf3d& Vt, gf3d& Vp, double factor, gf3d& u) {
        int ierr = 0;
        std::vector<int> MyGlobalElements(NumMyElements);
        Map->MyGlobalElements(&MyGlobalElements[0]);
        // assume equidistant grid
        const double d_r = (*r)[1] - (*r)[0];
        const double d_phi = (*phi)[1] - (*phi)[0];
        const double d_theta = (*theta)[1] - (*theta)[0];
        // define factors for second derivatives in 4th order differencing
        const double ddr_fact = 1.0 / (12.0 * d_r * d_r);
        const double dr_fact = 1.0 / (12.0 * d_r);
        const double ddphi_fact = 1.0 / (12.0 * d_phi * d_phi);
        const double dphi_fact = 1.0 / (12.0 * d_phi);
        const double ddtheta_fact = 1.0 / (12.0 * d_theta * d_theta);
        const double dtheta_fact = 1.0 / (12.0 * d_theta);
        // ... and mixed derivatives in 2nd order differencing
        const double drdtheta_fact = 0.25 / (d_r * d_theta);
        const double drdphi_fact = 0.25 / (d_r * d_phi);
        const double dthetadphi_fact = 0.25 / (d_theta * d_phi);
        //
         // set up vectors for indices and corresponding matrix entries 
         // of off-diagonal entries
         //
        std::vector<int> Indices(25);
        std::vector<double> Values(25);
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
            const double rl = (*r)[i];
            // const double r2 = rl*rl;
            const double sintheta = sin((*theta)[j]);
            const double sin2theta = sintheta * sintheta;
            const double costheta = cos((*theta)[j]);
            const double cottheta = costheta / sintheta;
            //
            // compute Gamma^i's from rescaled DG's and background
            //
            const double Gam_r = DG_r(i, j, k) - gup_tt(i, j, k) * rl - gup_pp(i, j, k) * rl * sin2theta;
            const double Gam_t = DG_t(i, j, k) / rl + 2.0 * gup_rt(i, j, k) / rl - gup_pp(i, j, k) * sintheta * costheta;
            const double Gam_p = DG_p(i, j, k) / (rl * sintheta) + 2.0 * gup_rp(i, j, k) / rl + 2.0 * gup_tp(i, j, k) * cottheta;
            //
           // prepare diagonal element
           //
            double diag = -30.0 * (gup_rr(i, j, k) * ddr_fact +
                gup_tt(i, j, k) * ddtheta_fact +
                gup_pp(i, j, k) * ddphi_fact) + factor * u(i, j, k);
            // k - 2:
            Values[1] = -gup_pp(i, j, k) * ddphi_fact + (Vp(i, j, k) - Gam_p) * dphi_fact;
            Indices[1] = II_ind(i, j, k - 2);
            // k - 1:
            Values[2] = 16.0 * gup_pp(i, j, k) * ddphi_fact - 8.0 * (Vp(i, j, k) - Gam_p) * dphi_fact;
            Indices[2] = II_ind(i, j, k - 1);
            // j - 2
            Values[3] = -gup_tt(i, j, k) * ddtheta_fact + (Vt(i, j, k) - Gam_t) * dtheta_fact;
            Indices[3] = II_ind(i, j - 2, k);
            // j - 1
            Values[4] = 16.0 * gup_tt(i, j, k) * ddtheta_fact - 8.0 * (Vt(i, j, k) - Gam_t) * dtheta_fact;
            Indices[4] = II_ind(i, j - 1, k);
            // i - 2:
            Values[5] = -gup_rr(i, j, k) * ddr_fact + (Vr(i, j, k) - Gam_r) * dr_fact;
            Indices[5] = II_ind(i - 2, j, k);
            // i - 1:
            Values[6] = 16.0 * gup_rr(i, j, k) * ddr_fact - 8.0 * (Vr(i, j, k) - Gam_r) * dr_fact;
            Indices[6] = II_ind(i - 1, j, k);
            // j + 1
            Values[7] = 16.0 * gup_tt(i, j, k) * ddtheta_fact + 8.0 * (Vt(i, j, k) - Gam_t) * dtheta_fact;
            Indices[7] = II_ind(i, j + 1, k);
            // j + 2
            Values[8] = -gup_tt(i, j, k) * ddtheta_fact - (Vt(i, j, k) - Gam_t) * dtheta_fact;
            Indices[8] = II_ind(i, j + 2, k);
            // k + 1:
            Values[9] = 16.0 * gup_pp(i, j, k) * ddphi_fact + 8.0 * (Vp(i, j, k) - Gam_p) * dphi_fact;
            Indices[9] = II_ind(i, j, k + 1);
            // k + 2:
            Values[10] = -gup_pp(i, j, k) * ddphi_fact - (Vp(i, j, k) - Gam_p) * dphi_fact;
            Indices[10] = II_ind(i, j, k + 2);
            //
            // now add mixed derivatives, treated only to second order
            //
            // i-1, j-1, k:
            Values[11] = 2.0 * gup_rt(i, j, k) * drdtheta_fact;
            Indices[11] = II_ind(i - 1, j - 1, k);
            // i-1, j+1, k:
            Values[12] = -2.0 * gup_rt(i, j, k) * drdtheta_fact;
            Indices[12] = II_ind(i - 1, j + 1, k);
            // i-1, j, k-1:
            Values[13] = 2.0 * gup_rp(i, j, k) * drdphi_fact;
            Indices[13] = II_ind(i - 1, j, k - 1);
            // i-1, j, k+1:
            Values[14] = -2.0 * gup_rp(i, j, k) * drdphi_fact;
            Indices[14] = II_ind(i - 1, j, k + 1);
            // i, j-1, k-1:
            Values[15] = 2.0 * gup_tp(i, j, k) * dthetadphi_fact;
            Indices[15] = II_ind(i, j - 1, k - 1);
            // i, j-1, k+1:
            Values[16] = -2.0 * gup_tp(i, j, k) * dthetadphi_fact;
            Indices[16] = II_ind(i, j - 1, k + 1);
            // i, j+1, k-1:
            Values[17] = -2.0 * gup_tp(i, j, k) * dthetadphi_fact;
            Indices[17] = II_ind(i, j + 1, k - 1);
            // i, j+1, k+1:
            Values[18] = 2.0 * gup_tp(i, j, k) * dthetadphi_fact;
            Indices[18] = II_ind(i, j + 1, k + 1);

            //
            // rest depends on where i is...
            //
            if (i < n_r - N_g - 2) {   // all neighbors are in interior of grid
                //
                // mixed derivatives
                //
                // i+1, j-1, k:
                Values[19] = -2.0 * gup_rt(i, j, k) * drdtheta_fact;
                Indices[19] = II_ind(i + 1, j - 1, k);
                // i+1, j+1, k:
                Values[20] = 2.0 * gup_rt(i, j, k) * drdtheta_fact;
                Indices[20] = II_ind(i + 1, j + 1, k);
                // i+1, j, k-1:
                Values[21] = -2.0 * gup_rp(i, j, k) * drdphi_fact;
                Indices[21] = II_ind(i + 1, j, k - 1);
                // i+1, j, k+1:
                Values[22] = 2.0 * gup_rp(i, j, k) * drdphi_fact;
                Indices[22] = II_ind(i + 1, j, k + 1);
                //
                // i + 1:
                Values[23] = 16.0 * gup_rr(i, j, k) * ddr_fact + 8.0 * (Vr(i, j, k) - Gam_r) * dr_fact;
                Indices[23] = II_ind(i + 1, j, k);
                // i + 2:
                Values[24] = -gup_rr(i, j, k) * ddr_fact - (Vr(i, j, k) - Gam_r) * dr_fact;
                Indices[24] = II_ind(i + 2, j, k);
                //
                // diagonal element
                //
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 25, &Values[0], &Indices[0]);
                assert(ierr == 0);
            } else if (i == n_r - N_g - 2) {  // point i+2 is in ghost zone
                //
                // mixed derivatives
                //
                // i+1, j-1, k:
                Values[19] = -2.0 * gup_rt(i, j, k) * drdtheta_fact;
                Indices[19] = II_ind(i + 1, j - 1, k);
                // i+1, j+1, k:
                Values[20] = 2.0 * gup_rt(i, j, k) * drdtheta_fact;
                Indices[20] = II_ind(i + 1, j + 1, k);
                // i+1, j, k-1:
                Values[21] = -2.0 * gup_rp(i, j, k) * drdphi_fact;
                Indices[21] = II_ind(i + 1, j, k - 1);
                // i+1, j, k+1:
                Values[22] = 2.0 * gup_rp(i, j, k) * drdphi_fact;
                Indices[22] = II_ind(i + 1, j, k + 1);
                // i + 1:
                Values[23] = 16.0 * gup_rr(i, j, k) * ddr_fact + 8.0 * (Vr(i, j, k) - Gam_r) * dr_fact;
                Indices[23] = II_ind(i + 1, j, k);
                // i + 2:
                double ip2_term = -gup_rr(i, j, k) * ddr_fact - (Vr(i, j, k) - Gam_r) * dr_fact;
                Values[23] += fall_off_factor((*r)[i + 1], (*r)[i + 2]) * ip2_term;
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 24, &Values[0], &Indices[0]);
                assert(ierr == 0);

            } else if (i == n_r - N_g - 1) { // both upper neighbors in ghost zone
                const double ip1_term = 16.0 * gup_rr(i, j, k) * ddr_fact + 8.0 * (Vr(i, j, k) - Gam_r) * dr_fact;
                const double ip2_term = -gup_rr(i, j, k) * ddr_fact - (Vr(i, j, k) - Gam_r) * dr_fact;
                // i+1, j, k and i+2, j, k correct diagonal term:
                diag += fall_off_factor((*r)[i], (*r)[i + 2]) * ip2_term + fall_off_factor((*r)[i], (*r)[i + 1]) * ip1_term;
                const double jmoterm = -2.0 * gup_rt(i, j, k) * drdtheta_fact;
                const double jpoterm = 2.0 * gup_rt(i, j, k) * drdtheta_fact;
                const double kmoterm = -2.0 * gup_rp(i, j, k) * drdphi_fact;
                const double kpoterm = 2.0 * gup_rp(i, j, k) * drdphi_fact;
                const double factor = fall_off_factor((*r)[i], (*r)[i + 1]);
                // i+1, j-1, k corrects i, j-1, k term [4]
                Values[4] += factor * jmoterm;
                // i+1, j+1, k corrects i, j+1, k term [7]
                Values[7] += factor * jpoterm;
                // i+1, j, k-1 corrects i, j, k-1 term [2]
                Values[2] += factor * kmoterm;
                // i+1, j, k+1 corrects i, j, k+1 term [9]
                Values[9] += factor * kpoterm;
                // diagonal element
                Values[0] = diag;
                Indices[0] = II;
                // 
                // now store matrix elements
                //
                ierr = A->InsertGlobalValues(MyGlobalElements[ii], 19, &Values[0], &Indices[0]);
                assert(ierr == 0);
            } else {
                cout << " ELLSOLVER3D: Ooops - should never have reached i = " << i << endl;
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
        // first use periodicity for i indices in inner ghost zones...
        if ((i == n_r - 1) || (i == n_r - 2))
            cerr << " ELLSOLVER3D: II_ind should never be called in fall-off region! " << endl;
        if (i == 0) {
            i = 3;
            j = n_theta - j - 1;                       // theta -> pi - theta
            k = 2 + (k + n_phi / 2 - 4) % (n_phi - 4);   // phi -> phi + pi
        }
        if (i == 1) {
            i = 2;
            j = n_theta - j - 1;                       // theta -> pi - theta
            k = 2 + (k + n_phi / 2 - 4) % (n_phi - 4);   // phi -> phi + pi
        }
        // then use periodicity for j indices in ghost zones...
        if (j == 0) {
            j = 3;
            k = 2 + (k + n_phi / 2 - 4) % (n_phi - 4);   // phi -> phi + pi
        }
        if (j == 1) {
            j = 2;
            k = 2 + (k + n_phi / 2 - 4) % (n_phi - 4);   // phi -> phi + pi
        }
        if (j == n_theta - 2) {
            j = n_theta - 3;
            k = 2 + (k + n_phi / 2 - 4) % (n_phi - 4);   // phi -> phi + pi
        }
        if (j == n_theta - 1) {
            j = n_theta - 4;
            k = 2 + (k + n_phi / 2 - 4) % (n_phi - 4);   // phi -> phi + pi
        }
        // ... and finally use periodicity for k indices in ghost zones...
        if (k == 0) k = n_phi - 4;
        if (k == 1) k = n_phi - 3;
        if (k == n_phi - 1) k = 3;
        if (k == n_phi - 2) k = 2;
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

#endif  /* _ELLSOLVER3D_H_ */

#endif  /* NoEllSolver */
