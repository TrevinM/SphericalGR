// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing flat, 2nd order 3D solver for vector Laplacian, 
// but assuming that phi-components are ONLY NONZERO components!
// In that case, the vector laplacian, applied to the geometrically 
// rescaled component w^phi = r sin(theta) W^phi, reduces to
//
// r^{-2} \partial_r ( r^2 \partial_r w^\phi ) 
//  + r{-2} \sin^{-1} \theta \partial_\theta ( \sin \partial \theta w^\phi )
//  - r{-2} \sin^{-2} w^phi 
//
// Can be solved in same way as Laplacian - we just modify powers of r and
// \sin \theta in the functions \Gamma^i below.
//
//================================================

#ifndef NoEllSolver

//================================================
//
// Uses Trilinos software to solve covariant 3D Laplace
// operator in spherical polar coordinates
//
//  \nabla_L^2 W^phi + factor * u * W^phi = RHS^phi
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
// maximum value, (n_r - N_g)*(n_theta - N_g)*(n_phi - N_g) - 1.
//
// N_g is hard-coded to three.
//
//================================================

#ifndef _VECLAPLACE_H_
#define _VECLAPLACE_H_


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

class VecLaplace {
private:
  int n_r;
  int n_theta;
  int n_phi;
  int n_global;
  double c;    // factor for logarithmic grid
  int N_g;      // number of ghost zones - hardcoded to two
  int NumMyElements;
  VecDoub *r, *Delta_r, *theta, *phi;
  int fall_off;   // default set to 1 - change with SetFallOff()
#ifdef EPETRA_MPI
  Epetra_MpiComm *Comm;
#else
  Epetra_SerialComm *Comm;
#endif
  Epetra_Map *Map;
  Epetra_Vector *rhs, *sol;
  Epetra_CrsMatrix *A;
  bool verbose;
  double Ep1, E0, Em1;                // coefficients for centered first derivative (second order)
  double Fp1, F0, Fm1;                // coefficients for centered second derivative (first order)
public:
  //================================================
  // Constructor
  //================================================
  VecLaplace(int n_r_i, int n_theta_i, int n_phi_i, double c_i, 
	     VecDoub * r_i, VecDoub * Delta_r_i, 
	     VecDoub * theta_i, VecDoub * phi_i, 
	     bool verbose_i = false) :
    n_r(n_r_i), n_theta(n_theta_i), n_phi(n_phi_i), c(c_i), N_g(3), 
    r(r_i), Delta_r(Delta_r_i), theta(theta_i), phi(phi_i),
    fall_off(1), verbose(verbose_i)
  {    

    cout << " CCCCAAAUUUTTTIIIOOONNN - dont use VecLaplace - currently not functional... " << endl;
    cout << " (need to implement new differencing, N_g = 3...) " << endl;
    if (verbose) cout << " Constructing VecLaplace";
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
    const int NumGlobalElements = (n_r - 2*N_g) * (n_theta - 2*N_g) * (n_phi - 2*N_g);
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
    // 
    // Finally: set up coefficients for logarithmic differencing in r-direction
    // 
    const Doub opc = 1.0 + c;
    // const Doub c2 = c*c;
    // Centered first derivative (second order) FIXED
    Ep1 = 1.0 / (c * opc);
    E0  = (c - 1.0) / c;
    Em1 = - c / opc;
    // Centered second derivative (first order) FIXED
    Fp1 = 2.0 / (c * opc);
    F0 = - 2.0 / c;
    Fm1 = 2.0 / opc;
    if (verbose) cout << " Ep1 = " << Ep1 << " E0 = " << E0 << " Em1 = " << Em1 << endl;
  };
  //================================================
  // Destructor
  //================================================
  ~VecLaplace() {
    if (verbose) cout << " Destructing VecLaplace... " << endl;
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
  int SetupSolver(double factor, gf3d & u) {
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
    for (int i = 0; i<NumMyElements; i++)
      NumNz[i] = 7;
    //
    // create matrix
    // 
    if ( A != NULL ) delete A;
    A = new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]);
    ierr = SetupCovLaplace(factor, u);
    ierr = A->FillComplete();
    //    cout << *A << endl;
    assert(ierr==0);
    return ierr;
  }
  //================================================
  // Set up Solver (max_it on input is maximum number of iterations, num_it
  //                number of iterations used; returns residual)
  //================================================
  double Solve(int max_it, int & num_it, double tol) {
    Epetra_LinearProblem linprob(A,sol,rhs);
    AztecOO Solver(linprob);
    ofstream aztec_output;
    aztec_output.open("output/Aztec_Output_CovEllSolver");
    Solver.SetOutputStream(aztec_output);
    //    Solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
    //    Solver.SetAztecOption(AZ_precond, AZ_ilut);
    //    Solver.SetAztecOption(AZ_output, AZ_summary);
    Solver.SetAztecOption(AZ_conv, AZ_rhs);
    Solver.Iterate(max_it,tol);    
    num_it = Solver.NumIters();
    aztec_output.close();
    return Solver.TrueResidual();
  }
  //================================================
  // Set right hand side
  //================================================
  int SetRHS(gf3d & rhs_grid) { return SetRHS(1.0,rhs_grid); };
  int SetRHS(double factor, gf3d & rhs_grid) {
    //    cout << " Setting up RHS..." << endl;
    for (int II = 0; II < n_global; II++) { // loop over superindex
      (*rhs)[II] = factor * rhs_grid(i_ind(II),j_ind(II),k_ind(II));
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
  int SetInitialGuess(gf3d & guess) {
    //    cout << " Setting up initial guess..." << endl;
    for (int II = 0; II < n_global; II++) { // loop over superindex
      (*sol)[II] = guess(i_ind(II),j_ind(II),k_ind(II));
    }
    return 1;
  };
  //================================================
  // Get solution
  //================================================
  void GetSolution(gf3d & sol_grid) {
    for (int i = N_g; i < n_r-N_g; i++) // loop over grid indices
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) 
	  sol_grid[i][j][k] =  (*sol)[II_ind(i,j,k)];
    //
    // now fill outer ghostzones with correct falloff
    //
    int i_in = n_r - 3;
    for (int j = N_g; j < n_theta-N_g; j++) 
      for (int k = N_g; k < n_phi-N_g; k++) {
	int i = n_r - 3;
	sol_grid[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid(i_in,j,k);
	i = n_r - 2;
	sol_grid[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid(i_in,j,k);
	i = n_r - 1;
	sol_grid[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid(i_in,j,k);
      }
    sol_grid.fill_ghosts();
  };
private:
  //================================================
  // Set up covariant Laplace operator
  //================================================
  int SetupCovLaplace(double factor, gf3d & u) {
    int ierr = 0;
    std::vector<int> MyGlobalElements(NumMyElements);
    Map->MyGlobalElements(&MyGlobalElements[0]);
    // assume equidistant grid in theta and phi directions, but 
    // logarithmic grid in r direction
    const double d_phi = (*phi)[1] - (*phi)[0];
    const double d_theta = (*theta)[1] - (*theta)[0];
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
      // define factors for derivatives in 2nd order differencing
      const double d_r = (*Delta_r)[i];
      const double ddr_fact     = 1.0/(d_r*d_r);
      const double ddtheta_fact = 1.0/(d_theta*d_theta); 
      const double ddphi_fact   = 1.0/(d_phi*d_phi); 
      const double dr_fact      = 1.0/(d_r);
      const double dtheta_fact  = 0.5/(d_theta);
      // const double dphi_fact    = 0.5/(d_phi);
      //      const double drdtheta_fact   = 0.25/(d_r*d_theta);
      //      const double drdphi_fact     = 0.25/(d_r*d_phi);
      //      const double dthetadphi_fact = 0.25/(d_theta*d_phi);
      //
      const double rl = (*r)[i];
      const double r2 = rl*rl;
      const double sintheta = sin((*theta)[j]); 
      const double sin2theta = sintheta*sintheta;
      const double costheta = cos((*theta)[j]);
      // const double cottheta = costheta/sintheta;
      //
      // HACK!!!  These are not really \Gamma^i, as they would appear
      // in the flat Laplace operator - instead they accommodate the 
      // powers of r and \sin \theta appropriate for vector Laplacian 
      // - see comments at top of file
      //
      const double Gam_r = - 2.0/rl;
      const double Gam_t = - costheta/(r2*sintheta);
      //
      // prepare diagonal element
      //
      double diag = F0 * ddr_fact - E0 * Gam_r* dr_fact - 2.0 * ddtheta_fact/r2 - 
	2.0 * ddphi_fact/(r2*sin2theta) 
	- 1.0/(r2*sin2theta)
	+ factor * u(i,j,k);
      //
      // now go through grid, but ignore i+1 terms for now
      //
      // i-1, j, k:
      int par = 1;
      Indices[1] = II_ind(i-1,j,k,par);
      Values[1] = par * ( Fm1 * ddr_fact - Em1*Gam_r*dr_fact ); 
      // i, j-1, k:
      Indices[2] = II_ind(i,j-1,k,par);
      Values[2] = par * ( ddtheta_fact/r2 + Gam_t*dtheta_fact ); 
      // i, j, k-1:
      Indices[3] = II_ind(i,j,k-1,par);
      Values[3] = par * ddphi_fact/(r2*sin2theta);
      // i, j+1, k:
      Indices[4] = II_ind(i,j+1,k,par);
      Values[4] = par * ( ddtheta_fact/r2 - Gam_t*dtheta_fact ); 
      // i, j, k+1:
      Indices[5] = II_ind(i,j,k+1,par);
      Values[5] = par * ddphi_fact/(r2*sin2theta);
      //
      // rest depends on where i is...
      //
      if (i < n_r - N_g - 1) {   // i+1 is in interior of grid...
	// i+1, j, k:
	Indices[6] = II_ind(i+1,j,k,par);
	Values[6] = par * ( Fp1*ddr_fact - Ep1*Gam_r*dr_fact ); 
	// diagonal element:
	Values[0] = diag; 
	Indices[0] = II;
	// 
	// now store matrix elements
	//
	ierr = A->InsertGlobalValues(MyGlobalElements[ii],7,&Values[0],&Indices[0]);
	assert(ierr==0);
      } else if (i == n_r - N_g - 1) { // i+1 is in ghost zone
	const double iterm =  Fp1 * ddr_fact - Ep1*Gam_r * dr_fact; 
	const double factor = fall_off_factor((*r)[i],(*r)[i+1]);
	// i+1, j, k corrects diagonal term:
	diag += factor * iterm;
	// diagonal element:
	Values[0] = diag; 
	Indices[0] = II;
	// 
	// now store matrix elements
	//
	ierr = A->InsertGlobalValues(MyGlobalElements[ii],6,&Values[0],&Indices[0]);
	assert(ierr==0);
      } else {
	cout << " VECLAPLACE: Ooops - should never have reached i = " << i << endl;
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
  // Also: return correct parity for phi-component of vector: -1, -1, +1
  //
  //================================================
  inline int II_ind(int i, int j, int k) { int par; return II_ind(i, j, k, par); }
  inline int II_ind(int i, int j, int k, int & par) { 
    int center_par = -1;
    int axis_par = -1;
    int eq_par = 1;
    par = 1;
    // first use periodicity for i indices in inner ghost zones...
    if ( (i == n_r-1) || (i == n_r-2) )
      cerr << " VECLAPLACE: II_ind should never be called in fall-off region! " << endl;
    if (i == 0) {
      cerr << " VECLAPLACE: why is II_ind called for i == 0??? " << endl;
      par *= center_par;
      i = 3;
#ifndef EQSYMMETRY
      j = n_theta - j - 1;                       // theta -> pi - theta
#else
      par *= eq_par;
#endif
#ifndef AXISYMMETRY
      k = 2 + (k + n_phi/2 - 4) % (n_phi - 4);   // phi -> phi + pi
#endif
    }
    if (i == 1) {
      par *= center_par;
      i = 2;
#ifndef EQSYMMETRY
      j = n_theta - j - 1;                       // theta -> pi - theta
#else
      par *= eq_par;
#endif
#ifndef AXISYMMETRY
      k = 2 + (k + n_phi/2 - 4) % (n_phi - 4);   // phi -> phi + pi
#endif
    }
    // then use periodicity for j indices in ghost zones...
    if (j == 0) {
      par *= axis_par;
      j = 3;
#ifndef AXISYMMETRY
      k = 2 + (k + n_phi/2 - 4) % (n_phi - 4);   // phi -> phi + pi
#endif
    }
    if (j == 1) {
      par *= axis_par;
      j = 2;
#ifndef AXISYMMETRY
      k = 2 + (k + n_phi/2 - 4) % (n_phi - 4);   // phi -> phi + pi
#endif
    }
    if (j == n_theta-2) {
#ifndef EQSYMMETRY   
      par *= axis_par;
#else
      par *= eq_par;
#endif
      j = n_theta-3;
#ifndef AXISYMMETRY
      k = 2 + (k + n_phi/2 - 4) % (n_phi - 4);   // phi -> phi + pi
#endif
    }
    if (j == n_theta-1) {
#ifndef EQSYMMETRY   
      par *= axis_par;
#else
      par *= eq_par;
#endif
      j = n_theta-4;
#ifndef AXISYMMETRY
      k = 2 + (k + n_phi/2 - 4) % (n_phi - 4);   // phi -> phi + pi
#endif
    }
    // ... and finally use periodicity for k indices in ghost zones...
#ifndef AXISYMMETRY
    if (k == 0) k = n_phi-4;
    if (k == 1) k = n_phi-3;
    if (k == n_phi-1) k = 3;
    if (k == n_phi-2) k = 2;
#else
    k = 2;
#endif
    // now return corresponding superindex
    return (i - N_g) + (j - N_g)*(n_r - 2*N_g) + (k - N_g)*(n_r - 2*N_g)*(n_theta - 2*N_g); 
  };
  inline int i_ind(int II) { return ( II % ((n_r - 2*N_g)*(n_theta - 2*N_g)) ) % (n_r - 2*N_g) + N_g; };
  inline int j_ind(int II) { return ( II % ((n_r - 2*N_g)*(n_theta - 2*N_g)) ) / (n_r - 2*N_g) + N_g; };
  inline int k_ind(int II) { return II / ((n_r - 2*N_g)*(n_theta - 2*N_g)) + N_g; };
  //===============================================================
  // provide fall_off factor
  //===============================================================
  double fall_off_factor(double r_in, double r_out) {
    const double ratio = r_in/r_out;
    if (fall_off == 1) return ratio;
    else if (fall_off == 2) return ratio*ratio;
    else return pow(ratio,fall_off);
  };
public:
  //===============================================================
  //
  // As a test, solve the problem
  //
  //  \nabla^2 f + \lambda f = (A + \lambda) Y^{11}
  //
  //===============================================================
  
  void Test(double AA, gf3d & rhs_grid, gf3d & sol_grid, gf3d & res) {
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
  
#endif  /* _VECLAPLACE_H_ */

#endif  /* NoEllSolver */
