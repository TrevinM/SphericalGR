// Tell emacs that this is -*-c++-*- mode
//================================================
//
// Class containing 3D vector Laplace solver
// for all three vector components, assuming that vector components
// are independent of phi
//
// (Not to be confused with the simpler version in VecLaplaceSolver.h,
// which solves for phi component only!)
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

#ifndef _VECLAPLACIAN_H_
#define _VECLAPLACIAN_H_


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

enum {R, THETA, PHI};

using namespace std;

class VecLaplacian {
private:
  Grid *grid;
  int fall_off;   // default set to 1 - change with SetFallOff()
  int n_r;
  int n_theta;
  int n_phi;
  int n_global, n_vec; 
  int N_g;      // number of ghost zones - hardcoded to three
  int NumMyElements;
  VecDoub *r;  // , *theta, *phi;
#ifdef EPETRA_MPI
  Epetra_MpiComm *Comm;
#else
  Epetra_SerialComm *Comm;
#endif
  Epetra_Map *Map;
  Epetra_Vector *rhs, *sol;
  Epetra_CrsMatrix *A;
  bool verbose;
public:
  //================================================
  // Constructor
  //================================================
  VecLaplacian(Grid * grid_i, bool verbose_i = false) :
    grid(grid_i), fall_off(1), verbose(verbose_i)
  {    
    N_g = grid->N_ghosts();
    n_r = grid->N_r_tot();
    n_theta = grid->N_theta_tot();
    n_phi = grid->N_phi_tot();
    r = grid->r();
    if (verbose) cout << " VECLAPLACIAN: Constructing VecLaplacian";
    // If Trilinos was built with MPI, initialize MPI, otherwise
    // initialize the serial "communicator" that stands in for MPI.
#ifdef EPETRA_MPI
    Comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
    Comm = new Epetra_SerialComm();
#endif
    //    cout << *Comm << endl;
    //
    // Create Epetra Map - note factor of three for three components...
    //
    n_vec = (n_r - 2*N_g) * (n_theta - 2*N_g) * (n_phi - 2*N_g);
    const int NumGlobalElements = 3 * n_vec;
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
  ~VecLaplacian() {
    if (verbose) cout << " VECLAPLACIAN: Destructing VecLaplacian... " << endl;
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
  int SetupSolver() {
    int ierr = 0;
    // even though rest of code is serial, we will allow for
    // parallel implementation here
    NumMyElements = Map->NumMyElements();  // # of elements on local proc.
    // create vector that contains (global) indices of local elements
    //
    std::vector<int> MyGlobalElements(NumMyElements);
    Map->MyGlobalElements(&MyGlobalElements[0]);
    //
    // in our application, need
    //     8 off-diagonal terms for both r and theta component
    //     4  off-diagonal terms for phi component
    //     3  diagonal terms (one each...)
    //
    std::vector<int> NumNz(NumMyElements);
    for (int i = 0; i<NumMyElements; i++)
      NumNz[i] = 23;
    //
    // create matrix
    // 
    if ( A != NULL ) delete A;
    A = new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]);
    ierr = SetupCovLaplace();
    ierr = A->FillComplete();
    //    cout << *A << endl;
    assert(ierr==0);
    //
    // Sanity check: make sure indices are handled correctly...
    //
    for (int ii = 0; ii < NumMyElements; ii++) { // loop over superindex
      int II = MyGlobalElements[ii];               // global superindex
      //
      // check index handling...
      int i = i_ind(II);
      int j = j_ind(II);
      int k = k_ind(II);
      int dim = dim_ind(II);
      double parity;
      int II_return = II_ind(i,j,k,dim,parity);
      if (II != II_return) {
	cout << " VECLAPLACIAN: found II = " << II
	     << " but II_return = " << II_return << endl;
	exit(0);
      }
    }      
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
    aztec_output.open("output/Aztec_Output_VecLaplacian");
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
  int SetRHS(gf3d & rhs_grid_r, gf3d & rhs_grid_t, gf3d & rhs_grid_p) {
    return SetRHS(1.0, rhs_grid_r, rhs_grid_t, rhs_grid_p);
  };
  int SetRHS(double factor,
	     gf3d & rhs_grid_r, gf3d & rhs_grid_t, gf3d & rhs_grid_p) {
    //    cout << " Setting up RHS..." << endl;
    for (int II = 0; II < n_global; II++) { // loop over superindex
      if (II < n_vec) {
	(*rhs)[II] = factor * rhs_grid_r(i_ind(II),j_ind(II),k_ind(II));
      } else if (II < 2*n_vec) {
	(*rhs)[II] = factor * rhs_grid_t(i_ind(II),j_ind(II),k_ind(II));
      } else {
	(*rhs)[II] = factor * rhs_grid_p(i_ind(II),j_ind(II),k_ind(II));
      }
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
  int SetInitialGuess(gf3d & guess_r, gf3d & guess_t, gf3d & guess_p) {
    //    cout << " Setting up initial guess..." << endl;
    for (int II = 0; II < n_global; II++) { // loop over superindex
      if (II < n_vec) {
	(*rhs)[II] = guess_r(i_ind(II),j_ind(II),k_ind(II));
      } else if (II < 2*n_vec) {
	(*rhs)[II] = guess_t(i_ind(II),j_ind(II),k_ind(II));
      } else {
	(*rhs)[II] = guess_p(i_ind(II),j_ind(II),k_ind(II));
      }
    }
    return 1;
  };
  //================================================
  // Get solution
  //================================================
  void GetSolution(gf3d & sol_grid_r, gf3d & sol_grid_t, gf3d & sol_grid_p) {
    for (int i = N_g; i < n_r-N_g; i++) // loop over grid indices
      for (int j = N_g; j < n_theta-N_g; j++) 
	for (int k = N_g; k < n_phi-N_g; k++) {
	  sol_grid_r[i][j][k] =  (*sol)[II_ind(i,j,k,R)];
	  sol_grid_t[i][j][k] =  (*sol)[II_ind(i,j,k,THETA)];
	  sol_grid_p[i][j][k] =  (*sol)[II_ind(i,j,k,PHI)];
	}
    //
    // now fill outer ghostzones with correct falloff
    //
    //    int i_in = n_r - 4;
    int i_in = n_r - N_g - 1;
    for (int j = N_g; j < n_theta-N_g; j++) 
      for (int k = N_g; k < n_phi-N_g; k++)
	for (int i = n_r - N_g; i < n_r; i++) {
	  sol_grid_r[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid_r(i_in,j,k);
	  sol_grid_t[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid_t(i_in,j,k);
	  sol_grid_p[i][j][k] = fall_off_factor((*r)[i_in],(*r)[i]) * sol_grid_p(i_in,j,k);
	}
    sol_grid_r.fill_ghosts();
    sol_grid_t.fill_ghosts();
    sol_grid_p.fill_ghosts();
  };
private:
  //================================================
  // Set up covariant Laplace operator
  //================================================
  int SetupCovLaplace() {
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
    std::vector<int> Indices(23);
    std::vector<double> Values(23);
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
      const double r2 = rl*rl;
      const double sintheta = grid->sintheta(j);
      const double sin2theta = sintheta*sintheta;
      const double costheta = grid->costheta(j);
      const double cottheta = costheta/sintheta;
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
      const double Gam_r = - 2.0/rl;
      const double Gam_t = - cottheta/r2;
      //
      // define factors for derivatives in 2nd order differencing
      //
      const double ddx_term     = 1.0/(d_x*d_x)*dxdr*dxdr;
      const double ddy_term     = 1.0/(d_y*d_y)*dydtheta*dydtheta/r2; 
      const double dx2_term     = 1.0/(2.0*d_x)*ddxdr;
      const double dy2_term     = 1.0/(2.0*d_y)*ddydtheta/r2;
      const double dx_term      = 1.0/(2.0*d_x)*dxdr/rl;
      const double dy_term      = 1.0/(2.0*d_y)*dydtheta/r2; 
      const double dxdy_term    = 1.0/(4.0*d_x*d_y)*dxdr*dydtheta/rl;
      //
      // now go through grid, but ignore i+1 terms for now
      //
      // R components:
      if (dim_ind(II) == R) {
	//
	// prepare diagonal element
	//
	double parity = 1.0;
	// i, j k, R: diagonal element - prepare only...
	double diag = - 2.0 * ((4./3.)*ddx_term + ddy_term) - (8./3.)/r2;
	// i-1, j, k, R: d^2_r b^r, d_r b^r 
	Indices[1] = II_ind(i-1, j, k, R, parity);
	Values[1] = parity * ( (4./3.)*(ddx_term - dx2_term) - (8./3.)*dx_term );
	// i, j-1, k, R: d^2_t b^r, d_t b^r
	Indices[2] = II_ind(i, j-1, k, R, parity);
	Values[2] = parity * ( ddy_term - dy2_term - cottheta*dy_term );
	// i, j+1, k, R: d^2_t b^r, d_t b^r
	Indices[3] = II_ind(i, j+1, k, R, parity);
	Values[3] = parity * ( ddy_term + dy2_term + cottheta*dy_term );
	// i, j, k, THETA: b^t prepare only
	double IJ = - (7./3.)*cottheta/r2;
	// i-1, j, k, THETA: d_r b^t
	Indices[4] = II_ind(i-1, j, k, THETA, parity);
	Values[4] = parity * ( - cottheta/3.*dx_term );
	// i, j-1, k, THETA: d_t b^t prepare only
	double IJm1 = 7./3.*dy_term;
	// i, j+1, k, THETA: d_t b^t prepare only
	double IJp1 = - 7./3.*dy_term;
	// i-1, j-1, k, THETA: d_r d_t b^t
	Indices[5] = II_ind(i-1, j-1, k, THETA, parity);
	Values[5] = parity * ( dxdy_term / 3.0 );
	// i-1, j+1, k, THETA: d_r d_t b^t
	Indices[6] = II_ind(i-1, j+1, k, THETA, parity);
	Values[6] = parity * ( - dxdy_term / 3.0 );
       
	//
	// rest depends on where i is...
	//
	if (i < n_r - N_g - 1) {   // i+1 is in interior of grid...
	  // i+1, j, k, R: d^2_r b^r, d_r b^r  
	  Indices[7] = II_ind(i+1, j, k, R, parity);
	  Values[7] = parity * ( (4./3.)*(ddx_term + dx2_term) + (8./3.)*dx_term );
	  // i, j, k, THETA: b^t 
	  Indices[8] = II_ind(i, j, k, THETA, parity);
	  Values[8] = parity * IJ;
	  // i+1, j, k, THETA: d_r b^t
	  Indices[9] = II_ind(i+1, j, k, THETA, parity);
	  Values[9] = parity * cottheta/3.*dx_term;
	  // i+1, j-1, k, THETA: d_r d_t b^t
	  Indices[10] = II_ind(i+1, j-1, k, THETA, parity);
	  Values[10] = parity * ( - dxdy_term / 3.0);
	  // i+1, j+1, k, THETA: d_r d_t b^t
	  Indices[11] = II_ind(i+1, j+1, k, THETA, parity);
	  Values[11] = parity * ( dxdy_term / 3.0);
	  // i, j-1, k, THETA: d_t b^t 
	  Indices[12] = II_ind(i, j-1, k, THETA, parity);
	  Values[12] = parity * IJm1;
	  // i, j+1, k, THETA: d_t b^t 
	  Indices[13] = II_ind(i, j+1, k, THETA, parity);
	  Values[13] =  parity * IJp1;
	  // diagonal element:
	  Indices[0] = II;	  
	  Values[0] = diag; 
	  // 
	  // now store matrix elements
	  //
	  ierr = A->InsertGlobalValues(MyGlobalElements[ii],14,&Values[0],&Indices[0]);
	  assert(ierr==0);
	} else if (i == n_r - N_g - 1) { // i+1 is in ghost zone
	  // i+1, j, k, R: corrects diagonal term:
	  double iterm = (4./3.)*(ddx_term + dx2_term) + (8./3.)*dx_term; 
	  const double factor = fall_off_factor((*r)[i],(*r)[i+1]);
	  diag += factor * iterm;
	  Indices[0] = II;	  
	  Values[0] = diag;
	  // i+1, j, k, THETA: corrects IJ term
	  iterm = cottheta/3.*dx_term;
	  IJ += factor * iterm;
	  Indices[7] = II_ind(i, j, k, THETA, parity);
	  Values[7] = parity * IJ; 
	  // i+1, j-1, k, THETA: corrects IJm1 term
	  iterm = - dxdy_term / 3.0;
	  IJm1 += factor * iterm;
	  Indices[8] = II_ind(i, j-1, k, THETA, parity);
	  Values[8] = parity * IJm1;
	  // i+1, j+1, k, THETA: corrects IJp1 term
	  iterm = dxdy_term / 3.0;
	  IJp1 += factor * iterm;
	  Indices[9] = II_ind(i, j+1, k, THETA, parity);
	  Values[9] = parity * IJp1;
	  // 
	  // now store matrix elements
	  //
	  ierr = A->InsertGlobalValues(MyGlobalElements[ii],10,&Values[0],&Indices[0]);
	  assert(ierr==0);
	} else {
	  cout << " VECLAPLACIAN: Ooops - should never have reached i = " << i << endl;
	}
      } else if (dim_ind(II) == THETA) {
	//
	// prepare diagonal element
	//
	double parity = 1.0;
	// i, j k, THETA: diagonal element - prepare only...
	double diag = - 2.0 * (ddx_term + (4./3.)*ddy_term) - (4./3.)/(r2*sin2theta);
	// i-1, j, k, THETA: d^2_r b^t, d_r b^t 
	Indices[1] = II_ind(i-1, j, k, THETA, parity);
	Values[1] = parity * ( ddx_term - dx2_term - 2.0*dx_term );
	// i, j-1, k, THETA: d^2_t b^t, d_t b^t
	Indices[2] = II_ind(i, j-1, k, THETA, parity);
	Values[2] = parity * (4./3.)*( ddy_term - dy2_term - cottheta*dy_term );
	// i, j+1, k, THETA: d^2_t b^t, d_t b^t
	Indices[3] = II_ind(i, j+1, k, THETA, parity);
	Values[3] = parity * (4./3.)*( ddy_term + dy2_term + cottheta*dy_term );
	// i, j-1, k, R: d_t b^r prepare only
	double IJm1 = - 8./3.*dy_term;
	// i, j+1, k, R: d_t b^r prepare only
	double IJp1 =   8./3.*dy_term;
	// i-1, j-1, k, R: d_r d_t b^r
	Indices[4] = II_ind(i-1, j-1, k, R, parity);
	Values[4] = parity * ( dxdy_term / 3.0 );
	// i-1, j+1, k, R: d_r d_t b^r
	Indices[5] = II_ind(i-1, j+1, k, R, parity);
	Values[5] = parity * ( - dxdy_term / 3.0 );
	//
	// rest depends on where i is...
	//
	if (i < n_r - N_g - 1) {   // i+1 is in interior of grid...
	  // i+1, j, k, THETA: d^2_r b^t, d_r b^t  
	  Indices[6] = II_ind(i+1, j, k, THETA, parity);
	  Values[6] = parity * ( ddx_term + dx2_term + 2.0*dx_term);
	  // i+1, j-1, k, R: d_r d_t b^r
	  Indices[7] = II_ind(i+1, j-1, k, R, parity);
	  Values[7] = parity * ( - dxdy_term / 3.0);
	  // i+1, j+1, k, R: d_r d_t b^r
	  Indices[8] = II_ind(i+1, j+1, k, R, parity);
	  Values[8] = parity * ( dxdy_term / 3.0);
	  // i, j-1, k, R: d_t b^r 
	  Indices[9] = II_ind(i, j-1, k, R, parity);
	  Values[9] = parity * IJm1;
	  // i, j+1, k, R: d_t b^r 
	  Indices[10] = II_ind(i, j+1, k, R, parity);
	  Values[10] = parity * IJp1;
	  // diagonal element:
	  Indices[0] = II;	  
	  Values[0] = diag; 
	  // 
	  // now store matrix elements
	  //
	  // Values[0] = 1.0;
	  // for (int ind = 1; ind < 11; ind++) Values[ind] = 0.0;
	  ierr = A->InsertGlobalValues(MyGlobalElements[ii],11,&Values[0],&Indices[0]);
	  assert(ierr==0);
	} else if (i == n_r - N_g - 1) { // i+1 is in ghost zone
	  // i+1, j, k, THETA: corrects diagonal term:
	  double iterm = ddx_term + dx2_term + 2.0*dx_term; 
	  const double factor = fall_off_factor((*r)[i],(*r)[i+1]);
	  diag += factor * iterm;
	  Values[0] = diag;
	  Indices[0] = II;
	  // i+1, j-1, k, R: corrects IJm1 term
	  iterm = - dxdy_term / 3.0;
	  IJm1 += factor * iterm;
	  Indices[6] = II_ind(i, j-1, k, R, parity);
	  Values[6] = parity * IJm1;
	  // i+1, j+1, k, R: corrects IJp1 term
	  iterm = dxdy_term / 3.0;
	  IJp1 += factor * iterm;
	  Indices[7] = II_ind(i, j+1, k, R, parity);
	  Values[7] = parity * IJp1;
	  // 
	  // now store matrix elements
	  //
	  // Values[6] = 0;
	  // Values[7] = 0;
	  // Values[0] = 1.0;
	  // for (int ind = 1; ind < 10; ind++) Values[ind] = 0.0;
	  ierr = A->InsertGlobalValues(MyGlobalElements[ii],8,&Values[0],&Indices[0]);
	  assert(ierr==0);
	} else {
	  cout << " VECLAPLACIAN: Ooops - should never have reached i = " << i << endl;
	}
      } else if (dim_ind(II) == PHI) {
	//
	// prepare diagonal element
	//
	double parity = 1.0;
	// i, j k, R: diagonal element - prepare only...
	double diag = - 2.0 * (ddx_term + ddy_term) - 1.0/(r2*sin2theta);
	// i-1, j, k, PHI: d^2_r b^p, d_r b^p 
	Indices[1] = II_ind(i-1, j, k, PHI, parity);
	Values[1] = parity * ( ddx_term - dx2_term - 2.0*dx_term );
	// i, j-1, k, PHI: d^2_t b^p, d_t b^p
	Indices[2] = II_ind(i, j-1, k, PHI, parity);
	Values[2] = parity * ( ddy_term - dy2_term - cottheta*dy_term );
	// i, j+1, k, PHI: d^2_t b^p, d_t b^p
	Indices[3] = II_ind(i, j+1, k, PHI, parity);
	Values[3] = parity * ( ddy_term + dy2_term + cottheta*dy_term );
	//
	// rest depends on where i is...
	//
	if (i < n_r - N_g - 1) {   // i+1 is in interior of grid...
	  // i+1, j, k, PHI: d^2_r b^p, d_r b^p  
	  Indices[4] = II_ind(i+1, j, k, PHI, parity);
	  Values[4] = parity * ( ddx_term + dx2_term + 2.0*dx_term);
	  // diagonal element:
	  Indices[0] = II;	  
	  Values[0] = diag; 
	  // 
	  // now store matrix elements
	  //
	  ierr = A->InsertGlobalValues(MyGlobalElements[ii],5,&Values[0],&Indices[0]);
	  assert(ierr==0);
	} else if (i == n_r - N_g - 1) { // i+1 is in ghost zone
	  // i+1, j, k, PHI: corrects diagonal term:
	  double iterm = ddx_term + dx2_term + 2.0*dx_term; 
	  const double factor = fall_off_factor((*r)[i],(*r)[i+1]);
	  diag += factor * iterm;
	  Indices[0] = II;	  
	  Values[0] = diag;
	  // 
	  // now store matrix elements
	  //
	  ierr = A->InsertGlobalValues(MyGlobalElements[ii],4,&Values[0],&Indices[0]);
	  assert(ierr==0);
	} else {
	  cout << " VECLAPLACIAN: Ooops - should never have reached i = " << i << endl;
	}
      } else {
	cout << " VECLAPLACIAN: index mess-up!! " << endl;
	exit(0);
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
  inline int II_ind(int i, int j, int k, int dim) {
    double parity;
    return II_ind(i,j,k,dim, parity);
  }
  inline int II_ind(int i, int j, int k, int dim, double & parity) {
    double center_par, ax_par, eq_par;
    parity = 1.0;
    if (dim == R) {
      center_par = -1;
      ax_par = +1;
      eq_par = +1;
    } else if (dim == THETA) {
      center_par = +1;
      ax_par = -1;
      eq_par = -1;
    } else if (dim == PHI) {
      center_par = -1;
      ax_par = -1;
      eq_par = +1;
    } else {
      cout << " VECLAPLACIAN: dim in II_ind not recognized! " << endl;
      exit(1);
    }
#ifdef AXISYMMETRY
    k = N_g;
#endif
    // first use periodicity for i indices in inner ghost zones...
    if ( i >= n_r - N_g )
      cerr << " VECLAPLACIAN: II_ind should never be called in fall-off region! " << endl;
    //
    // deal with i in inner ghost zone
    // 
    if (i < N_g) {
      parity *= center_par;
      i = 2*N_g - 1 - i;
#ifdef EQSYMMETRY
      parity *= eq_par;
#else
      j = n_theta - j - 1;                       // theta -> pi - theta
#endif
#ifndef AXISYMMETRY
      k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
#endif
    }
    //
    // then use periodicity for j indices in ghost zones...
    //
    if (j < N_g) {
      parity *= ax_par;
      j = 2*N_g - 1 - j;
#ifndef AXISYMMETRY
      k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
#endif
    }
    if (j >= n_theta-N_g) {
#ifdef EQSYMMETRY
      parity *= eq_par;
#else
      parity *= ax_par;
#endif /* EQSYMMETRY */ 
      int delta_j = j - (n_theta - N_g);
      j = n_theta - N_g - 1 - delta_j;
#ifndef AXISYMMETRY
      k = (k - N_g + (n_phi - 2*N_g)/2) % (n_phi - 2*N_g) + N_g;  // phi -> phi + pi
#endif
    }
    // ... and finally use periodicity for k indices in ghost zones...
    if (k < N_g) k = n_phi - 2*N_g + k;
    if (k >= n_phi - N_g) k = N_g + (k - (n_phi - N_g));
    // now return corresponding superindex
    return (i - N_g) + (j - N_g)*(n_r - 2*N_g)
      + (k - N_g)*(n_r - 2*N_g)*(n_theta - 2*N_g)
      + dim * n_vec;       // add this last line for r, theta, or phi component... 
  };
  inline int i_ind(int II) {
    const int II_loc = II % n_vec;
    return ( II_loc % ((n_r - 2*N_g)*(n_theta - 2*N_g)) ) % (n_r - 2*N_g) + N_g;
  };
  inline int j_ind(int II) {
    const int II_loc = II % n_vec;
    return ( II_loc % ((n_r - 2*N_g)*(n_theta - 2*N_g)) ) / (n_r - 2*N_g) + N_g;
  };
  inline int k_ind(int II) {
    const int II_loc = II % n_vec;
    return II_loc / ((n_r - 2*N_g)*(n_theta - 2*N_g)) + N_g;
  };
  inline int dim_ind(int II) { return II / n_vec; };
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
  
#endif  /* _VECLAPLACIAN_H_ */

#endif  /* NoEllSolver */
