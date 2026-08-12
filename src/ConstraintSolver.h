// Tell emacs that this is -*-c++-*- mode
//================================================
// Class containing routines for solution of 
// Hamiltonian constraint 
//================================================
//
#ifndef CONSTRAINTSOLVER_H
#define CONSTRAINTSOLVER_H
//
// #include "nr3.h"
#include "State.h"
#include "Curvature.h"
#include "Grid.h"
#include "Matter.h"
// #include "FlatEllSolver3D2ndOrder.h"

class ConstraintSolver {
private:
    Grid* grid;
    state* s;
    curvature* curve;
    int N_g, N_r, N_theta, N_phi;
#ifndef NoEllSolver
    FlatEllSolver3D* ellsolver;
#endif
    // Auxiliary functions:
    gf3d u, v, psi, A2, D2_psi, HamRes, delta_psi;
    double PI;
public:
    //
    // Constructor
    //
    ConstraintSolver(Grid* grid_i, state* s_i, curvature* curve_i,
        Matter* matter_i);
    //
    // Destructor
    //
    ~ConstraintSolver();
    //
    // interface for solving Hamiltonian constraint
    //
public:
    double SolveHamiltonian(double tol = 1.e-10, int itmax = 50,
        bool precollapsed_lapse = true);
    //
    // Compute Residual
    // 
private:
    double HamiltonianResidual();
    //
    // Compute auxiliary functions u and v
    // 
    int Compute_uv();
};


#endif /* CONSTRAINTSOLVER_H */
