// Tell emacs that this is -*-c++-*- mode
// 
//================================================
//
// Routines that transform tensors from cartesian coordinates and back
//
//================================================
//
#ifndef TRANS
#define TRANS

#include "tensors.h"

//
//================================================
//
// Rank 1 lower indices: 
//
//================================================
//
vect Cartesian_to_Spherical_lower(vect V, Doub x, Doub y, Doub z) {
    Doub r = sqrt(x * x + y * y + z * z);
    Doub theta = acos(z / r);
    Doub phi = atan(y / x);
    Doub st = sin(theta);
    Doub ct = cos(theta);
    Doub sp = sin(phi);
    Doub cp = cos(phi);
    //
    // compute transformation matrix \partial x^i / partial x^j'
    // 

    tensor Lambda(x / r, y / r, z / r,
        r * ct * cp, r * ct * sp, -r * st,
        -x, y, 0.0);
    vect V_sc;
    for (int i = 0; i < 3; i++) {
        V_sc[i] = 0.0;
        for (int k = 0; k < 3; k++)
            V_sc[i] += Lambda[i][k] * V[k];
    }
    return V_sc;
};

//
//================================================
//
// Rank 1 upper indices: 
//
//================================================
//
vect Cartesian_to_Spherical_upper(vect V, Doub r, Doub theta, Doub phi) {
    const Doub st = sin(theta);
    const Doub ct = cos(theta);
    const Doub sp = sin(phi);
    const Doub cp = cos(phi);

    //
    // compute transformation matrix \partial x^i / partial x^j'
    // 

    tensor Lambda(st * cp, st * sp, ct,
        ct * cp / r, ct * sp / r, -st / r,
        -sp / (r * st), cp / (r * st), 0.0);

    vect V_sc;
    for (int i = 0; i < 3; i++) {
        V_sc[i] = 0.0;
        for (int k = 0; k < 3; k++)
            V_sc[i] += Lambda[i][k] * V[k];
    }
    return V_sc;
};



//
//================================================
//
// Rank 2 lower indices: 
//
//================================================
//
tensor Cartesian_to_Spherical(tensor g, Doub r, Doub theta, Doub phi) {
    Doub st = sin(theta);
    Doub ct = cos(theta);
    Doub sp = sin(phi);
    Doub cp = cos(phi);
    //
    // compute transformation matrix \partial x^i / partial x^j'
    // 

    tensor Lambda(st * cp, st * sp, ct,
        r * ct * cp, r * ct * sp, -r * st,
        -r * st * cp, r * st * sp, 0.0);
    tensor g_sc;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            g_sc[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    g_sc[i][j] += Lambda[i][k] * Lambda[j][l] * g[k][l];
        }
    return g_sc;
};


#endif  /* TRANS */
