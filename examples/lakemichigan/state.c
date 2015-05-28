/*
 * File: state.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains functions that define and implement the equation of state for
 * the density from the salinity, temperature, and pressure.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "state.h"

// Local function
REAL pden(REAL s, REAL t, REAL p, REAL pr);
REAL dens(REAL s, REAL t, REAL p);
REAL dens0(REAL s, REAL t);
REAL smow(REAL t);
REAL ptmp(REAL s, REAL t, REAL p, REAL pr);
REAL adtg( REAL s, REAL t, REAL p);
REAL T90conv(REAL t);
REAL T68conv(REAL t);
REAL seck(REAL s, REAL t, REAL p);

/*
 * Function: StateEquation
 * Usage: rho = StateEquation(prop,s,T,p);
 * ---------------------------------------
 * Returns the density as a function of temperature, salinity, and
 * pressure, where pressure is the hydrostatic pressure p=RHO0*prop->grav*z,
 * and RHO0 and prop->grav are defined in suntans.h.  Note that rho should
 * always normalized by RHO0 so that this function returns a dimensionless
 * quantity.
 *
 */
REAL StateEquation(const propT *prop, const REAL s, const REAL T, const REAL p) {
  REAL sigma;
  // If temperature is a passive scalar
  //  return prop->beta*s;
  // Otherwise temperature is not passive
  // return prop->beta*s - prop->alpha*T;

  sigma=pden(s,T,p*1e-4,0)-RHO0;

  return sigma/RHO0;
}

REAL pden(REAL s, REAL t, REAL p, REAL pr){
    REAL pt;

    pt = ptmp(s, t, p, pr);
    return dens(s, pt, pr);
}
REAL dens(REAL s, REAL t, REAL p){

    REAL densP0, K;

    densP0 = dens0(s, t);
    K = seck(s, t, p);
    p = p / 10.;  // Convert from db to atm pressure units.
    return densP0 / (1 - p / K);
}

REAL dens0(REAL s, REAL t){
    REAL T68;
    REAL b0 = 8.24493e-1;
    REAL b1 = -4.0899e-3; 
    REAL b2 = 7.6438e-5;
    REAL b3 = -8.2467e-7;
    REAL b4 = 5.3875e-9;
    REAL c0 = -5.72466e-3;
    REAL c1 = 1.0227e-4;
    REAL c2 = -1.6546e-6;
    REAL d = 4.8314e-4;


    /*
    # UNESCO 1983 Eqn.(13) p17.
    b = (8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9)
    c = (-5.72466e-3, 1.0227e-4, -1.6546e-6)
    d = 4.8314e-4
    */

    T68 = T68conv(t);
    return (smow(t) + (b0 + (b1 + (b2 + (b3 + b4 * T68) * T68) *
            T68) * T68) * s + (c0 + (c1 + c2 * T68) * T68) * s *
            sqrt(s) + d * pow(s,2));
}

REAL smow(REAL t){
    REAL T68;
    REAL a0 =999.842594;
    REAL a1 = 6.793952e-2;
    REAL a2 = -9.095290e-3;
    REAL a3 = 1.001685e-4;
    REAL a4 = -1.120083e-6;
    REAL a5 =6.536332e-9;

    /*
    a = (999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6,
         6.536332e-9)
    */

    T68 = T68conv(t);
    return (a0 + (a1 + (a2 + (a3 + (a4 + a5 * T68) * T68) * T68) *
            T68) * T68);
}

REAL ptmp(REAL s, REAL t, REAL p, REAL pr){
    REAL del_P, del_th, q, th;
    REAL sqrt2 = sqrt(2);

    // Theta1.
    del_P = pr - p;
    del_th = del_P * adtg(s, t, p);
    th = T68conv(t) + 0.5 * del_th;
    q = del_th;

    // Theta2.
    del_th = del_P * adtg(s, T90conv(th), p + 0.5 * del_P);
    th = th + (1 - 1 / sqrt2 ) * (del_th - q);
    q = (2 - sqrt2 ) * del_th + (-2 + 3 / sqrt2 ) * q;

    // Theta3.
    del_th = del_P * adtg(s, T90conv(th), p + 0.5 * del_P);
    th = th + (1 + 1 / sqrt2) * (del_th - q);
    q = (2 + sqrt2) * del_th + (-2 - 3 / sqrt2) * q;

    // Theta4.
    del_th = del_P * adtg(s, T90conv(th), p + del_P);
    return T90conv(th + (del_th - 2 * q) / 6);
}

REAL adtg( REAL s, REAL t, REAL p){
    REAL T68;
    
    REAL a0 = 3.5803e-5;
    REAL a1 = 8.5258e-6;
    REAL a2 = -6.836e-8;
    REAL a3 = 6.6228e-10;
    REAL b0 = 1.8932e-6;
    REAL b1 = -4.2393e-8;
    REAL c0 = 1.8741e-8;
    REAL c1 = -6.7795e-10;
    REAL c2 = 8.733e-12;
    REAL c3 = -5.4481e-14;
    REAL d0 = -1.1351e-10;
    REAL d1 = 2.7759e-12;
    REAL e0 = -4.6206e-13;
    REAL e1 = 1.8676e-14;
    REAL e2 = -2.1687e-16;

    /*
    a = [3.5803e-5, 8.5258e-6, -6.836e-8, 6.6228e-10]
    b = [1.8932e-6, -4.2393e-8]
    c = [1.8741e-8, -6.7795e-10, 8.733e-12, -5.4481e-14]
    d = [-1.1351e-10, 2.7759e-12]
    e = [-4.6206e-13, 1.8676e-14, -2.1687e-16]
    */

    T68 = T68conv(t);
    return (a0 + (a1 + (a2 + a3 * T68) * T68) * T68 +
            (b0 + b1 * T68) * (s - 35) +
            ((c0 + (c1 + (c2 + c3 * T68) * T68) * T68) +
             (d0 + d1 * T68) * (s - 35)) * p +
            (e0 + (e1 + e2 * T68) * T68) * p * p);
}

REAL T90conv(REAL t){
    return t/1.00024;
}
REAL T68conv(REAL t){
    return t*1.00024;
}

REAL seck(REAL s, REAL t, REAL p){
    REAL T68, AW, BW, KW, A, B, K0;

    REAL h[] = {3.239908, 1.43713e-3, 1.16092e-4, -5.77905e-7};
    REAL k[] = {8.50935e-5, -6.12293e-6, 5.2787e-8};
    REAL e[] = {19652.21, 148.4206, -2.327105, 1.360477e-2, -5.155288e-5};
    REAL j0 = 1.91075e-4;
    REAL i[] = {2.2838e-3, -1.0981e-5, -1.6078e-6};
    REAL m[] = {-9.9348e-7, 2.0816e-8, 9.1697e-10};
    REAL f[] = {54.6746, -0.603459, 1.09987e-2, -6.1670e-5};
    REAL g[] = {7.944e-2, 1.6483e-2, -5.3009e-4};

    // Compute compression terms.
    p = p / 10.0;  // Convert from db to atmospheric pressure units.
    T68 = T68conv(t);

    //# Pure water terms of the secant bulk modulus at atmos pressure.
    //# UNESCO Eqn 19 p 18.
    //# h0 = -0.1194975
    AW = h[0] + (h[1] + (h[2] + h[3] * T68) * T68) * T68;

    //# k0 = 3.47718e-5
    BW = k[0] + (k[1] + k[2] * T68) * T68;

    //# e0 = -1930.06
    KW = e[0] + (e[1] + (e[2] + (e[3] + e[4] * T68) * T68) * T68) * T68;

    //# Sea water terms of secant bulk modulus at atmos. pressure.
    A = AW + (i[0] + (i[1] + i[2] * T68) * T68 + j0 * sqrt(s)) * s;

    B = BW + (m[0] + (m[1] + m[2] * T68) * T68) * s; // # Eqn 18.

    K0 = (KW + (f[0] + (f[1] + (f[2] + f[3] * T68) * T68) * T68 +
                (g[0] + (g[1] + g[2] * T68) * T68) * sqrt(s)) * s);  //# Eqn 16.
    return K0 + (A + B * p) * p; // # Eqn 15.

}
