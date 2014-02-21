/*
 * File: turbulence.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for turbulence.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _turbulence_h
#define _turbulence_h

// Von Karman's Constant
#define KAPPA_VK 0.42

// Background length scale
#define LBACKGROUND 1e-12

void my25(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **q, REAL **l, REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc);
//void my25(gridT *grid, physT *phys, propT *prop, REAL **q, REAL **l, REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc) ;
#endif
