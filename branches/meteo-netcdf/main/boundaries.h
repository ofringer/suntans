/*
 * File: boundaries.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for boundaries.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _boundaries_h
#define _boundaries_h

#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "met.h"

// Enumerated type for open/specified bc specification
enum {
  specified, open
};

// Structure array to store netcdf boundary condition data
typedef struct _boundT{
  
  // Dimension sizes
  size_t Ntype2;
  size_t Ntype3;
  size_t Nseg;
  size_t Nt;
  size_t Nk;
  
  // boolean operators
  int hasType2;
  int hasType3;
  int hasSeg;

  // Grid cell indices
  int *edgep;
  int *localedgep;
  int *cellp;
  int *segedgep;
  int *segp;

  // Indices that point the grid to the cell in the file
  int *ind2;
  int *ind3;
  int *ind3edge;

  // Boundary coordinates
  REAL *xe;
  REAL *ye;
  REAL *xv;
  REAL *yv;
  REAL *z;
  REAL *time;	
  REAL *segarea;
  REAL *localsegarea;

  // Time record locators
  int t0;
  int t1;

  // Data arrays at forward (_f) and backward (_b) timestep
  // Type-2 (edge centred) boundaries
  REAL **boundary_u_f;
  REAL **boundary_v_f;
  REAL **boundary_w_f;
  REAL **boundary_T_f;
  REAL **boundary_S_f;
  REAL *boundary_Q_f;

  REAL **boundary_u_b;
  REAL **boundary_v_b;
  REAL **boundary_w_b;
  REAL **boundary_T_b;
  REAL **boundary_S_b;
  REAL *boundary_Q_b;

  REAL **boundary_u;
  REAL **boundary_v;
  REAL **boundary_w;
  REAL **boundary_T;
  REAL **boundary_S;
  REAL *boundary_Q;

  // Type-3 (cell centred) boundaries
  REAL **uc_f;
  REAL **vc_f;
  REAL **wc_f;
  REAL **T_f;
  REAL **S_f;
  REAL * h_f;

  REAL **uc_b;
  REAL **vc_b;
  REAL **wc_b;
  REAL **T_b;
  REAL **S_b;
  REAL *h_b;

  REAL **uc;
  REAL **vc;
  REAL **wc;
  REAL **T;
  REAL **S;
  REAL *h;
} boundT;

// Declare the boundary structure global
boundT *bound;

void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop);
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void BoundaryScalars(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void WindStress(gridT *grid, physT *phys, propT *prop, metT *met, int myproc);

FILE *windFID;

#ifdef USENETCDF
    void InitBoundaryData(propT *prop, gridT *grid, int myproc);
    void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc);
    void UpdateBdyNC(propT *prop, gridT *grid, int myproc, MPI_Comm comm);
#endif

#endif
