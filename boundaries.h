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

// Enumerated type for open/specified bc specification
enum {
  specified, open
};

void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop);
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc);
void BoundaryScalars(gridT *grid, physT *phys, propT *prop);
void WindStress(gridT *grid, physT *phys, propT *prop, int myproc);

FILE *windFID;

#endif
