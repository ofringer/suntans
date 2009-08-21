/*
 * File: tvd.h
 * Author: Zhonghua Zhang and Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for tvd.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _tvd_h
#define _tvd_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"
#include "fileio.h"

// TVD Method.  0=No TVD (Uses method in original scheme, which is First-order upwind)
// Otherwise, TVD is implemented, with TVDMACRO=1: First-order upwind, 2: Lax-Wendroff, 3: Superbee, 4: Van Leer
// Note that this TVD implementation does not work with wetting and drying and is 
// not strictly tvd.
#define TVDMACRO 0

void HorizontalFaceScalars(gridT *grid, physT *phys, propT *prop, REAL **boundary_scal, int TVD,
			   MPI_Comm comm, int myproc); 
void GetApAm(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	     REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD);

#endif


