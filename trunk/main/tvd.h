/*
 * File: tvd.h
 * Author: Zhonghua Zhang
 * Institution: Stanford University
 * Date: 03/06/06
 * ----------------------------------------
 * Header file for tvd.c.
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
#define TVDMACRO 0

void HorizontalFaceScalars(gridT *grid, physT *phys, propT *prop, REAL **boundary_scal, int TVD);
void GetApAm(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	     REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD);

#endif


