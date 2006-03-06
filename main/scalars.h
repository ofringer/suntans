/*
 * File: scalars.h
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 03/06/06 
 * ----------------------------------------
 * Header file for scalars.c
 *
 */
#ifndef _scalars_h
#define _scalars_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

void UpdateScalars(gridT *grid, physT *phys, propT *prop, REAL **scal, REAL **boundary_scal, REAL **Cn, 
		   REAL kappa, REAL kappaH, REAL **kappa_tv, REAL theta,
		   REAL **src1, REAL **src2, REAL *Ftop, REAL *Fbot, int alpha_top, int alpha_bot);

#endif
