/*
 * Header file for boundaries.c
 *
 * $Id: boundaries.h,v 1.3 2004-07-27 20:36:28 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2004/05/29 20:25:02  fringer
 * Revision before converting to CVS.
 *
 * Revision 1.1  2004/03/12 06:22:13  fringer
 * Initial revision
 *
 *
 */
#ifndef _boundaries_h
#define _boundaries_h

#include "suntans.h"
#include "phys.h"
#include "grid.h"

void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop);
void SetBoundaryScalars(REAL **boundary_scal, gridT *grid, physT *phys, propT *prop, char *type);

#endif
