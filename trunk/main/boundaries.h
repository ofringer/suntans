/*
 * Header file for boundaries.c
 *
 * $Id: boundaries.h,v 1.4 2005-04-01 22:39:35 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2004/07/27 20:36:28  fringer
 * Added SetBoundaryScalars function prototype.
 *
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

// Enumerated type for open/specified bc specification
enum {
  specified, open
};

void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop);
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop);
void BoundaryScalars(gridT *grid, physT *phys, propT *prop);
void WindStress(gridT *grid, physT *phys, propT *prop);

#endif
