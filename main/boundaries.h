/*
 * Header file for boundaries.c
 *
 * $Id: boundaries.h,v 1.1 2004-03-12 06:22:13 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#ifndef _boundaries_h
#define _boundaries_h

#include "suntans.h"
#include "phys.h"
#include "grid.h"

void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop);

#endif
