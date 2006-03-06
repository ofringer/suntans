/*
 * File: check.h
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 01/27/06 (Mozart's 250th birthday)
 * ----------------------------------------
 * Header file for check.c.
 *
 */
#ifndef _check_h
#define _check_h

#include "grid.h"
#include "phys.h"
#include "mympi.h"

int Check(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
int CheckDZ(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void Progress(propT *prop, int myproc);

#endif
