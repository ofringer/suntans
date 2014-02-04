/*
 * File: merge.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for merge.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _merge_h
#define _merge_h

#include "suntans.h"
#include "grid.h"
#include "mympi.h"

gridT *mergedGrid;
int *Nc_all, Nc_max;
int **mnptr_all, *send3DSize;
REAL *localTempMergeArray, *merged2DArray, **merged3DArray;

void InitializeMerging(gridT *grid, int numprocs, int myproc, MPI_Comm comm);
void MergeCellCentered2DArray(REAL *localArray, gridT *grid, int numprocs, int myproc, MPI_Comm comm);
void MergeCellCentered3DArray(REAL **localArray, gridT *grid, int numprocs, int myproc, MPI_Comm comm);
void FreeMergingArrays(gridT *grid, int myproc);

#endif
