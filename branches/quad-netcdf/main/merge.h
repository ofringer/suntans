#ifndef _merge_h
#define _merge_h

#include "suntans.h"
#include "grid.h"
#include "mympi.h"

gridT *mergedGrid;
int *Nc_all, Nc_max;
int **mnptr_all, *send3DSize;
REAL *localTempMergeArray, **mergedArray;

void InitializeMerging(gridT *grid, int numprocs, int myproc, MPI_Comm comm);
void MergeCellCenteredArray(REAL **localArray, gridT *grid, int dims, int numprocs, int myproc, MPI_Comm comm);
void FreeMergingArrays(gridT *grid, int myproc);

#endif
