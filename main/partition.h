/*
 * Header file for partition.c
 *
 */
#ifndef _partition_h
#define _partition_h

#include "grid.h"

void GetPartitioning(gridT *maingrid, gridT **localgrid, int myproc, int numprocs, MPI_Comm comm);

#endif
