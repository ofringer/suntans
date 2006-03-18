/*
 * File: partition.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for partition.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _partition_h
#define _partition_h

#include "grid.h"

void GetPartitioning(gridT *maingrid, gridT **localgrid, int myproc, int numprocs, MPI_Comm comm);

#endif
