/*
 * File: physio.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for physio.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _physio_h
#define _physio_h

#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "mympi.h"

void Write2DData(REAL *array, int merge, FILE *fid, char *error_message, 
		 gridT *grid, int numprocs, int myproc, MPI_Comm comm);
void Write3DData(REAL **array, REAL *temp_array, int merge, FILE *fid, char *error_message, 
		 gridT *grid, int numprocs, int myproc, MPI_Comm comm);
void OpenFiles(propT *prop, int myproc);
void ReadPhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void OutputPhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
void OutputRestartVariables(gridT *grid, physT *phys, propT *prop,int myproc, int numprocs, int blowup, MPI_Comm comm);

#endif
