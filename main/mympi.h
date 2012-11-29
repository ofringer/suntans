/*
 * File: mympi.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for mympi.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _mympi_h
#define _mympi_h

#ifndef NOMPI
#include "mpi.h"
#else
#include "no-mpi.h"
#endif

#include "suntans.h"
#include "fileio.h"

void StartMpi(int *argc, char **argv[], MPI_Comm *comm, int *myproc, int *numprocs);
void EndMpi(MPI_Comm *comm);
REAL MPI_GetValue(char *file, char *str, char *call, int myproc);
void MPI_GetString(char *string, char *file, char *str, char *call, int myproc);
void MPI_GetFile(char *string, char *file, char *str, char *call, int myproc);
FILE *MPI_FOpen(char *file, char *perms, char *caller, int myproc);
int MPI_GetSize(char *file, char *caller, int myproc);
int MPI_NCOpen(char *file, int perms, char *caller, int myproc);

#endif
