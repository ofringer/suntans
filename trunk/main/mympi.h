/*
 * Header file for mympi.c
 *
 * $Id: mympi.h,v 1.2 2002-11-05 01:31:17 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:21:54  fringer
 * Initial revision
 *
 *
 */
#ifndef _mympi_h
#define _mympi_h

#include<mpi.h>
#include "suntans.h"
#include "fileio.h"

void StartMpi(int *argc, char **argv[], MPI_Comm *comm, int *myproc, int *numprocs);
void EndMpi(MPI_Comm *comm);
REAL MPI_GetValue(char *file, char *str, char *call, int myproc);
void MPI_GetString(char *string, char *file, char *str, char *call, int myproc);

#endif
