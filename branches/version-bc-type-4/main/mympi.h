/*
 * Header file for mympi.c
 *
 * $Id: mympi.h,v 1.4 2004-05-29 20:25:02 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2003/05/01 00:31:58  fringer
 * Added MPI_FOpend and MPI_GetSize function prototypes.
 *
 * Revision 1.2  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
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
void MPI_GetFile(char *string, char *file, char *str, char *call, int myproc);
FILE *MPI_FOpen(char *file, char *perms, char *caller, int myproc);
int MPI_GetSize(char *file, char *caller, int myproc);

#endif
