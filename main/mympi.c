/*
 * File: mympi.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 09/30/02
 * --------------------------------
 * This file contains mpi-based functions.
 *
 * $Id: mympi.c,v 1.2 2003-04-29 16:37:23 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:21:16  fringer
 * Initial revision
 *
 *
 */
#include <errno.h>
#include "mympi.h"

/*
 * Function: StartMpi
 * Usage: StartMpi(argc,argv,&comm,&numprocs,&myproc);
 * ----------------------------------------------------
 * Start up the mpi communicator and determine the number of processors
 * and the id of the current processor.  Also, parse the command line
 * and exit if it is incorrect.
 *
 */
void StartMpi(int *argc, char **argv[], MPI_Comm *comm, int *myproc, int *numprocs)
{
  MPI_Init(argc,argv);
  MPI_Comm_dup(MPI_COMM_WORLD,comm);
  MPI_Comm_size(*comm, numprocs);
  MPI_Comm_rank(*comm, myproc);
}

/*
 * Function: DndMpi
 * Usage: EndMpi(comm);
 * ---------------------
 * Free the communicator and kill mpi.
 *
 */
void EndMpi(MPI_Comm *comm)
{
  MPI_Comm_free(comm);
  MPI_Finalize();
}

/*
 * Function: MPI_GetValue
 * Usage: x = MPI_GetValue("file.dat","xval",myproc);
 * --------------------------------------------------
 * This is a wrapper function for the GetValue function
 * in fileio.c.  It is identical, except here it checks
 * for an error and quits MPI if there is one.
 *
 */
REAL MPI_GetValue(char *file, char *str, char *call, int myproc)
{
  int status;
  REAL val = GetValue(file,str,&status);

  if(!status) {
    if(myproc==0) {
      printf("Error in MPI_GetValue (called from %s) while reading %s.\n",call,file);
      printf("Couldn't find value for %s!\n",str);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  return val;
}

void MPI_GetString(char *string, char *file, char *str, char *call, int myproc)
{
  int status;
  GetString(string,file,str,&status);

  if(!status) {
    if(myproc==0) {
      printf("Error in MPI_GetString (called from %s) while reading %s.\n",call,file);
      printf("Couldn't find string for %s!\n",str);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
}

/*
 * Function: MPI_FOpen
 * Usage: fid = MPI_FOpen(string,"r","GetValue",myproc);
 * -----------------------------------------------------
 * Exits if the requested file does not exist and closes
 * MPI cleanly. The third string is useful for determining which
 * function the function was called from.  When two processes
 * are trying to read the same file at the same time, an error
 * code of EAGAIN results.  This is ommitted as a possible error
 * code.
 *
 */
FILE *MPI_FOpen(char *file, char *perms, char *caller, int myproc) {
  extern int errno;
  char str[BUFFERLENGTH];
  FILE *fid = fopen(file,perms);

  if(errno && errno!=EAGAIN) {
    if(myproc==0) {
      sprintf(str,"Error in Function %s while trying to open %s",caller,file);
      perror(str);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else
    return fid;
}

/*
 * Function: MPI_GetSize
 * Usage: N = MPI_GetSize(filename,"ReadMainGrid",myproc);
 * -------------------------------------------------------
 * Uses the getsize function defined in fileio.c but checks
 * for file existence first.
 *
 */
int MPI_GetSize(char *file, char *caller, int myproc) {
  FILE *fid = MPI_FOpen(file,"r",caller,myproc);
  fclose(fid);

  return getsize(file);
}
