/*
 * File: partition-noparmetis.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains functions that are defined when there are no parmetis
 * libraries present.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "partition.h"

/*
 * Function: GetPartitioning
 * Usage: GetPartitioning(maingrid,localgrid,myproc,numprocs,comm);
 * ----------------------------------------------------------------
 * This is a dummy function that is used when the parmetis libraries are not defined.
 * When the parmetis libraries are defined, see the file partition.c.
 *
 */
void GetPartitioning(gridT *maingrid, gridT **localgrid, int myproc, int numprocs, MPI_Comm comm) {
  int j, proc, Nclocal;
  if(numprocs>1) {
    if(myproc==0) printf("Error!  Cannot create grid for multiple processor simulation without ParMETIS!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else
    for(j=0;j<maingrid->Nc;j++)
      maingrid->part[j]=0;
}
