/*
 * File: partition.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 09/30/02
 * --------------------------------
 * This file reads in and partitions the unstructured grid and 
 * writes files that contain grid information for each processor.
 *
 * $Id: suntans.c,v 1.2 2002-11-05 01:31:17 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:16:10  fringer
 * Initial revision
 *
 *
 */
#include "suntans.h"
#include "mympi.h"
#include "grid.h"
#include "phys.h"
#include "report.h"

main(int argc, char *argv[])
{
  int myproc, numprocs, j;
  MPI_Comm comm;
  gridT *grid;
  physT *phys;

  StartMpi(&argc,&argv,&comm,&myproc,&numprocs);

  ParseFlags(argc,argv,myproc);
  
  if(GRID)
    GetGrid(&grid,myproc,numprocs,comm);

  if(SOLVE) {
    InitializeVerticalGrid(&grid);
    AllocatePhysicalVariables(grid,&phys);
    InitializePhysicalVariables(grid,phys);
    Solve(grid,phys,myproc,numprocs,comm);
    FreePhysicalVariables(grid,phys);
  }

  EndMpi(&comm);
}





