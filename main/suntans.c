/*
 * File: suntans.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * SUNTANS - Stanford Unstructured Nonhydrostatic Terrain-following Adaptive
 * Navier-Stokes Simulator
 *
 * http://suntans.stanford.edu
 * 
 * Written by Oliver B. Fringer
 * Dept. of Civil and Environmental Engineering
 * Stanford University
 * 
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
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
  propT *prop;

  StartMpi(&argc,&argv,&comm,&myproc,&numprocs);

  ParseFlags(argc,argv,myproc);
  
  if(GRID)
    GetGrid(&grid,myproc,numprocs,comm);
  else
    ReadGrid(&grid,myproc,numprocs,comm);
  

  if(SOLVE) {
    ReadProperties(&prop,myproc);
    InitializeVerticalGrid(&grid,myproc);
    AllocatePhysicalVariables(grid,&phys,prop);
    AllocateTransferArrays(&grid,myproc,numprocs,comm);
    OpenFiles(prop,myproc);
    if(RESTART)
      ReadPhysicalVariables(grid,phys,prop,myproc,comm);
    else
      InitializePhysicalVariables(grid,phys,prop,myproc,comm);

    Solve(grid,phys,prop,myproc,numprocs,comm);
    //    FreePhysicalVariables(grid,phys,prop);
    //    FreeTransferArrays(grid,myproc,numprocs,comm);
  }

  EndMpi(&comm);
}





