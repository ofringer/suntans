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
 * By using this software, you understand that:
 * 
 * 1) Because of its complexity and the expertise required to use
 *    SUNTANS, and because it is under constant development, we make no
 *    guarantees as to the stability of the code nor the accuracy of the
 *    results you obtain, and while you will be a member of an active
 *    user group and an associated email list and will have access to
 *    extensive documentation, we cannot guarantee timely support or
 *    advice for its use. We will, however, try to answer your specific
 *    problems when they are posted to the SUNTANS user email list.
 * 
 * 2) You are using the source code for academic purposes alone and are
 *    not using SUNTANS for any contractual or commercial work
 *    whatsoever.
 * 
 * 3) While you are free to edit and alter the source code, you are
 *    forbidden to distribute any portion of the source to anyone or any
 *    other organization.
 * 
 * 4) You will site the use of SUNTANS in any publication that comes out
 *    of the results obtained from SUNTANS as:
 * 
 *    Fringer, O.B., Stanford University. 2006. "Source code for SUNTANS:
 *    Stanford Unstructured Nonhydrostatic Terrain-following Adaptive
 *    Navier-Stokes Simulator", URL: http://suntans.stanford.edu. Copyright
 *    (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior
 *    University. All Rights Reserved.
 * 
 *    where (C) is the copyright symbol.
 * 
 * 5) Upon our request, you will make all changes and/or additions to the
 *    source available to us, and therefore all changes may become part
 *    of the downloadable source that is available to registered SUNTANS
 *    users.
 * 
 * 6) We reserve the right to deny a request for downloading the SUNTANS
 *    source code if, based on the information you provide, we feel you
 *    do not possess the expertise required for a successful
 *    implementation or if we feel that the application is not
 *    appropriate for the use of SUNTANS. In that case we will recommend
 *    an alternate tool that may be more appropriate for your needs and
 *    abilities.
 * 
 * 7) You will use SUNTANS at your own risk, and Stanford does not
 *    represent that it is accurate or up-to-date. Stanford University
 *    will have no liability to you or to any third party as a result of
 *    its use of SUNTANS, and you will indemnify and hold Stanford
 *    University harmless from any claims related to your use of
 *    SUNTANS. By clicking on this link you are agreeing to the terms of
 *   downloading the source code for SUNTANS.
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
    InitializeVerticalGrid(&grid);
    AllocatePhysicalVariables(grid,&phys,prop);
    AllocateTransferArrays(&grid,myproc,numprocs,comm);
    OpenFiles(prop,myproc);
    if(RESTART)
      ReadPhysicalVariables(grid,phys,prop,myproc,comm);
    else
      InitializePhysicalVariables(grid,phys,prop);
    Solve(grid,phys,prop,myproc,numprocs,comm);
    //    FreePhysicalVariables(grid,phys,prop);
    //    FreeTransferArrays(grid,myproc,numprocs,comm);
  }

  EndMpi(&comm);
}





