/*
 * File: partition.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 09/30/02
 * --------------------------------
 * This file reads in and partitions the unstructured grid and 
 * writes files that contain grid information for each processor.
 *
 * $Id: suntans.c,v 1.5 2004-09-22 06:29:54 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.4  2004/05/29 20:25:02  fringer
 *  Revision before converting to CVS.
 *
 * Revision 1.3  2003/05/05 01:28:14  fringer
 * Added ReadGrid line to add functionality for use of -s flag without
 * -g flag .
 *
 * Revision 1.2  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
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
  propT *prop;

  StartMpi(&argc,&argv,&comm,&myproc,&numprocs);

  ParseFlags(argc,argv,myproc);
  
  if(GRID)
    GetGrid(&grid,myproc,numprocs,comm);
  else
    ReadGrid(&grid,myproc,numprocs,comm);

  if(SOLVE) {
    InitializeVerticalGrid(&grid);
    AllocatePhysicalVariables(grid,&phys,prop);
    ReadProperties(&prop,myproc);
    OpenFiles(prop,myproc);
    if(RESTART)
      ReadPhysicalVariables(grid,phys,prop,myproc);
    else
      InitializePhysicalVariables(grid,phys,prop);
    SetDragCoefficients(grid,phys,prop);
    Solve(grid,phys,prop,myproc,numprocs,comm);
    //    FreePhysicalVariables(grid,phys,prop);
  }

  EndMpi(&comm);
}





