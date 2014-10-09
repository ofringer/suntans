/*
 * File: triangulate-notriangle.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Returns an error if the -t option is used when the triangle libraries are 
 * not installed.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "mympi.h"
#include "grid.h"

/*
 * Function: GetTriangulation
 * Usage: GetTriangulation(grid);
 * ------------------------------
 * Dummy function to return an error if the triangle libraries are not installed.
 *
 */
int GetTriangulation(gridT **grid, int myproc) {
  if(myproc==0) {
    printf("Error! Cannot run with the -t option without the triangle libraries!\n");
    printf("Create a grid using a grid generation tool and use the -g option.\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  return -1;
}
