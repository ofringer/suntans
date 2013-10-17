/*
 * File: suntans.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Main header file
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _suntans_h
#define _suntans_h

#include "math.h"
#include "time.h"

// number of faces (for triangle)
#define REAL double
#define BUFFERLENGTH 256
#define NUMEDGECOLUMNS 3
#define BREAK printf("%d\n",*((int *)0));
#define PI 3.141592654
#define RHO0 1000.0
#define GRAV 9.81
#define INFTY 1e20
#define CONSERVED 1e-5
#define EMPTY 999999
#define SMALL 1e-15
#define CHECKCONSISTENCY 0
#define DRYCELLHEIGHT 1e-10
#define BUFFERHEIGHT 1e-2

// Error/Exit codes
#define EXIT_WRITING 1

#define DEFAULTDATAFILE "suntans.dat"

// define global variables for filenames
char DATADIR[BUFFERLENGTH],
  DATAFILE[BUFFERLENGTH],
  PSLGFILE[BUFFERLENGTH], 
  POINTSFILE[BUFFERLENGTH], 
  EDGEFILE[BUFFERLENGTH], 
  CELLSFILE[BUFFERLENGTH], 
  NODEFILE[BUFFERLENGTH], 
  INPUTDEPTHFILE[BUFFERLENGTH],
  CELLCENTEREDFILE[BUFFERLENGTH], 
  EDGECENTEREDFILE[BUFFERLENGTH], 
  VERTSPACEFILE[BUFFERLENGTH], 
  TOPOLOGYFILE[BUFFERLENGTH];
// define global variables
int TRIANGULATE, GRID, SOLVE, VERBOSE, WARNING, ASCII, RESTART, NUMPROCS, STEPSPERFILE;

#endif
