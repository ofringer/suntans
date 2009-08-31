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

#define NFACES 3
#define REAL double
#define BUFFERLENGTH 256
#define NUMEDGECOLUMNS 3
#define BREAK printf("%d\n",*((int *)0));
#define PI 3.141592654
#define GRAV 9.81
#define RHO0 1000.0
#define INFTY 1e20
#define CONSERVED 1e-10
#define EMPTY 999999
#define SMALL 1e-15
#define CHECKCONSISTENCY 1
#define DEFAULT_hprecond 1
#define DRYCELLHEIGHT 1e-10

// Error/Exit codes
#define EXIT_WRITING 1

#define DEFAULTDATAFILE "suntans.dat"

char DATADIR[BUFFERLENGTH],
  DATAFILE[BUFFERLENGTH],
  PSLGFILE[BUFFERLENGTH], 
  POINTSFILE[BUFFERLENGTH], 
  EDGEFILE[BUFFERLENGTH], 
  CELLSFILE[BUFFERLENGTH], 
  INPUTDEPTHFILE[BUFFERLENGTH],
  CELLCENTEREDFILE[BUFFERLENGTH], 
  EDGECENTEREDFILE[BUFFERLENGTH], 
  VERTSPACEFILE[BUFFERLENGTH], 
  TOPOLOGYFILE[BUFFERLENGTH];
int TRIANGULATE, GRID, SOLVE, VERBOSE, WARNING, ASCII, RESTART;

#endif
