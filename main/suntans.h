/*
 * Main header file
 *
 * $Id: suntans.h,v 1.2 2002-11-03 00:52:35 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:20:49  fringer
 * Initial revision
 *
 *
 */
#ifndef _suntans_h
#define _suntans_h

#define NFACES 3
#define REAL double
#define BUFFERLENGTH 256
#define NUMEDGECOLUMNS 3
#define BREAK printf("%d\n",*((int *)0));
#define PI 3.141592654
#define GRAV 9.81
#define INFTY 10000
#define CONSERVED 1e-5

#define DATADIR "/home/fringer/research/SUNTANS/data"
#define GRIDDATAFILELIST DATADIR"/suntans_grid_files.dat"
#define DATAFILE "suntans.dat"

char POINTSFILE[BUFFERLENGTH], EDGEFILE[BUFFERLENGTH], CELLSFILE[BUFFERLENGTH], INPUTDEPTHFILE[BUFFERLENGTH],
  CELLCENTEREDFILE[BUFFERLENGTH], EDGECENTEREDFILE[BUFFERLENGTH], VERTSPACEFILE[BUFFERLENGTH], 
  TOPOLOGYFILE[BUFFERLENGTH];
int GRID, SOLVE, VERBOSE, WARNING;

#endif
