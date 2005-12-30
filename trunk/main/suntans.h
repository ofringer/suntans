/*
 * Main header file
 *
 * $Id: suntans.h,v 1.9 2005-12-30 23:27:28 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.8  2005/07/12 01:07:29  fringer
 * Decreased the conserved cutoff to 1e-10 instead of 1e-5.
 *
 * Revision 1.7  2004/09/22 06:31:25  fringer
 * Added a definition of SMALL which is used in turbulence.c to prevent
 * denominators from vanishing.
 *
 * Revision 1.6  2004/05/29 20:25:02  fringer
 *  Revision before converting to CVS.
 *
 * Revision 1.5  2003/04/29 00:20:32  fringer
 * Changed INFTY to 1e20, added EMPTY 999999, and added TRIANGULATE and ASCII.
 *
 * Revision 1.4  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
 * Revision 1.3  2002/11/03 00:54:57  fringer
 * Also removed the DATADIR line since the main file is the only
 * file that contains any information about directories, etc...
 *
 * Revision 1.2  2002/11/03 00:52:35  fringer
 * Removed the requirement of specifying a separate file that contains the
 * grid data files.  They are specified in the main suntans data file now.
 *
 * Revision 1.1  2002/11/03 00:20:49  fringer
 * Initial revision
 *
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
