/*
 * File: initialization.c
 * Description:  This file contains the functions that are used
 * to initialize the depth, free-surface, and salinity.
 *
 * $Id: initialization.c,v 1.1 2003-04-26 14:16:59 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#include "math.h"
#include "suntans.h"
#include "initialization.h"

REAL ReturnDepth(REAL x, REAL y) {
  return 1.0-0.25*exp(-pow(x-5,2)-pow(y-5,2));
}

REAL ReturnFreeSurface(REAL x, REAL y) {
  return 0.05*cos(PI*x/10);0.25*exp(-pow(x-5,2)-pow(y-5,2));
}

REAL ReturnSalinity(REAL x, REAL y, REAL z) {
  return 0;
}

//      grid->dv[n] = 3000 - 0*500*exp(-pow((grid->xv[n]-50000)/10000,2)
//				   -pow((grid->yv[n]-25000)/10000,2));

//2500*.5*(1+tanh((grid->xv[n]-60000)/5000))
//	- 250*exp(-pow((grid->xv[n]-30000)/500,2));

/*if(grid->xv[n]<=70000)
     grid->dv[n]=3000;
     else if(grid->xv[n]>70000 && grid->xv[n]<80000)
     grid->dv[n]=3000-2500*(grid->xv[n]-70000)/10000;
     else
     grid->dv[n]=500;
*/
//      grid->dv[n]=1-.75*exp(-pow((grid->xv[n]-0)/2,2)
//      			    -pow((grid->yv[n]-5)/2,2));
//      grid->dv[n]=1-.25*0.5*(1+tanh(-(sqrt(pow(grid->xv[n]-0,2)
//			    +pow(grid->yv[n]-7,2))-1)/.0001));
//grid->dv[n]=1;
//grid->dv[n] = 3000 - 1000*exp(-pow((grid->xv[n]-50000)/5000,2));
//      grid->dv[n] = 5;
//      if(grid->xv[n]<500) grid->dv[n]=20;
//      grid->dv[n]=1;
