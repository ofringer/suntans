/*
 * File: initialization.c
 * Description:  This file contains the functions that are used
 * to initialize the depth, free-surface, and salinity.
 *
 * $Id: initialization.c,v 1.2 2003-04-29 00:19:22 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2003/04/26 14:16:59  fringer
 * Initial revision
 *
 *
 */
#include "math.h"
#include "suntans.h"
#include "initialization.h"

/*
 * Function: ReturnDepth
 * Usage: grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
 * --------------------------------------------------------
 * Helper function to create a bottom bathymetry.  Used in
 * grid.c in the GetDepth function when IntDepth is 0.
 *
 */
REAL ReturnDepth(REAL x, REAL y) {
  REAL length;
  /*
  return 1.0-0.25*exp(-pow(x-5,2)-pow(y-5,2));
  */

  // Deep grid with a Gaussian hill.
  //return 3000 - 500*exp(-pow((x-50000)/10000,2));
  //			-pow((y-25000)/10000,2));

  // Tanh slope
  /*
  return 2500*.5*(1+tanh((x-60000)/5000))
    - 250*exp(-pow(x-30000)/500,2));
  */

  // Simple continental shelf
  length = 5000;
  if(x<=70000)
    return 3000;
  else if(x>70000 && x<70000+length)
    return 3000-2500*(x-70000)/length;
  else
    return 500;

  // Shallow grid with a Gaussian hill
  /*
  return 1-.75*exp(-pow((x-0)/2,2)
		   -pow((y-5)/2,2));
  */

  // Constant depth
  /*
  return 1;
  */
}

/*
 * Function: ReturnFreeSurface
 * Usage: grid->h[n]=ReturnFreeSurface(grid->xv[n],grid->yv[n]);
 * -------------------------------------------------------------
 * Helper function to create an initial free-surface. Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnFreeSurface(REAL x, REAL y) {
  return 0.0*cos(PI*x/10);0.25*exp(-pow(x-5,2)-pow(y-5,2));
}

/*
 * Function: ReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnSalinity(REAL x, REAL y, REAL z) {
  REAL delta_s = 0.0147;
  REAL dmax = 3000;

  return -delta_s*(z/dmax+.5);
  //phys->s[i][k]=-1.07*0*(z+1500+0*cos(PI*(grid->xv[i]-10000)/10000))/1500;
  //phys->s[i][k]=-7.333e-3*(.5*(z+1500)/1500);
  //phys->s0[i][k]=-7.333e-3*(.5*(z+1500)/1500);
  //if(grid->xv[i]>70000) phys->s[i][k]=.01;
  //if(grid->xv[i]>500 && grid->dzz[i][k]>0)
  //if(grid->xv[i]>5) phys->s[i][k]=0.01;
  //phys->s[i][k]=-(1-tanh((grid->xv[i]-5)/.25));
  //phys->s[i][k]=-tanh((z+0.5*grid->dv[i]-.1*cos(PI*grid->xv[i]/10))/.01);
  //phys->s[i][k]=-0*tanh((z+0.5+.3*exp(-pow(grid->xv[i]/2,2)))/.01);
  //if(grid->xv[i]>500 && grid->dzz[i][k]>0)
  //  phys->s[i][k]=0;
  //  r=sqrt(pow(grid->xv[i]-11.25,2)+pow(grid->yv[i]-7.5,2));
  //if(r<2)
  //  phys->s[i][k]=1;
  //else
  //  phys->s[i][k]=0;
  //phys->s[i][k]=1-tanh((grid->xv[i]-7.5)/2);
}

/*
 * Function: ReturnHorizontalVelocity
 * Usage: grid->u[n]=ReturnHorizontalVelocity(grid->xv[n],grid->yv[n],
 *                                            grid->n1[n],grid->n2[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial velocity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z) {
  REAL u, v, umag=0;
  
  // Irrotational vortex.
  u = -umag*(y-5)/5;
  v = umag*(x-5)/5;
  return u*n1+v*n2;
}

