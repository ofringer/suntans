#include "math.h"
#include "fileio.h"
#include "suntans.h"
#include "initialization.h"

#define sech 1/cosh
/*
 * Function: GetDZ
 * Usage: GetDZ(dz,depth,Nkmax,myproc);
 * ------------------------------------
 * Returns the vertical grid spacing in the array dz.
 *
 */
int GetDZ(REAL *dz, REAL depth, REAL localdepth, int Nkmax, int myproc) {
  int k, status, ktop=15;
  REAL z=0, dz0, r = GetValue(DATAFILE,"rstretch",&status), dztop=2, d0;

  if(Nkmax==1) {
    dz[0]=depth;
  } else {
    dz[0]=5;
    for(k=1;k<ktop+1;k++)
      dz[k]=dztop;
    if(r==1)
      dz[ktop+1] = (depth-ktop*dztop-dz[0])/(Nkmax-ktop-1);
    else
      dz[ktop+1] = (depth-ktop*dztop-dz[0])*(r-1)/(pow(r,(Nkmax-ktop-1))-1);
    for(k=ktop+2;k<Nkmax;k++) 
      dz[k]=r*dz[k-1];
    
    d0=0;
    for(k=0;k<Nkmax;k++) {
      //printf("dz[%d]=%f\n",k,dz[k]);
      d0+=dz[k];
    }
    //printf("dz sum = %f\n",d0);
  }

  return 0;
}
  
/*
 * Function: ReturnDepth
 * Usage: grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
 * --------------------------------------------------------
 * Helper function to create a bottom bathymetry.  Used in
 * grid.c in the GetDepth function when IntDepth is 0.
 *
 */
REAL ReturnDepth(REAL x, REAL y) {
  REAL length, xmid, shelfdepth, depth;

  return 10;
}

 /*
  * Function: ReturnFreeSurface
  * Usage: grid->h[n]=ReturnFreeSurface(grid->xv[n],grid->yv[n]);
  * -------------------------------------------------------------
  * Helper function to create an initial free-surface. Used
  * in phys.c in the InitializePhysicalVariables function.
  *
  */
REAL ReturnFreeSurface(REAL x, REAL y, REAL d) {
  return -5;
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
  REAL x0 = 5.7e5;
  if(y>4.2e6)
    return 32*0.5*(1+tanh(-(x-x0)/10000));
  return 32;
}

/*
 * Function: ReturnTemperature
 * Usage: grid->T[n]=ReturnTemperaturegrid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial temperature field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnTemperature(REAL x, REAL y, REAL z, REAL depth) {
  REAL x0 = 545691, y0 = 4185552, R=1000;
  if(pow(x-x0,2)+pow(y-y0,2)<R*R)
    return 1;
  return 0;
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
  return 0;
}

