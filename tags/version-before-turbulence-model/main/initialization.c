/*
 * File: initialization.c
 * Description:  This file contains the functions that are used
 * to initialize the depth, free-surface, and salinity.
 *
 * $Id: initialization.c,v 1.10 2004-07-22 20:24:52 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.9  2004/06/23 05:24:44  fringer
 * Committing as a test since fringer is not a writer.
 *
 * Revision 1.8  2004/06/17 02:54:53  fringer
 * Added new Monterey density field which uses the swstate m-file to obtain
 * the density from salinity and temperature.
 *
 * Revision 1.7  2004/06/15 18:24:40  fringer
 * Cleaned up and removed extraneous code for testing.
 *
 * Revision 1.6  2004/06/14 07:56:41  fringer
 * Changes made during debugging.
 *
 * Revision 1.5  2004/06/01 04:31:29  fringer
 * Initialization for entrainment problem:
 * D = 2 m when x>4 & x<7, 1 m otherwise
 * s -> -.005*tanh((x-3)/.01), +.01 for s<-1
 *
 * Revision 1.4  2004/05/31 06:59:33  fringer
 * Added vertical diffusion of scalar.  Left in Ftop and Fbot
 * as commented out terms.  These are where top and bottom diffusive
 * fluxes are included.
 *
 * Revision 1.3  2004/05/29 20:25:02  fringer
 * Revision before converting to CVS.
 *
 * Revision 1.2  2003/04/29 00:19:22  fringer
 * Added all initialization functions.
 *
 * Revision 1.1  2003/04/26 14:16:59  fringer
 * Initial revision
 *
 *
 */
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
  int k, status;
  REAL z=0, dz0, r = GetValue(DATAFILE,"rstretch",&status);

  if(dz!=NULL) {
    if(r==1) 
      for(k=0;k<Nkmax;k++)
	dz[k]=depth/Nkmax;
    else if(r>1 && r<=1.1) {    
      dz[0] = depth*(r-1)/(pow(r,Nkmax)-1);
      if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[0]);
      for(k=1;k<Nkmax;k++) 
	dz[k]=r*dz[k-1];
    } else {
      printf("Error in GetDZ when trying to create vertical grid:\n");
      printf("Stretching parameter rstretch must  be in the range (1,1.1).\n");
      exit(1);
    }
  } else {
    if(r!=1)
      dz0 = depth*(r-1)/(pow(r,Nkmax)-1);
    else
      dz0 = depth/Nkmax;
    z = dz0;
    for(k=1;k<Nkmax;k++) {
      dz0*=r;
      z+=dz0;
      if(z>=localdepth) {
	return k;
      }
    }
  }
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
  int nc, np, Nc = 13;
  REAL length, xmid, *xc, *yc, R=0.025, shelfdepth=200;

  return 1000-0*975*exp(-pow((x-70000)/10000,2)-pow((y-60000)/10000,2));

  // Huntington beach continental shelf
  /*
  if(x<6418)
   return 590-145*exp(x/5000);
  else
   return 62-55*(x-7000)/7000;
  */
  /*
   length = 2000;
   if(x<=6500-length/2)
     return 120;
   else if(x>6500-length/2 && x<=6500+length/2 && length>0)
     return 120-60*(x-6500+length/2)/length;
   else
     return 60;
   */
   // Monterey Bay continental shelf
   length = 10000;
   xmid = 75000;
   shelfdepth = 500;
   if(x<=xmid-length/2)
     return 3000;
   else if(x>xmid-length/2 && x<=xmid+length/2 && length>0)
     return 3000-(3000-shelfdepth)*(x-xmid+length/2)/length;
   else
     return shelfdepth;
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
  REAL r, rp, zeta, zetav, s=tan(PI/3), x0, y0, sig;
  
  return 0*cos(PI*x/100000);
  sig = 3000;
  x0 = 40000;
  y0 = 65000;
  r = s*(x-x0)-(y-y0);
  rp = (x-x0)+s*(y-y0);
  zeta = .5*(-tanh((rp-40000)/sig)+tanh((rp+40000)/sig));
  zeta = pow(sech(-(r/sig)),2)*zeta;
  //  zeta = exp(-fabs(r/sig))*zeta;

  x0 = 50000;
  y0 = 40000;
  zetav = .5*(-tanh((y-y0-20000)/sig)+tanh((y-y0+20000)/sig));
  zetav = zetav*pow(sech(-s*(x-x0)/sig),2);
  //  zetav = zetav*exp(-s*fabs(x-50000)/sig);

  return 1000*(zeta+zetav);
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
  REAL delta_s, thermocline_depth=20, z0, zeta, r, rp, zetav, x0, y0, sig, s=tan(PI/3);
  REAL power, factor,dmax=120,dshelf=60;//dmax=3000,dshelf=500;//dmax=120,dshelf=60;//dmax = 3000;
  thermocline_depth=35;
  dmax = 120;
  dshelf = 60;

  sig = 3000;
  x0 = 40000;
  y0 = 65000;
  r = s*(x-x0)-(y-y0);
  rp = (x-x0)+s*(y-y0);
  zeta = .5*(-tanh((rp-40000)/sig)+tanh((rp+40000)/sig));
  zeta = pow(sech(-(r/sig)),2)*zeta;
  //  zeta = exp(-fabs(r/sig))*zeta;

  x0 = 50000;
  y0 = 40000;
  zetav = .5*(-tanh((y-y0-20000)/sig)+tanh((y-y0+20000)/sig));
  zetav = zetav*pow(sech(-s*(x-x0)/sig),2);
  //  zetav = zetav*exp(-s*fabs(x-50000)/sig);

  return -.5*.04286*tanh((z+50+100*(zeta+zetav))/10);
  //  return -.0728*tanh(-(x-.4)/.001);

  thermocline_depth=25;
  dmax = 3000;
  
  if(z>-thermocline_depth)
    return 3.4286*pow(fabs(thermocline_depth),0.0187)-3.6;
  return 3.4286*pow(fabs(z),0.0187)-3.6;

  // Monterey (critical when length=10000 slope=2500/10000)
  delta_s = 0.0147;
  // Huntington (critical when length=2000 slope=60/2000)
  //delta_s = 0.038377;

  power = .35;
  factor = 1/power*pow((dmax-thermocline_depth)/(dshelf-thermocline_depth),power-1);

  if(z>-thermocline_depth) {
    //    printf("%e %e\n",z,0.09*fabs(z/dmax));
    return 0.09*fabs(z/dmax);
  } else {
    //    printf("%e %e\n",z,0.09*thermocline_depth/dmax+delta_s*pow(fabs((z+thermocline_depth)/dmax),0.35));
    return 0.09*thermocline_depth/dmax+delta_s*pow(fabs((z+thermocline_depth)/dmax),0.35);
  }
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

  REAL x0=7500, x1=5500, z0=-60, z1=-120, w=500, h=40,
    x2=7000, z2=-90;

  if(z<-depth+100)
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
  REAL u, v, umag=0;
  
  // Irrotational vortex.
  u = -umag*(y-5)/5;
  v = umag*(x-5)/5;
  return u*n1+v*n2;
}

