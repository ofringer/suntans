/*
 * File: initialization.c
 * Description:  This file contains the functions that are used
 * to initialize the depth, free-surface, and salinity.
 *
 * $Id: initialization.c,v 1.4 2004-05-31 06:59:33 fringer Exp $
 * $Log: not supported by cvs2svn $
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

  return 1-.5*exp(-pow(x/4,2));
  if(x<6)
    return .5;
  return 2;
  /*
  if((x>4.556 && x<6 && y<1.73*(x-4.556)) | (x>6 && y<2.5))
    return 5;
  return 10;
  */
  //  if(x>2 && x<3) return .25;
  //  return .5;
  //  return 10.0+0*8*exp(-pow((x-100)/40,2)-0*pow((y-500)/50,2));
  /*
  if(x>3 & x<6 & y>2)
    return 1;
  return 2;
  */
   xc = (REAL *)malloc(Nc*sizeof(REAL));
   yc = (REAL *)malloc(Nc*sizeof(REAL));

   xc[0]=0.300000;yc[0]=0.150000;
   xc[1]=0.300000;yc[1]=0.250000;
   xc[2]=0.300000;yc[2]=0.350000;
   xc[3]=0.500000;yc[3]=0.150000;
   xc[4]=0.500000;yc[4]=0.250000;
   xc[5]=0.500000;yc[5]=0.350000;
   xc[6]=0.700000;yc[6]=0.150000;
   xc[7]=0.700000;yc[7]=0.250000;
   xc[8]=0.700000;yc[8]=0.350000;
   xc[9]=0.400000;yc[9]=0.200000;
   xc[10]=0.400000;yc[10]=0.300000;
   xc[11]=0.600000;yc[11]=0.200000;
   xc[12]=0.600000;yc[12]=0.300000;
   /*
   for(nc=0;nc<Nc;nc++) 
     if(sqrt(pow(x-xc[nc],2)+pow(y-yc[nc],2))<R)
       return .4;
   return .5;
   */

   // Constant depth
   //  return 200;

   // Gaussian slope
   //  return 1000-750*exp(-pow((x-5000.0)/1500,2)-0*pow((y-50000)/20000,2));

   // Step
   /*
   if((x>350 && x<450 && y>100 && y<210) || (x>550 && x<650 && y>100 && y<210))
     return 100;
   else
     return 200;
   */

   // Deep grid with a Gaussian hill.
   //  return 1000 - 250*exp(-pow((x-2500)/100,2)
   //    			 -0*pow((y-7000)/10000,2));

   // Tanh slope
   //return 3000-2500*.5*(1+tanh((x-60000)/5000));
   //    - 250*exp(-pow(x-30000)/500,2);

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
   // Simple continental shelf
   length = 10000;
   xmid = 75000;
   shelfdepth = 500;
   if(x<=xmid-length/2)
     return 3000-0*exp(-pow((x-50000)/2000,2));
   else if(x>xmid-length/2 && x<=xmid+length/2 && length>0)
     return 3000-(3000-shelfdepth)*(x-xmid+length/2)/length;
   else
     return shelfdepth;//-450*(x-xmid-length/2)/(100000-xmid-length/2);

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
REAL ReturnFreeSurface(REAL x, REAL y, REAL d) {
   return -d+.25;-.5+.6*cos(PI*x/10);
   return 0;
   if(x>4 && x<6)
     return 0;
   return
     -1;
   return .02+.01*cos(PI*x/10);
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
  REAL delta_s, thermocline_depth=20, z0;
  REAL power, factor,dmax=120,dshelf=60;//dmax=3000,dshelf=500;//dmax=120,dshelf=60;//dmax = 3000;
  thermocline_depth=35;
  dmax = 120;
  dshelf = 60;

  //  thermocline_depth=180;
  //  dmax = 3000;

  if(z>-1)
    return -0.005*tanh((x-5-0*(y-0*1)/tan(PI/3))/.01);
  return .01;

  /*
  if(x>10)
    return -0.005*tanh((z+5)/.01);
  else
    return 0;
  */
  // Monterey (critical when length=10000 slope=2500/10000)
  delta_s = 0.0147;
  // Huntington (critical when length=2000 slope=60/2000)
  //delta_s = 0.038377;

  //  return -3.43e-8/.007*z;
  //return -delta_s*z/dmax;

  /*
  z0=z;
  if(z>-thermocline_depth)
    z=-thermocline_depth;
  */

  //  return -1.1594894e-4*z;  // Huntington critical stratification for beta=.002
  //  return -2.3189788e-7*z;  // Huntington critical stratification for beta=1

  /*
  if(x>450 && x<550 && z<-100)
    return delta_s;
  else 
    return -delta_s;
  */
  //  return -delta_s/2*tanh((z+625-400*cos(PI*x/1250))/50);
  //  return -delta_s/2*tanh((z+500+250*exp(-pow(x/1500,2)))/50);
  //  return -delta_s/2*tanh((x-500)/50);

  power = .35;
  factor = 1/power*pow((dmax-thermocline_depth)/(dshelf-thermocline_depth),power-1);

  if(z>-thermocline_depth) {
    //    printf("%f %f\n",z,0.09*fabs(z/dmax));
    //return delta_s*pow(fabs(thermocline_depth/dmax),0.35);
    return 0.09*fabs(z/dmax);
  } else {
    //    printf("%f %f\n",z,delta_s*pow(fabs(z/dmax),0.35));
    return delta_s*pow(fabs(z/dmax),0.35);
  }
    //    return factor*delta_s*pow(fabs((z+thermocline_depth)/(dmax-thermocline_depth)),power);
  //phys->s[i][k]=-1.07*0*(z+1500+0*cos(PI*(grid->xv[i]-10000)/10000))/1500;
  //phys->s[i][k]=-7.333e-3*(.5*(z+1500)/1500);
  //phys->s0[i][k]=-7.333e-3*(.5*(z+1500)/1500);
  //if(grid->xv[i]>70000) phys->s[i][k]=.01;
  //if(grid->xv[i]>500 && grid->dzz[i][k]>0)
  //if(grid->xv[i]>5) phys->s[i][k]=0.01;
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
  
  if(x>5)
    return 1;
  return 0;

  if(z<-depth+100)
    return 1;
  return 0;

  z0 = -depth;
  /*
  if(z<z0+h/2) return 1;
  return 0;
  */

  if((z<-53.4+h/2 && x>x0 && x<x0+w))
     //     (z<-264.0+h/2 && x>3200 && x<3700) || 
     //     (z<-104.0+h/2 && x>5350 && x<5850)) 
    return 1;
  return 0;

  if((z<z0+h/2 && z>z0-h/2 && x>x0 && x<x0+w) || 
     (z<z0+h/2 && z>z0-h/2 && x>x0-1000 && x<x0+w-1000) || 
     (z<z0+h/2 && z>z0-h/2 && x>x0-2000 && x<x0+w-2000) || 
     (z<z0+h/2 && z>z0-h/2 && x>x0-3000 && x<x0+w-3000))
    return 1;
  //  else if(z<z1+h/2 && z>z1-h/2 && x>x1-w && x<x1)
  //    return 1;
  return 0;


  if(depth>400 && depth<600)
    return 1;
  else 
    return 0;

  if(x>50)
    return 1;
  //  if(x>6500 && x<7500)
  //    return 1;

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

