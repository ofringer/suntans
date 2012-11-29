/*
 * File: diffusion.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Compute eddy-viscosity via Smagorinsky or for central-differencing as dictated by
 * the lax-wendroff scheme.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "diffusion.h"
#include "util.h"

/*
 * Function: LaxWendroff
 * Usage: LaxWendroff(grid,phys,prop,myproc,comm);
 * -----------------------------------------------
 * Compute the numerical diffusion coefficients required to stablize momentum advection
 * when central-differencing is used (nonlinear=2).
 *
 */
void LaxWendroff(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int i, j, k, nc1, nc2;
  REAL numax=-INFTY, numin=INFTY, numinall, numaxall, laxVertical;

  laxVertical=1;
  if(grid->Nkmax==1)
    laxVertical=0;

  if(prop->nonlinear==2) {
    for(i=0;i<grid->Nc;i++) 
      for(k=0;k<grid->Nk[i];k++) {
	phys->nu_lax[i][k]=
	  0.5*(pow(laxVertical*0.5*(phys->wtmp2[i][k]+phys->wtmp2[i][k+1]),2)+
	       pow(phys->uc[i][k],2)+pow(phys->vc[i][k],2))*prop->dt;
	if(phys->nu_lax[i][k]<numin) numin=phys->nu_lax[i][k];
	if(phys->nu_lax[i][k]>numax) numax=phys->nu_lax[i][k];
      }

    if(VERBOSE>2) {
      MPI_Reduce(&numin,&numinall,1,MPI_DOUBLE,MPI_MAX,0,comm);
      MPI_Reduce(&numax,&numaxall,1,MPI_DOUBLE,MPI_MAX,0,comm);
      
      if(myproc==0) printf("Lax-Wendroff diffusion coefficients: numin = %.3e, numax = %.3e\n",numin,numax);
    }
  }
}

