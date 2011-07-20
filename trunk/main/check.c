/*
 * File: check.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains functions that check for instability or
 * wetting and drying problems.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "check.h"
#include "suntans.h"
#include "stdio.h"
#include "timer.h"
#include "memory.h"

#define DASHES "----------------------------------------------------------------------\n"
#define CMAXSUGGEST 0.5

/*
 * Function: Check
 * Usage: Check(grid,phys,prop,myproc,numprocs,comm);
 * --------------------------------------------------
 * Check to make sure the run isn't blowing up.
 *
 */    
int Check(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, k, icu, kcu, icw, kcw, Nc=grid->Nc, Ne=grid->Ne, ih, is, ks, iu, ku, iw, kw, nc1, nc2;
  int uflag=1, wflag=1, sflag=1, hflag=1, myalldone, alldone, progout;
  REAL C, CmaxU, CmaxW, allCmaxU, allCmaxW, dtsuggestU, dtsuggestW;

  icu=kcu=icw=kcw=ih=is=ks=iu=ku=iw=kw=0;

  for(i=0;i<Nc;i++) 
    if(phys->h[i]!=phys->h[i]) {
	hflag=0;
	ih=i;
	break;
    }

  for(i=0;i<Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      if(phys->s[i][k]!=phys->s[i][k]) {
	sflag=0;
	is=i;
	ks=k;
	break;
      }
    if(!sflag)
      break;
  }

  for(i=0;i<Ne;i++) {
    for(k=0;k<grid->Nke[i];k++)
      if(phys->u[i][k]!=phys->u[i][k]) {
	uflag=0;
	iu=i;
	ku=k;
	break;
      }
    if(!uflag)
      break;
  }

  for(i=0;i<Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      if(phys->w[i][k]!=phys->w[i][k]) {
	wflag=0;
	iw=i;
	kw=k;
	break;
      }
    if(!wflag)
      break;
  }

  CmaxU=0;
  for(i=0;i<Ne;i++) 
    for(k=grid->etop[i];k<grid->Nke[i];k++) {
      C = fabs(phys->u[i][k])*prop->dt/grid->dg[i];
      if(C>CmaxU) {
	icu = i;
	kcu = k;
	CmaxU = C;
      }
    }

  CmaxW=0;
  if(prop->nonlinear!=0) {
    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	C = 0.5*fabs(phys->w[i][k]+phys->w[i][k+1])*prop->dt/grid->dzz[i][k];
	if(C>CmaxW) {
	  icw = i;
	  kcw = k;
	  CmaxW = C;
	}
      }
  }

  progout = (int)(prop->nsteps*(double)prop->ntprog/100);
  if(progout>0 && !(prop->n%progout)) {
    MPI_Reduce(&CmaxU,&allCmaxU,1,MPI_DOUBLE,MPI_MAX,0,comm);
    MPI_Reduce(&CmaxW,&allCmaxW,1,MPI_DOUBLE,MPI_MAX,0,comm);
  }

  if(myproc==0) {
    prop->CmaxU = allCmaxU;
    prop->CmaxW = allCmaxW;
  }

  myalldone=0;
  if(!uflag || !wflag || !sflag || !hflag || CmaxU>prop->Cmax || CmaxW>prop->Cmax) {
    printf(DASHES);
    printf("Time step %d: Processor %d, Run is blowing up!\n",prop->n,myproc);

    if(CmaxU>prop->Cmax) {
      nc1 = grid->grad[2*icu];
      nc2 = grid->grad[2*icu+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      dtsuggestU = CMAXSUGGEST*grid->dg[icu]/fabs(phys->u[icu][kcu]);

      printf("Horizontal Courant number problems:\n");
      printf("  Grid indices: j=%d k=%d (Nke=%d)\n", icu, kcu, grid->Nke[icu]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xe[icu],grid->ye[icu],
	     0.5*(DepthFromDZ(grid,phys,grid->grad[2*icu],kcu)+
		  DepthFromDZ(grid,phys,grid->grad[2*icu+1],kcu)));
      printf("  Free-surface heights (on either side): %.3e, %.3e\n",phys->h[nc1],phys->h[nc2]);
      printf("  Depths (on either side): %.3e, %.3e\n",grid->dv[nc1],grid->dv[nc2]);
      printf("  Flux-face height = %.3e, CdB = %.3e\n",grid->dzf[icu][kcu],phys->CdB[icu]);
      printf("  Umax = %.3e\n",phys->u[icu][kcu]);
      printf("  Horizontal grid spacing grid->dg[%d] = %.3e\n",icu,grid->dg[icu]);
      printf("  Horizontal Courant number is CmaxU = %.2f.\n",CmaxU);
      printf("  You specified a maximum of %.2f in suntans.dat\n",prop->Cmax);
      printf("  Your time step size is %.2f.\n",prop->dt);
      printf("  Reducing it to at most %.2f (For C=0.5) might solve this problem.\n",dtsuggestU);
    }
    
    if(!uflag) {
	printf("Problem with U (U=NaN):\n");
	printf("  Grid indices: j=%d k=%d (Nke=%d)\n", iu, ku, grid->Nke[iu]);
	printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xe[iu],grid->ye[iu],
	       0.5*(DepthFromDZ(grid,phys,grid->grad[2*iu],ku)+
		    DepthFromDZ(grid,phys,grid->grad[2*iu+1],ku)));	
    }

    if(CmaxW>prop->Cmax) {
      dtsuggestW = CMAXSUGGEST*grid->dzz[icw][kcw]/fabs(0.5*(phys->w[icw][kcw]+phys->w[icw][kcw]));

      printf("Vertical Courant number problems:\n");
      printf("  Grid indices: i=%d k=%d (Nkc=%d)\n", icw, kcw, grid->Nkc[icw]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[icw],grid->yv[icw],DepthFromDZ(grid,phys,icw,kcw));
      printf("  Free-surface height: %.3e\n",phys->h[icw]);
      printf("  Depth: %.3e\n",grid->dv[icw]);
      printf("  Wmax = %.3e (located half-way between faces)\n",0.5*(phys->w[icw][kcw]+phys->w[icw][kcw+1]));
      printf("  Vertical grid spacing dz = %.3e\n",grid->dzz[icw][kcw]);
      printf("  Vertical Courant number is CmaxW = %.2f.\n",CmaxW);
      printf("  You specified a maximum of %.2f in suntans.dat\n",prop->Cmax);
      printf("  Your time step size is %.2f.\n",prop->dt);
      printf("  Reducing it to at most %.2f (For C=0.5) might solve this problem.\n",dtsuggestW);
    }

    if(!wflag) {
      printf("Problem with W (W=NaN):\n");
      printf("  Grid indices: i=%d k=%d (Nkc=%d)\n", iw, kw, grid->Nkc[iw]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[iw],grid->yv[iw],DepthFromDZ(grid,phys,iw,kw));
    }

    if(!sflag) {
      printf("Problem with the scalar s (s=NaN):\n");
      printf("  Grid indices: i=%d k=%d (Nkc=%d)\n", is, ks, grid->Nkc[is]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[is],grid->yv[is],DepthFromDZ(grid,phys,is,ks));
    }

    if(!hflag) {
      printf("Problem with the free surface (h=NaN):\n");
      printf("  Grid index: i=%d\n", ih);
      printf("  Location: x=%.3e, y=%.3e\n",grid->xv[ih],grid->yv[ih]);
    }
    printf(DASHES);

    myalldone=1;
  }

  MPI_Reduce(&myalldone,&alldone,1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&alldone,1,MPI_INT,0,comm);

  return alldone;
}

/*
 * Function: CheckDZ
 * Usage: CheckDZ(grid,phys,prop,myproc,numprocs,comm);
 * ----------------------------------------------------
 * Check to make sure the vertical grid spacing is >= 0 and that the free surface is
 * not crossing through cells when wetdry != 0.
 *
 */    
int CheckDZ(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, k, iz, kz, iw, kw, Nc=grid->Nc, Ne=grid->Ne;
  int zflag=1, wflag=1, myalldone, alldone;

  iw=kw=iz=kz=0;

  for(i=0;i<Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      if(grid->dzz[i][k]<=0) {
	zflag=0;
	iz=i;
	kz=k;
	break;
      }
    if(!zflag)
      break;
  }

  if(!prop->wetdry)
    for(i=0;i<Nc;i++) {
      if(grid->ctop[i]!=grid->ctopold[i]) {
	wflag=0;
	iw=i;
	kw=k;
	break;
      }
  }

  myalldone=0;
  if(!zflag || !wflag) {
    printf(DASHES);
    printf("Time step %d: Processor %d, Wetting and drying problems!\n",prop->n,myproc);
    
    if(!zflag) {
      printf("Problems with the vertical grid spacing:\n");
      printf("  Grid indices: j=%d k=%d\n", iz, kz);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[iz],grid->yv[iz],
	     0.5*(DepthFromDZ(grid,phys,grid->grad[2*iz],kz)+
		  DepthFromDZ(grid,phys,grid->grad[2*iz+1],kz)));
      printf("  Vertical grid spacing = %.3e <= 0.\n",grid->dzz[iz][kz]);
      printf("This problem may be occuring because the horizontal Courant number\n");
      printf("is exceeding the maximum allowable (Cmax in suntans.dat).\n");
    }

    if(!wflag) {
      printf("Cells are wetting and drying although this is not allowed.\n");
      printf("Wetting and drying is only allowed if wetdry=1 in suntans.dat.\n");
      printf("  Grid indices: j=%d k=%d\n", iw, kw);
      printf("  Location: x=%.3e, y=%.3e\n",grid->xv[iw],grid->yv[iw]);
    }
    printf(DASHES);

    myalldone=1;
  }

  MPI_Reduce(&myalldone,&alldone,1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&alldone,1,MPI_INT,0,comm);

  return alldone;
}

/*
 * Function: Progress
 * Usage: Progress(prop,myproc);
 * -----------------------------
 * Output the progress of the calculation to the terminal.
 *
 */
void Progress(propT *prop, int myproc, int numprocs) 
{
  int progout, prog;
  char filename[BUFFERLENGTH];
  FILE *fid;
  REAL timeperstep = (Timer()-t_start)/(prop->n-prop->nstart);
  REAL t_sim, t_rem;
  
  MPI_GetFile(filename,DATAFILE,"ProgressFile","Progress",myproc);

  if(myproc==0) {
    fid = fopen(filename,"w");
    fprintf(fid,"On %d of %d, t=%.2f (%d%% Complete, %d output)",
	    prop->n,prop->nstart+prop->nsteps,prop->rtime,100*(prop->n-prop->nstart)/prop->nsteps,
	    1+(prop->n-prop->nstart)/prop->ntout);      
    fclose(fid);
  }
  
  if(myproc==0 && prop->ntprog>0 && VERBOSE>0) {
    progout = (int)(prop->nsteps*(double)prop->ntprog/100);
    prog=(int)(100.0*(double)(prop->n-prop->nstart)/(double)prop->nsteps);
    if(progout>0)
      if(!(prop->n%progout)) {
	if(prop->nonlinear) {
	  printf("%d%% Complete. CmaxU=%.2e, CmaxW=%.2e, %.2e s/step; %.2f s remaining.\n",
		 prog,prop->CmaxU,prop->CmaxW,timeperstep,timeperstep*(prop->nsteps+prop->nstart-prop->n));
	} else {
	  printf("%d%% Complete. CmaxU=%.2e, %.2e s/step; %.2f s remaining.\n",
		 prog,prop->CmaxU,timeperstep,timeperstep*(prop->nsteps+prop->nstart-prop->n));	  
	}
      }
    if(prop->n==prop->nsteps+prop->nstart) {
      t_sim = Timer()-t_start;
      t_rem = t_sim-t_nonhydro-t_predictor-t_source-t_transport-t_turb-t_io-t_check;

      printf("Total simulation time: %.2f s\n",t_sim);
      printf("Average per time step: %.2e s\n",t_sim/prop->nsteps);
      printf("Timing Summary:\n");
      if(prop->nonhydrostatic)
	printf("  Nonhydrostatic pressure: %.2f s (%.2f%)\n",t_nonhydro,
	       100*t_nonhydro/t_sim);
      printf("  Free surface and vertical friction: %.2f s (%.2f%)\n",t_predictor,
	     100*t_predictor/t_sim);
      printf("  Explicit terms: %.2f s (%.2f%)\n",t_source,
	     100*t_source/t_sim);
      if(prop->beta || prop->gamma)
	printf("  Scalar transport: %.2f s (%.2f%)\n",t_transport,
	       100*t_transport/t_sim);
      if(prop->turbmodel)
	printf("  Turbulence: %.2f s (%.2f%)\n", t_turb,
	       100*t_turb/t_sim);
      printf("  Bounds checking: %.2f s (%.2f%)\n", t_check,
	     100*t_check/t_sim);
      printf("  I/O: %.2f s (%.2f%)\n", t_io,
	     100*t_io/t_sim);
      printf("  Remainder: %.2f s (%.2f%)\n", t_rem,
	     100*t_rem/t_sim);
      if(numprocs>1) {
	printf("  Communication time: %.2e s\n",t_comm);
	printf("  Computation/Communication: %.2e\n",(t_sim-t_comm)/t_comm);
      }
    }
  }
}

/*
 * Function: MemoryStats
 * Usage: MemoryStats(myproc,numprocs,comm);
 * -----------------------------------------
 * Print out statistics on total memory and grid points.
 *
 */
void MemoryStats(gridT *grid, int myproc, int numprocs, MPI_Comm comm) {
  int i, ncells, allncells;
  unsigned AllSpace, TotSpacekb = TotSpace>>10;
  
  ncells=0;
  for(i=0;i<grid->Nc;i++)
    ncells+=grid->Nk[i];

  MPI_Reduce(&TotSpacekb,&(AllSpace),1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&AllSpace,1,MPI_INT,0,comm);
  MPI_Reduce(&ncells,&(allncells),1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&allncells,1,MPI_INT,0,comm);

  if(numprocs>0)
    printf("Processor %d,  Total memory: %u Mb, %d cells\n",
	   myproc,TotSpacekb>>10,ncells);
  if(myproc==0) 
    printf("All processors: %u Mb, %d cells (%d bytes/cell)\n",
	   AllSpace>>10,allncells,
	   (int)(1024.0*(REAL)AllSpace/(REAL)allncells));
}

