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
  int i, k, icu, kcu, icw, kcw, Nc=grid->Nc, Ne=grid->Ne, ih, is, ks, iu, ku, iw, kw;
  int uflag=1, wflag=1, sflag=1, hflag=1, myalldone, alldone;
  REAL C, CmaxU, CmaxW, dtsuggestU, dtsuggestW;

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
    for(k=grid->etop[i]+1;k<grid->Nke[i];k++) {
      C = fabs(phys->u[i][k])*prop->dt/grid->dg[i];
      if(C>CmaxU) {
	icu = i;
	kcu = k;
	CmaxU = C;
      }
    }

  CmaxW=0;
  if(!prop->wetdry)
    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	C = 0.5*fabs(phys->w[i][k]+phys->w[i][k+1])*prop->dt/grid->dzz[i][k];
	if(C>CmaxW) {
	  icw = i;
	  kcw = k;
	  CmaxW = C;
	}
      }
  
  myalldone=0;
  if(!uflag || !wflag || !sflag || !hflag || CmaxU>prop->Cmax || CmaxW>prop->Cmax) {
    printf(DASHES);
    printf("Time step %d: Processor %d, Run is blowing up!\n",prop->n,myproc);

    if(CmaxU>prop->Cmax) {
      dtsuggestU = CMAXSUGGEST*grid->dg[icu]/fabs(phys->u[icu][kcu]);

      printf("Horizontal Courant number problems:\n");
      printf("  Grid indices: j=%d k=%d\n", icu, kcu);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xe[icu],grid->ye[icu],
	     0.5*(DepthFromDZ(grid,phys,grid->grad[2*icu],kcu)+
		  DepthFromDZ(grid,phys,grid->grad[2*icu+1],kcu)));
      printf("  Umax = %.3e\n",phys->u[icu][kcu]);
      printf("  Horizontal grid spacing grid->dg[%d] = %.3e\n",icu,grid->dg[icu]);
      printf("  Horizontal Courant number is CmaxU = %.2f.\n",CmaxU);
      printf("  You specified a maximum of %.2f in suntans.dat\n",prop->Cmax);
      printf("  Your time step size is %.2f.\n",prop->dt);
      printf("  Reducing it to at most %.2f (For C=0.5) might solve this problem.\n",dtsuggestU);
    }
    
    if(!uflag) {
	printf("Problem with U (U=NaN):\n");
	printf("  Grid indices: j=%d k=%d\n", iu, ku);
	printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xe[iu],grid->ye[iu],
	       0.5*(DepthFromDZ(grid,phys,grid->grad[2*iu],ku)+
		    DepthFromDZ(grid,phys,grid->grad[2*iu+1],ku)));	
    }

    if(CmaxW>prop->Cmax) {
      dtsuggestW = CMAXSUGGEST*grid->dzz[icw][kcw]/fabs(0.5*(phys->w[icw][kcw]+phys->w[icw][kcw]));

      printf("Vertical Courant number problems:\n");
      printf("  Grid indices: i=%d k=%d\n", icw, kcw);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[icw],grid->yv[icw],DepthFromDZ(grid,phys,icw,kcw));
      printf("  Wmax = %.3e (located half-way between faces)\n",0.5*(phys->w[icw][kcw]+phys->w[icw][kcw+1]));
      printf("  Vertical grid spacing dz = %.3e\n",icu,grid->dzz[icw][kcw]);
      printf("  Vertical Courant number is CmaxW = %.2f.\n",CmaxW);
      printf("  You specified a maximum of %.2f in suntans.dat\n",prop->Cmax);
      printf("  Your time step size is %.2f.\n",prop->dt);
      printf("  Reducing it to at most %.2f (For C=0.5) might solve this problem.\n",dtsuggestW);
    }

    if(!wflag) {
      printf("Problem with W (W=NaN):\n");
      printf("  Grid indices: i=%d k=%d\n", iw, kw);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[iw],grid->yv[iw],DepthFromDZ(grid,phys,iw,kw));
    }

    if(!sflag) {
      printf("Problem with the scalar s (s=NaN):\n");
      printf("  Grid indices: i=%d k=%d\n", is, ks);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[is],grid->yv[is],DepthFromDZ(grid,phys,is,ks));
    }

    if(!hflag) {
      printf("Problem with the free surface (h=NaN):\n");
      printf("  Grid indices: i=%d k=%d\n", ih);
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
void Progress(propT *prop, int myproc) 
{
  int progout, prog;
  char filename[BUFFERLENGTH];
  FILE *fid;
  
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
      if(!(prop->n%progout))
	printf("%d%% Complete.\n",prog);
  }
}


