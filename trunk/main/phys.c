/*
 * File: phys.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains physically-based functions.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "util.h"
#include "initialization.h"
#include "memory.h"
#include "turbulence.h"
#include "boundaries.h"
#include "check.h"
#include "scalars.h"
#include "timer.h"
#include "profiles.h"
#include "state.h"

/*
 * Private Function declarations.
 *
 */
static void UpdateDZ(gridT *grid, physT *phys, int option);
static void UPredictor(gridT *grid, physT *phys, 
		       propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void Corrector(REAL **qc, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs);
static void CGSolve(gridT *grid, physT *phys, propT *prop, 
		    int myproc, int numprocs, MPI_Comm comm);
static void CGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys, propT *prop, 
		     int myproc, int numprocs, MPI_Comm comm);
static void ConditionQ(REAL **x, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
static void Preconditioner(REAL **x, REAL **xc, REAL **coef, gridT *grid, physT *phys, propT *prop);
static void GuessQ(REAL **q, REAL **wold, REAL **w, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void GSSolve(gridT *grid, physT *phys, propT *prop, 
		    int myproc, int numprocs, MPI_Comm comm);
static REAL InnerProduct(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static REAL InnerProduct3(REAL **x, REAL **y, gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static void OperatorH(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop);
static void OperatorQC(REAL **coef, REAL **fcoef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop);
static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, physT *phys, propT *prop);
static void OperatorQ(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop);
static void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop);
static void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs,
			  MPI_Comm comm);
static void EddyViscosity(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);
static void HorizontalSource(gridT *grid, physT *phys, propT *prop,
			     int myproc, int numprocs, MPI_Comm comm);
static void StoreVariables(gridT *grid, physT *phys);
static void NewCells(gridT *grid, physT *phys, propT *prop);
static void WPredictor(gridT *grid, physT *phys, propT *prop,
		       int myproc, int numprocs, MPI_Comm comm);
static void ComputeVelocityVector(REAL **u, REAL **uc, REAL **vc, gridT *grid);
static void OutputData(gridT *grid, physT *phys, propT *prop,
		int myproc, int numprocs, int blowup, MPI_Comm comm);
static REAL InterpToFace(int j, int k, REAL **phi, REAL **u, gridT *grid);
static void SetDensity(gridT *grid, physT *phys, propT *prop);

/*
 * Function: AllocatePhysicalVariables
 * Usage: AllocatePhysicalVariables(grid,phys,prop);
 * -------------------------------------------------
 * This function allocates space for the physical arrays but does not
 * allocate space for the grid as this has already been allocated.
 *
 */
void AllocatePhysicalVariables(gridT *grid, physT **phys, propT *prop)
{
  int flag=0, i, j, jptr, ib, Nc=grid->Nc, Ne=grid->Ne, nf;

  *phys = (physT *)SunMalloc(sizeof(physT),"AllocatePhysicalVariables");

  (*phys)->u = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->uc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->uold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->D = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->utmp = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->utmp2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->ut = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_U = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");

  (*phys)->wf = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");

  for(j=0;j<Ne;j++) {
    if(grid->Nkc[j]<grid->Nke[j]) {
      printf("Error!  Nkc(=%d)<Nke(=%d) at edge %d\n",grid->Nkc[j],grid->Nke[j],j);
      flag = 1;
    }
    (*phys)->u[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->ut[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_U[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wf[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
  }
  if(flag) {
    MPI_Finalize();
    exit(0);
  }

  (*phys)->h = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->hold = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  
  (*phys)->w = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wtmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wtmp2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_W = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->q = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->qc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->qtmp = (REAL **)SunMalloc(NFACES*Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s0 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->rho = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_R = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp3 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->nu_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->kappa_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  if(prop->turbmodel) {
    (*phys)->qT = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->lT = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->Cn_q = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->Cn_l = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  }
  (*phys)->tau_T = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->tau_B = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdT = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdB = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  
  for(i=0;i<Nc;i++) {
    (*phys)->uc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->uold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp2[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_W[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->q[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->qc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    for(nf=0;nf<NFACES;nf++)
      (*phys)->qtmp[i*NFACES+nf] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s0[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->rho[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_R[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    if(prop->turbmodel) {
      (*phys)->Cn_q[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->Cn_l[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->qT[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->lT[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    }
    (*phys)->stmp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp2[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp3[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->nu_tv[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->kappa_tv[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  (*phys)->boundary_u = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_v = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_w = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_s = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_T = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_rho = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_tmp = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_h = (REAL *)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->boundary_flag = (REAL *)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL),"AllocatePhysicalVariables");
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
    j=grid->edgep[jptr];

    (*phys)->boundary_u[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_v[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_w[jptr-grid->edgedist[2]] = (REAL *)SunMalloc((grid->Nke[j]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_s[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_T[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_tmp[jptr-grid->edgedist[2]] = (REAL *)SunMalloc((grid->Nke[j]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_rho[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  (*phys)->ap = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->am = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->bp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->bm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->a = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->b = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->c = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->d = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");

  // Allocate for the face scalar
  (*phys)->SfHp = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->SfHm = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  for(j=0;j<Ne;j++) {
    (*phys)->SfHp[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->SfHm[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  // Allocate for TVD schemes
  (*phys)->Cp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->Cm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->rp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->rm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");

  (*phys)->wp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->wm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");

  (*phys)->gradSx = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->gradSy = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  for(i=0;i<Nc;i++) {
    (*phys)->gradSx[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->gradSy[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }
}

/*
 * Function: FreePhysicalVariables
 * Usage: FreePhysicalVariables(grid,phys,prop);
 * ---------------------------------------------
 * This function frees all space allocated in AllocatePhysicalVariables
 *
 */
void FreePhysicalVariables(gridT *grid, physT *phys, propT *prop)
{
  int i, j, Nc=grid->Nc, Ne=grid->Ne, nf;
  
  for(j=0;j<Ne;j++) {
    free(phys->u[j]);
    free(phys->utmp[j]);
    free(phys->utmp2[j]);
    free(phys->ut[j]);
    free(phys->Cn_U[j]);
    free(phys->wf[j]);
  }

  for(i=0;i<Nc;i++) {
    free(phys->uc[i]);
    free(phys->vc[i]);
    free(phys->uold[i]);
    free(phys->vold[i]);
    free(phys->w[i]);
    free(phys->wtmp[i]);
    free(phys->wtmp2[i]);
    free(phys->Cn_W[i]);
    free(phys->q[i]);
    free(phys->qc[i]);
    for(nf=0;nf<NFACES;nf++)
      free(phys->qtmp[i*NFACES+nf]);
    free(phys->s[i]);
    free(phys->T[i]);
    free(phys->s0[i]);
    free(phys->rho[i]);
    free(phys->Cn_R[i]);
    free(phys->Cn_T[i]);
    if(prop->turbmodel) {
      free(phys->Cn_q[i]);
      free(phys->Cn_l[i]);
      free(phys->qT[i]);
      free(phys->lT[i]);
    }
    free(phys->stmp[i]);
    free(phys->stmp2[i]);
    free(phys->stmp3[i]);
    free(phys->nu_tv[i]);
    free(phys->kappa_tv[i]);
  }

  free(phys->h);
  free(phys->htmp);
  free(phys->uc);
  free(phys->vc);
  free(phys->w);
  free(phys->wtmp);
  free(phys->wtmp2);
  free(phys->Cn_W);
  free(phys->wf);
  free(phys->q);
  free(phys->qtmp);
  free(phys->s);
  free(phys->T);
  free(phys->s0);
  free(phys->rho);
  free(phys->Cn_R);
  free(phys->Cn_T);
  if(prop->turbmodel) {
    free(phys->Cn_q);
    free(phys->Cn_l);
    free(phys->qT);
    free(phys->lT);
  }  
  free(phys->stmp);
  free(phys->stmp2);
  free(phys->stmp3);
  free(phys->nu_tv);
  free(phys->kappa_tv);
  free(phys->tau_T);
  free(phys->tau_B);
  free(phys->CdT);
  free(phys->CdB);
  free(phys->u);
  free(phys->D);
  free(phys->utmp);
  free(phys->ut);
  free(phys->Cn_U);

  free(phys->ap);
  free(phys->am);
  free(phys->bp);
  free(phys->bm);
  free(phys->a);
  free(phys->b);
  free(phys->c);
  free(phys->d);

  // Free the horizontal facial scalar  
  for(j=0;j<Ne;j++) {
    free( phys->SfHp[j] );
    free( phys->SfHm[j] );
  }
  free(phys->SfHp);
  free(phys->SfHm);

  // Free the variables for TVD scheme
  free(phys->Cp);
  free(phys->Cm);
  free(phys->rp);
  free(phys->rm);
  free(phys->wp);
  free(phys->wm);

  free(phys->gradSx);
  free(phys->gradSy);

  free(phys);
}
    
/*
 * Function: ReadPhysicalVariables
 * Usage: ReadPhysicalVariables(grid,phys,prop,myproc);
 * ----------------------------------------------------
 * This function reads in physical variables for a restart run
 * from the restart file defined by prop->StartFID.
 *
 */
void ReadPhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc) {

  int i, j;

  if(VERBOSE>1 && myproc==0) printf("Reading from rstore...\n");
    
  fread(&(prop->nstart),sizeof(int),1,prop->StartFID);

  fread(phys->h,sizeof(REAL),grid->Nc,prop->StartFID);
  for(j=0;j<grid->Ne;j++) 
    fread(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StartFID);

  if(prop->turbmodel) {
    for(i=0;i<grid->Nc;i++) 
      fread(phys->Cn_q[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
    for(i=0;i<grid->Nc;i++) 
      fread(phys->Cn_l[i],sizeof(REAL),grid->Nk[i],prop->StartFID);

    for(i=0;i<grid->Nc;i++) 
      fread(phys->qT[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
    for(i=0;i<grid->Nc;i++) 
      fread(phys->lT[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  }
  for(i=0;i<grid->Nc;i++) 
    fread(phys->nu_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->kappa_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID);

  for(j=0;j<grid->Ne;j++) 
    fread(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StartFID);

  for(i=0;i<grid->Nc;i++) 
    fread(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  fclose(prop->StartFID);

  UpdateDZ(grid,phys,0);
  ComputeVelocityVector(phys->u,phys->uc,phys->vc,grid);
}

/*
 * Function: InitializePhyiscalVariables
 * Usage: InitializePhyiscalVariables(grid,phys,prop);
 * ---------------------------------------------------
 * This function initializes the physical variables by calling
 * the routines defined in the file initialize.c
 *
 */
void InitializePhysicalVariables(gridT *grid, physT *phys, propT *prop)
{
  int i, j, k, Nc=grid->Nc;
  REAL z, *stmp;

  prop->nstart=0;

  // Initialize the free surface
  for(i=0;i<Nc;i++) {
    phys->h[i]=ReturnFreeSurface(grid->xv[i],grid->yv[i],grid->dv[i]);
    if(phys->h[i]<-grid->dv[i])
      phys->h[i]=-grid->dv[i] + 1e-10*grid->dz[grid->Nk[i]-1];
  }

  // Need to update the vertical grid after updating the free surface.
  // The 1 indicates that this is the first call to UpdateDZ
  UpdateDZ(grid,phys,1);

  for(i=0;i<Nc;i++) {
    phys->w[i][grid->Nk[i]]=0;
    for(k=0;k<grid->Nk[i];k++) {
      phys->w[i][k]=0;
      phys->q[i][k]=0;
      phys->s[i][k]=0;
      phys->T[i][k]=0;
      phys->s0[i][k]=0;
    }
  }

  // Initialize the temperature, salinity, and background salinity
  // distributions.  Since z is not stored, need to use dz[k] to get
  // z[k].
  if(prop->readSalinity) {
    stmp = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializePhysicalVariables");
    fread(stmp,sizeof(REAL),grid->Nkmax,prop->InitSalinityFID);
    fclose(prop->InitSalinityFID);

    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	phys->s[i][k]=stmp[k];
	phys->s0[i][k]=stmp[k];
      }
    SunFree(stmp,grid->Nkmax,"InitializePhysicalVariables");
  } else {
    for(i=0;i<Nc;i++) {
      z = 0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	z-=grid->dz[k]/2;
	phys->s[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z);
	phys->s0[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z);
	z-=grid->dz[k]/2;
      }
    }
  }

  if(prop->readTemperature) {
    stmp = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializePhysicalVariables");
    fread(stmp,sizeof(REAL),grid->Nkmax,prop->InitTemperatureFID);
    fclose(prop->InitTemperatureFID);    

    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	phys->T[i][k]=stmp[k];

    SunFree(stmp,grid->Nkmax,"InitializePhysicalVariables");
  } else {
    for(i=0;i<Nc;i++) {
      z = 0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	z-=grid->dz[k]/2;
	phys->T[i][k]=ReturnTemperature(grid->xv[i],grid->yv[i],z,grid->dv[i]);
	z-=grid->dz[k]/2;
      }
    }
  }

  // Initialize the velocity field 
  for(j=0;j<grid->Ne;j++) {
    z = 0;
    for(k=0;k<grid->Nkc[j];k++) {
      z-=grid->dz[k]/2;
      phys->u[j][k]=ReturnHorizontalVelocity(grid->xe[j],grid->ye[j],grid->n1[j],grid->n2[j],z);
      z-=grid->dz[k]/2;
    }
    for(k=0;k<grid->Nke[j];k++)
      phys->wf[j][k]=0;
  }

  // Need to compute the velocity vectors at the cell centers based
  // on the initialized velocities at the faces.
  ComputeVelocityVector(phys->u,phys->uc,phys->vc,grid);
  ComputeVelocityVector(phys->u,phys->uold,phys->vold,grid);

  // Determine minimum and maximum salinity
  phys->smin=phys->s[0][0];
  phys->smax=phys->s[0][0];
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) {
      if(phys->s[i][k]<phys->smin) phys->smin=phys->s[i][k];      
      if(phys->s[i][k]>phys->smax) phys->smax=phys->s[i][k];      
    }

  // Set the density from s and T using the equation of state 
  SetDensity(grid,phys,prop);

  // Initialize the eddy-viscosity and scalar diffusivity
  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) {
      phys->nu_tv[i][k]=0;
      phys->kappa_tv[i][k]=0;
    }

  if(prop->turbmodel) 
    for(i=0;i<grid->Nc;i++) 
      for(k=0;k<grid->Nk[i];k++) {
	phys->qT[i][k]=0;
	phys->lT[i][k]=0;
      }
}

/*
 * Function: SetDragCoefficients
 * Usage: SetDragCoefficents(grid,phys,prop);
 * ------------------------------------------
 * Set the drag coefficients based on the log law as well as the applied shear stress.
 *
 */
void SetDragCoefficients(gridT *grid, physT *phys, propT *prop) {
  int i, j, jptr, nc1, nc2;

  if(prop->z0T==0) 
    for(j=0;j<grid->Ne;j++) 
      phys->CdT[j]=prop->CdT;
  else
    for(j=0;j<grid->Ne;j++) {
      nc1=grid->grad[2*j];
      nc2=grid->grad[2*j+1];
      if(nc1==-1) nc1=nc2; 
      if(nc2==-1) nc2=nc1; 
      if(grid->Nk[nc2]>grid->Nk[nc1]) nc1=nc2;
      if(grid->Nk[nc1]>grid->Nk[nc2]) nc2=nc1;
      
      phys->CdT[j]=pow(log(0.25*(grid->dzz[nc1][grid->ctop[nc1]]+grid->dzz[nc2][grid->ctop[nc2]])/prop->z0T)/KAPPA_VK,-2);
    }

  if(prop->z0B==0) 
    for(j=0;j<grid->Ne;j++) 
      phys->CdB[j]=prop->CdB;
  else
    for(j=0;j<grid->Ne;j++) {
      nc1=grid->grad[2*j];
      nc2=grid->grad[2*j+1];
      if(nc1==-1) nc1=nc2; 
      if(nc2==-1) nc2=nc1; 
      if(grid->Nk[nc2]>grid->Nk[nc1]) nc1=nc2;
      if(grid->Nk[nc1]>grid->Nk[nc2]) nc2=nc1;
      
      phys->CdB[j]=pow(log(0.25*(grid->dzz[nc1][grid->Nke[j]-1]+grid->dzz[nc2][grid->Nke[j]-1])/prop->z0B)/KAPPA_VK,-2);
    }
}

/*
 * Function: InitializeVerticalGrid
 * Usage: InitializeVerticalGrid(grid);
 * ------------------------------------
 * Initialize the vertical grid by allocating space for grid->dzz and grid->dzzold
 * This just sets dzz and dzzold to dz since upon initialization dzz and dzzold
 * do not vary in the horizontal.
 *
 */
void InitializeVerticalGrid(gridT **grid)
{
  int i, k, Nc=(*grid)->Nc;

  (*grid)->dzz = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->dzzold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");

  for(i=0;i<Nc;i++) {
    (*grid)->dzz[i]=(REAL *)SunMalloc(((*grid)->Nk[i])*sizeof(REAL),"InitializeVerticalGrid");
    (*grid)->dzzold[i]=(REAL *)SunMalloc(((*grid)->Nk[i])*sizeof(REAL),"InitializeVerticalGrid");
    for(k=0;k<(*grid)->Nk[i];k++) {
      (*grid)->dzz[i][k]=(*grid)->dz[k];  
      (*grid)->dzzold[i][k]=(*grid)->dz[k];  
    }
  }
}

/*
 * Function: UpdateDZ
 * Usage: UpdateDZ(grid,phys,0);
 * -----------------------------
 * This function updates the vertical grid spacings based on the free surface and
 * the bottom bathymetry.  That is, if the free surface cuts through cells and leaves
 * any cells dry (or wet), then this function will set the vertical grid spacing 
 * accordingly.
 *
 * If option==1, then it assumes this is the first call and sets dzzold to dzz at
 *   the end of this function.
 * Otherwise it sets dzzold to dzz at the beginning of the function and updates dzz
 *   thereafter.
 *
 */
static void UpdateDZ(gridT *grid, physT *phys, int option)
{
  int i, j, k, ne1, ne2, Nc=grid->Nc, Ne=grid->Ne, flag;
  REAL z;

  // If this is not an initial call then set dzzold to store the old value of dzz
  // and also set the etopold and ctopold pointers to store the top indices of
  // the grid.
  if(!option) {
    for(j=0;j<Ne;j++)
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->Nk[i];k++)
	grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }

  // First set the thickness of the bottom grid layer.  If this is a partial-step
  // grid then the dzz will vary over the horizontal at the bottom layer.  Otherwise,
  // the dzz at the bottom will be equal to dz at the bottom.
  for(i=0;i<Nc;i++) {
    z = 0;
    for(k=0;k<grid->Nk[i];k++)
      z-=grid->dz[k];
    grid->dzz[i][grid->Nk[i]-1]=grid->dz[grid->Nk[i]-1]+grid->dv[i]+z;
  }
  
  // Loop through and set the vertical grid thickness when the free surface cuts through 
  // a particular cell.
  if(grid->Nkmax>1) {
    for(i=0;i<Nc;i++) {
      z = 0;
      flag = 0;

      for(k=0;k<grid->Nk[i];k++) {
	z-=grid->dz[k];
	if(phys->h[i]>=z) 
	  if(!flag) {
	    if(k==grid->Nk[i]-1) {
	      grid->dzz[i][k]=phys->h[i]+grid->dv[i];
	      grid->ctop[i]=k;
	    } else if(phys->h[i]==z) {
	      grid->dzz[i][k]=0;
	      grid->ctop[i]=k+1;
	    } else {
	      grid->dzz[i][k]=phys->h[i]-z;
	      grid->ctop[i]=k;
	    }
	    flag=1;
	  } else {
	    if(k==grid->Nk[i]-1) 
	      grid->dzz[i][k]=grid->dz[k]+grid->dv[i]+z;
	    else 
	      if(z<-grid->dv[i])
		grid->dzz[i][k]=0;
	      else 
		grid->dzz[i][k]=grid->dz[k];
	  } 
	else 
	  grid->dzz[i][k]=0;
      }
    }
  } else 
    for(i=0;i<Nc;i++) 
      grid->dzz[i][0]=grid->dv[i]+phys->h[i];

  // Now set grid->etop and ctop which store the index of the top cell  
  for(j=0;j<grid->Ne;j++) {
    ne1 = grid->grad[2*j];
    ne2 = grid->grad[2*j+1];
    if(ne1 == -1)
      grid->etop[j]=grid->ctop[ne2];
    else if(ne2 == -1)
      grid->etop[j]=grid->ctop[ne1];
    else if(grid->ctop[ne1]<grid->ctop[ne2])
      grid->etop[j]=grid->ctop[ne1];
    else
      grid->etop[j]=grid->ctop[ne2];
  }

  // If this is an initial call set the old values to the new values.
  if(option) {
    for(j=0;j<Ne;j++) 
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->Nk[i];k++)
	grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }
}

/*
 * Function: DepthFromDZ
 * Usage: z = DepthFromDZ(grid,phys,i,k);
 * --------------------------------------
 * Return the depth beneath the free surface at location i, k.
 *
 */
REAL DepthFromDZ(gridT *grid, physT *phys, int i, int kind) {
  int k;
  REAL z = phys->h[i]-0.5*grid->dzz[i][grid->ctop[i]];
  for(k=grid->ctop[i];k<kind;k++) {
    z-=0.5*grid->dzz[i][k-1];
    z-=0.5*grid->dzz[i][k];
  }
  return z;
}

/*
 * Function: Solve
 * Usage: Solve(grid,phys,prop,myproc,numprocs,comm);
 * --------------------------------------------------
 * This is the main solving routine and is called from suntans.c.
 *
 */
void Solve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, k, n, blowup=0;
  REAL t0;

  // Compute the initial quantities for comparison to determine conservative properties
  prop->n=0;
  ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);

  // Print out memory usage per processor and total memory if this is the first time step
  if(VERBOSE>1) MemoryStats(grid,myproc,numprocs,comm);

  prop->theta0=prop->theta;

  t_start=Timer();
  t_source=t_predictor=t_nonhydro=t_turb=t_transport=t_io=t_comm=t_check=0;
  for(n=prop->nstart+1;n<=prop->nsteps+prop->nstart;n++) {
    prop->n = n;
    prop->rtime = (n-1)*prop->dt;
    if(prop->nsteps>0) {

      // Set boundary values
      BoundaryVelocities(grid,phys,prop);
      BoundaryScalars(grid,phys,prop);
      WindStress(grid,phys,prop);

      // Ramp down theta from 1 to the value specified in suntans.dat over
      // the time thetaramptime specified in suntans.dat to damp out transient
      // oscillations
      if(prop->thetaramptime!=0)
	prop->theta=(1-exp(-prop->rtime/prop->thetaramptime))*prop->theta0+
	  exp(-prop->rtime/prop->thetaramptime);

      // Store the old velocity and scalar fields
      StoreVariables(grid,phys);

      // Compute the horizontal source term phys->utmp which contains the explicit part
      // or the right hand side of the free-surface equation. 
      t0=Timer();
      HorizontalSource(grid,phys,prop,myproc,numprocs,comm);
      t_source+=Timer()-t0;

      // Use the explicit part created in HorizontalSource and solve for the free-surface
      // and hence compute the predicted or hydrostatic horizontal velocity field.  Then
      // send and receive the free surface interprocessor boundary data to the neighboring processors.
      // The predicted horizontal velocity is now in phys->u
      t0=Timer();
      UPredictor(grid,phys,prop,myproc,numprocs,comm);
      ISendRecvCellData2D(phys->h,grid,myproc,comm);

      // Adjust the velocity field in the new cells if the newcells variable is set to 1 in
      // suntans.dat.  Once this is done, send the interprocessor u-velocities to the neighboring
      // processors.
      if(prop->newcells) {
	NewCells(grid,phys,prop);
	ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
      }
      t_predictor+=Timer()-t0;
      
      t0=Timer();
      blowup = CheckDZ(grid,phys,prop,myproc,numprocs,comm);
      t_check+=Timer()-t0;

      // Compute vertical momentum and the nonhydrostatic pressure
      t0=Timer();
      if(prop->nonhydrostatic && !blowup) {

	// Predicted vertical velocity field is in phys->w
	WPredictor(grid,phys,prop,myproc,numprocs,comm);

	// Source term for the pressure-Poisson equation is in phys->stmp
	ComputeQSource(phys->stmp,grid,phys,prop,myproc,numprocs);

	// Solve for the nonhydrostatic pressure.  
	// phys->stmp2 contains the initial guess
	// phys->stmp contains the source term
	// phys->stmp3 is used for temporary storage
	CGSolveQ(phys->qc,phys->stmp,phys->stmp3,grid,phys,prop,myproc,numprocs,comm);

	// Correct the nonhydrostatic velocity field with the nonhydrostatic pressure
	// correction field phys->stmp2.  This will correct phys->u so that it is now
	// the volume-conserving horizontal velocity field.  phys->w is not corrected since
	// it is obtained via continuity.  Also, update the total nonhydrostatic pressure
	// with the pressure correction. 
	Corrector(phys->qc,grid,phys,prop,myproc,numprocs,comm);

	// Send/recv the horizontal velocity data after it has been corrected.
	ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
	// Send q to the boundary cells now that it has been updated
	ISendRecvCellData3D(phys->q,grid,myproc,comm);
      }	else {
	// Compute the vertical velocity based on continuity and then send/recv to
	// neighboring processors.
	Continuity(phys->w,grid,phys,prop);
      }
      // Send/recv the vertical velocity data 
      ISendRecvWData(phys->w,grid,myproc,comm);
      t_nonhydro+=Timer()-t0;

      // Compute the eddy viscosity
      t0=Timer();
      EddyViscosity(grid,phys,prop,comm,myproc);
      t_turb+=Timer()-t0;

      // Update the salinity only if beta is nonzero in suntans.dat
      if(prop->beta) {
	t0=Timer();
	UpdateScalars(grid,phys,prop,phys->s,phys->boundary_s,phys->Cn_R,prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->thetaS,
		      NULL,NULL,NULL,NULL,0,0,comm,myproc);
	ISendRecvCellData3D(phys->s,grid,myproc,comm);
	t_transport+=Timer()-t0;
      }

      // Update the temperature only if gamma is nonzero in suntans.dat
      if(prop->gamma) {
	t0=Timer();
	UpdateScalars(grid,phys,prop,phys->T,phys->boundary_T,phys->Cn_T,prop->kappa_T,prop->kappa_TH,phys->kappa_tv,prop->thetaS,
		      NULL,NULL,NULL,NULL,0,0,comm,myproc);
	ISendRecvCellData3D(phys->T,grid,myproc,comm);
	t_transport+=Timer()-t0;
      }

      if(prop->beta || prop->gamma)
	SetDensity(grid,phys,prop);

      // utmp2 contains the velocity field at time step n, u contains
      // it at time step n+1.  This is so that at the next time step
      // phys->uold contains velocity at time step n-1 and phys->uc contains
      // that at time step n.
      ComputeVelocityVector(phys->u,phys->uc,phys->vc,grid);
      ISendRecvCellData3D(phys->uc,grid,myproc,comm);
      ISendRecvCellData3D(phys->vc,grid,myproc,comm);
    }

    // Check whether or not run is blowing up
    t0=Timer();
    blowup=(Check(grid,phys,prop,myproc,numprocs,comm) || blowup);
    t_check+=Timer()-t0;
    // Output data based on ntout specified in suntans.dat
    t0=Timer();
    OutputData(grid,phys,prop,myproc,numprocs,blowup,comm);
    InterpData(grid,phys,prop,comm,numprocs,myproc);
    t_io+=Timer()-t0;
    // Output progress
    Progress(prop,myproc,numprocs);

    if(blowup)
      break;
  }
}

/*
 * Function: StoreVariables
 * Usage: StoreVariables(grid,phys);
 * ---------------------------------
 * Store the old values of s, u, and w into stmp3, utmp2, and wtmp2,
 * respectively.
 *
 */
static void StoreVariables(gridT *grid, physT *phys) {
  int i, j, k, iptr, jptr;

  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) {
      phys->stmp3[i][k]=phys->s[i][k];
      phys->wtmp2[i][k]=phys->w[i][k];
    }

  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=0;
    for(k=0;k<grid->Nke[j];k++)
      phys->utmp2[j][k]=phys->u[j][k];
  }
}

/*
 * Function: HorizontalSource
 * Usage: HorizontalSource(grid,phys,prop,myproc,numprocs);
 * --------------------------------------------------------
 * Compute the horizontal source term that is used to obtain the free surface.
 *
 * This function adds the following to the horizontal source term:
 *
 * 1) Old nonhydrostatic pressure gradient with theta method
 * 2) Coriolis terms with AB2
 * 3) Baroclinic term with AB2
 * 4) Horizontal and vertical advection of horizontal momentum with AB2
 * 5) Horizontal laminar+turbulent diffusion of horizontal momentum
 *
 * Cn_U contains the Adams-Bashforth terms at time step n-1.
 * If wetting and drying is employed, no advection is computed in 
 * the upper cell.
 *
 */
static void HorizontalSource(gridT *grid, physT *phys, propT *prop,
			     int myproc, int numprocs, MPI_Comm comm) {
  int i, iptr, nf, j, jptr, k, nc, nc1, nc2, ne, k0, kmin, kmax;
  REAL *a, *b, *c, fab, sum, def1, def2, dgf;

  a = phys->a;
  b = phys->b;
  c = phys->c;

  // fab is 1 for a forward Euler calculation on the first time step,
  // for which Cn_U is 0.  Otherwise, fab=3/2 and Cn_U contains the
  // Adams-Bashforth terms at time step n-1
  if(prop->n==1) {
    fab=1;
    for(j=0;j<grid->Ne;j++)
      for(k=0;k<grid->Nke[j];k++)
	phys->Cn_U[j][k]=0;
  } else
    fab=1.5;

  // Set utmp and ut to zero since utmp will store the source term of the
  // horizontal momentum equation
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nke[j];k++) {
      phys->utmp[j][k]=0;
      phys->ut[j][k]=0;
    }
  }

  // Sponge layer at x=0 that decays exponentially with distance sponge_distance
  // over a timescale given by sponge_decay.  Both defined in suntans.dat
  // First use Cn_U from step n-1 then set it to 0.
  if(prop->sponge_distance==0) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->utmp[j][k]=(1-fab)*phys->Cn_U[j][k]+phys->u[j][k]
	  -prop->dt/grid->dg[j]*(phys->q[nc1][k]-phys->q[nc2][k]);
	
	phys->Cn_U[j][k]=0;
      }
    }
    // Add on explicit term to boundary edges
    for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr]; 
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->utmp[j][k]=(1-fab)*phys->Cn_U[j][k]+phys->u[j][k];
	
	phys->Cn_U[j][k]=0;
      }
    }
  } else {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->utmp[j][k]=(1-fab)*phys->Cn_U[j][k]
	  -prop->dt/grid->dg[j]*(phys->q[nc1][k]-phys->q[nc2][k])
	  +(1.0-prop->dt*exp(-0.5*(grid->xv[nc1]+grid->xv[nc2])/prop->sponge_distance)/
	    prop->sponge_decay)*phys->u[j][k];
	
	phys->Cn_U[j][k]=0;
      }
    }
    for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->utmp[j][k]=(1-fab)*phys->Cn_U[j][k]
	  +(1.0-prop->dt*exp(-grid->xv[nc1]/prop->sponge_distance)/
	    prop->sponge_decay)*phys->u[j][k];
	
	phys->Cn_U[j][k]=0;
      }
    }
  }
  
  // Coriolis terms
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->Cn_U[j][k]+=prop->dt*prop->Coriolis_f*(InterpToFace(j,k,phys->vc,phys->u,grid)*grid->n1[j]-
						   InterpToFace(j,k,phys->uc,phys->u,grid)*grid->n2[j]);
  }
  
  // Baroclinic term
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(grid->etop[j]<grid->Nke[j]-1) 
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	k0=grid->etop[j];
	
	for(k0=grid->etop[j];k0<k;k0++)
	  phys->Cn_U[j][k]-=0.5*GRAV*prop->dt*(phys->rho[nc1][k0]-phys->rho[nc2][k0])*
	    (grid->dzz[nc1][k0]+grid->dzz[nc2][k0])/grid->dg[j];
	phys->Cn_U[j][k]-=0.25*GRAV*prop->dt*(phys->rho[nc1][k]-phys->rho[nc2][k])*
	  (grid->dzz[nc1][k]+grid->dzz[nc2][k])/grid->dg[j];
      }
  }
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      k0=grid->etop[j];
      
      for(k0=grid->etop[j];k0<k;k0++)
	phys->Cn_U[j][k]-=GRAV*prop->dt*(phys->rho[nc1][k0]-phys->boundary_rho[jptr-grid->edgedist[2]][k0])*
	  grid->dzz[nc1][k0]/grid->dg[j];
      phys->Cn_U[j][k]-=0.5*GRAV*prop->dt*(phys->rho[nc1][k]-phys->boundary_rho[jptr-grid->edgedist[2]][k])*
	grid->dzz[nc1][k]/grid->dg[j];
    }
  }

  // Set stmp and stmp2 to zero since these are used as temporary variables for advection and
  // diffusion.
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;
  
  // Compute Eulerian advection of momentum (nonlinear!=0)
  if(prop->nonlinear) {
    
    // U-fluxes at boundary cells
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->ut[j][k]=phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->dzz[grid->grad[2*j]][k];
    }
    for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->ut[j][k]=phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->dzz[grid->grad[2*j]][k];
    }
    
    // Compute the u-component fluxes at the faces
    if(prop->nonlinear==1)  // Upwind
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	
	for(k=kmin;k<grid->Nke[j];k++) {
	  if(phys->u[j][k]>0)
	    nc=nc2;
	  else
	    nc=nc1;
	  phys->ut[j][k]=phys->uc[nc][k]*grid->dzz[nc][k];
	}
      }
    else // Central
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	
	for(k=kmin;k<grid->Nke[j];k++) 
	  phys->ut[j][k]=0.5*InterpToFace(j,k,phys->uc,phys->u,grid)*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
      }
    
    // Now compute the cell-centered source terms and put them into stmp
    for(i=0;i<grid->Nc;i++) {
      
      for(nf=0;nf<NFACES;nf++) {
	
	ne = grid->face[i*NFACES+nf];
	
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	  phys->stmp[i][k]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][k]);
	
	// Top cell is filled with momentum from neighboring cells
	for(k=grid->etop[ne];k<=grid->ctop[i];k++) 
	  phys->stmp[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][grid->ctop[i]]);
      }
    }
    
    // V-fluxes at boundary cells
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->ut[j][k]=phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->dzz[grid->grad[2*j]][k];
    }
    for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->ut[j][k]=phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->dzz[grid->grad[2*j]][k];
    }
    
    // Compute the v-component fluxes at the faces
    if(prop->nonlinear==1)  // Upwind
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	for(k=kmin;k<grid->Nke[j];k++) {
	  if(phys->u[j][k]>0)
	    nc=nc2;
	  else
	    nc=nc1;
	  phys->ut[j][k]=phys->vc[nc][k]*grid->dzz[nc][k];
	}
      }
    else // Central
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	
	for(k=kmin;k<grid->Nke[j];k++) 
	  phys->ut[j][k]=0.5*InterpToFace(j,k,phys->vc,phys->u,grid)*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
      }
    
    // Now compute the cell-centered source terms and put them into stmp.
    for(i=0;i<grid->Nc;i++) {
      
      for(k=0;k<grid->Nk[i];k++) 
	phys->stmp2[i][k]=0;
      
      for(nf=0;nf<NFACES;nf++) {
	
	ne = grid->face[i*NFACES+nf];
	
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	  phys->stmp2[i][k]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][k]);
	
	// Top cell is filled with momentum from neighboring cells
	for(k=grid->etop[ne];k<=grid->ctop[i];k++) 
	  phys->stmp2[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][grid->ctop[i]]);
      }
    }
    
    // Now do vertical advection
    for(i=0;i<grid->Nc;i++) {
      
      if(prop->nonlinear==1)  // Upwind
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	  a[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
		      (phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
	  b[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
		      (phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
	}
      else  // Central
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	  a[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
				 grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1]));
	  b[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
				 grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1]));
	}
      
      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) {
	phys->stmp[i][k]+=(a[k]-a[k+1])/grid->dzz[i][k];
	phys->stmp2[i][k]+=(b[k]-b[k+1])/grid->dzz[i][k];
      }
      
      if(grid->ctop[i]!=grid->Nk[i]-1) {
	// Top - advection only comes in through bottom of cell.
	phys->stmp[i][grid->ctop[i]]-=a[grid->ctop[i]+1]/grid->dzz[i][grid->ctop[i]];
	phys->stmp2[i][grid->ctop[i]]-=b[grid->ctop[i]+1]/grid->dzz[i][grid->ctop[i]];
	// Bottom - advection only comes in through top of cell.
	phys->stmp[i][grid->Nk[i]-1]+=a[grid->Nk[i]-1]/grid->dzz[i][grid->Nk[i]-1];
	phys->stmp2[i][grid->Nk[i]-1]+=b[grid->Nk[i]-1]/grid->dzz[i][grid->Nk[i]-1];
      }
    }
  }

  // Now add on horizontal diffusion to stmp and stmp2
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];
    
    for(k=kmin;k<grid->Nke[j];k++) {
      a[k]=prop->nu_H*(phys->uc[nc2][k]-phys->uc[nc1][k])*grid->df[j]/grid->dg[j];
      b[k]=prop->nu_H*(phys->vc[nc2][k]-phys->vc[nc1][k])*grid->df[j]/grid->dg[j];
      phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
      phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
      phys->stmp2[nc1][k]-=b[k]/grid->Ac[nc1];
      phys->stmp2[nc2][k]+=b[k]/grid->Ac[nc2];
    }
    
    for(k=grid->Nke[j];k<grid->Nk[nc1];k++) {
      phys->stmp[nc1][k]+=prop->CdW*fabs(phys->uc[nc1][k])*phys->uc[nc1][k]*grid->df[j]/grid->Ac[nc1];
      phys->stmp2[nc1][k]+=prop->CdW*fabs(phys->vc[nc1][k])*phys->vc[nc1][k]*grid->df[j]/grid->Ac[nc1];
    }
    for(k=grid->Nke[j];k<grid->Nk[nc2];k++) {
      phys->stmp[nc2][k]+=prop->CdW*fabs(phys->uc[nc2][k])*phys->uc[nc2][k]*grid->df[j]/grid->Ac[nc2];
      phys->stmp2[nc2][k]+=prop->CdW*fabs(phys->vc[nc2][k])*phys->vc[nc2][k]*grid->df[j]/grid->Ac[nc2];
    }
  }
  
  // Check to make sure integrated fluxes are 0 for conservation
  // This will not be conservative if CdW or nu_H are nonzero!
  if(WARNING && prop->CdW==0 && prop->nu_H==0) {
    sum=0;
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	sum+=grid->Ac[i]*phys->stmp[i][k]*grid->dzz[i][k];
    }
    if(fabs(sum)>CONSERVED) printf("Warning, not U-momentum conservative!\n");

    sum=0;
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	sum+=grid->Ac[i]*phys->stmp2[i][k]*grid->dzz[i][k];
    }
    if(fabs(sum)>CONSERVED) printf("Warning, not V-momentum conservative!\n");
  }
  
  // Send/recv stmp and stmp2 to account for advective fluxes in ghost cells at
  // interproc boundaries.
  ISendRecvCellData3D(phys->stmp,grid,myproc,comm);
  ISendRecvCellData3D(phys->stmp2,grid,myproc,comm);

  for(jptr=grid->edgedist[2];jptr<0*grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];
    
    i = grid->grad[2*j];
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

    sum=0;
    for(nf=0;nf<NFACES;nf++) {
      if((nc=grid->neigh[i*NFACES+nf])!=-1) {
	sum+=grid->Ac[nc];
	for(k=grid->ctop[nc2];k<grid->Nk[nc2];k++) {
	  phys->stmp[i][k]+=grid->Ac[nc]*phys->stmp[nc][k];
	  phys->stmp2[i][k]+=grid->Ac[nc]*phys->stmp2[nc][k];
	}
      }
    }
    sum=1/sum;
    for(k=grid->ctop[nc2];k<grid->Nk[nc2];k++) {
      phys->stmp[i][k]*=sum;
      phys->stmp2[i][k]*=sum;
    }
  }
      
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    // Note that dgf==dg only when the cells are orthogonal!
    def1 = grid->def[nc1*NFACES+grid->gradf[2*j]];
    def2 = grid->def[nc2*NFACES+grid->gradf[2*j+1]];
    dgf = def1+def2;
    
    if(grid->ctop[nc1]>grid->ctop[nc2])
      k0=grid->ctop[nc1];
    else
      k0=grid->ctop[nc2];
    
    for(k=k0;k<grid->Nk[nc1];k++) 
      phys->Cn_U[j][k]-=def1/dgf
	*prop->dt*(phys->stmp[nc1][k]*grid->n1[j]+phys->stmp2[nc1][k]*grid->n2[j]);
    
    for(k=k0;k<grid->Nk[nc2];k++) 
      phys->Cn_U[j][k]-=def2/dgf
	*prop->dt*(phys->stmp[nc2][k]*grid->n1[j]+phys->stmp2[nc2][k]*grid->n2[j]);
  }
  
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 
    
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      phys->utmp[j][k]+=fab*phys->Cn_U[j][k];
  }

  // Now add on to the open boundaries
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr]; 
    
    nc1 = grid->grad[2*j];
    k0=grid->ctop[nc1];

    if(phys->boundary_flag[jptr-grid->edgedist[2]]==open)
      for(k=k0;k<grid->Nk[nc1];k++) 
	phys->Cn_U[j][k]-=0.5*prop->dt*(phys->stmp[nc1][k]*grid->n1[j]+phys->stmp2[nc1][k]*grid->n2[j]);
  }
  
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr]; 
    
    if(phys->boundary_flag[jptr-grid->edgedist[2]]==open)
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->utmp[j][k]+=fab*phys->Cn_U[j][k];
  }
}

/*
 * Function: NewCells
 * Usage: NewCells(grid,phys,prop);
 * --------------------------------
 * Adjust the velocity in the new cells that were previously dry.
 * This function is required for an Eulerian advection scheme because
 * it is difficult to compute the finite volume form of advection in the
 * upper cells when wetting and drying is employed.  In this function
 * the velocity in the new cells is set such that the quantity u*dz is
 * conserved from one time step to the next.  This works well without
 * wetting and drying.  When wetting and drying is employed, it is best
 * to extrapolate from the lower cells to obtain the velocity in the new
 * cells.
 *
 */
static void NewCells(gridT *grid, physT *phys, propT *prop) {
 
  int j, jptr, k, nc1, nc2;
  REAL dz;
 
  // Correct the velocity resulting from changes in volume
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
 
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
 
    if(grid->etop[j]<=grid->etopold[j]) {
      dz = 0;
      for(k=grid->etop[j];k<=grid->etopold[j];k++)
        dz+=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
 
      for(k=grid->etop[j];k<=grid->etopold[j];k++)
        phys->u[j][k]=phys->u[j][grid->etopold[j]]/dz*
	  (0.5*(grid->dzzold[nc1][k]+grid->dzzold[nc2][k]));
    }
  }
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
 
    nc1 = grid->grad[2*j];
 
    if(grid->etop[j]<=grid->etopold[j]) {
      dz = 0;
      for(k=grid->etop[j];k<=grid->etopold[j];k++)
        dz+=grid->dzz[nc1][k];
 
      for(k=grid->etop[j];k<=grid->etopold[j];k++)
        phys->u[j][k]=phys->u[j][grid->etopold[j]]/dz*
	  grid->dzzold[nc1][k];
    }
  }
}

/*
 * Function: WPredictor
 * Usage: WPredictor(grid,phys,prop,myproc,numprocs,comm);
 * -------------------------------------------------------
 * This function updates the vertical predicted velocity field with:
 *
 * 1) Old nonhydrostatic pressure gradient with theta method
 * 2) Horizontal and vertical advection of vertical momentum with AB2
 * 3) Horizontal laminar+turbulent diffusion of vertical momentum with AB2
 * 4) Vertical laminar+turbulent diffusion of vertical momentum with theta method
 *
 * Cn_W contains the Adams-Bashforth terms at time step n-1.
 * If wetting and drying is employed, no advection is computed in 
 * the upper cell.
 *
 */
static void WPredictor(gridT *grid, physT *phys, propT *prop,
		       int myproc, int numprocs, MPI_Comm comm) {
  int i, iptr, j, jptr, k, ne, nf, nc, nc1, nc2, kmin;
  REAL fab, sum, *a, *b, *c;

  a = phys->a;
  b = phys->b;
  c = phys->c;

  if(prop->n==1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
	phys->Cn_W[i][k]=0;
  } else
    fab=1.5;

  // Add on the nonhydrostatic pressure gradient from the previous time
  // step to compute the source term for the tridiagonal inversion.
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->wtmp[i][k]=phys->w[i][k]+(1-fab)*phys->Cn_W[i][k];
      phys->Cn_W[i][k]=0;
    }
    
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      phys->wtmp[i][k]-=2.0*prop->dt/(grid->dzz[i][k-1]+grid->dzz[i][k])*
	(phys->q[i][k-1]-phys->q[i][k]);
    phys->wtmp[i][grid->ctop[i]]+=2.0*prop->dt/grid->dzz[i][grid->ctop[i]]*
      phys->q[i][grid->ctop[i]];
  }

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=0;

  // Compute Eulerian advection (nonlinear!=0)
  if(prop->nonlinear) {
    // Compute the w-component fluxes at the faces

    // Fluxes at boundary faces
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->ut[j][k]=phys->boundary_w[jptr-grid->edgedist[2]][k]*grid->dzz[grid->grad[2*j]][k];
    }
    // Fluxes at boundary faces
    for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->ut[j][k]=phys->boundary_w[jptr-grid->edgedist[2]][k]*grid->dzz[grid->grad[2*j]][k];
    }

    if(prop->nonlinear==1) // Upwind
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	
	for(k=kmin;k<grid->Nke[j];k++) {
	  if(phys->u[j][k]>0)
	    nc = nc2;
	  else
	    nc = nc1;

	  phys->ut[j][k]=0.5*(phys->w[nc][k]+phys->w[nc][k+1])*grid->dzz[nc][k];
	}
      }
    else // Central 
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	for(k=kmin;k<grid->Nke[j];k++) 
	  phys->ut[j][k]=0.25*(InterpToFace(j,k,phys->w,phys->u,grid)+
			       InterpToFace(j,k+1,phys->w,phys->u,grid))*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
      }
    
    for(i=0;i<grid->Nc;i++) {
      
      for(nf=0;nf<NFACES;nf++) {
	
	ne = grid->face[i*NFACES+nf];
	
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	  phys->stmp[i][k]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][k]);
	
	// Top cell is filled with momentum from neighboring cells
	for(k=grid->etop[ne];k<=grid->ctop[i];k++) 
	  phys->stmp[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][grid->ctop[i]]);
      }
      
      // Vertical advection
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	phys->stmp[i][k]+=(pow(phys->w[i][k],2)-pow(phys->w[i][k+1],2))/grid->dzz[i][k];
      }
      // Top cell
      phys->stmp[i][grid->ctop[i]]-=pow(phys->w[i][grid->ctop[i]+1],2)/grid->dzz[i][grid->ctop[i]];
    }
    
    // Check to make sure integrated fluxes are 0 for conservation
    if(WARNING && prop->CdW==0 && prop->nu_H==0) {
      sum=0;
      for(i=0;i<grid->Nc;i++) {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  sum+=grid->Ac[i]*phys->stmp[i][k]*grid->dzz[i][k];
      }
      if(fabs(sum)>CONSERVED) printf("Warning, not W-momentum conservative!\n");
    }
  }

  // Add horizontal diffusion to stmp
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];

    for(k=kmin;k<grid->Nke[j];k++) {
      a[k]=.5*prop->nu_H*(phys->w[nc2][k]-phys->w[nc1][k]+phys->w[nc2][k+1]-phys->w[nc1][k+1])*grid->df[j]/grid->dg[j];
      phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
      phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
      //      phys->stmp[nc1][k]-=prop->nu_H*(phys->w[nc2][k]-phys->w[nc1][k])*grid->df[j]/grid->dg[j]/grid->Ac[nc1];
      //      phys->stmp[nc2][k]-=prop->nu_H*(phys->w[nc1][k]-phys->w[nc2][k])*grid->df[j]/grid->dg[j]/grid->Ac[nc2];
    }
    for(k=grid->Nke[j];k<grid->Nk[nc1];k++) 
      phys->stmp[nc1][k]+=0.25*prop->CdW*fabs(phys->w[nc1][k]+phys->w[nc1][k+1])*
	(phys->w[nc1][k]+phys->w[nc1][k+1])*grid->df[j]/grid->Ac[nc1];
    for(k=grid->Nke[j];k<grid->Nk[nc2];k++) 
      phys->stmp[nc2][k]+=0.25*prop->CdW*fabs(phys->w[nc2][k]+phys->w[nc2][k+1])*
	(phys->w[nc2][k]+phys->w[nc2][k+1])*grid->df[j]/grid->Ac[nc2];
  }

  //Now use the cell-centered advection terms to update the advection at the faces
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      phys->Cn_W[i][k]-=prop->dt*(grid->dzz[i][k-1]*phys->stmp[i][k-1]+grid->dzz[i][k]*phys->stmp[i][k])/
	(grid->dzz[i][k-1]+grid->dzz[i][k]);
    
    // Top flux advection consists only of top cell
    k=grid->ctop[i];
    phys->Cn_W[i][k]-=prop->dt*phys->stmp[i][k];
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->wtmp[i][k]+=fab*phys->Cn_W[i][k];
  }

  // wtmp now contains the right hand side without the vertical diffusion terms.  Now we
  // add the vertical diffusion terms to the explicit side and invert the tridiagonal for
  // vertical diffusion (only if grid->Nk[i]-grid->ctop[i]>=2)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    if(grid->Nk[i]-grid->ctop[i]>1) {
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	a[k] = 2*(prop->nu+phys->nu_tv[i][k-1])/grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
	b[k] = 2*(prop->nu+phys->nu_tv[i][k])/grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
      }
      b[grid->ctop[i]]=(prop->nu+phys->nu_tv[i][grid->ctop[i]])/pow(grid->dzz[i][grid->ctop[i]],2);
      a[grid->ctop[i]]=b[grid->ctop[i]];
      
      // Add on the explicit part of the vertical diffusion term
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
	phys->wtmp[i][k]+=prop->dt*(1-prop->theta)*(a[k]*phys->w[i][k-1]
						    -(a[k]+b[k])*phys->w[i][k]
						    +b[k]*phys->w[i][k+1]);
      phys->wtmp[i][grid->ctop[i]]+=prop->dt*(1-prop->theta)*(-(a[grid->ctop[i]]+b[grid->ctop[i]])*phys->w[i][grid->ctop[i]]
							      +(a[grid->ctop[i]]+b[grid->ctop[i]])*phys->w[i][grid->ctop[i]+1]);
      
      // Now formulate the components of the tridiagonal inversion.
      // c is the diagonal entry, a is the lower diagonal, and b is the upper diagonal.
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	c[k]=1+prop->dt*prop->theta*(a[k]+b[k]);
	a[k]*=(-prop->dt*prop->theta);
	b[k]*=(-prop->dt*prop->theta);
      }
      b[grid->ctop[i]]+=a[grid->ctop[i]];

      TriSolve(&(a[grid->ctop[i]]),&(c[grid->ctop[i]]),&(b[grid->ctop[i]]),
	       &(phys->wtmp[i][grid->ctop[i]]),&(phys->w[i][grid->ctop[i]]),grid->Nk[i]-grid->ctop[i]);
    } else {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	phys->w[i][k]=phys->wtmp[i][k];
    }
  }
}

/*
 * Function: Corrector
 * Usage: Corrector(qc,grid,phys,prop,myproc,numprocs,comm);
 * ---------------------------------------------------------
 * Correct the horizontal velocity field with the pressure correction.
 * Do not correct velocities for which D[j]==0 since these are boundary
 * cells.
 *
 * After correcting the horizontal velocity, update the total nonhydrostatic
 * pressure with the pressure correction.
 *
 */
static void Corrector(REAL **qc, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr, j, jptr, k;

  // Correct the horizontal velocity only if this is not a boundary point.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 
    if(phys->D[j]!=0 && grid->etop[j]<grid->Nke[j]-1)
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->u[j][k]-=prop->dt/grid->dg[j]*
	  (qc[grid->grad[2*j]][k]-qc[grid->grad[2*j+1]][k]);
  }

  // Correct the vertical velocity
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
      phys->w[i][k]-=2.0*prop->dt/(grid->dzz[i][k-1]+grid->dzz[i][k])*
	(qc[i][k-1]-qc[i][k]);
    phys->w[i][grid->ctop[i]]+=2.0*prop->dt/grid->dzz[i][grid->ctop[i]]*
      qc[i][grid->ctop[i]];
  }

  // Update the pressure since qc is a pressure correction
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    if(grid->ctop[i]<grid->Nk[i]-1)
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	phys->q[i][k]+=qc[i][k];
  }
}

/*
 * Function: ComputeQSource
 * Usage: ComputeQSource(src,grid,phys,prop,myproc,numprocs);
 * ----------------------------------------------------------
 * Compute the source term for the nonhydrostatic pressure by computing
 * the divergence of the predicted velocity field, which is in phys->u and
 * phys->w.  The upwind flux face heights are used in order to ensure
 * consistency with continuity.
 *
 */
static void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs) {
  
  int i, iptr, j, jptr, k, nf, ne, nc1, nc2;
  REAL *ap=phys->a, *am=phys->b;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k] = grid->Ac[i]*(phys->w[i][k]-phys->w[i][k+1]);

    for(nf=0;nf<NFACES;nf++) {

      ne = grid->face[i*NFACES+nf];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      for(k=grid->ctop[i];k<grid->Nke[ne];k++) {
	ap[k] = 0.5*(phys->u[ne][k]+fabs(phys->u[ne][k]));
	am[k] = 0.5*(phys->u[ne][k]-fabs(phys->u[ne][k]));
      }
      for(k=grid->Nke[ne];k<grid->Nkc[ne];k++) {
	ap[k] = 0;
	am[k] = 0;
      }

      for(k=grid->ctop[i];k<grid->Nke[ne];k++) 
	src[i][k]+=(grid->dzz[nc2][k]*ap[k]+grid->dzz[nc1][k]*am[k])*
	  grid->normal[i*NFACES+nf]*grid->df[ne];
    }

    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k]/=prop->dt;
  }

  // D[j] is used in OperatorQ, and it must be zero to ensure no gradient
  // at the hydrostatic faces.
  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=grid->df[j]/grid->dg[j];
  }

  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(nf=0;nf<NFACES;nf++)
      phys->D[grid->face[i*NFACES+nf]]=0;
  }
}

/*
 * Function: CGSolveQ
 * Usage: CGSolveQ(q,src,c,grid,phys,prop,myproc,numprocs,comm);
 * -------------------------------------------------------------
 * Solve for the nonhydrostatic pressure with the preconditioned
 * conjugate gradient algorithm.
 *
 * The preconditioner stores the diagonal preconditioning elements in 
 * the temporary c array.
 *
 * This function replaces q with x and src with p.  phys->uc and phys->vc
 * are used as temporary arrays as well to store z and r.
 *
 */
static void CGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr, k, n, niters;

  REAL **x, **r, **rtmp, **p, **z, mu, nu, alpha, alpha0, eps, eps0;

  z = phys->stmp2;
  x = q;
  r = phys->stmp3;
  rtmp = phys->uold;
  p = src;

  // Compute the preconditioner and the preconditioned solution
  // and send it to neighboring processors
  if(prop->qprecond==1) {
    ConditionQ(c,grid,phys,prop,myproc,comm);
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	p[i][k]/=c[i][k];
	x[i][k]*=c[i][k];
      }
    }
  }
  ISendRecvCellData3D(x,grid,myproc,comm);

  niters = prop->qmaxiters;

  // Create the coefficients for the operator
  QCoefficients(phys->wtmp,phys->qtmp,c,grid,phys,prop);

  // Initialization for CG
  if(prop->qprecond==1) OperatorQC(phys->wtmp,phys->qtmp,x,z,c,grid,phys,prop);
  else OperatorQ(phys->wtmp,x,z,c,grid,phys,prop);
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      r[i][k] = p[i][k]-z[i][k];
  }    
  if(prop->qprecond==2) {
    Preconditioner(r,rtmp,phys->wtmp,grid,phys,prop);
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	p[i][k] = rtmp[i][k];
    }
    alpha = alpha0 = InnerProduct3(r,rtmp,grid,myproc,numprocs,comm);
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	p[i][k] = r[i][k];
    }
    alpha = alpha0 = InnerProduct3(r,r,grid,myproc,numprocs,comm);
  }
  if(!prop->resnorm) alpha0 = 1;

  if(prop->qprecond==2)
    eps=eps0=InnerProduct3(r,r,grid,myproc,numprocs,comm);
  else
    eps=eps0=alpha0;

  // Iterate until residual is less than prop->qepsilon
  for(n=0;n<niters && eps!=0;n++) {

    ISendRecvCellData3D(p,grid,myproc,comm);
    if(prop->qprecond==1) OperatorQC(phys->wtmp,phys->qtmp,p,z,c,grid,phys,prop);
    else OperatorQ(phys->wtmp,p,z,c,grid,phys,prop);

    mu = 1/alpha;
    nu = alpha/InnerProduct3(p,z,grid,myproc,numprocs,comm);

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	x[i][k] += nu*p[i][k];
	r[i][k] -= nu*z[i][k];
      }
    }
    if(prop->qprecond==2) {
      Preconditioner(r,rtmp,phys->wtmp,grid,phys,prop);
      alpha = InnerProduct3(r,rtmp,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	i = grid->cellp[iptr];

	for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	  p[i][k] = rtmp[i][k] + mu*p[i][k];
      }
    } else {
      alpha = InnerProduct3(r,r,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	i = grid->cellp[iptr];

	for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	  p[i][k] = r[i][k] + mu*p[i][k];
      }
    }

    if(prop->qprecond==2)
      eps=InnerProduct3(r,r,grid,myproc,numprocs,comm);
    else
      eps=alpha;

    //    if(myproc==0) printf("%d %e %e %e %e\n",myproc,mu,nu,eps0,sqrt(eps/eps0));
    if(VERBOSE>3) printf("CGSolve Pressure Iteration: %d, resid=%e, proc=%d\n",n,sqrt(eps/eps0),myproc);
    if(sqrt(eps/eps0)<prop->qepsilon) 
      break;
  }
  if(myproc==0 && VERBOSE>2) 
    if(eps==0)
      printf("Warning...Time step %d, norm of pressure source is 0.\n",prop->n);
    else
      if(n==niters)  printf("Warning... Time step %d, Pressure iteration not converging after %d steps! RES=%e > %.2e\n",
			    prop->n,n,sqrt(eps/eps0),prop->qepsilon);
      else printf("Time step %d, CGSolve pressure converged after %d iterations, res=%e < %.2e\n",
		  prop->n,n,sqrt(eps/eps0),prop->qepsilon);

  // Rescale the preconditioned solution 
  if(prop->qprecond==1) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	x[i][k]/=c[i][k];
    }
  }
  
  // Send the solution to the neighboring processors
  ISendRecvCellData3D(x,grid,myproc,comm);
}

/*
 * Function: EddyViscosity
 * Usage: EddyViscosity(grid,phys,prop);
 * -------------------------------------
 * This function is used to compute the eddy viscosity, the
 * shear stresses, and the drag coefficients at the upper and lower
 * boundaries.
 *
 */
static void EddyViscosity(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  if(prop->turbmodel) 
    my25(grid,phys,prop,phys->qT,phys->lT,phys->Cn_q,phys->Cn_l,phys->nu_tv,phys->kappa_tv,comm,myproc);
}

/*
 * Function: UPredictor 
 * Usage: UPredictor(grid,phys,prop,myproc,numprocs,comm);
 * -------------------------------------------------------
 * Predictor step for the horizontal velocity field.  This function
 * computes the free surface using the theta method and then uses it to update the predicted
 * velocity field in the absence of the nonhydrostatic pressure.
 *
 * Upon entry, phys->utmp contains the right hand side of the u-momentum equation
 *
 */
static void UPredictor(gridT *grid, physT *phys, 
		      propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, j, jptr, ne, nf, nf1, normal, nc1, nc2, k;
  REAL sum, dt=prop->dt, theta=prop->theta, fluxheight, h0, boundary_flag;
  REAL *a, *b, *c, *d, *e1, **E, *a0, *b0, *c0, *d0;

  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;
  e1 = phys->ap;
  E = phys->ut;

  a0 = phys->am;
  b0 = phys->bp;
  c0 = phys->bm;

  // Set D[j] = 0 
  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i]+1;k++) 
      phys->wtmp2[i][k]=phys->w[i][k];

  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=0;
    for(k=0;k<grid->Nke[j];k++)
      phys->utmp2[j][k]=phys->u[j][k];
  }

  // phys->u contains the velocity specified at the open boundaries
  // It is also the velocity at time step n.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]=phys->u[j][k];
  }

  // Update the velocity in the interior nodes with the old free-surface gradient
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // Add the explicit part of the free-surface to create U**.
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]-=GRAV*(1-theta)*dt*(phys->h[nc1]-phys->h[nc2])/grid->dg[j];
  }

  // Update the boundary faces with the linearized free-surface gradient at the boundary
  // i.e. using the radiative condition, dh/dn = 1/c dh/dt
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    boundary_flag=phys->boundary_flag[jptr-grid->edgedist[2]];
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]+=(-2.0*GRAV*(1-boundary_flag)*(1-theta)*dt/grid->dg[j]
			 +boundary_flag*sqrt(GRAV/(grid->dv[i]+phys->h[i])))*phys->h[i]
	+2.0*GRAV*(1-boundary_flag)*dt/grid->dg[j]*phys->boundary_h[jptr-grid->edgedist[2]];
  }

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;

    // Add the shear stress from the top cell
    phys->utmp[j][grid->etop[j]]+=2.0*dt*phys->tau_T[j]/
      (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]]);

    // Create the tridiagonal entries and formulate U***
    if(!(grid->dzz[nc1][grid->etop[j]]==0 && grid->dzz[nc2][grid->etop[j]]==0)) {

      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	a[k]=0;
	b[k]=0;
	c[k]=0;
	d[k]=0;
      }

      // Vertical eddy-viscosity interpolated to faces since it is stored
      // at cell-centers.
      for(k=grid->etop[j]+1;k<grid->Nke[j];k++) 
	c[k]=0.25*(phys->nu_tv[nc1][k-1]+phys->nu_tv[nc2][k-1]+
		   phys->nu_tv[nc1][k]+phys->nu_tv[nc2][k]);
      
      // Coefficients for the viscous terms.  Face heights are taken as
      // the average of the face heights on either side of the face (not upwinded).
      for(k=grid->etop[j]+1;k<grid->Nke[j];k++) 
	a[k]=2.0*(prop->nu+c[k])/(0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k])*
				  (grid->dzz[nc1][k-1]+grid->dzz[nc2][k-1]+
				   grid->dzz[nc1][k]+grid->dzz[nc2][k]));
      
      for(k=grid->etop[j];k<grid->Nke[j]-1;k++) {
	b[k]=2.0*(prop->nu+c[k+1])/(0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k])*
				    (grid->dzz[nc1][k]+grid->dzz[nc2][k]+
				     grid->dzz[nc1][k+1]+grid->dzz[nc2][k+1]));
      }

      if(grid->Nke[j]-grid->etop[j]>1) {

	// Explicit part of the viscous term
	for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++)
	  phys->utmp[j][k]+=dt*(1-theta)*(a[k]*phys->u[j][k-1]-
					  (a[k]+b[k])*phys->u[j][k]+
					  b[k]*phys->u[j][k+1]);
	
	// Top cell
	phys->utmp[j][grid->etop[j]]+=dt*(1-theta)*(-(b[grid->etop[j]]+2.0*phys->CdT[j]*
						      fabs(phys->u[j][grid->etop[j]])/
						      (grid->dzz[nc1][grid->etop[j]]+
						       grid->dzz[nc2][grid->etop[j]]))*
						    phys->u[j][grid->etop[j]]
						    +b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);

	// Bottom cell
	phys->utmp[j][grid->Nke[j]-1]+=dt*(1-theta)*(a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2]-
						     (a[grid->Nke[j]-1]+2.0*phys->CdB[j]*
						      fabs(phys->u[j][grid->Nke[j]-1])/
						      (grid->dzz[nc1][grid->Nke[j]-1]+
						       grid->dzz[nc2][grid->Nke[j]-1]))*
						     phys->u[j][grid->Nke[j]-1]);
      
      } else 
	phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdB[j]+phys->CdT[j])/
	  (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
	  fabs(phys->u[j][grid->etop[j]])*phys->u[j][grid->etop[j]];

      // Now set up the coefficients for the tridiagonal inversion for the
      // implicit part.  These are given from the arrays above in the discrete operator
      // d^2U/dz^2 = -theta dt a_k U_{k-1} + (1+theta dt (a_k+b_k)) U_k - theta dt b_k U_{k+1}

      // Right hand side U** is given by d[k] here.
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	e1[k]=1.0;
	d[k]=phys->utmp[j][k];
      }
      
      if(grid->Nke[j]-grid->etop[j]>1) {
	// Top cells
	c[grid->etop[j]]=-theta*dt*b[grid->etop[j]];
	b[grid->etop[j]]=1.0+theta*dt*(b[grid->etop[j]]+
				       2.0*phys->CdT[j]*fabs(phys->u[j][grid->etop[j]])/
				       (grid->dzz[nc1][grid->etop[j]]+
					grid->dzz[nc2][grid->etop[j]]));
	a[grid->etop[j]]=0;
	
	// Bottom cell
	c[grid->Nke[j]-1]=0;
	b[grid->Nke[j]-1]=1.0+theta*dt*(a[grid->Nke[j]-1]+
					2.0*phys->CdB[j]*fabs(phys->u[j][grid->Nke[j]-1])/
					(grid->dzz[nc1][grid->Nke[j]-1]+
					 grid->dzz[nc2][grid->Nke[j]-1]));
	a[grid->Nke[j]-1]=-theta*dt*a[grid->Nke[j]-1];
      
	// Interior cells
	for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
	  c[k]=-theta*dt*b[k];
	  b[k]=1.0+theta*dt*(a[k]+b[k]);
	  a[k]=-theta*dt*a[k];
	}
      } else {
	b[grid->etop[j]]=1.0+2.0*theta*dt*fabs(phys->u[j][grid->etop[j]])/
	  (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
	  (phys->CdB[j]+phys->CdT[j]);
      }	  

      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	if(grid->dzz[nc1][k]==0 && grid->dzz[nc2][k]==0) {
	  printf("Exiting because dzz[%d][%d]=%f or dzz[%d][%d]=%f\n",
		 nc1,k,grid->dzz[nc1][k],nc2,k,grid->dzz[nc2][k]);
	  exit(0);
	}
	if(a[k]!=a[k]) printf("a[%d] problems, dzz[%d][%d]=%f\n",k,j,k,grid->dzz[j][k]);
	if(b[k]!=b[k] || b[k]==0) printf("b[%d] problems\n",k);
	if(c[k]!=c[k]) printf("c[%d] problems\n",k);
      }

      // Now utmp will have U*** in it, which is given by A^{-1}U**, and E will have
      // A^{-1}e1, where e1 = [1,1,1,1,1,...,1]^T 
      // Store the tridiagonals so they can be used twice (TriSolve alters the values
      // of the elements in the diagonals!!!
      for(k=0;k<grid->Nke[j];k++) {
	a0[k]=a[k];
	b0[k]=b[k];
	c0[k]=c[k];
      }
      if(grid->Nke[j]-grid->etop[j]>1) {
	TriSolve(&(a[grid->etop[j]]),&(b[grid->etop[j]]),&(c[grid->etop[j]]),
		 &(d[grid->etop[j]]),&(phys->utmp[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);
	TriSolve(&(a0[grid->etop[j]]),&(b0[grid->etop[j]]),&(c0[grid->etop[j]]),
		 &(e1[grid->etop[j]]),&(E[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);	
      } else {
	phys->utmp[j][grid->etop[j]]/=b[grid->etop[j]];
	E[j][grid->etop[j]]=1.0/b[grid->etop[j]];
      }

      // Now vertically integrate E to create the vertically integrated flux-face
      // values that comprise the coefficients of the free-surface solver.  This
      // will create the D vector, where D=DZ^T E (which should be given by the
      // depth when there is no viscosity.
      phys->D[j]=0;
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	fluxheight=UpWind(phys->u[j][k],grid->dzz[nc1][k],grid->dzz[nc2][k]);

	phys->D[j]+=fluxheight*E[j][k];
      }
    }
  }

  for(j=0;j<grid->Ne;j++) 
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      if(phys->utmp[j][k]!=phys->utmp[j][k]) {
	printf("Error in function Predictor at j=%d k=%d (U***=nan)\n",j,k);
	exit(1);
      }

  // So far we have U*** and D.  Now we need to create h* in htmp.   This
  // will comprise the source term for the free-surface solver.  Before we
  // do this we need to set the new velocity at the open boundary faces and
  // place them into utmp.  These need the velocity from the old time step uold
  // As well as the current vectors.
  OpenBoundaryFluxes(phys->q,phys->utmp,phys->utmp2,grid,phys,prop);

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    sum = 0;
    for(nf=0;nf<NFACES;nf++) {
      
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      for(k=grid->etop[ne];k<grid->Nke[ne];k++) {
	fluxheight=UpWind(phys->u[ne][k],grid->dzz[nc1][k],grid->dzz[nc2][k]);

	sum+=((1-theta)*phys->u[ne][k]+theta*phys->utmp[ne][k])*
	  fluxheight*grid->df[ne]*normal;
      }
    }
    phys->htmp[i]=phys->h[i]-dt/grid->Ac[i]*sum;
  }

  // Now we have the required components for the CG solver for the free-surface:
  //
  // h^{n+1} - g*(theta*dt)^2/Ac * Sum_{faces} D_{face} dh^{n+1}/dn df N = htmp
  //
  // L(h) = b
  //
  // L(h) = h + 1/Ac * Sum_{faces} D_{face} dh^{n+1}/dn N
  // b = htmp
  //
  // As the initial guess let h^{n+1} = h^n, so just leave it as it is to
  // begin the solver.
  if(prop->cgsolver==0)
    GSSolve(grid,phys,prop,myproc,numprocs,comm);
  else if(prop->cgsolver==1)
    CGSolve(grid,phys,prop,myproc,numprocs,comm);

  // Add on the implicit barotropic term to obtain the hydrostatic horizontal velocity field.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k]=phys->utmp[j][k]-GRAV*theta*dt*E[j][k]*
	(phys->h[nc1]-phys->h[nc2])/grid->dg[j];

    if(grid->etop[j]==grid->Nke[j]-1 && grid->dzz[nc1][grid->etop[j]]==0 &&
       grid->dzz[nc2][grid->etop[j]]==0) 
      phys->u[j][grid->etop[j]]=0;
  }

  // Set the flux values at boundary cells if specified (marker=4)
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    boundary_flag = phys->boundary_flag[jptr-grid->edgedist[2]];
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k]=phys->utmp[j][k]-(2.0*GRAV*(1-boundary_flag)*theta*dt/grid->dg[j]+
				      boundary_flag*sqrt(GRAV/(grid->dv[i]+phys->h[i])))*E[j][k]*phys->h[i];
  }

  // Now update the vertical grid spacing with the new free surface.
  UpdateDZ(grid,phys,0);

  // Use the new free surface to add the implicit part of the free-surface
  // pressure gradient to the horizontal momentum.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(grid->etop[j]>grid->etopold[j]) 
      for(k=0;k<grid->etop[j];k++)
	phys->u[j][k]=0;
    else 
      for(k=grid->etop[j];k<grid->etopold[j];k++)
	phys->u[j][k]=phys->utmp[j][k]-GRAV*theta*dt*
	  (phys->h[nc1]-phys->h[nc2])/grid->dg[j];
  }

  // Set the flux values at boundary cells if specified (marker=4)
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    boundary_flag = phys->boundary_flag[jptr-grid->edgedist[2]];
    if(grid->etop[j]>grid->etopold[j]) 
      for(k=0;k<grid->etop[j];k++)
	phys->u[j][k]=0;
    else 
      for(k=grid->etop[j];k<grid->etopold[j];k++) 
	phys->u[j][k]=phys->utmp[j][k]-(2.0*GRAV*(1-boundary_flag)*theta*dt/grid->dg[j]+
					boundary_flag*sqrt(GRAV/(grid->dv[i]+phys->h[i])))*phys->h[i];
  }

  // Set the flux values at the open boundary (marker=2).  These
  // were set to utmp previously in OpenBoundaryFluxes.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k] = phys->utmp[j][k];
  }

  // Now set the fluxes at the free-surface boundary by assuming dw/dz = 0
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      if(grid->mark[ne]==3) {
	for(k=grid->etop[ne];k<grid->Nke[ne];k++) {
	  phys->u[ne][k] = 0;
	  sum=0;
	  for(nf1=0;nf1<NFACES;nf1++)
	    sum+=phys->u[grid->face[i*NFACES+nf1]][k]*grid->df[grid->face[i*NFACES+nf1]]*grid->normal[i*NFACES+nf1];
	  phys->u[ne][k]=-sum/grid->df[ne]/grid->normal[i*NFACES+nf];
	}
      }
    }
  }
}

/*
 * Function: CGSolve
 * Usage: CGSolve(grid,phys,prop,myproc,numprocs,comm);
 * ----------------------------------------------------
 * Solve the free surface equation using the conjugate gradient algorithm.
 *
 * The source term upon entry is in phys->htmp, which is placed into p, and
 * the free surface upon entry is in phys->h, which is placed into x.
 *
 */
static void CGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, n, niters;
  REAL *x, *r, *D, *p, *z, mu, nu, eps, eps0;

  z = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"CGSolve");
  x = phys->h;
  r = phys->hold;
  D = phys->D;
  p = phys->htmp;

  niters = prop->maxiters;

  // For the boundary term (marker of type 3):
  // 1) Need to set x to zero in the interior points, but
  //    leave it as is for the boundary points.
  // 2) Then set z=Ax and substract b = b-z so that
  //    the new problem is Ax=b with the boundary values
  //    on the right hand side acting as forcing terms.
  // 3) After b=b-z for the interior points, then need to
  //    set b=0 for the boundary points.
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    x[i]=0;
  }
  ISendRecvCellData2D(x,grid,myproc,comm);
  OperatorH(x,z,grid,phys,prop);

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    p[i] = p[i] - z[i];
    r[i] = p[i];
    x[i] = 0;
  }    
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    p[i] = 0;
  }    
  eps0 = eps = InnerProduct(r,r,grid,myproc,numprocs,comm);
  if(!prop->resnorm) eps0 = 1;

  for(n=0;n<niters && eps!=0;n++) {

    ISendRecvCellData2D(p,grid,myproc,comm);
    OperatorH(p,z,grid,phys,prop);

    mu = 1/eps;
    nu = eps/InnerProduct(p,z,grid,myproc,numprocs,comm);

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      x[i] += nu*p[i];
      r[i] -= nu*z[i];
    }
    eps = InnerProduct(r,r,grid,myproc,numprocs,comm);
    mu*=eps;

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      p[i] = r[i] + mu*p[i];
    }

    if(VERBOSE>3) printf("CGSolve Iteration: %d, resid=%e, proc=%d\n",n,sqrt(eps/eps0),myproc);
    if(sqrt(eps/eps0)<prop->epsilon) 
      break;
  }
  if(myproc==0 && VERBOSE>2) 
    if(eps==0)
      printf("Warning...Time step %d, norm of free-surface source is 0.\n",prop->n);
    else
      if(n==niters)  printf("Warning... Time step %d, Free-surface iteration not converging after %d steps! RES=%e > %.2e\n",
			    prop->n,n,sqrt(eps/eps0),prop->epsilon);
      else printf("Time step %d, CGSolve surface converged after %d iterations, res=%e < %.2e\n",
		  prop->n,n,sqrt(eps/eps0),prop->epsilon);

  ISendRecvCellData2D(x,grid,myproc,comm);
  SunFree(z,grid->Nc*sizeof(REAL),"CGSolve");
}

/*
 * Function: InnerProduct
 * Usage: InnerProduct(x,y,grid,myproc,numprocs,comm);
 * ---------------------------------------------------
 * Compute the inner product of two one-dimensional arrays x and y.
 * Used for the CG method to solve for the free surface.
 *
 */
static REAL InnerProduct(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, MPI_Comm comm) {
  
  int i, iptr;
  REAL sum, mysum=0;
  
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    mysum+=x[i]*y[i];
  }
  MPI_Reduce(&mysum,&(sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&sum,1,MPI_DOUBLE,0,comm);

  return sum;
}  

/*
 * Function: InnerProduct3
 * Usage: InnerProduct3(x,y,grid,myproc,numprocs,comm);
 * ---------------------------------------------------
 * Compute the inner product of two two-dimensional arrays x and y.
 * Used for the CG method to solve for the nonhydrostatic pressure.
 *
 */
static REAL InnerProduct3(REAL **x, REAL **y, gridT *grid, int myproc, int numprocs, MPI_Comm comm) {
  
  int i, k, iptr;
  REAL sum, mysum=0;
  
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      mysum+=x[i][k]*y[i][k];
  }
  MPI_Reduce(&mysum,&(sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&sum,1,MPI_DOUBLE,0,comm);

  return sum;
}  

/*
 * Function: OperatorH
 * Usage: OperatorH(x,y,grid,phys,prop);
 * -------------------------------------
 * Given a vector x, computes the left hand side of the free surface 
 * Poisson equation and places it into y with y = L(x).
 *
 */
static void OperatorH(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop) {
  
  int i, j, iptr, jptr, ne, nf;
  REAL tmp = GRAV*pow(prop->theta*prop->dt,2), h0, boundary_flag;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      y[i] = x[i];
      for(nf=0;nf<NFACES;nf++) 
	if(grid->neigh[i*NFACES+nf]!=-1) {
	  ne = grid->face[i*NFACES+nf];

	  y[i]+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
	    (x[i]-x[grid->neigh[i*NFACES+nf]])/grid->Ac[i];
	}
  }

  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    boundary_flag = phys->boundary_flag[jptr-grid->edgedist[2]];
    y[i] += prop->dt*prop->theta*(2.0*GRAV*(1-boundary_flag)*prop->theta*prop->dt/grid->dg[j]
				  +boundary_flag*sqrt(GRAV/(grid->dv[i]+phys->boundary_h[jptr-grid->edgedist[2]])))*
      (grid->dv[i]+phys->boundary_h[jptr-grid->edgedist[2]])*grid->df[j]/grid->Ac[i]*x[i];
  }
}

/*
 * Function: OperatorQC
 * Usage: OperatorQC(coef,fcoef,x,y,c,grid,phys,prop);
 * ---------------------------------------------------
 * Given a vector x, computes the left hand side of the nonhydrostatic pressure
 * Poisson equation and places it into y with y = L(x) for the preconditioned
 * solver.
 *
 * The coef array contains coefficients for the vertical derivative terms in the operator
 * while the fcoef array contains coefficients for the horizontal derivative terms.  These
 * are computed before the iteration in QCoefficients. The array c stores the preconditioner.
 *
 */
static void OperatorQC(REAL **coef, REAL **fcoef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop) {

  int i, iptr, k, ne, nf, nc, kmin, kmax;
  REAL *a = phys->a;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	y[i][k]=-x[i][k];

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {

	    ne = grid->face[i*NFACES+nf];

	    if(grid->ctop[nc]>grid->ctop[i])
	      kmin = grid->ctop[nc];
	    else
	      kmin = grid->ctop[i];
	    
	    for(k=kmin;k<grid->Nke[ne];k++) 
	      y[i][k]+=x[nc][k]*fcoef[i*NFACES+nf][k];
	  }

      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	y[i][k]+=coef[i][k]*x[i][k-1]+coef[i][k+1]*x[i][k+1];

      if(grid->ctop[i]<grid->Nk[i]-1) {
	// Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
	k=grid->ctop[i];
	y[i][k]+=coef[i][k+1]*x[i][k+1];

	// Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
	k=grid->Nk[i]-1;
	y[i][k]+=coef[i][k]*x[i][k-1];
      }
  }
}

/*
 * Function: QCoefficients
 * Usage: QCoefficients(coef,fcoef,c,grid,phys,prop);
 * --------------------------------------------------
 * Compute coefficients for the pressure-Poisson equation.  fcoef stores
 * coefficients at the vertical flux faces while coef stores coefficients
 * at the horizontal faces for vertical derivatives of q.
 *
 */
static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, physT *phys, propT *prop) {

  int i, iptr, k, kmin, nf, nc, ne;

  if(prop->qprecond==1) 
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]]/c[i][grid->ctop[i]];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
	coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1])/(c[i][k]*c[i][k-1]);

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {

	    ne = grid->face[i*NFACES+nf];

	    if(grid->ctop[nc]>grid->ctop[i])
	      kmin = grid->ctop[nc];
	    else
	      kmin = grid->ctop[i];
	    
	    for(k=kmin;k<grid->Nke[ne];k++) 
	      fcoef[i*NFACES+nf][k]=grid->dzz[i][k]*phys->D[ne]/(c[i][k]*c[nc][k]);
	}
    }
  else
    for(i=0;i<grid->Nc;i++) {//iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      //      i = grid->cellp[iptr];
      
      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
	coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
  }
}

/*
 * Function: OperatorQ
 * Usage: OperatorQ(coef,x,y,c,grid,phys,prop);
 * --------------------------------------------
 * Given a vector x, computes the left hand side of the nonhydrostatic pressure
 * Poisson equation and places it into y with y = L(x) for the non-preconditioned
 * solver.
 *
 * The coef array contains coefficients for the vertical derivative terms in the operator.
 * This is computed before the iteration in QCoefficients. The array c stores the preconditioner.
 * The preconditioner stored in c is not used.
 *
 */
static void OperatorQ(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop) {

  int i, iptr, k, ne, nf, nc, kmin, kmax;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	y[i][k]=0;

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {
	  
	    ne = grid->face[i*NFACES+nf];

	    if(grid->ctop[nc]>grid->ctop[i])
	      kmin = grid->ctop[nc];
	    else
	      kmin = grid->ctop[i];

	    for(k=kmin;k<grid->Nke[ne];k++) 
	      y[i][k]+=(x[nc][k]-x[i][k])*grid->dzz[i][k]*phys->D[ne];
	  }

      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	y[i][k]+=coef[i][k]*x[i][k-1]-(coef[i][k]+coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

      if(grid->ctop[i]<grid->Nk[i]-1) {
	// Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
	k=grid->ctop[i];
	y[i][k]+=(-2*coef[i][k]-coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

	// Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
	k=grid->Nk[i]-1;
	y[i][k]+=coef[i][k]*x[i][k-1]-coef[i][k]*x[i][k];
      } else
	y[i][grid->ctop[i]]-=2.0*coef[i][grid->ctop[i]]*x[i][grid->ctop[i]];
  }
}

/*
 * Function: GuessQ
 * Usage: Guessq(q,wold,w,grid,phys,prop,myproc,numprocs,comm);
 * ------------------------------------------------------------
 * Guess a pressure correction field that will enforce the hydrostatic velocity
 * field to speed up the convergence of the pressure Poisson equation.
 *
 */
static void GuessQ(REAL **q, REAL **wold, REAL **w, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {
  
  int i, iptr, k;
  REAL qerror;

  // First compute the vertical velocity field that would satisfy continuity
  Continuity(w,grid,phys,prop);

  // Then use this velocity field to back out the required pressure field by
  // integrating w^{n+1}=w*-dt dq/dz and imposing the q=0 boundary condition at
  // the free-surface.
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
  
    q[i][grid->ctop[i]]=grid->dzz[i][grid->ctop[i]]/2/prop->dt/prop->theta*
      (w[i][grid->ctop[i]]-wold[i][grid->ctop[i]]);
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      q[i][k]=q[i][k-1]+(grid->dzz[i][k]+grid->dzz[i][k-1])/(2.0*prop->dt*prop->theta)*
	(w[i][k]-wold[i][k]);
    }
  }
}

/*
 * Function: Preconditioner
 * Usage: Preconditioner(x,xc,coef,grid,phys,prop);
 * ------------------------------------------------
 * Multiply the vector x by the inverse of the preconditioner M with
 * xc = M^{-1} x
 *
 */
static void Preconditioner(REAL **x, REAL **xc, REAL **coef, gridT *grid, physT *phys, propT *prop) {
  int i, iptr, k, nf, ne, nc, kmin;
  REAL *a = phys->a, *b = phys->b, *c = phys->c, *d = phys->d;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i=grid->cellp[iptr];

    if(grid->ctop[i]<grid->Nk[i]-1) {
      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) {
	a[k]=coef[i][k];
	b[k]=-coef[i][k]-coef[i][k+1];
	c[k]=coef[i][k+1];
	d[k]=x[i][k];
      }

      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      b[k]=-2*coef[i][k]-coef[i][k+1];
      c[k]=coef[i][k+1];
      d[k]=x[i][k];
      
      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      a[k]=coef[i][k];
      b[k]=-coef[i][k];
      d[k]=x[i][k];

      TriSolve(&(a[grid->ctop[i]]),&(b[grid->ctop[i]]),&(c[grid->ctop[i]]),
      	       &(d[grid->ctop[i]]),&(xc[i][grid->ctop[i]]),grid->Nk[i]-grid->ctop[i]);
      //      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      //	xc[i][k]=x[i][k];

      /*
      if(i==1) {
	printf("i = 1\n");
	for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) 
	  printf("%f %f %f %f %f\n",a[i],b[i],c[i],d[i],xc[i][k]);
      }
      */
    } else 
      xc[i][grid->ctop[i]]=-0.5*x[i][grid->ctop[i]]/coef[i][grid->ctop[i]];
  }
}

/*
 * Function: ConditionQ
 * Usage: ConditionQ(x,grid,phys,prop,myproc,comm);
 * ------------------------------------------------
 * Compute the magnitude of the diagonal elements of the coefficient matrix
 * for the pressure-Poisson equation and place it into x after taking its square root.
 *
 */
static void ConditionQ(REAL **x, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {

  int i, iptr, k, ne, nf, nc, kmin, warn=0;
  REAL *a = phys->a;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	x[i][k]=0;

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {

	  ne = grid->face[i*NFACES+nf];
	  
	  if(grid->ctop[nc]>grid->ctop[i])
	    kmin = grid->ctop[nc];
	  else
	    kmin = grid->ctop[i];
	  
	  for(k=kmin;k<grid->Nke[ne];k++) 
	    x[i][k]+=grid->dzz[i][k]*phys->D[ne];
	}
      
      a[grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
	a[k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1]);

      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	x[i][k]+=(a[k]+a[k+1]);

      if(grid->ctop[i]<grid->Nk[i]-1) {
	// Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
	k=grid->ctop[i];
	x[i][k]+=2*a[k]+a[k+1];

	// Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
	k=grid->Nk[i]-1;
	x[i][k]+=a[k];
      }
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	if(x[i][k]<=0) {
	  x[i][k]=1;
	  warn=1;
	}
	x[i][k]=sqrt(x[i][k]);
      }
  }
  if(WARNING && warn) printf("Warning...invalid preconditioner!\n");

  // Send the preconditioner to the neighboring processors.
  ISendRecvCellData3D(x,grid,myproc,comm);
}

/*
 * Function: GSSolve
 * Usage: GSSolve(grid,phys,prop,myproc,numprocs,comm);
 * ----------------------------------------------------
 * Solve the free surface equation with a Gauss-Seidell relaxation.
 * This function is used for debugging only.
 *
 */
static void GSSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, nf, ne, n, niters, *N;
  REAL *h, *hold, *D, *hsrc, myresid, resid, residold, tmp, relax, myNsrc, Nsrc, coef;

  h = phys->h;
  hold = phys->hold;
  D = phys->D;
  hsrc = phys->htmp;
  N = grid->normal;

  tmp = GRAV*pow(prop->theta*prop->dt,2);

  ISendRecvCellData2D(h,grid,myproc,comm);

  relax = prop->relax;
  niters = prop->maxiters;
  resid=0;
  myresid=0;

  myNsrc=0;
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    myNsrc+=pow(hsrc[i],2);
  }
  MPI_Reduce(&myNsrc,&(Nsrc),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&Nsrc,1,MPI_DOUBLE,0,comm);
  Nsrc=sqrt(Nsrc);

  for(n=0;n<niters;n++) {

    for(i=0;i<grid->Nc;i++) {
      hold[i] = h[i];
    }

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      h[i] = hsrc[i];

      coef=1;
      for(nf=0;nf<NFACES;nf++) 
	if(grid->neigh[i*NFACES+nf]!=-1) {
	  ne = grid->face[i*NFACES+nf];

	  coef+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	  h[i]+=relax*tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
	    phys->h[grid->neigh[i*NFACES+nf]]/grid->Ac[i];
	}
      h[i]/=coef;
    }

    residold=resid;
    myresid=0;
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      hold[i] = hsrc[i];

      coef=1;
      for(nf=0;nf<NFACES;nf++) 
	if(grid->neigh[i*NFACES+nf]!=-1) {
	  ne = grid->face[i*NFACES+nf];
	  coef+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	  hold[i]+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
	    phys->h[grid->neigh[i*NFACES+nf]]/grid->Ac[i];
	}
      myresid+=pow(hold[i]/coef-h[i],2);
    }
    MPI_Reduce(&myresid,&(resid),1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Bcast(&resid,1,MPI_DOUBLE,0,comm);
    resid=sqrt(resid);

    ISendRecvCellData2D(h,grid,myproc,comm);
    MPI_Barrier(comm);

    if(fabs(resid)<prop->epsilon)
      break;
  }
  if(n==niters && myproc==0 && WARNING) 
    printf("Warning... Iteration not converging after %d steps! RES=%e\n",n,resid);
  
  for(i=0;i<grid->Nc;i++)
    if(h[i]!=h[i]) 
      printf("NaN h[%d] in cgsolve!\n",i);
}

/*
 * Function: Continuity
 * Usage: Continuity(w,grid,phys,prop);
 * ------------------------------------
 * Compute the vertical velocity field that satisfies continuity.  Use
 * the upwind flux face heights to ensure consistency with continuity.
 *
 */
static void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop)
{
  int i, k, nf, iptr, ne, nc1, nc2, j, jptr;
  REAL ap, am, d;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->ctop[i];k++) 
      w[i][k] = 0;

    w[i][grid->Nk[i]] = 0;
    d=grid->dzz[i][grid->Nk[i]-1];
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      w[i][k] = w[i][k+1];
      for(nf=0;nf<NFACES;nf++) {
	ne = grid->face[i*NFACES+nf];
	nc1 = grid->grad[2*ne];
	nc2 = grid->grad[2*ne+1];
	if(nc1==-1) nc1=nc2;
	if(nc2==-1) nc2=nc1;

	ap=0;
	if(k<grid->Nk[nc2])
	  ap=0.5*(phys->u[ne][k]+fabs(phys->u[ne][k]))*grid->dzz[nc2][k];

	am=0;
	if(k<grid->Nk[nc1])
	  am=0.5*(phys->u[ne][k]-fabs(phys->u[ne][k]))*grid->dzz[nc1][k];

	w[i][k]-=(ap+am)*grid->df[ne]*grid->normal[i*NFACES+nf]/grid->Ac[i];
      }
    }
  }
}

/*
 * Function: ComputeConservatives
 * Usage: ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);
 * -----------------------------------------------------------------
 * Compute the total mass, volume, and potential energy within the entire
 * domain and return a warning if the mass and volume are not conserved to within
 * the tolerance CONSERVED specified in suntans.h 
 *
 */
static void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs,
			  MPI_Comm comm)
{
  int i, iptr, k;
  REAL mass, volume, volh, height, Ep;

  if(myproc==0) phys->mass=0;
  if(myproc==0) phys->volume=0;
  if(myproc==0) phys->Ep=0;

  // volh is the horizontal integral of h+d, whereas vol is the
  // 3-d integral of dzz
  mass=0;
  volume=0;
  volh=0;
  Ep=0;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    height = 0;
    volh+=grid->Ac[i]*(grid->dv[i]+phys->h[i]);
    Ep+=0.5*GRAV*grid->Ac[i]*(phys->h[i]+grid->dv[i])*(phys->h[i]-grid->dv[i]);
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      height += grid->dzz[i][k];
      volume+=grid->Ac[i]*grid->dzz[i][k];
      mass+=phys->s[i][k]*grid->Ac[i]*grid->dzz[i][k];
    }
  }

  // Comment out the volh reduce if that integral is desired.  The
  // volume integral is used since the volh integral is useful
  // only for debugging.
  MPI_Reduce(&mass,&(phys->mass),1,MPI_DOUBLE,MPI_SUM,0,comm);
  //MPI_Reduce(&volh,&(phys->volume),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&volume,&(phys->volume),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&Ep,&(phys->Ep),1,MPI_DOUBLE,MPI_SUM,0,comm);

  // Compare the quantities to the original values at the beginning of the
  // computation.  If prop->n==0 (beginning of simulation), then store the
  // starting values for comparison.
  if(myproc==0) {
    if(prop->n==0) {
      phys->volume0 = phys->volume;
      phys->mass0 = phys->mass;
      phys->Ep0 = phys->Ep;
    } else {
      if(fabs((phys->volume-phys->volume0)/phys->volume0)>CONSERVED && prop->volcheck)
	printf("Warning! Not volume conservative! V(0)=%e, V(t)=%e\n",
	       phys->volume0,phys->volume);
      if(fabs((phys->mass-phys->mass0)/phys->volume0)>CONSERVED && prop->masscheck) 
	printf("Warning! Not mass conservative! M(0)=%e, M(t)=%e\n",
	       phys->mass0,phys->mass);
    }
  }
}

/*
 * Function: ComputeVelocityVector
 * Usage: ComputeVelocityVector(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  This function estimates the velocity vector with
 *
 * u = 1/Area * Sum_{faces} u_{face} normal_{face} df_{face}/d_{ef,face}
 *
 */
static void ComputeVelocityVector(REAL **u, REAL **uc, REAL **vc, gridT *grid) {
  
  int k, n, ne, nf;
  REAL sum;

  for(n=0;n<grid->Nc;n++) {
    for(k=0;k<grid->Nk[n];k++) {
      uc[n][k]=0;
      vc[n][k]=0;
    }
    for(k=grid->ctop[n];k<grid->Nk[n];k++) {
      for(nf=0;nf<NFACES;nf++) {
	ne = grid->face[n*NFACES+nf];
	uc[n][k]+=u[ne][k]*grid->n1[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
	vc[n][k]+=u[ne][k]*grid->n2[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
      }
      uc[n][k]/=grid->Ac[n];
      vc[n][k]/=grid->Ac[n];
    }
  }
}

/*
 * Function: OutputData
 * Usage: OutputData(grid,phys,prop,myproc,numprocs,blowup,comm);
 * --------------------------------------------------------------
 * Output the data every ntout steps as specified in suntans.dat
 * If this is the last time step or if the run is blowing up (blowup==1),
 * then output the data to the restart file specified by the file pointer
 * prop->StoreFID.
 *
 * If ASCII is specified then the data is output in ascii format, otherwise
 * it is output in binary format.  This variable is specified in suntans.h.
 *
 */
static void OutputData(gridT *grid, physT *phys, propT *prop,
		int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int i, j, k, nwritten;
  REAL *tmp = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"OutputData");

  if(!(prop->n%prop->ntconserve) && !blowup) {
    ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);
    if(myproc==0)
      fprintf(prop->ConserveFID,"%e %e %e %e %e %e %e %e\n",prop->rtime,phys->mass,phys->volume,
	      phys->Ep-phys->Ep0,phys->Eflux1,phys->Eflux2,phys->Eflux3,phys->Eflux4);
  }

  if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) {

    if(myproc==0 && VERBOSE>1) 
      if(!blowup) 
	printf("Outputting data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      else
	printf("Outputting blowup data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);


    if(ASCII) 
      for(i=0;i<grid->Nc;i++)
	fprintf(prop->FreeSurfaceFID,"%f\n",phys->h[i]);
    else {
      nwritten=fwrite(phys->h,sizeof(REAL),grid->Nc,prop->FreeSurfaceFID);
      if(nwritten!=grid->Nc) {
	printf("Error outputting free-surface data!\n");
	exit(EXIT_WRITING);
      }
    }
    fflush(prop->FreeSurfaceFID);

    // ut stores the tangential component of velocity on the faces.
    if(ASCII) 
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++) 
	  fprintf(prop->HorizontalVelocityFID,"%e %e %e\n",
		  phys->uc[i][k],phys->vc[i][k],0.5*(phys->w[i][k]+phys->w[i][k+1]));
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->HorizontalVelocityFID,"0.0 0.0 0.0\n");
      }
    else 
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i]) 
	    tmp[i]=phys->uc[i][k];
	  else
	    tmp[i]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Nc,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting Horizontal Velocity data!\n");
	  exit(EXIT_WRITING);
	}
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    tmp[i]=phys->vc[i][k];
	  else
	    tmp[i]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Nc,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting Horizontal Velocity data!\n");
	  exit(EXIT_WRITING);
	}
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    tmp[i]=0.5*(phys->w[i][k]+phys->w[i][k+1]);
	  else
	    tmp[i]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Nc,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting Face Velocity data!\n");
	  exit(EXIT_WRITING);
	}
      }
    fflush(prop->HorizontalVelocityFID);

    if(ASCII)
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i]+1;k++)
	  fprintf(prop->VerticalVelocityFID,"%e\n",phys->w[i][k]);
	for(k=grid->Nk[i]+1;k<grid->Nkmax+1;k++)
	  fprintf(prop->VerticalVelocityFID,"0.0\n");
      }
    else {
      for(k=0;k<grid->Nkmax+1;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i]+1)
	    phys->htmp[i]=phys->w[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->VerticalVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting vertical velocity data!\n");
	  exit(EXIT_WRITING);
	}
      }
    }
    fflush(prop->VerticalVelocityFID);

    if(ASCII) {
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++)
	  fprintf(prop->SalinityFID,"%e\n",phys->s[i][k]);
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->SalinityFID,"0.0\n");
      }
      if(prop->n==1+prop->nstart) {
	for(i=0;i<grid->Nc;i++) {
	  for(k=0;k<grid->Nk[i];k++)
	    fprintf(prop->BGSalinityFID,"%e\n",phys->s0[i][k]);
	  for(k=grid->Nk[i];k<grid->Nkmax;k++)
	    fprintf(prop->BGSalinityFID,"0.0\n");
	}
      }
    } else {
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    phys->htmp[i]=phys->s[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->SalinityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting salinity data!\n");
	  exit(EXIT_WRITING);
	}
      }
      if(prop->n==1+prop->nstart) {
	for(k=0;k<grid->Nkmax;k++) {
	  for(i=0;i<grid->Nc;i++) {
	    if(k<grid->Nk[i])
	      phys->htmp[i]=phys->s0[i][k];
	    else
	      phys->htmp[i]=EMPTY;
	  }
	  nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->BGSalinityFID);
	  if(nwritten!=grid->Nc) {
	    printf("Error outputting background salinity data!\n");
	    exit(EXIT_WRITING);
	  }
	}
      }
    }
    fflush(prop->SalinityFID);

    if(ASCII) 
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++)
	  fprintf(prop->TemperatureFID,"%e\n",phys->T[i][k]);
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->TemperatureFID,"0.0\n");
      }
    else 
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    phys->htmp[i]=phys->T[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->TemperatureFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting temperature data!\n");
	  exit(EXIT_WRITING);
	}
      }
    fflush(prop->TemperatureFID);

    if(ASCII) 
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++)
	  fprintf(prop->PressureFID,"%e\n",phys->q[i][k]);
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->PressureFID,"0.0\n");
      }
    else 
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    phys->htmp[i]=phys->q[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->PressureFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting pressure data!\n");
	  exit(EXIT_WRITING);
	}
      }
    fflush(prop->PressureFID);

    if(prop->turbmodel) {
      if(ASCII) 
	for(i=0;i<grid->Nc;i++) {
	  for(k=0;k<grid->Nk[i];k++)
	    fprintf(prop->EddyViscosityFID,"%e\n",phys->nu_tv[i][k]);
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->EddyViscosityFID,"0.0\n");
	}
      else 
	for(k=0;k<grid->Nkmax;k++) {
	  for(i=0;i<grid->Nc;i++) {
	    if(k<grid->Nk[i])
	      phys->htmp[i]=phys->nu_tv[i][k];
	    else
	      phys->htmp[i]=EMPTY;
	  }
	  nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->EddyViscosityFID);
	  if(nwritten!=grid->Nc) {
	    printf("Error outputting eddy viscosity data!\n");
	    exit(EXIT_WRITING);
	  }
	}
      fflush(prop->EddyViscosityFID);
      
      if(ASCII) 
	for(i=0;i<grid->Nc;i++) {
	  for(k=0;k<grid->Nk[i];k++)
	    fprintf(prop->ScalarDiffusivityFID,"%e\n",phys->kappa_tv[i][k]);
	  for(k=grid->Nk[i];k<grid->Nkmax;k++)
	    fprintf(prop->ScalarDiffusivityFID,"0.0\n");
	}
      else 
	for(k=0;k<grid->Nkmax;k++) {
	  for(i=0;i<grid->Nc;i++) {
	    if(k<grid->Nk[i])
	      phys->htmp[i]=phys->kappa_tv[i][k];
	    else
	      phys->htmp[i]=EMPTY;
	  }
	  nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->ScalarDiffusivityFID);
	  if(nwritten!=grid->Nc) {
	    printf("Error outputting scalar diffusivity data!\n");
	    exit(EXIT_WRITING);
	  }
	}
      fflush(prop->ScalarDiffusivityFID);
    }
    
    for(i=0;i<grid->Nc;i++) {
      for(k=0;k<grid->Nk[i];k++)
	fprintf(prop->VerticalGridFID,"%e\n",grid->dzz[i][k]);
      for(k=grid->Nk[i];k<grid->Nkmax;k++)
	fprintf(prop->VerticalGridFID,"0.0\n");
    }
  }
  fflush(prop->VerticalGridFID);

  if(prop->n==1)
    fclose(prop->BGSalinityFID);

  if(prop->n==prop->nsteps+prop->nstart) {
    fclose(prop->FreeSurfaceFID);
    fclose(prop->HorizontalVelocityFID);
    fclose(prop->VerticalVelocityFID);
    fclose(prop->SalinityFID);
    fclose(prop->VerticalGridFID);
    if(myproc==0) fclose(prop->ConserveFID);
  }

  if(prop->n==prop->nsteps+prop->nstart || blowup) {
    if(VERBOSE>1 && myproc==0) printf("Writing to rstore...\n");
    
    nwritten=fwrite(&(prop->n),sizeof(int),1,prop->StoreFID);

    fwrite(phys->h,sizeof(REAL),grid->Nc,prop->StoreFID);
    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    if(prop->turbmodel) {
      for(i=0;i<grid->Nc;i++) 
	fwrite(phys->Cn_q[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
      for(i=0;i<grid->Nc;i++) 
	fwrite(phys->Cn_l[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

      for(i=0;i<grid->Nc;i++) 
	fwrite(phys->qT[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
      for(i=0;i<grid->Nc;i++) 
	fwrite(phys->lT[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    }
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->nu_tv[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->kappa_tv[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    fclose(prop->StoreFID);
  }

  SunFree(tmp,grid->Ne*sizeof(REAL),"OutputData");
}

/*
 * Function: ReadProperties
 * Usage: ReadProperties(prop,myproc);
 * -----------------------------------
 * This function reads in the properties specified in the suntans.dat
 * data file.
 *
 */
void ReadProperties(propT **prop, int myproc)
{
  *prop = (propT *)SunMalloc(sizeof(propT),"ReadProperties");
  
  (*prop)->thetaramptime = MPI_GetValue(DATAFILE,"thetaramptime","ReadProperties",myproc);
  (*prop)->theta = MPI_GetValue(DATAFILE,"theta","ReadProperties",myproc);
  (*prop)->thetaS = MPI_GetValue(DATAFILE,"thetaS","ReadProperties",myproc);
  (*prop)->thetaB = MPI_GetValue(DATAFILE,"thetaB","ReadProperties",myproc);
  (*prop)->beta = MPI_GetValue(DATAFILE,"beta","ReadProperties",myproc);
  (*prop)->kappa_s = MPI_GetValue(DATAFILE,"kappa_s","ReadProperties",myproc);
  (*prop)->kappa_sH = MPI_GetValue(DATAFILE,"kappa_sH","ReadProperties",myproc);
  (*prop)->gamma = MPI_GetValue(DATAFILE,"gamma","ReadProperties",myproc);
  (*prop)->kappa_T = MPI_GetValue(DATAFILE,"kappa_T","ReadProperties",myproc);
  (*prop)->kappa_TH = MPI_GetValue(DATAFILE,"kappa_TH","ReadProperties",myproc);
  (*prop)->nu = MPI_GetValue(DATAFILE,"nu","ReadProperties",myproc);
  (*prop)->nu_H = MPI_GetValue(DATAFILE,"nu_H","ReadProperties",myproc);
  (*prop)->tau_T = MPI_GetValue(DATAFILE,"tau_T","ReadProperties",myproc);
  (*prop)->z0T = MPI_GetValue(DATAFILE,"z0T","ReadProperties",myproc);
  (*prop)->z0B = MPI_GetValue(DATAFILE,"z0B","ReadProperties",myproc);
  (*prop)->CdT = MPI_GetValue(DATAFILE,"CdT","ReadProperties",myproc);
  (*prop)->CdB = MPI_GetValue(DATAFILE,"CdB","ReadProperties",myproc);
  (*prop)->CdW = MPI_GetValue(DATAFILE,"CdW","ReadProperties",myproc);
  (*prop)->turbmodel = (int)MPI_GetValue(DATAFILE,"turbmodel","ReadProperties",myproc);
  (*prop)->dt = MPI_GetValue(DATAFILE,"dt","ReadProperties",myproc);
  (*prop)->Cmax = MPI_GetValue(DATAFILE,"Cmax","ReadProperties",myproc);
  (*prop)->nsteps = (int)MPI_GetValue(DATAFILE,"nsteps","ReadProperties",myproc);
  (*prop)->ntout = (int)MPI_GetValue(DATAFILE,"ntout","ReadProperties",myproc);
  (*prop)->ntprog = (int)MPI_GetValue(DATAFILE,"ntprog","ReadProperties",myproc);
  (*prop)->ntconserve = (int)MPI_GetValue(DATAFILE,"ntconserve","ReadProperties",myproc);
  (*prop)->nonhydrostatic = (int)MPI_GetValue(DATAFILE,"nonhydrostatic","ReadProperties",myproc);
  (*prop)->cgsolver = (int)MPI_GetValue(DATAFILE,"cgsolver","ReadProperties",myproc);
  (*prop)->maxiters = (int)MPI_GetValue(DATAFILE,"maxiters","ReadProperties",myproc);
  (*prop)->qmaxiters = (int)MPI_GetValue(DATAFILE,"qmaxiters","ReadProperties",myproc);
  (*prop)->qprecond = (int)MPI_GetValue(DATAFILE,"qprecond","ReadProperties",myproc);
  (*prop)->epsilon = MPI_GetValue(DATAFILE,"epsilon","ReadProperties",myproc);
  (*prop)->qepsilon = MPI_GetValue(DATAFILE,"qepsilon","ReadProperties",myproc);
  (*prop)->resnorm = MPI_GetValue(DATAFILE,"resnorm","ReadProperties",myproc);
  (*prop)->relax = MPI_GetValue(DATAFILE,"relax","ReadProperties",myproc);
  (*prop)->amp = MPI_GetValue(DATAFILE,"amp","ReadProperties",myproc);
  (*prop)->omega = MPI_GetValue(DATAFILE,"omega","ReadProperties",myproc);
  (*prop)->timescale = MPI_GetValue(DATAFILE,"timescale","ReadProperties",myproc);
  (*prop)->flux = MPI_GetValue(DATAFILE,"flux","ReadProperties",myproc);
  (*prop)->volcheck = MPI_GetValue(DATAFILE,"volcheck","ReadProperties",myproc);
  (*prop)->masscheck = MPI_GetValue(DATAFILE,"masscheck","ReadProperties",myproc);
  (*prop)->nonlinear = MPI_GetValue(DATAFILE,"nonlinear","ReadProperties",myproc);
  (*prop)->newcells = MPI_GetValue(DATAFILE,"newcells","ReadProperties",myproc);
  (*prop)->wetdry = MPI_GetValue(DATAFILE,"wetdry","ReadProperties",myproc);
  (*prop)->Coriolis_f = MPI_GetValue(DATAFILE,"Coriolis_f","ReadProperties",myproc);
  (*prop)->sponge_distance = MPI_GetValue(DATAFILE,"sponge_distance","ReadProperties",myproc);
  (*prop)->sponge_decay = MPI_GetValue(DATAFILE,"sponge_decay","ReadProperties",myproc);
  (*prop)->readSalinity = MPI_GetValue(DATAFILE,"readSalinity","ReadProperties",myproc);
  (*prop)->readTemperature = MPI_GetValue(DATAFILE,"readTemperature","ReadProperties",myproc);
}

/* 
 * Function: OpenFiles
 * Usage: OpenFiles(prop,myproc);
 * ------------------------------
 * Open all of the files used for i/o to store the file pointers.
 *
 */
void OpenFiles(propT *prop, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  if(prop->readSalinity) {
    MPI_GetFile(filename,DATAFILE,"InitSalinityFile","OpenFiles",myproc);
    prop->InitSalinityFID = MPI_FOpen(filename,"r","OpenFiles",myproc);
  }
  if(prop->readTemperature) {
    MPI_GetFile(filename,DATAFILE,"InitTemperatureFile","OpenFiles",myproc);
    prop->InitTemperatureFID = MPI_FOpen(filename,"r","OpenFiles",myproc);
  }

  MPI_GetFile(filename,DATAFILE,"FreeSurfaceFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->FreeSurfaceFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"HorizontalVelocityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->HorizontalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"VerticalVelocityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->VerticalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"SalinityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->SalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"BGSalinityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->BGSalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"TemperatureFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->TemperatureFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"PressureFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->PressureFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"EddyViscosityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->EddyViscosityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"ScalarDiffusivityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->ScalarDiffusivityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"VerticalGridFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->VerticalGridFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"StoreFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->StoreFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  if(RESTART) {
    MPI_GetFile(filename,DATAFILE,"StartFile","OpenFiles",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    prop->StartFID = MPI_FOpen(str,"r","OpenFiles",myproc);
  }

  if(myproc==0) {
    MPI_GetFile(filename,DATAFILE,"ConserveFile","OpenFiles",myproc);
    sprintf(str,"%s",filename);
    prop->ConserveFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  }
}

/*
 * Function: InterpToFace
 * Usage: uface = InterpToFace(j,k,phys->uc,u,grid);
 * -------------------------------------------------
 * Linear interpolation of a Voronoi-centered value to the face, using the equation
 * 
 * uface = 1/Dj*(def1*u2 + def2*u1);
 *
 * Note that def1 and def2 are not the same as grid->def[] unless the
 * triangles have not been corrected.
 *
 */
static REAL InterpToFace(int j, int k, REAL **phi, REAL **u, gridT *grid) {
  int nc1, nc2;
  REAL def1, def2, Dj;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  Dj = grid->dg[j];
  def1 = sqrt(pow(grid->xv[nc1]-grid->xe[j],2)+
	      pow(grid->yv[nc1]-grid->ye[j],2));
  def2 = Dj-def1;

  if(def1==0 || def2==0)
    return UpWind(u[j][k],phi[nc1][k],phi[nc2][k]);
  else    
    return (phi[nc1][k]*def2+phi[nc2][k]*def1)/(def1+def2);
}

/*
 * Function: SetDensity
 * Usage: SetDensity(grid,phys,prop);
 * ----------------------------------
 * Sets the values of the density in the density array rho and
 * at the boundaries.
 *
 */
static void SetDensity(gridT *grid, physT *phys, propT *prop) {
  int i, j, k, jptr, ib;
  REAL z, p;

  for(i=0;i<grid->Nc;i++) {
    z=phys->h[i];
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      z+=0.5*grid->dzz[i][k];
      p=RHO0*GRAV*z;
      phys->rho[i][k]=StateEquation(prop,phys->s[i][k],phys->T[i][k],p);
      z+=0.5*grid->dzz[i][k];
    }
  }

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      z=phys->h[ib];
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
	z+=0.5*grid->dzz[ib][k];
	p=RHO0*GRAV*z;
	phys->boundary_rho[jptr-grid->edgedist[2]][k]=
	  StateEquation(prop,phys->boundary_s[jptr-grid->edgedist[2]][k],
			phys->boundary_T[jptr-grid->edgedist[2]][k],p);
	z+=0.5*grid->dzz[ib][k];
      }
  }
}
