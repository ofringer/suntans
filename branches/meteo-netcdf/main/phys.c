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
#include "diffusion.h"
#include "sources.h"
/*
 * Private Function declarations.
 *
 */
static void UpdateDZ(gridT *grid, physT *phys, propT *prop, int option);
static void UPredictor(gridT *grid, physT *phys, 
    propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void Corrector(REAL **qc, gridT *grid, physT *phys, propT *prop, int myproc, 
    int numprocs, MPI_Comm comm);
static void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, 
    int myproc, int numprocs);
static void CGSolve(gridT *grid, physT *phys, propT *prop, 
    int myproc, int numprocs, MPI_Comm comm);
static void HPreconditioner(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop);
static void HCoefficients(REAL *coef, REAL *fcoef, gridT *grid, physT *phys, 
    propT *prop);
static void CGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys, 
    propT *prop, 
    int myproc, int numprocs, MPI_Comm comm);
static void ConditionQ(REAL **x, gridT *grid, physT *phys, propT *prop, int myproc, 
    MPI_Comm comm);
static void Preconditioner(REAL **x, REAL **xc, REAL **coef, gridT *grid, physT *phys, 
    propT *prop);
static void GuessQ(REAL **q, REAL **wold, REAL **w, gridT *grid, physT *phys, 
    propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void GSSolve(gridT *grid, physT *phys, propT *prop, 
    int myproc, int numprocs, MPI_Comm comm);
static REAL InnerProduct(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, 
    MPI_Comm comm);
static REAL InnerProduct3(REAL **x, REAL **y, gridT *grid, int myproc, int numprocs, 
    MPI_Comm comm);
static void OperatorH(REAL *x, REAL *y, REAL *coef, REAL *fcoef, gridT *grid, 
    physT *phys, propT *prop);
static void OperatorQC(REAL **coef, REAL **fcoef, REAL **x, REAL **y, REAL **c, 
    gridT *grid, physT *phys, propT *prop);
static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, 
    physT *phys, propT *prop);
static void OperatorQ(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, 
    physT *phys, propT *prop);
static void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop);
void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop);
static void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, 
    int numprocs, MPI_Comm comm);
static void EddyViscosity(gridT *grid, physT *phys, propT *prop, REAL **wnew, 
    MPI_Comm comm, int myproc);
static void HorizontalSource(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, MPI_Comm comm);
static void StoreVariables(gridT *grid, physT *phys);
static void NewCells(gridT *grid, physT *phys, propT *prop);
static void WPredictor(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, MPI_Comm comm);
inline void ComputeUC(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc, interpolation interp);
static void ComputeUCPerot(REAL **u, REAL **uc, REAL **vc, gridT *grid);
inline static void ComputeUCRT(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc);
static void ComputeNodalVelocity(physT *phys, gridT *grid, interpolation interp, int myproc);
static void  ComputeTangentialVelocity(physT *phys, gridT *grid, interpolation ninterp, interpolation tinterp,int myproc);
static void  ComputeQuadraticInterp(REAL x, REAL y, int ic, int ik, REAL **uc, 
    REAL **vc, physT *phys, gridT *grid, interpolation ninterp, 
    interpolation tinterp, int myproc);
inline static void ComputeRT0Velocity(REAL* tempu, REAL* tempv, REAL e1n1, REAL e1n2, 
    REAL e2n1, REAL e2n2, REAL Uj1, REAL Uj2);
static void BarycentricCoordsFromCartesian(gridT *grid, int cell, 
    REAL x, REAL y, REAL* lambda);
static void BarycentricCoordsFromCartesianEdge(gridT *grid, int cell, 
    REAL x, REAL y, REAL* lambda);
static void OutputData(gridT *grid, physT *phys, propT *prop,
//void OutputData(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, int blowup, MPI_Comm comm);
static REAL UFaceFlux(int j, int k, REAL **phi, REAL **u, gridT *grid, REAL dt, 
    int method);
static REAL HFaceFlux(int j, int k, REAL *phi, REAL **u, gridT *grid, REAL dt, 
    int method);
static void SetDensity(gridT *grid, physT *phys, propT *prop);
//static void SetFluxHeight(gridT *grid, physT *phys, propT *prop);
void SetFluxHeight(gridT *grid, physT *phys, propT *prop);
static void GetMomentumFaceValues(REAL **uface, REAL **ui, REAL **boundary_ui, REAL **U, gridT *grid, physT *phys, propT *prop,
				  MPI_Comm comm, int myproc, int nonlinear);

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
  int flag=0, i, j, jptr, ib, Nc=grid->Nc, Ne=grid->Ne, Np=grid->Np, nf, k;

  // allocate physical structure
  *phys = (physT *)SunMalloc(sizeof(physT),"AllocatePhysicalVariables");

  // allocate  variables in plan
  (*phys)->u = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->uc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");

  // new variables for higher-order interpolation following Wang et al 2011
  (*phys)->nRT1u = (REAL ***)SunMalloc(Np*sizeof(REAL **),"AllocatePhysicalVariables");
  (*phys)->nRT1v = (REAL ***)SunMalloc(Np*sizeof(REAL **),"AllocatePhysicalVariables");
  (*phys)->nRT2u = (REAL **)SunMalloc(Np*sizeof(REAL*),"AllocatePhysicalVariables");
  (*phys)->nRT2v = (REAL **)SunMalloc(Np*sizeof(REAL*),"AllocatePhysicalVariables");
  (*phys)->tRT1 = (REAL **)SunMalloc(Ne*sizeof(REAL*),"AllocatePhysicalVariables");
  (*phys)->tRT2 = (REAL **)SunMalloc(Ne*sizeof(REAL*),"AllocatePhysicalVariables");

  // allocate rest of variables in plan
  (*phys)->uold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->D = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->utmp = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->utmp2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->ut = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_U = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_U2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables"); //AB3


  // for each variable in plan consider the number of layers it affects
  for(j=0;j<Ne;j++) {
    // the following line seems somewhat dubious...
    if(grid->Nkc[j] < grid->Nke[j]) {
      printf("Error!  Nkc(=%d)<Nke(=%d) at edge %d\n",grid->Nkc[j],grid->Nke[j],j);
      flag = 1;
    }
    // allocate memory for the max (cell-centered) quanity on the edge (from definition above)
    (*phys)->u[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->ut[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_U[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_U2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");//AB3
    /* new interpolation variables */
    // loop over the edges (Nkc vs Nke since for cells Nkc < ik < Nke there should be 0 velocity
    // on face to prevent mass from leaving the system)
    (*phys)->tRT1[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->tRT2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
  }
  // if we have an error quit MPI
  if(flag) {
    MPI_Finalize();
    exit(0);
  }

  // cell-centered physical variables in plan (no vertical direction)
  (*phys)->h = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->hcorr = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->active = (unsigned char *)SunMalloc(Nc*sizeof(char),"AllocatePhysicalVariables");
  (*phys)->hold = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp2 = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp3 = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->hcoef = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->hfcoef = (REAL *)SunMalloc(NFACES*Nc*sizeof(REAL),"AllocatePhysicalVariables");

  // cell-centered values that are also depth-varying
  (*phys)->w = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wnew = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wtmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wtmp2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_W = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_W2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables"); //AB3
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
  (*phys)->nu_lax = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  if(prop->turbmodel==1) {
    (*phys)->qT = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->lT = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->Cn_q = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->Cn_l = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  }
  (*phys)->tau_T = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->tau_B = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdT = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdB = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");

  /* new interpolation variables */
  // loop over the nodes
  for(i=0; i < Np; i++) {
    // most complex one...
    (*phys)->nRT1u[i] = (REAL **)SunMalloc(grid->Nkp[i]*sizeof(REAL *),
        "AllocatePhysicalVariables");
    (*phys)->nRT1v[i] = (REAL **)SunMalloc(grid->Nkp[i]*sizeof(REAL *),
        "AllocatePhysicalVariables");
    // loop over all the cell over depth and allocate (note that for rapidly varying
    // bathymetery this will result in too much memory used in some places)
    for(k=0; k < grid->Nkp[i]; k++) {
      (*phys)->nRT1u[i][k] = (REAL *)SunMalloc(grid->numpcneighs[i]*sizeof(REAL),
          "AllocatePhysicalVariables");
      (*phys)->nRT1v[i][k] = (REAL *)SunMalloc(grid->numpcneighs[i]*sizeof(REAL),
          "AllocatePhysicalVariables");
    }
    // simpler one
    (*phys)->nRT2u[i] = (REAL *)SunMalloc(grid->Nkp[i]*sizeof(REAL),
        "AllocatePhysicalVariables");
    (*phys)->nRT2v[i] = (REAL *)SunMalloc(grid->Nkp[i]*sizeof(REAL),
        "AllocatePhysicalVariables");
  }

  // for each cell allocate memory for the number of layers at that location
  for(i=0;i<Nc;i++) {
    (*phys)->uc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->uold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wnew[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp2[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_W[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_W2[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables"); //AB3
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
    if(prop->turbmodel==1) {
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
    (*phys)->nu_lax[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  // allocate boundary value memory
  (*phys)->boundary_u = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_v = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_w = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_s = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_T = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_rho = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_tmp = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_h = (REAL *)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->boundary_flag = (REAL *)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL),"AllocatePhysicalVariables");
  // allocate over vertical layers
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

  // allocate coefficients
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
  int i, j, Nc=grid->Nc, Ne=grid->Ne, Np=grid->Np, nf;

  /* free variables for higher-order interpolation */
  // note that this isn't even currently called!
  // over each node
  for(i=0; i < Np; i++) {
    free(phys->nRT2u[i]);
    free(phys->nRT2v[i]);
    // over each layer
    for(j=0; j < grid->Nkp[i]; j++) {
      // free the memory
      free(phys->nRT1u[i][j]);
      free(phys->nRT1v[i][j]);
    }
  }
  // over each edge
  for(i=0; i < Ne; i++) {
    free(phys->tRT1[i]);
    free(phys->tRT2[i]);
  }

  // free all the arrays over depth for edge-oriented
  for(j=0;j<Ne;j++) {
    free(phys->u[j]);
    free(phys->utmp[j]);
    free(phys->utmp2[j]);
    free(phys->ut[j]);
    free(phys->Cn_U[j]);
    free(phys->Cn_U2[j]); //AB3
  }

  // free all the arrays over depth for cell-oriented
  for(i=0;i<Nc;i++) {
    free(phys->uc[i]);
    free(phys->vc[i]);
    free(phys->wc[i]);
    free(phys->uold[i]);
    free(phys->vold[i]);
    free(phys->w[i]);
    free(phys->wnew[i]);
    free(phys->wtmp[i]);
    free(phys->wtmp2[i]);
    free(phys->Cn_W[i]);
    free(phys->Cn_W2[i]); //AB3
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
    if(prop->turbmodel==1) {
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
    free(phys->nu_lax[i]);
  }

  free(phys->h);
  free(phys->hcorr);
  free(phys->htmp);
  free(phys->htmp2);
  free(phys->htmp3);
  free(phys->hcoef);
  free(phys->hfcoef);
  free(phys->uc);
  free(phys->vc);
  free(phys->wc);
  free(phys->w);
  free(phys->wnew);
  free(phys->wtmp);
  free(phys->wtmp2);
  free(phys->Cn_W);
  free(phys->Cn_W2); //AB3
  free(phys->q);
  free(phys->qc);
  free(phys->qtmp);
  free(phys->s);
  free(phys->T);
  free(phys->s0);
  free(phys->rho);
  free(phys->Cn_R);
  free(phys->Cn_T);
  if(prop->turbmodel==1) {
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
  free(phys->nu_lax);
  free(phys->tau_T);
  free(phys->tau_B);
  free(phys->CdT);
  free(phys->CdB);
  free(phys->u);
  free(phys->D);
  free(phys->utmp);
  free(phys->ut);
  free(phys->Cn_U);
  free(phys->Cn_U2);

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
void ReadPhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {

  int i, j;

  if(VERBOSE>1 && myproc==0) printf("Reading from rstore...\n");
  //fixdzz
  UpdateDZ(grid,phys,prop,-1); 

  if(fread(&(prop->nstart),sizeof(int),1,prop->StartFID) != 1)
    printf("Error reading prop->nstart\n");

  if(fread(phys->h,sizeof(REAL),grid->Nc,prop->StartFID) != grid->Nc)
    printf("Error reading phys->h\n");
  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j])
      printf("Error reading phys->Cn_U[j]\n");
  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->Cn_U2[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j]) //AB3
      printf("Error reading phys->Cn_U2[j]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_W[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_W2[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_W[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_R[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_T[i]\n");

  if(prop->turbmodel==1) {
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->Cn_q[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->Cn_q[i]\n");
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->Cn_l[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->Cn_l[i]\n");

    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->qT[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->qT[i]\n");
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->lT[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->lT[i]\n");
  }
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->nu_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->nu_tv[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->kappa_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->kappa_tv[i]\n");

  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j])
      printf("Error reading phys->u[j]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StartFID) != grid->Nk[i]+1)
      printf("Error reading phys->w[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->q[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->qc[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->qc[i]\n");

  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->s[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->T[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->s0[i]\n");
  fclose(prop->StartFID);

  UpdateDZ(grid,phys,prop, 0);

  // cell centered velocity computed so that this does not 
  // need to be reconsidered 
  ComputeUC(phys->uc, phys->vc, phys,grid, myproc, prop->interp);

  ISendRecvCellData3D(phys->uc,grid,myproc,comm);
  ISendRecvCellData3D(phys->vc,grid,myproc,comm);

  // Set the density from s and T using the equation of state
  SetDensity(grid,phys,prop);
}

/*
 * Function: InitializePhyiscalVariables
 * Usage: InitializePhyiscalVariables(grid,phys,prop,myproc,comm);
 * ---------------------------------------------------------------
 * This function initializes the physical variables by calling
 * the routines defined in the file initialize.c
 *
 */
void InitializePhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm)
{
  int i, j, k, Nc=grid->Nc;
  REAL z, *stmp;

  prop->nstart=0;

  // Need to update the vertical grid and fix any cells in which
  // dzz is too small when h=0.
   UpdateDZ(grid,phys,prop, -1);
 
/*  for(i=0;i<Nc;i++) {
    grid->dv[i]=0;
    for(k=0;k<grid->Nk[i];k++)
      grid->dv[i]+=grid->dzz[i][k];
  }

  REAL mindepth=INFTY;
  for(i=0;i<Nc;i++)
    if(grid->dv[i]<mindepth)
      mindepth=grid->dv[i];
  printf("MINDEPTH=%f\n",mindepth);
  */

  // Initialize the free surface
  for(i=0;i<Nc;i++) {
    phys->h[i]=ReturnFreeSurface(grid->xv[i],grid->yv[i],grid->dv[i]);
    if(phys->h[i]<-grid->dv[i] + DRYCELLHEIGHT) 
      phys->h[i]=-grid->dv[i] + DRYCELLHEIGHT;
  }

  // Need to update the vertical grid after updating the free surface.
  // The 1 indicates that this is the first call to UpdateDZ
  UpdateDZ(grid,phys,prop, 1);

  // initailize variables to 0 (except for filter "pressure")
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
    if(fread(stmp,sizeof(REAL),grid->Nkmax,prop->InitSalinityFID) != grid->Nkmax)
      printf("Error reading stmp first\n");
    fclose(prop->InitSalinityFID);

    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        phys->s[i][k]=stmp[k];
        phys->s0[i][k]=stmp[k];
      }
    SunFree(stmp,grid->Nkmax,"InitializePhysicalVariables");
  } 
  else {
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
    if(fread(stmp,sizeof(REAL),grid->Nkmax,prop->InitTemperatureFID) != grid->Nkmax)
      printf("Error reading stmp second\n");
    fclose(prop->InitTemperatureFID);    

    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        phys->T[i][k]=stmp[k];

    SunFree(stmp,grid->Nkmax,"InitializePhysicalVariables");
  } 
  else {
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
      phys->u[j][k]=ReturnHorizontalVelocity(
          grid->xe[j],grid->ye[j],grid->n1[j],grid->n2[j],z);
      z-=grid->dz[k]/2;
    }
  }

  // Need to compute the velocity vectors at the cell centers based
  // on the initialized velocities at the faces.
  ComputeUC(phys->uc, phys->vc, phys, grid, myproc, prop->interp);
  ComputeUC(phys->uold, phys->vold, phys, grid, myproc, prop->interp);

  // send and receive interprocessor data
  ISendRecvCellData3D(phys->uc,grid,myproc,comm);
  ISendRecvCellData3D(phys->vc,grid,myproc,comm);
  ISendRecvCellData3D(phys->uold,grid,myproc,comm);
  ISendRecvCellData3D(phys->vold,grid,myproc,comm);

  // Determine minimum and maximum salinity
  phys->smin=phys->s[0][0];
  phys->smax=phys->s[0][0];
  // overall cells in plan
  for(i=0;i<grid->Nc;i++) {
    // and in the vertical
    for(k=0;k<grid->Nk[i];k++) {
      if(phys->s[i][k]<phys->smin) phys->smin=phys->s[i][k];      
      if(phys->s[i][k]>phys->smax) phys->smax=phys->s[i][k];      
    }
  }

  // Set the density from s and T using the equation of state 
  SetDensity(grid,phys,prop);

  // Initialize the eddy-viscosity and scalar diffusivity
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++) {
      phys->nu_tv[i][k]=0;
      phys->kappa_tv[i][k]=0;
      phys->nu_lax[i][k]=0;
    }
  }

  if(prop->turbmodel==1) {
    for(i=0;i<grid->Nc;i++) {
      for(k=0;k<grid->Nk[i];k++) {
        phys->qT[i][k]=0;
        phys->lT[i][k]=0;
      }
    }
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
  int i, j;

  if(prop->z0T==0) 
    for(j=0;j<grid->Ne;j++) 
      phys->CdT[j]=prop->CdT;
  else
    for(j=0;j<grid->Ne;j++) 
      phys->CdT[j]=pow(log(0.5*grid->dzf[j][grid->etop[j]]/prop->z0T)/KAPPA_VK,-2);

  if(prop->z0B==0) 
    for(j=0;j<grid->Ne;j++) 
      phys->CdB[j]=prop->CdB;
  else
    for(j=0;j<grid->Ne;j++) 
      phys->CdB[j]=pow(log(0.5*grid->dzf[j][grid->Nke[j]-1]/prop->z0B)/KAPPA_VK,-2);

  for(j=0;j<grid->Ne;j++)
    if(grid->dzf[j][grid->Nke[j]-1]<BUFFERHEIGHT && grid->etop[j]==grid->Nke[j]-1){
      phys->CdB[j]=100;
      //printf("Making CdB a large value due to small cell!\n");
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
 * This function is necessary so that gridding in veritcal can be redone each
 * simulation to account for changes in suntans.dat for Nkmax
 *
 */
void InitializeVerticalGrid(gridT **grid,int myproc)
{
  int i, j, k, Nc=(*grid)->Nc, Ne=(*grid)->Ne;

  // initialize in plan for the face (dzf dzfB) and 
  // cell centered (dzz dzzold) quantities
  (*grid)->stairstep = MPI_GetValue(DATAFILE,"stairstep","InitializeVerticalGrid",myproc);
  (*grid)->fixdzz = MPI_GetValue(DATAFILE,"fixdzz","InitializeVerticalGrid",myproc);
  (*grid)->dzsmall = (REAL)MPI_GetValue(DATAFILE,"dzsmall","InitializeVerticalGrid",myproc);
  (*grid)->smoothbot = (REAL)MPI_GetValue(DATAFILE,"smoothbot","InitializeVerticalGrid",myproc);  

 

  (*grid)->dzf = (REAL **)SunMalloc(Ne*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->dzfB = (REAL *)SunMalloc(Ne*sizeof(REAL),"InitializeVerticalGrid");
  (*grid)->dzz = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->dzzold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->dzbot = (REAL *)SunMalloc(Nc*sizeof(REAL),"InitializeVerticalGrid");


  // initialize over depth for edge-oriented quantities
  for(j=0;j<Ne;j++) 
    (*grid)->dzf[j]=(REAL *)SunMalloc(((*grid)->Nkc[j])*sizeof(REAL),"InitializeVerticalGrid");

  // initialize over depth for cell-centered quantities
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
static void UpdateDZ(gridT *grid, physT *phys, propT *prop, int option)
{
  int i, j, k, ne1, ne2, Nc=grid->Nc, Ne=grid->Ne, flag, nc1, nc2;
  REAL z, dzz1, dzz2;

  // don't need to recompute for linearized FS
  if(prop->linearFS) {
    return;
  }

  // If this is not an initial call then set dzzold to store the old value of dzz
  // and also set the etopold and ctopold pointers to store the top indices of
  // the grid.
  if(!option) {
    for(j=0;j<Ne;j++)
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->ctop[i];k++)
        grid->dzzold[i][k]=0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }
  //fixdzz
  if(option==-1)
    for(i=0;i<Nc;i++)
      phys->h[i]=.0; 



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




  /*
  for(i=0;i<Nc;i++)
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(grid->dzz[i][k]<=grid->dz[k]/50)
	grid->dzz[i][k]=grid->dz[k]/50;
    }
  */

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




  // If this is an initial call set the old values to the new values and
  // Determine the bottom-most flux-face height in the absence of h
  // Check the smallest dzz and set to minimum  
  if(option==-1) {
    for(i=0;i<Nc;i++){
      k=grid->Nk[i]-1;      
      grid->dzbot[i]=grid->dzz[i][k];
      //printf("cell %d dzsmall=%e dzz=%e\n",i,grid->dzsmall*grid->dz[k],grid->dzz[i][k]);
      if(!grid->stairstep && grid->fixdzz )   
        if(grid->dzz[i][k]<grid->dz[k]*grid->dzsmall) {
          //printf("cell %d dzsmall=%e dzz=%e\n",i,grid->dzsmall*grid->dz[k],grid->dzz[i][k]);
          grid->dv[i]+= (grid->dz[k]*grid->dzsmall-grid->dzz[i][k]);	
          //grid->dzz[i][k]=grid->dz[k]*grid->dzsmall;
        }
    }
  }

  if(option==1) {    
    for(j=0;j<Ne;j++) 
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->Nk[i];k++)
        grid->dzzold[i][k]=grid->dzz[i][k];
    }

    for(j=0;j<Ne;j++) {
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      dzz1 = grid->dzz[nc1][grid->Nk[nc1]-1];
      dzz2 = grid->dzz[nc2][grid->Nk[nc2]-1];
      z=0;
      for(k=0;k<grid->Nke[j]-1;k++)
        z-=grid->dz[k];
      if(phys->h[nc1]<z) 
        dzz1 = dzz1-z+phys->h[nc1];
      if(phys->h[nc2]<z) 
        dzz2 = dzz2-z+phys->h[nc2];
      grid->dzfB[j] = Min(dzz1,dzz2);
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
//  printf("DepthfromDZ\n");
//  printf("h[%d] = %f\n",i,phys->h[i]);
//  printf("dzz = %f\n",grid->dzz[i]);
//  printf("ctop= %f\n",grid->ctop[i]);
  if(i==-1) {
    printf("!!Error with pointer => h[-1]!!\n");
    //return NAN;  // not consistent with all C compilers
    return -1;
  }
  else {
    int k;
    REAL z = phys->h[i]-0.5*grid->dzz[i][grid->ctop[i]];
    for(k=grid->ctop[i];k<kind;k++) {
      z-=0.5*grid->dzz[i][k-1];
      z-=0.5*grid->dzz[i][k];
    }
    //  printf("DepthfromDZ done\n");
    return z;
  }
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
  // this make sure that we aren't loosing mass/energy
  ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);

  // Print out memory usage per processor and total memory if this is the first time step
  if(VERBOSE>2) MemoryStats(grid,myproc,numprocs,comm);

  // initialize theta0
  prop->theta0=prop->theta;

  // initialize the timers
  t_start=Timer();
  t_source=t_predictor=t_nonhydro=t_turb=t_transport=t_io=t_comm=t_check=0;

  // Set all boundary values at time t=nstart*dt;
  prop->n=prop->nstart;
  // initialize the time (often used for boundary/initial conditions)
  prop->rtime=prop->nstart*prop->dt;

  // get the boundary velocities (boundaries.c)
  BoundaryVelocities(grid,phys,prop,myproc); 
  // get the openboundary flux (boundaries.c)
  OpenBoundaryFluxes(NULL,phys->u,NULL,grid,phys,prop);
  // get the boundary scalars (boundaries.c)
  BoundaryScalars(grid,phys,prop);
  // get the windstress (boundaries.c)
  WindStress(grid,phys,prop,myproc);
  // set the height of the face bewteen cells to compute the flux
  SetFluxHeight(grid,phys,prop);
  // set the drag coefficients for bottom friction
  SetDragCoefficients(grid,phys,prop);
  // for laxWendroff and central differencing- compute the numerical diffusion 
  // coefficients required for stability
  if(prop->laxWendroff && prop->nonlinear==2) LaxWendroff(grid,phys,prop,myproc,comm);

  // Initialize the Sponge Layer
  InitSponge(grid,myproc);

  // main time loop
  for(n=prop->nstart+1;n<=prop->nsteps+prop->nstart;n++) {

    prop->n = n;
    // compute the runtime 
    prop->rtime = n*prop->dt;

    if(prop->nsteps>0) {

      // Ramp down theta from 1 to the value specified in suntans.dat over
      // the time thetaramptime specified in suntans.dat to damp out transient
      // oscillations
      if(prop->thetaramptime!=0)
        prop->theta=(1-exp(-prop->rtime/prop->thetaramptime))*prop->theta0+
          exp(-prop->rtime/prop->thetaramptime);

      // Store the old velocity and scalar fields
      // Store the old values of s, u, and w into stmp3, utmp2, and wtmp2
      StoreVariables(grid,phys);

      // Compute the horizontal source term phys->utmp which contains the explicit part
      // or the right hand side of the free-surface equation. 
      // begin the timer
      t0=Timer();
      // get the flux height (since free surface is changing) which is stored in dzf
      SetFluxHeight(grid,phys,prop);
      // laxWendroff central differencing
      if(prop->laxWendroff && prop->nonlinear==2) 
        LaxWendroff(grid,phys,prop,myproc,comm);
      // compute the horizontal source terms (like 
      /* 
       * 1) Old nonhydrostatic pressure gradient with theta method
       * 2) Coriolis terms with AB2
       * 3) Baroclinic term with AB2
       * 4) Horizontal and vertical advection of horizontal momentum with AB2
       * 5) Horizontal laminar+turbulent diffusion of horizontal momentum
       */
      HorizontalSource(grid,phys,prop,myproc,numprocs,comm);
      // compute the time required for the source
      t_source+=Timer()-t0;

      // Use the explicit part created in HorizontalSource and solve for the 
      // free-surface
      // and hence compute the predicted or hydrostatic horizontal velocity field.  Then
      // send and receive the free surface interprocessor boundary data 
      // to the neighboring processors.
      // The predicted horizontal velocity is now in phys->u
      t0=Timer();
      // compute U^* and h^* (Eqn 40 and Eqn 31)
      UPredictor(grid,phys,prop,myproc,numprocs,comm);
      ISendRecvCellData2D(phys->h,grid,myproc,comm);
      t_predictor+=Timer()-t0;

      t0=Timer();
      blowup = CheckDZ(grid,phys,prop,myproc,numprocs,comm);
      t_check+=Timer()-t0;

      // apply continuity via Eqn 82
      Continuity(phys->wnew,grid,phys,prop);
      ISendRecvWData(phys->wnew,grid,myproc,comm);

      // Compute the eddy viscosity
      t0=Timer();
      EddyViscosity(grid,phys,prop,phys->wnew,comm,myproc);
      t_turb+=Timer()-t0;

      // Update the salinity only if beta is nonzero in suntans.dat
      if(prop->beta) {
        t0=Timer();
        UpdateScalars(grid,phys,prop,phys->wnew,phys->s,phys->boundary_s,phys->Cn_R,
            prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,
            NULL,NULL,NULL,NULL,0,0,comm,myproc,1,prop->TVDsalt);
        ISendRecvCellData3D(phys->s,grid,myproc,comm);
        t_transport+=Timer()-t0;
      }

      // Update the temperature only if gamma is nonzero in suntans.dat
      if(prop->gamma) {
        t0=Timer();

        HeatSource(phys->wtmp,phys->uold,grid,phys,prop);

        UpdateScalars(grid,phys,prop,phys->wnew,phys->T,phys->boundary_T,phys->Cn_T,
            prop->kappa_T,prop->kappa_TH,phys->kappa_tv,prop->theta,
            phys->uold,phys->wtmp,NULL,NULL,0,0,comm,myproc,0,prop->TVDtemp);
        ISendRecvCellData3D(phys->T,grid,myproc,comm);
        t_transport+=Timer()-t0;
      }

      // Compute vertical momentum and the nonhydrostatic pressure
      t0=Timer();
      if(prop->nonhydrostatic && !blowup) {

        // Predicted vertical velocity field is in phys->w
        WPredictor(grid,phys,prop,myproc,numprocs,comm);

        // Source term for the pressure-Poisson equation is in phys->stmp
        ComputeQSource(phys->stmp,grid,phys,prop,myproc,numprocs);

        // Solve for the nonhydrostatic pressure.  
        // phys->stmp2/qc contains the initial guess
        // phys->stmp contains the source term
        // phys->stmp3 is used for temporary storage
        CGSolveQ(phys->qc,phys->stmp,phys->stmp3,grid,phys,prop,myproc,numprocs,comm);

        // Correct the nonhydrostatic velocity field with the nonhydrostatic pressure
        // correction field phys->stmp2/qc.  This will correct phys->u so that it is now
        // the volume-conserving horizontal velocity field.  
        // phys->w is not corrected since
        // it is obtained via continuity.  
        // Also, update the total nonhydrostatic pressure
        // with the pressure correction. 
        Corrector(phys->qc,grid,phys,prop,myproc,numprocs,comm);

        // Send/recv the horizontal velocity data after it has been corrected.
        ISendRecvEdgeData3D(phys->u,grid,myproc,comm);

        // Send q to the boundary cells now that it has been updated
        ISendRecvCellData3D(phys->q,grid,myproc,comm);
      }
      else if(!(prop->interp == PEROT)) {
        // support quadratic interpolation work
        // Send/recv the horizontal velocity data for use with more complex interpolation
        ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
      }
      t_nonhydro+=Timer()-t0;


      // Send/recv the vertical velocity data 
      Continuity(phys->w,grid,phys,prop);
      ISendRecvWData(phys->w,grid,myproc,comm);

      // Set scalar and wind stress boundary values at new time step n; 
      // (n-1) is old time step.
      // BoundaryVelocities and OpenBoundaryFluxes were called in UPredictor to set the
      // boundary velocities to the new time step values for use in the 
      // free surface calculation.
      BoundaryScalars(grid,phys,prop);
      WindStress(grid,phys,prop,myproc);
      SetDragCoefficients(grid,phys,prop);

      if(prop->beta || prop->gamma)
        SetDensity(grid,phys,prop);

      // u now contains velocity on all edges at the new time step
      ComputeUC(phys->uc, phys->vc, phys,grid, myproc, prop->interp);
      
      // now send interprocessor data
      ISendRecvCellData3D(phys->uc,grid,myproc,comm);
      ISendRecvCellData3D(phys->vc,grid,myproc,comm);
    }

    // Adjust the velocity field in the new cells if the newcells variable is set 
    // to 1 in suntans.dat.  Once this is done, send the interprocessor 
    // u-velocities to the neighboring processors.
    if(prop->newcells) {
      NewCells(grid,phys,prop);
      ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
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
      phys->utmp[j][k]=phys->utmp2[j][k]=phys->u[j][k];
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
  int i, ib, iptr, boundary_index, nf, j, jptr, k, nc, nc1, nc2, ne, 
  k0, kmin, kmax;
  REAL *a, *b, *c, fab1, fab2, fab3, sum, def1, def2, dgf, Cz, tempu; //AB3
  // additions to test divergence averaging for w in momentum calc
  REAL wedge[3], lambda[3], wik;
  int aneigh;

  a = phys->a;
  b = phys->b;
  c = phys->c;

  // get the Adams-Bashforth multi-step integration started
  // fab is 1 for a forward Euler calculation on the first time step,
  // for which Cn_U is 0.  Otherwise, fab=3/2 and Cn_U contains the
  // Adams-Bashforth terms at time step n-1

 // Adams Bashforth coefficients
  if(prop->n==1 || prop->wetdry) {
    fab1=1;
    fab2=fab3=0;

    for(j=0;j<grid->Ne;j++)
      for(k=0;k<grid->Nke[j];k++)
        phys->Cn_U[j][k]=phys->Cn_U2[j][k]=0;
  } else if(prop->n==2) {
    fab1=3.0/2.0;
    fab2=-1.0/2.0;
    fab3=0;
  } else {
    if(prop->AB==2) {
      fab1=3.0/2.0;
      fab2=-1.0/2.0;
      fab3=0;
    } else {
      // AB3:
      fab1=23.0/12.0;
      fab2=-4.0/3.0;
      fab3=5.0/12.0;
    }
  }

  // Set utmp and ut to zero since utmp will store the source term of the
  // horizontal momentum equation
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nke[j];k++) {
      phys->utmp[j][k]=0;
      phys->ut[j][k]=0;
    }
  }

  // Update with old AB term
  // correct velocity based on non-hydrostatic pressure
  // over all computational edges
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    //AB3
    // for each edge over depth
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      // Equation 39: U_j,k^n+1 = U_j,k^* - dt*(qc_G2j,k - qc_G1j,k)/Dj
      // note that EE is used for the first time step
      phys->utmp[j][k]=fab2*phys->Cn_U[j][k]+fab3*phys->Cn_U2[j][k]+phys->u[j][k]
        -prop->dt/grid->dg[j]*(phys->q[nc1][k]-phys->q[nc2][k]);

      phys->Cn_U2[j][k]=phys->Cn_U[j][k];
//      phys->utmp[j][k]=(1-fab)*phys->Cn_U[j][k]+phys->u[j][k]
//        -prop->dt/grid->dg[j]*(phys->q[nc1][k]-phys->q[nc2][k]);
//
      phys->Cn_U[j][k]=0;


      // for(k=grid->etop[j];k<grid->Nke[j];k++) {
      //phys->utmp[j][k]=(1-fab)*phys->Cn_U[j][k]+phys->u[j][k]
      //-prop->dt/grid->dg[j]*(phys->q[nc1][k]-phys->q[nc2][k]);
      //phys->Cn_U[j][k]=0;
    }
  }
  // Add on explicit term to boundary edges (type 4 BCs)
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr]; 

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->utmp[j][k]=(1-fab1)*phys->Cn_U[j][k]+phys->u[j][k];

      phys->Cn_U[j][k]=0;
    }
  }

  // note that the above lines appear to allow use to "flush" Cn_U so
  // that it's values are utilizes and then it is free for additional
  // computations

  // Add on a momentum source to momentum equation
  // currently covers the sponge layer and can be used for Coriolis for the 
  // 2D problem
  MomentumSource(phys->utmp,grid,phys,prop);

  // 3D Coriolis terms
  // note that this uses linear interpolation to the faces from the cell centers
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->Cn_U[j][k]+=prop->dt*prop->Coriolis_f*(
          InterpToFace(j,k,phys->vc,phys->u,grid)*grid->n1[j]-
          InterpToFace(j,k,phys->uc,phys->u,grid)*grid->n2[j]);
  }

  // Baroclinic term
  // over computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(grid->etop[j]<grid->Nke[j]-1) 
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        // this next line seems completely unncessary and computationally wasteful
        k0=grid->etop[j];
        // from the top of the water surface to the bottom edge layer for the given step
        for(k0=Max(grid->ctop[nc1],grid->ctop[nc2]);k0<k;k0++) {
          // for the Cn_U at the particular layer, integrate density gradient over the depth
          // using standard integration to the half cell (extra factor of 1/2 in last term)
          // this is probably only 1st order accurate for the integration
          phys->Cn_U[j][k]-=0.5*prop->grav*prop->dt*
            (phys->rho[nc1][k0]-phys->rho[nc2][k0])*
            (grid->dzz[nc1][k0]+grid->dzz[nc2][k0])/grid->dg[j];
        }
        phys->Cn_U[j][k]-=0.25*prop->grav*prop->dt*
          (phys->rho[nc1][k]-phys->rho[nc2][k])*
          (grid->dzz[nc1][k]+grid->dzz[nc2][k])/grid->dg[j];
      }
  }

  // Set stmp and stmp2 to zero since these are used as temporary variables for advection and
  // diffusion.
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

  // Compute Eulerian advection of momentum (nonlinear!=0)
  if(prop->nonlinear) {

    // Interpolate uc to faces and place into ut
    GetMomentumFaceValues(phys->ut,phys->uc,phys->boundary_u,phys->u,grid,phys,prop,comm,myproc,prop->nonlinear);

    // Conservative method assumes ut is a flux
    if(prop->conserveMomentum)
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
	j=grid->edgep[jptr];

	for(k=grid->etop[j];k<grid->Nke[j];k++)
	  phys->ut[j][k]*=grid->dzf[j][k];
      }

    // Now compute the cell-centered source terms and put them into stmp
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];

      // Store dzz in a since for conservative scheme need to divide by depth (since ut is a flux)
      if(prop->conserveMomentum) {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  a[k]=grid->dzz[i][k];
      } else {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  a[k]=1.0;
      }

      // for each face
      for(nf=0;nf<NFACES;nf++) {

        // get the edge pointer
        ne = grid->face[i*NFACES+nf];

        // for all but the top cell layer
        for(k=grid->ctop[i];k<grid->Nk[i];k++) 
          // this is basically Eqn 50, u-component
          phys->stmp[i][k]+=
            phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/(a[k]*grid->Ac[i]);

        // Top cell is filled with momentum from neighboring cells
	if(prop->conserveMomentum)
	  for(k=grid->etop[ne];k<grid->ctop[i];k++) 
	    phys->stmp[i][grid->ctop[i]]+=
	      phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/(a[grid->ctop[i]]*grid->Ac[i]);
      }
    }

    // Interpolate vc to faces and place into ut
    GetMomentumFaceValues(phys->ut,phys->vc,phys->boundary_v,phys->u,grid,phys,prop,comm,myproc,prop->nonlinear);

    // Conservative method assumes ut is a flux
    if(prop->conserveMomentum)
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
	j=grid->edgep[jptr];

	for(k=grid->etop[j];k<grid->Nke[j];k++)
	  phys->ut[j][k]*=grid->dzf[j][k];
      }

    // Now compute the cell-centered source terms and put them into stmp.
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];

      for(k=0;k<grid->Nk[i];k++) 
        phys->stmp2[i][k]=0;

      // Store dzz in a since for conservative scheme need to divide by depth (since ut is a flux)
      if(prop->conserveMomentum) {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  a[k]=grid->dzz[i][k];
      } else {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  a[k]=1.0;
      }

      for(nf=0;nf<NFACES;nf++) {

        ne = grid->face[i*NFACES+nf];

        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          // Eqn 50, v-component
          phys->stmp2[i][k]+=
            phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/(a[k]*grid->Ac[i]);

        // Top cell is filled with momentum from neighboring cells
	if(prop->conserveMomentum)
	  for(k=grid->etop[ne];k<grid->ctop[i];k++) 
	    phys->stmp2[i][grid->ctop[i]]+=
	      phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/(a[k]*grid->Ac[i]);
      }
    }

    // stmp2 holds the v-component of advection of momentum and
    // note that this could also be sped up by pulling out common
    // factors to reduce redundant calculations

    /* Vertical advection of momentum calc */
    // Only if thetaM<0, otherwise use implicit scheme in UPredictor()

    if(prop->thetaM<0) {
      // Now do vertical advection of momentum
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	i=grid->cellp[iptr];
	switch(prop->nonlinear) {
        case 1:
          for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
            a[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
			(phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
            b[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
			(phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
          }
          break;
        case 2: case 5:
          for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
            a[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
				   grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1]));
            b[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
				   grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1]));
          }
          break;
        case 4:
          for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
            Cz = 2.0*phys->w[i][k]*prop->dt/(grid->dzz[i][k]+grid->dzz[i][k-1]);
            a[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
				   grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1])
				  -0.5*Cz*(phys->uc[i][k-1]-phys->uc[i][k]));
            b[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
				   grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1])
				  -0.5*Cz*(phys->vc[i][k-1]-phys->vc[i][k]));
          }
          break;
        default:
          for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
            a[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
			(phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
            b[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
			(phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
          }
          break;
	}
	
	// Always do first-order upwind in bottom cell if partial stepping is on
	if(prop->stairstep==0) {
	  k = grid->Nk[i]-1;
	  a[k] = 0.5*(
		      (phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
		      (phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
	  b[k] = 0.5*(
		      (phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
		      (phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
	}
	
	a[grid->ctop[i]]=phys->w[i][grid->ctop[i]]*phys->uc[i][grid->ctop[i]];
	b[grid->ctop[i]]=phys->w[i][grid->ctop[i]]*phys->vc[i][grid->ctop[i]];
	a[grid->Nk[i]]=0;
	b[grid->Nk[i]]=0;
	
	for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	  phys->stmp[i][k]+=(a[k]-a[k+1])/grid->dzz[i][k];
	  phys->stmp2[i][k]+=(b[k]-b[k+1])/grid->dzz[i][k];
	}
      }
    } // end of nonlinear computation
  }

  // stmp and stmp2 just store the summed C_H and C_V values for horizontal
  // advection (prior to utilization via Eqn 41 or Eqn 47)

  /* Horizontal diffusion calculations */

  // now compute for no slip regions for type 4 boundary conditions
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++)
  {
    // get index for edge pointers
    j = grid->edgep[jptr];
    //    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth using a ghost cell for the ghost-neighbor
    // on type 4 boundary condition
    for (k=grid->ctop[nc]; k<grid->Nk[nc]; k++){
      phys->stmp[nc][k]  += -2.0*prop->nu_H*(
          phys->boundary_u[boundary_index][k] - phys->uc[nc][k])/grid->dg[j]*
        grid->df[j]/grid->Ac[nc];
      phys->stmp2[nc][k] += -2.0*prop->nu_H*(
          phys->boundary_v[boundary_index][k] - phys->vc[nc][k])/grid->dg[j]*
        grid->df[j]/grid->Ac[nc];
    }
  }

  // Now add on horizontal diffusion to stmp and stmp2
  // for the computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];

    for(k=kmin;k<grid->Nke[j];k++) {
      // Eqn 58 and Eqn 59
      // seems like nu_lax should be distance weighted as in Eqn 60
      a[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->uc[nc2][k]-phys->uc[nc1][k])*grid->df[j]/grid->dg[j];
      b[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->vc[nc2][k]-phys->vc[nc1][k])*grid->df[j]/grid->dg[j];
      phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
      phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
      phys->stmp2[nc1][k]-=b[k]/grid->Ac[nc1];
      phys->stmp2[nc2][k]+=b[k]/grid->Ac[nc2];
    }

    // compute wall-drag for diffusion of momentum BC (only on side walls)
    // Eqn 64 and Eqn 65 and Eqn 58
    for(k=grid->Nke[j];k<grid->Nk[nc1];k++) {
      phys->stmp[nc1][k]+=
        prop->CdW*fabs(phys->uc[nc1][k])*phys->uc[nc1][k]*grid->df[j]/grid->Ac[nc1];
      phys->stmp2[nc1][k]+=
        prop->CdW*fabs(phys->vc[nc1][k])*phys->vc[nc1][k]*grid->df[j]/grid->Ac[nc1];
    }
    for(k=grid->Nke[j];k<grid->Nk[nc2];k++) {
      phys->stmp[nc2][k]+=
        prop->CdW*fabs(phys->uc[nc2][k])*phys->uc[nc2][k]*grid->df[j]/grid->Ac[nc2];
      phys->stmp2[nc2][k]+=
        prop->CdW*fabs(phys->vc[nc2][k])*phys->vc[nc2][k]*grid->df[j]/grid->Ac[nc2];
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

  // type 2 boundary condition (specified flux in)
  for(jptr=grid->edgedist[2];jptr<0*grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    // zero existing calculations for type 2 edge
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

    sum=0;
    for(nf=0;nf<NFACES;nf++) {
      if((nc=grid->neigh[i*NFACES+nf])!=-1) {
        sum+=grid->Ac[nc];
        for(k=grid->ctop[nc2];k<grid->Nk[nc2];k++) {
          // get fluxes from other non-boundary cells
          phys->stmp[i][k]+=grid->Ac[nc]*phys->stmp[nc][k];
          phys->stmp2[i][k]+=grid->Ac[nc]*phys->stmp2[nc][k];
        }
      }
    }
    sum=1/sum;
    for(k=grid->ctop[nc2];k<grid->Nk[nc2];k++) {
      // area-averaged calc
      phys->stmp[i][k]*=sum;
      phys->stmp2[i][k]*=sum;
    }
  }

  // computational cells
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


    // compute momentum advection and diffusion contributions to Cn_U, note
    // the minus sign (why we needed it for the no-slip boundary condition)
    // for each face compute Cn_U performing averaging operation such as in
    // Eqn 41 and Eqn 42 and Eqn 55 etc.
    // the two equations correspond to the two adjacent cells to the edge and 
    // their contributions
    for(k=k0;k<grid->Nk[nc1];k++) 
      phys->Cn_U[j][k]-=def1/dgf
        *prop->dt*(phys->stmp[nc1][k]*grid->n1[j]+phys->stmp2[nc1][k]*grid->n2[j]);

    for(k=k0;k<grid->Nk[nc2];k++) 
      phys->Cn_U[j][k]-=def2/dgf
        *prop->dt*(phys->stmp[nc2][k]*grid->n1[j]+phys->stmp2[nc2][k]*grid->n2[j]);
  }

  // Now add on stmp and stmp2 from the boundaries 
  // for type 3 boundary condition
  for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    k0=grid->ctop[nc1];

    for(nf=0;nf<NFACES;nf++) {
      if((nc2=grid->neigh[nc1*NFACES+nf])!=-1) {
        ne=grid->face[nc1*NFACES+nf];
        for(k=k0;k<grid->Nk[nc1];k++) {
          phys->Cn_U[ne][k]-=
            grid->def[nc1*NFACES+nf]/grid->dg[ne]*
            prop->dt*(
                phys->stmp[nc2][k]*grid->n1[ne]+phys->stmp2[nc2][k]*grid->n2[ne]);
        }
      }
    }
  }

  // note that we now basically have the term dt*F_j,k in Equation 33

  // update utmp 
  // this will complete the adams-bashforth time stepping
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 

    for(k=grid->etop[j];k<grid->Nke[j];k++)
      phys->utmp[j][k]+=fab1*phys->Cn_U[j][k];
  }
//  // check to make sure we don't have a blow-up
//  for(j=0;j<grid->Ne;j++) 
//    for(k=grid->etop[j];k<grid->Nke[j];k++) 
//      if(phys->utmp[j][k]!=phys->utmp[j][k]) {
//        printf("0Error in utmp at j=%d k=%d (U***=nan)\n",j,k);
//        exit(1);
//      }
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

    if(grid->etop[j]<grid->Nke[j]-1 && grid->etop[j]<=grid->etopold[j]) 
      for(k=grid->etop[j];k<=grid->etopold[j];k++)
	phys->u[j][k]=phys->u[j][grid->etopold[j]+1];
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
  int i, ib, iptr, j, jptr, k, ne, nf, nc, nc1, nc2, kmin, boundary_index;
  REAL fab1,fab2,fab3, sum, *a, *b, *c, Cz;

  a = phys->a;
  b = phys->b;
  c = phys->c;

  /* if(prop->n==1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
	phys->Cn_W[i][k]=0;
  } else
    fab=1.5;
  */

 // AB3
  if(prop->n==1) {
    fab1=1;
    fab2=fab3=0;

    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
	phys->Cn_W[i][k]=phys->Cn_W2[i][k]=0;
  } else if(prop->n==2) {
    fab1=3.0/2.0;
    fab2=-1.0/2.0;
    fab3=0;
  } else {
    if(prop->AB==2) {
      fab1=3.0/2.0;
      fab2=-1.0/2.0;
      fab3=0;
    } 
    else {
      fab1=23.0/12.0;
      fab2=-4.0/3.0;
      fab3=5.0/12.0;
    }
  }

  // Add on the nonhydrostatic pressure gradient from the previous time
  // step to compute the source term for the tridiagonal inversion.
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      // phys->wtmp[i][k]=phys->w[i][k]+(1-fab)*phys->Cn_W[i][k];
      //AB3
      phys->wtmp[i][k]=phys->w[i][k] + fab2*phys->Cn_W[i][k] + fab3*phys->Cn_W2[i][k];

      phys->Cn_W2[i][k]=phys->Cn_W[i][k];
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
    
    // First compute w at the cell centers (since w is defined at the faces)
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	phys->wc[i][k]=0.5*(phys->w[i][k]+phys->w[i][k+1]);
    }

    // Interpolate wc to faces and place into ut
    GetMomentumFaceValues(phys->ut,phys->wc,phys->boundary_w,phys->utmp2,grid,phys,prop,comm,myproc,prop->nonlinear);

    if(prop->conserveMomentum)
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
	j=grid->edgep[jptr];

	for(k=grid->etop[j];k<grid->Nke[j];k++)
	  phys->ut[j][k]*=grid->dzf[j][k];
      }

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];

      // For conservative scheme need to divide by depth (since ut is a flux)
      if(prop->conserveMomentum) {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  a[k]=grid->dzz[i][k];
      } else {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  a[k]=1.0;
      }

      for(nf=0;nf<NFACES;nf++) {

        ne = grid->face[i*NFACES+nf];

        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          // this is basically Eqn 50, w-component
          phys->stmp[i][k]+=phys->ut[ne][k]*phys->utmp2[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/(a[k]*grid->Ac[i]);

        // Top cell is filled with momentum from neighboring cells
	if(prop->conserveMomentum)
	  for(k=grid->etop[ne];k<grid->ctop[i];k++) 
	    phys->stmp[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->utmp2[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/(a[k]*grid->Ac[i]);
      }

      // Vertical advection; note that in this formulation first-order upwinding is not implemented.
      if(prop->nonlinear==1 || prop->nonlinear==2 ||prop->nonlinear==5) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          phys->stmp[i][k]+=(pow(phys->w[i][k],2)-pow(phys->w[i][k+1],2))/grid->dzz[i][k];
        }
      }
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
      a[k]=.5*(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->w[nc2][k]-phys->w[nc1][k]+phys->w[nc2][k+1]-phys->w[nc1][k+1])*grid->df[j]/grid->dg[j];
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

  // do the same for type 4 boundary conditions but utilize no-slip boundary condition
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++){
    // get index for edge pointers
    j = grid->edgep[jptr];
    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth
    for (k=grid->ctop[nc]; k<grid->Nke[nc]; k++){
      phys->stmp[nc][k]  += -2.0*prop->nu_H*
        (phys->boundary_w[boundary_index][k] - 0.5*(phys->w[nc][k] + phys->w[nc][k+1]))/grid->dg[j]*
        grid->df[j]/grid->Ac[nc];
    }
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

  // Vertical advection using Lax-Wendroff
  if(prop->nonlinear==4) 
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr]; 

      for(k=grid->ctop[i]+1;k<grid->Nk[i]+1;k++) {
        Cz = 0.5*(phys->w[i][k-1]+phys->w[i][k])*prop->dt/grid->dzz[i][k-1];
        a[k]=0.5*(phys->w[i][k-1]+phys->w[i][k])*(0.5*(phys->w[i][k-1]+phys->w[i][k])-0.5*Cz*(phys->w[i][k-1]-phys->w[i][k]));
      }
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
        phys->Cn_W[i][k]-=2.0*prop->dt*(a[k]-a[k+1])/(grid->dzz[i][k]+grid->dzz[i][k+1]);
      }
    }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->wtmp[i][k]+=fab1*phys->Cn_W[i][k];  //AB3
  }

  // wtmp now contains the right hand side without the vertical diffusion terms.  Now we
  // add the vertical diffusion terms to the explicit side and invert the tridiagonal for
  // vertical diffusion (only if grid->Nk[i]-grid->ctop[i]>=2)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    if(grid->Nk[i]-grid->ctop[i]>1) {
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) { // multiple layers
        a[k] = 2*(prop->nu+prop->laxWendroff_Vertical*phys->nu_lax[i][k-1]+
            phys->nu_tv[i][k-1])/grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
        b[k] = 2*(prop->nu+prop->laxWendroff_Vertical*phys->nu_lax[i][k]+
            phys->nu_tv[i][k])/grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
      }
      b[grid->ctop[i]]=(prop->nu+prop->laxWendroff_Vertical*phys->nu_lax[i][grid->ctop[i]]+
          phys->nu_tv[i][grid->ctop[i]])/pow(grid->dzz[i][grid->ctop[i]],2);
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
    } else { // one layer
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
  // no boundary points are corrected for non-hydrostatic pressure!!!
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
static void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs) 
{

  int i, iptr, j, jptr, k, nf, ne, nc1, nc2;
  REAL *ap=phys->a, *am=phys->b, thetafactor=(1-prop->theta)/prop->theta;

  // for each cell
  for(i=0;i<grid->Nc;i++)
    // initialize to zero
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k] = 0;

  // for each computational cell
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get the cell pointer
    i = grid->cellp[iptr];

    /* VERTICAL CONTRIBUTION */
    // over all cells that are defined to a depth
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      // compute the vertical contributions to the source term
      src[i][k] = grid->Ac[i]*(phys->w[i][k]-phys->w[i][k+1]) 
        + thetafactor*grid->Ac[i]*(phys->wtmp2[i][k]-phys->wtmp2[i][k+1]);
    }

    /* HORIZONTAL CONTRIBUTION */
    // over each face to get the horizontal contributions to the source term
    for(nf=0;nf<NFACES;nf++) {

      // get the edge pointer
      ne = grid->face[i*NFACES+nf];

      // for each of the defined edges over depth
      for(k=grid->ctop[i];k<grid->Nke[ne];k++) 
        // compute the horizontal source term via the (D_H)(u^*)
        src[i][k]+=(phys->u[ne][k]+thetafactor*phys->utmp2[ne][k])*grid->dzf[ne][k]*
          grid->normal[i*NFACES+nf]*grid->df[ne];
    }

    // divide final result by dt
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k]/=prop->dt;
  }

  // D[j] is used in OperatorQ, and it must be zero to ensure no gradient
  // at the hydrostatic faces.
  // this can artificially be used to control the boundary condition
  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=grid->df[j]/grid->dg[j];
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
    if(VERBOSE>3 && myproc==0) printf("CGSolve Pressure Iteration: %d, resid=%e\n",n,sqrt(eps/eps0));
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
 * Usage: EddyViscosity(grid,phys,prop,w,comm,myproc);
 * ---------------------------------------------------
 * This function is used to compute the eddy viscosity, the
 * shear stresses, and the drag coefficients at the upper and lower
 * boundaries.
 *
 */
static void EddyViscosity(gridT *grid, physT *phys, propT *prop, REAL **wnew, MPI_Comm comm, int myproc)
{
  if(prop->turbmodel==1) 
    my25(grid,phys,prop,wnew,phys->qT,phys->lT,phys->Cn_q,phys->Cn_l,phys->nu_tv,phys->kappa_tv,comm,myproc);
}

/*
 * Function: UPredictor 
 * Usage: UPredictor(grid,phys,prop,myproc,numprocs,comm);
 * -------------------------------------------------------
 * Predictor step for the horizontal velocity field.  This function
 * computes the free surface using the theta method and then uses 
 * it to update the predicted
 * velocity field in the absence of the nonhydrostatic pressure.
 *
 * Upon entry, phys->utmp contains the right hand side of the u-momentum equation
 *
 */
static void UPredictor(gridT *grid, physT *phys, 
    propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, j, jptr, ne, nf, nf1, normal, nc1, nc2, k, n0, n1;
  REAL sum, dt=prop->dt, theta=prop->theta, h0, boundary_flag;
  REAL *a, *b, *c, *d, *e1, **E, *a0, *b0, *c0, *d0, theta0, alpha;

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
  // note that phys->utmp is the horizontalsource term computed 
  // in HorizontalSource for 1/2(3F_j,k^n -F_j,k^n-1).

  // phys->u contains the velocity specified at the open boundaries
  // It is also the velocity at time step n. (type 2)
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    // transfer boundary flux velocities onto the horizontal source
    // term (exact as specified)
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]=phys->u[j][k];
  }


  // Update the velocity in the interior nodes with the old free-surface gradient
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // Add the explicit part of the free-surface to create U**.
    // 5th term of Eqn 31
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]-=
        prop->grav*(1-theta)*dt*(phys->h[nc1]-phys->h[nc2])/grid->dg[j];
  }

  // Drag term must be fully implicit
  theta0=theta;
  theta=1;

  // Advection term for vertical momentum.  When alpha=1, first-order upwind,
  // alpha=0 is second-order central.  Always do first-order upwind when doing
  // vertically-implicit momentum advection
  alpha=1;

  // for each of the computational edges
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;

    // Add the wind shear stress from the top cell
    phys->utmp[j][grid->etop[j]]+=2.0*dt*phys->tau_T[j]/
      (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]]);

    // Create the tridiagonal entries and formulate U***
    // provided that we don't have a zero-depth top cell
    if(!(grid->dzz[nc1][grid->etop[j]]==0 && grid->dzz[nc2][grid->etop[j]]==0)) {

      // initialize coefficients
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
            phys->nu_tv[nc1][k]+phys->nu_tv[nc2][k]+
            prop->laxWendroff_Vertical*(phys->nu_lax[nc1][k-1]+phys->nu_lax[nc2][k-1]+
              phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]));

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

      // Coefficients for vertical momentum advection terms
      // d[] stores vertical velocity interpolated to faces vertically half-way between U locations
      // So d[k] contains w defined at the vertical w-location of cell k
      if(prop->nonlinear && prop->thetaM>=0 && grid->Nke[j]-grid->etop[j]>1) {
	if(grid->ctop[nc1]>grid->ctop[nc2]) {
	  n0=nc2;
	  n1=nc1;
	} else {
	  n0=nc1;
	  n1=nc2;
	}
	// Don't do advection on vertical faces without water on both sides.
	for(k=0;k<grid->ctop[n1];k++)
	  d[k]=0;
	for(k=grid->ctop[n1];k<grid->Nke[j];k++)
	  d[k] = 0.5*(phys->w[n0][k]+phys->w[n1][k]);
	d[grid->Nke[j]]=0; // Assume w=0 at a corners (even if w is nonzero on one side of the face)

	for(k=grid->etop[j];k<grid->Nke[j];k++) {
	  a0[k] = (alpha*0.5*(d[k]-fabs(d[k])) + 0.5*(1-alpha)*d[k])/(0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]));
	  b0[k] = (alpha*0.5*(d[k]+fabs(d[k])-d[k+1]+fabs(d[k+1]))+0.5*(1-alpha)*(d[k]-d[k+1]))/(0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]));
	  c0[k] = -(alpha*0.5*(d[k+1]+fabs(d[k+1])) + 0.5*(1-alpha)*d[k+1])/(0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]));
	}
      }

      // add on explicit diffusion to RHS (utmp)
      if(grid->Nke[j]-grid->etop[j]>1) { // more than one vertical layer on edge

        // Explicit part of the viscous term over wetted parts of edge
        // for the interior cells
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++)
          phys->utmp[j][k]+=dt*(1-theta)*(a[k]*phys->u[j][k-1]-
              (a[k]+b[k])*phys->u[j][k]+
              b[k]*phys->u[j][k+1]);

        // Top cell
        // account for no slip conditions which are assumed if CdT = -1 
        if(phys->CdT[j] == -1){ // no slip on top
          phys->utmp[j][grid->etop[j]]+=
            dt*(1-theta)*(a[grid->etop[j]]*-phys->u[j][grid->etop[j]]-
              (a[grid->etop[j]]+b[grid->etop[j]])*phys->u[j][grid->etop[j]]+
              b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);
        }
        else{ // standard drag law code
          phys->utmp[j][grid->etop[j]]+=dt*(1-theta)*(-(b[grid->etop[j]]+2.0*phys->CdT[j]*
                fabs(phys->u[j][grid->etop[j]])/
                (grid->dzz[nc1][grid->etop[j]]+
                 grid->dzz[nc2][grid->etop[j]]))*
              phys->u[j][grid->etop[j]]
              +b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);
        }

        // Bottom cell
        // account for no slip conditions which are assumed if CdB = -1
        if(phys->CdB[j] == -1){ // no slip on bottom
         // some sort of strange error here... in previous code, now fixed 
//          phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdB[j]+phys->CdT[j])/
//            (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
//            fabs(phys->u[j][grid->etop[j]])*phys->u[j][grid->etop[j]];
          phys->utmp[j][grid->Nke[j]-1]+=dt*(1-theta)*(a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2]-
              (a[grid->Nke[j]-1]+b[grid->Nke[j]-1])*phys->u[j][grid->Nke[j]-1]+
              b[grid->Nke[j]-1]*-phys->u[j][grid->Nke[j]-1]);

        }
        else{ // standard drag law code
          phys->utmp[j][grid->Nke[j]-1]+=dt*(1-theta)*(
              a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2] -
              (a[grid->Nke[j]-1] 
               + 2.0*phys->CdB[j]*fabs(phys->u[j][grid->Nke[j]-1])/
               (grid->dzz[nc1][grid->Nke[j]-1]+
                grid->dzz[nc2][grid->Nke[j]-1]))*
              phys->u[j][grid->Nke[j]-1]);
        }
      } 
      else{  // one layer for edge
        // drag on bottom boundary
        if(phys->CdB[j] == -1){ // no slip on bottom
          phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(
              2.0*(2.0*(prop->nu + c[k]))*phys->u[j][grid->etop[j]]/
              ((grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
               (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])));
        }
        else{ // standard drag law formation on bottom
          phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdB[j])/
            (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
            fabs(phys->u[j][grid->etop[j]])*phys->u[j][grid->etop[j]];
        }
        // drag on top boundary
        if(phys->CdT[j] == -1){ // no slip on top
          phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(
              2.0*(2.0*(prop->nu + c[k]))*phys->u[j][grid->etop[j]]/
              ((grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
               (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])));
        }
        else{ // standard drag law formulation on top
          phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdT[j])/
            (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
            fabs(phys->u[j][grid->etop[j]])*phys->u[j][grid->etop[j]];
        }
      }

      // add on explicit vertical momentum advection only if there is more than one vertical layer edge.
      if(prop->nonlinear && prop->thetaM>=0 && grid->Nke[j]-grid->etop[j]>1) {
	for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++)
	  phys->utmp[j][k]-=prop->dt*(1-prop->thetaM)*(a0[k]*phys->u[j][k-1]+b0[k]*phys->u[j][k]+c0[k]*phys->u[j][k+1]);
	
	// Top boundary
	phys->utmp[j][grid->etop[j]]-=prop->dt*(1-prop->thetaM)*((a0[grid->etop[j]]+b0[grid->etop[j]])*phys->u[j][grid->etop[j]]
								 +c0[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);
	
	// Bottom boundary
	phys->utmp[j][grid->Nke[j]-1]-=prop->dt*(1-prop->thetaM)*(a0[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2]
								  +(b0[grid->Nke[j]-1]+c0[grid->Nke[j]-1])*phys->u[j][grid->Nke[j]-1]);
      }

      // Now set up the coefficients for the tridiagonal inversion for the
      // implicit part.  These are given from the arrays above in the discrete operator
      // d^2U/dz^2 = -theta dt a_k U_{k-1} + (1+theta dt (a_k+b_k)) U_k - theta dt b_k U_{k+1}
      // = RHS of utmp

      // Right hand side U** is given by d[k] here.
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        e1[k]=1.0;
        d[k]=phys->utmp[j][k];
      }

      if(grid->Nke[j]-grid->etop[j]>1) { // for more than one vertical layer
        // Top cells
        c[grid->etop[j]]=-theta*dt*b[grid->etop[j]];
        // account for no slip conditions which are assumed if CdT = -1 
        if(phys->CdT[j] == -1){ // no slip
          b[grid->etop[j]]=1.0+theta*dt*(a[grid->etop[j]]+a[grid->etop[j]+1]+b[grid->etop[j]]);
        }
        else{ // standard drag law
          b[grid->etop[j]]=1.0+theta*dt*(b[grid->etop[j]]+
              2.0*phys->CdT[j]*fabs(phys->u[j][grid->etop[j]])/
              (grid->dzz[nc1][grid->etop[j]]+
               grid->dzz[nc2][grid->etop[j]]));
        }
        a[grid->etop[j]]=0;     // set a_1=0 (not used in tridiag solve)

        // Bottom cell
        c[grid->Nke[j]-1]=0;   // set c_N=0 (not used in tridiag solve)
        // account for no slip conditions which are assumed if CdB = -1  
        if(phys->CdB[j] == -1){ // no slip
          b[grid->Nke[j]-1]=1.0+theta*dt*(a[grid->Nke[j]-1]+b[grid->Nke[j]-1]+b[grid->Nke[j]-2]);
        }
        else{ // standard drag law
          b[grid->Nke[j]-1]=1.0+theta*dt*(a[grid->Nke[j]-1]+
              2.0*phys->CdB[j]*fabs(phys->u[j][grid->Nke[j]-1])/
              (grid->dzz[nc1][grid->Nke[j]-1]+
               grid->dzz[nc2][grid->Nke[j]-1]));
        }
        a[grid->Nke[j]-1]=-theta*dt*a[grid->Nke[j]-1];

        // Interior cells
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
          c[k]=-theta*dt*b[k];
          b[k]=1.0+theta*dt*(a[k]+b[k]);
          a[k]=-theta*dt*a[k];
        }
      } else { // for a single vertical layer
        b[grid->etop[j]] = 1.0;
        // account for no slip conditions which are assumed if CdB = -1  
        if(phys->CdB[j] == -1){ // no slip
          b[grid->etop[j]]+=4.0*theta*dt*2.0*(prop->nu+c[k])/
            ((grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
             (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]]));
        }
        else{
          b[grid->etop[j]]+=2.0*theta*dt*fabs(phys->utmp[j][grid->etop[j]])/
            (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
            (phys->CdB[j]);
        }
        // account for no slip conditions which are assumed if CdT = -1 
        if(phys->CdT[j] == -1){
          b[grid->etop[j]]+=4.0*theta*dt*2.0*(prop->nu+c[k])/
            ((grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
             (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]]));
        }
        else{
          b[grid->etop[j]]+=2.0*theta*dt*fabs(phys->utmp[j][grid->etop[j]])/
            (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
            (phys->CdT[j]);
        }

      }	  

      // Now add on implicit terms for vertical momentum advection, only if there is more than one layer
      if(prop->nonlinear && prop->thetaM>=0 && grid->Nke[j]-grid->etop[j]>1) {
	for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
	  a[k]+=prop->dt*prop->thetaM*a0[k];
	  b[k]+=prop->dt*prop->thetaM*b0[k];
	  c[k]+=prop->dt*prop->thetaM*c0[k];
	}

	// Top boundary
	b[grid->etop[j]]+=prop->dt*prop->thetaM*(a0[grid->etop[j]]+b0[grid->etop[j]]);
	c[grid->etop[j]]+=prop->dt*prop->thetaM*c0[grid->etop[j]];
	
	// Bottom boundary 
	a[grid->Nke[j]-1]+=prop->dt*prop->thetaM*a0[grid->Nke[j]-1];
	b[grid->Nke[j]-1]+=prop->dt*prop->thetaM*(b0[grid->Nke[j]-1]+c0[grid->Nke[j]-1]);
      }

      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        if(grid->dzz[nc1][k]==0 && grid->dzz[nc2][k]==0) {
          printf("Exiting because dzz[%d][%d]=%f or dzz[%d][%d]=%f\n",
              nc1,k,grid->dzz[nc1][k],nc2,k,grid->dzz[nc2][k]);
          exit(0);
        }
        if(a[k]!=a[k]) printf("a[%d] problems, dzz[%d][%d]=%f\n",k,j,k,grid->dzz[j][k]);
        if(b[k]!=b[k] || b[k]==0) printf("b[%d] problems, b=%f\n",k,b[k]);
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
      if(grid->Nke[j]-grid->etop[j]>1) { // more than one layer (z level)
        TriSolve(&(a[grid->etop[j]]),&(b[grid->etop[j]]),&(c[grid->etop[j]]),
            &(d[grid->etop[j]]),&(phys->utmp[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);
        TriSolve(&(a0[grid->etop[j]]),&(b0[grid->etop[j]]),&(c0[grid->etop[j]]),
            &(e1[grid->etop[j]]),&(E[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);	
      } else {  // one layer (z level)
        phys->utmp[j][grid->etop[j]]/=b[grid->etop[j]];
        E[j][grid->etop[j]]=1.0/b[grid->etop[j]];
      }

      // Now vertically integrate E to create the vertically integrated flux-face
      // values that comprise the coefficients of the free-surface solver.  This
      // will create the D vector, where D=DZ^T E (which should be given by the
      // depth when there is no viscosity.
      phys->D[j]=0;
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
        phys->D[j]+=E[j][k]*grid->dzf[j][k];
    }
  }
  theta=theta0;

  for(j=0;j<grid->Ne;j++) 
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      if(phys->utmp[j][k]!=phys->utmp[j][k]) {
        printf("Error in function Predictor at j=%d k=%d (U***=nan)\n",j,k);
        exit(1);
      }

  // So far we have U*** and D.  Now we need to create h* in htmp.   This
  // will comprise the source term for the free-surface solver.  Before we
  // do this we need to set the new velocity at the open boundary faces and
  // place them into utmp.  
  BoundaryVelocities(grid,phys,prop,myproc);
  OpenBoundaryFluxes(NULL,phys->utmp,NULL,grid,phys,prop);

  // for computational cells
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    sum = 0;
    for(nf=0;nf<NFACES;nf++) {

      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];

      for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
        sum+=((1-theta)*phys->u[ne][k]+theta*phys->utmp[ne][k])*
          grid->df[ne]*normal*grid->dzf[ne][k];
    }
    phys->htmp[i]=grid->Ac[i]*phys->h[i]-dt*sum;
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
    // Gauss-Siedel solver
    GSSolve(grid,phys,prop,myproc,numprocs,comm);
  else if(prop->cgsolver==1)
    // Conjugate-gradient solver
    CGSolve(grid,phys,prop,myproc,numprocs,comm);

  // correct cells drying below DRYCELLHEIGHT above the 
  // bathymetry
  for(i=0;i<grid->Nc;i++)
    if(phys->h[i]<=-grid->dv[i]+DRYCELLHEIGHT) {
      //phys->hcorr[i]=-grid->dv[i]+DRYCELLHEIGHT-phys->h[i];
      phys->h[i]=-grid->dv[i]+DRYCELLHEIGHT;
      phys->active[i]=0;
      phys->s[i][grid->Nk[i]-1]=0;
    } else {
      phys->hcorr[i]=0;
      phys->active[i]=1;
    }

  // Add back the implicit barotropic term to obtain the 
  // hydrostatic horizontal velocity field.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k]=phys->utmp[j][k]-prop->grav*theta*dt*E[j][k]*
        (phys->h[nc1]-phys->h[nc2])/grid->dg[j];

    // set dry cells (with zero height) to have zero velocity
    if(grid->etop[j]==grid->Nke[j]-1 && grid->dzz[nc1][grid->etop[j]]==0 &&
        grid->dzz[nc2][grid->etop[j]]==0) 
      phys->u[j][grid->etop[j]]=0;
  }

  // Now update the vertical grid spacing with the new free surface.
  // can comment this out to linearize the free surface 
  UpdateDZ(grid,phys,prop, 0);

  // Use the new free surface to add the implicit part of the free-surface
  // pressure gradient to the horizontal momentum.
  //
  // This was removed because it violates the assumption of linearity in that
  // the discretization only knows about grid cells below grid->ctopold[].
  /*
     for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
     j = grid->edgep[jptr];

     nc1 = grid->grad[2*j];
     nc2 = grid->grad[2*j+1];

     if(grid->etop[j]>grid->etopold[j]) 
     for(k=0;k<grid->etop[j];k++)
     phys->u[j][k]=0;
     else 
     for(k=grid->etop[j];k<grid->etopold[j];k++)
     phys->u[j][k]=phys->utmp[j][k]-prop->grav*theta*dt*
     (phys->h[nc1]-phys->h[nc2])/grid->dg[j];
     }
     */

  // Set the flux values at the open boundary (marker=2).  These
  // were set to utmp previously in OpenBoundaryFluxes.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k] = phys->utmp[j][k];
  }

  // Now set the fluxes at the free-surface boundary by assuming dw/dz = 0
  // this is for type 3 boundary conditions
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
static void CGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr, n, niters;
  REAL *x, *r, *rtmp, *p, *z, mu, nu, eps, eps0, alpha, alpha0;

  x = phys->h;
  r = phys->hold;
  rtmp = phys->htmp2;
  z = phys->htmp3;
  p = phys->htmp;

  niters = prop->maxiters;

  // Create the coefficients for the operator
  HCoefficients(phys->hcoef,phys->hfcoef,grid,phys,prop);

  // For the boundary term (marker of type 3):
  // 1) Need to set x to zero in the interior points, but
  //    leave it as is for the boundary points.
  // 2) Then set z=Ax and substract b = b-z so that
  //    the new problem is Ax=b with the boundary values
  //    on the right hand side acting as forcing terms.
  // 3) After b=b-z for the interior points, then need to
  //    set b=0 for the boundary points.

  /* Fix to account for boundary cells (type 3) in h */
  // 1) x=0 interior cells 
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    x[i]=0;
  }
  ISendRecvCellData2D(x,grid,myproc,comm);
  OperatorH(x,z,phys->hcoef,phys->hfcoef,grid,phys,prop);

  // 2) b = b-z
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    p[i] = p[i] - z[i];    
    r[i] = p[i];
    x[i] = 0;
  }    
  // 3) b=0 for the boundary cells
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) { 
    i = grid->cellp[iptr]; 

    p[i] = 0; 
  }     

  // continue with CG as expected now that boundaries are handled
  if(prop->hprecond==1) {
    HPreconditioner(r,rtmp,grid,phys,prop);
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      p[i] = rtmp[i];
    }
    alpha = alpha0 = InnerProduct(r,rtmp,grid,myproc,numprocs,comm);
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      p[i] = r[i];
    }
    alpha = alpha0 = InnerProduct(r,r,grid,myproc,numprocs,comm);
  }
  if(!prop->resnorm) alpha0 = 1;

  if(prop->hprecond==1)
    eps=eps0=InnerProduct(r,r,grid,myproc,numprocs,comm);
  else
    eps=eps0=alpha0;

  // Iterate until residual is less than prop->epsilon
  for(n=0;n<niters && eps!=0 && alpha!=0;n++) {

    ISendRecvCellData2D(p,grid,myproc,comm);
    OperatorH(p,z,phys->hcoef,phys->hfcoef,grid,phys,prop);

    mu = 1/alpha;
    nu = alpha/InnerProduct(p,z,grid,myproc,numprocs,comm);

    //    printf("mu = %f, nu = %f, IP = %f\n",mu,nu,InnerProduct(p,z,grid,myproc,numprocs,comm));
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      x[i] += nu*p[i];
      r[i] -= nu*z[i];
    }
    if(prop->hprecond==1) {
      HPreconditioner(r,rtmp,grid,phys,prop);
      alpha = InnerProduct(r,rtmp,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];

        p[i] = rtmp[i] + mu*p[i];
      }
    } else {
      alpha = InnerProduct(r,r,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];

        p[i] = r[i] + mu*p[i];
      }
    }

    if(prop->hprecond==1)
      eps=InnerProduct(r,r,grid,myproc,numprocs,comm);
    else
      eps=alpha;

    if(VERBOSE>3 && myproc==0) printf("CGSolve free-surface Iteration: %d, resid=%e\n",n,sqrt(eps/eps0));
    if(sqrt(eps/eps0)<prop->epsilon) 
      break;
  }
  if(myproc==0 && VERBOSE>2) 
    if(eps==0)
      printf("Warning...Time step %d, norm of free-surface source is 0.\n",prop->n);
    else
      if(n==niters)  printf("Warning... Time step %d, Free-surface iteration not converging after %d steps! RES=%e > %.2e\n",
          prop->n,n,sqrt(eps/eps0),prop->qepsilon);
      else printf("Time step %d, CGSolve free-surface converged after %d iterations, res=%e < %.2e\n",
          prop->n,n,sqrt(eps/eps0),prop->epsilon);

  // Send the solution to the neighboring processors
  ISendRecvCellData2D(x,grid,myproc,comm);
}

/*
 * Function: HPreconditioner
 * Usage: HPreconditioner(r,rtmp,grid,phys,prop);
 * ----------------------------------------------
 * Multiply the vector x by the inverse of the preconditioner M with
 * xc = M^{-1} x
 *
 */
static void HPreconditioner(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop) {
  int i, iptr;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    y[i]=x[i]/phys->hcoef[i];
  }
}

/*
 * Function: HCoefficients
 * Usage: HCoefficients(coef,fcoef,grid,phys,prop);
 * --------------------------------------------------
 * Compute coefficients for the free-surface solver.  fcoef stores
 * coefficients at the flux faces while coef stores coefficients
 * at the cell center.  If L is the linear operator on x, then
 *
 * L(x(i)) = coef(i)*x(i) - sum(m=1:3) fcoef(ne)*x(neigh)
 * coef(i) = (Ac(i) + sum(m=1:3) tmp*D(ne)*df(ne)/dg(ne))
 * fcoef(ne) = tmp*D(ne)*df(ne)/dg(ne)
 *
 * where tmp = prop->grav*(theta*dt)^2
 *
 */
static void HCoefficients(REAL *coef, REAL *fcoef, gridT *grid, physT *phys, propT *prop) {

  int i, j, iptr, jptr, ne, nf;
  REAL tmp = prop->grav*pow(prop->theta*prop->dt,2), h0, boundary_flag;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    coef[i] = grid->Ac[i];
    for(nf=0;nf<NFACES;nf++) 
      if(grid->neigh[i*NFACES+nf]!=-1) {
        ne = grid->face[i*NFACES+nf];

        fcoef[i*NFACES+nf]=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne];
        coef[i]+=fcoef[i*NFACES+nf];
      }
  }
}

/*
 *
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
 * Poisson equation and places it into y with y = L(x), where
 *
 * L(x(i)) = coef(i)*x(i) + sum(m=1:3) fcoef(ne)*x(neigh)
 * coef(i) = (Ac(i) + sum(m=1:3) tmp*D(ne)*df(ne)/dg(ne))
 * fcoef(ne) = tmp*D(ne)*df(ne)/dg(ne)
 *
 * where tmp = prop->grav*(theta*dt)^2
 *
 */
static void OperatorH(REAL *x, REAL *y, REAL *coef, REAL *fcoef, gridT *grid, physT *phys, propT *prop) {

  int i, j, iptr, jptr, ne, nf;
  REAL tmp = prop->grav*pow(prop->theta*prop->dt,2), h0, boundary_flag;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    y[i] = coef[i]*x[i];
    for(nf=0;nf<NFACES;nf++) 
      if(grid->neigh[i*NFACES+nf]!=-1)
        y[i]-=fcoef[i*NFACES+nf]*x[grid->neigh[i*NFACES+nf]];
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

  // sum over all computational cells
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    // over the depth of defined cells
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      y[i][k]=-x[i][k];

    // over each face
    for(nf=0;nf<NFACES;nf++) {
      if((nc=grid->neigh[i*NFACES+nf])!=-1) {

        ne = grid->face[i*NFACES+nf];

        if(grid->ctop[nc]>grid->ctop[i])
          kmin = grid->ctop[nc];
        else
          kmin = grid->ctop[i];

        // now sum over depth too
        for(k=kmin;k<grid->Nke[ne];k++) 
          y[i][k]+=x[nc][k]*fcoef[i*NFACES+nf][k];
      }
    }

    // over depth now
    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
      y[i][k]+=coef[i][k]*x[i][k-1]+coef[i][k+1]*x[i][k+1];

    // apply top/bottom bcs in vertical since we have no flux on sides
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
static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, 
    physT *phys, propT *prop) {

  int i, iptr, k, kmin, nf, nc, ne;

  // if we want to use the preconditioner
  if(prop->qprecond==1) {
    // over all computational cells
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]]/c[i][grid->ctop[i]];
      // over all the compuational cell depths
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
        // compute the cell coefficients
        coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1])/(c[i][k]*c[i][k-1]);

      // over all the faces 
      for(nf=0;nf<NFACES;nf++) 
        if((nc=grid->neigh[i*NFACES+nf])!=-1) {

          ne = grid->face[i*NFACES+nf];

          if(grid->ctop[nc]>grid->ctop[i])
            kmin = grid->ctop[nc];
          else
            kmin = grid->ctop[i];

          // over the edge depths
          for(k=kmin;k<grid->Nke[ne];k++) 
            // compute the face coefficients
            fcoef[i*NFACES+nf][k]=grid->dzz[i][k]*phys->D[ne]/(c[i][k]*c[nc][k]);
        }
    }
  }
  // no preconditioner
  else {
    /*for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++)  */
    // over each and every cell
    for(i=0;i<grid->Nc;i++) {
      //      i = grid->cellp[iptr];

      // compute the cell coeficient for the top cell
      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]];
      // compute the coefficients for the rest of the cells (towards bottom)
      // where dzz is averaged to the face 
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
        coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
    }
  }
  // where is fcoef solved for???
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

  // over each computational cell
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    // over cells that exist and aren't cut off by bathymetry
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      y[i][k]=0;

    // over each face
    for(nf=0;nf<NFACES;nf++) 
      // we only apply this to non-boundary cells 
      // (meaning phys->D implicitly 0 at this point)
      if((nc=grid->neigh[i*NFACES+nf])!=-1) {

        ne = grid->face[i*NFACES+nf];

        // determine the minimal k based on each cell over edge
        if(grid->ctop[nc]>grid->ctop[i])
          kmin = grid->ctop[nc];
        else
          kmin = grid->ctop[i];

        // create summed contributions for cell along the vertical
        for(k=kmin;k<grid->Nke[ne];k++) 
          y[i][k]+=(x[nc][k]-x[i][k])*grid->dzf[ne][k]*phys->D[ne];
      }
    // otherwise y[i][...] = 0 (boundary cell with no gradient in horiz.)

    // over all the depth (z dependence on Laplacian)
    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
      // we are adding on the vertical diffusional part 
      y[i][k]+=coef[i][k]*x[i][k-1]-(coef[i][k]+coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

    // apply top and bottom BCs in vertical, no flux on sides
    if(grid->ctop[i]<grid->Nk[i]-1) {
      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      y[i][k]+=(-2*coef[i][k]-coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      y[i][k]+=coef[i][k]*x[i][k-1]-coef[i][k]*x[i][k];
    } 
    else
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
// classic GS following nearly directly from Oliver's notes
static void GSSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, nf, ne, n, niters, *N;
  REAL *h, *hold, *D, *hsrc, myresid, resid, residold, tmp, relax, myNsrc, Nsrc, coef;

  h = phys->h;
  hold = phys->hold;
  D = phys->D;
  hsrc = phys->htmp;
  N = grid->normal;

  tmp = prop->grav*pow(prop->theta*prop->dt,2);

  // each processor must have all the boundary information 
  // for the starting values of h (as in Lh)
  ISendRecvCellData2D(h,grid,myproc,comm);

  relax = prop->relax;
  niters = prop->maxiters;
  resid=0;
  myresid=0;

  // debugging code since it's not used below
  //    myNsrc=0;
  //
  //    // for interior calculation cells 
  //    // this seems like debugging code
  //    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
  //      i = grid->cellp[iptr];
  //
  //      myNsrc+=pow(hsrc[i],2);
  //    }
  //    MPI_Reduce(&myNsrc,&(Nsrc),1,MPI_DOUBLE,MPI_SUM,0,comm);
  //    MPI_Bcast(&Nsrc,1,MPI_DOUBLE,0,comm);
  //    Nsrc=sqrt(Nsrc);

  for(n=0;n<niters;n++) {

    // hold = h;
    for(i=0;i<grid->Nc;i++) {
      hold[i] = h[i];
    }

    // for all the computational cells (since the boundary cells are 
    // already set in celldist[1] to celldist[...]
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      // get teh cell pointer
      i = grid->cellp[iptr];

      // right hand side (h_src)
      h[i] = hsrc[i];

      coef=1;
      for(nf=0;nf<NFACES;nf++) 
        if(grid->neigh[i*NFACES+nf]!=-1) {
          ne = grid->face[i*NFACES+nf];

          // coef is the diagonal coefficient term
          coef+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
          h[i]+=relax*tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
            phys->h[grid->neigh[i*NFACES+nf]]/grid->Ac[i];
        }
      // divide by diagonal coefficient term
      h[i]/=coef;
    }

    // now need to compare against the residual term to 
    // determine when GS has converged
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
      // compute the residual
      myresid+=pow(hold[i]/coef-h[i],2);
    }
    MPI_Reduce(&myresid,&(resid),1,MPI_DOUBLE,MPI_SUM,0,comm);
    // - is this line necessary?
    MPI_Bcast(&resid,1,MPI_DOUBLE,0,comm);
    resid=sqrt(resid);

    // send all final results out to each processor now
    ISendRecvCellData2D(h,grid,myproc,comm);
    // wait for communication to finish
    MPI_Barrier(comm);

    // if we have met the tolerance criteria
    if(fabs(resid)<prop->epsilon)
      break;
  }
  if(n==niters && myproc==0 && WARNING) 
    printf("Warning... Iteration not converging after %d steps! RES=%e\n",n,resid);

  for(i=0;i<grid->Nc;i++)
    if(h[i]!=h[i]) 
      printf("NaN h[%d] in gssolve!\n",i);
}

/*
 * Function: Continuity
 * Usage: Continuity(w,grid,phys,prop);
 * ------------------------------------
 * Compute the vertical velocity field that satisfies continuity.  Use
 * the upwind flux face heights to ensure consistency with continuity.
 *
 */
void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop)
//static void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop)
{
  int i, k, nf, iptr, ne, nc1, nc2, j, jptr;
  REAL ap, am, dzfnew, theta=prop->theta;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->Nk[i]+1;k++)
      w[i][k] = 0;

    // Continuity is written in terms of flux-face heights at time level n
    // Even though w is updated to grid->ctop[i], it only uses dzf which is
    // 0 for new cells (i.e. if grid->ctop[i]<grid->ctopold[i]).
    // no flux into bottom cells
    w[i][grid->Nk[i]] = 0;
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      // wtmp2 is the old value for w (basically at timestep n) where
      // w in this context is for n+1
      w[i][k] = w[i][k+1] - (1-theta)/theta*(phys->wtmp2[i][k]-phys->wtmp2[i][k+1]);
      for(nf=0;nf<NFACES;nf++) {
        ne = grid->face[i*NFACES+nf];

        if(k<grid->Nke[ne])
          w[i][k]-=(theta*phys->u[ne][k]+(1-theta)*phys->utmp2[ne][k])*
            grid->df[ne]*grid->normal[i*NFACES+nf]/grid->Ac[i]/theta*grid->dzf[ne][k];
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
    Ep+=0.5*prop->grav*grid->Ac[i]*(phys->h[i]+grid->dv[i])*(phys->h[i]-grid->dv[i]);
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
 * Function: ComputeUCPerot
 * Usage: ComputeUCPerot(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  This function estimates the velocity vector with
 *
 * u = 1/Area * Sum_{faces} u_{face} normal_{face} df_{face}*d_{ef,face}
 *
 */
static void ComputeUCPerot(REAL **u, REAL **uc, REAL **vc, gridT *grid) {

  int k, n, ne, nf, iptr;
  REAL sum;

  // for each computational cell (non-stage defined)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get cell pointer transfering from boundary coordinates 
    // to grid coordinates
    n=grid->cellp[iptr];

    // initialize over all depths
    for(k=0;k<grid->Nk[n];k++) {
      uc[n][k]=0;
      vc[n][k]=0;
    }
    // over the entire depth (cell depth)
    for(k=grid->ctop[n];k<grid->Nk[n];k++) {
      // over each face
      for(nf=0;nf<NFACES;nf++) {
        ne = grid->face[n*NFACES+nf];
        if(!(grid->smoothbot) || k<grid->Nke[ne]){
          uc[n][k]+=u[ne][k]*grid->n1[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
          vc[n][k]+=u[ne][k]*grid->n2[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
        }
        else{	
          uc[n][k]+=u[ne][grid->Nke[ne]-1]*grid->n1[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
          vc[n][k]+=u[ne][grid->Nke[ne]-1]*grid->n2[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
        }
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


// this was done for debugging purposes and should go back to the way 
// it was for the final release to prevent aliasing of the function
// with grid.c 
//static void OutputData(gridT *grid, physT *phys, propT *prop,
static void OutputData(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int i, j, jptr, k, nwritten;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
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

    // compute quadratic interpolated estimates for velocity using pretty plot and don't redo work
    if(prop->prettyplot==1 && prop->interp != QUAD) {
      ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
      ComputeUC(phys->uc, phys->vc, phys,grid, myproc, QUAD);
    }

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

  // probably should change to make a distinction between blowup and restarts
  if(!(prop->n%prop->ntoutStore) || blowup) {
    if(VERBOSE>1 && myproc==0) 
      printf("Outputting restart data at step %d\n",prop->n);

    MPI_GetFile(filename,DATAFILE,"StoreFile","OutputData",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    prop->StoreFID = MPI_FOpen(str,"w","OpenFiles",myproc);

    nwritten=fwrite(&(prop->n),sizeof(int),1,prop->StoreFID);

    fwrite(phys->h,sizeof(REAL),grid->Nc,prop->StoreFID);
    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->Cn_U2[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_W2[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    if(prop->turbmodel==1) {
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
      fwrite(phys->qc[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

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
 * data file.  Note that if an entry does not exist, a default can be used.
 *
 */
void ReadProperties(propT **prop, int myproc)
{
  // allocate memory
  *prop = (propT *)SunMalloc(sizeof(propT),"ReadProperties");

  // set values from suntans.dat file (DATAFILE)
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
  (*prop)->grav= MPI_GetValue(DATAFILE,"grav","ReadProperties",myproc);
  (*prop)->turbmodel = (int)MPI_GetValue(DATAFILE,"turbmodel","ReadProperties",myproc);
  (*prop)->dt = MPI_GetValue(DATAFILE,"dt","ReadProperties",myproc);
  (*prop)->Cmax = MPI_GetValue(DATAFILE,"Cmax","ReadProperties",myproc);
  (*prop)->nsteps = (int)MPI_GetValue(DATAFILE,"nsteps","ReadProperties",myproc);
  (*prop)->ntout = (int)MPI_GetValue(DATAFILE,"ntout","ReadProperties",myproc);
  (*prop)->ntoutStore = (int)MPI_GetValue(DATAFILE,"ntoutStore","ReadProperties",myproc);

  if((*prop)->ntoutStore==0)
    (*prop)->ntoutStore=(*prop)->nsteps;

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
  (*prop)->wetdry = MPI_GetValue(DATAFILE,"wetdry","ReadProperties",myproc);
  (*prop)->Coriolis_f = MPI_GetValue(DATAFILE,"Coriolis_f","ReadProperties",myproc);
  (*prop)->sponge_distance = MPI_GetValue(DATAFILE,"sponge_distance","ReadProperties",myproc);
  (*prop)->sponge_decay = MPI_GetValue(DATAFILE,"sponge_decay","ReadProperties",myproc);
  (*prop)->readSalinity = MPI_GetValue(DATAFILE,"readSalinity","ReadProperties",myproc);
  (*prop)->readTemperature = MPI_GetValue(DATAFILE,"readTemperature","ReadProperties",myproc);
  (*prop)->TVDsalt = MPI_GetValue(DATAFILE,"TVDsalt","ReadProperties",myproc);
  (*prop)->TVDtemp = MPI_GetValue(DATAFILE,"TVDtemp","ReadProperties",myproc);
  (*prop)->TVDturb = MPI_GetValue(DATAFILE,"TVDturb","ReadProperties",myproc);
  (*prop)->stairstep = MPI_GetValue(DATAFILE,"stairstep","ReadProperties",myproc);
  (*prop)->AB = MPI_GetValue(DATAFILE,"AB","ReadProperties",myproc); //AB3
  (*prop)->TVDmomentum = MPI_GetValue(DATAFILE,"TVDmomentum","ReadProperties",myproc); 
  (*prop)->conserveMomentum = MPI_GetValue(DATAFILE,"conserveMomentum","ReadProperties",myproc); 
  (*prop)->thetaM = MPI_GetValue(DATAFILE,"thetaM","ReadProperties",myproc); 
  (*prop)->wetdry = MPI_GetValue(DATAFILE,"wetdry","ReadProperties",myproc); 

  // When wetting and drying is desired:
  // -Do nonconservative momentum advection (conserveMomentum=0)
  // -Use backward Euler for vertical advection of horizontal momentum (thetaM=1)
  // -Update new cells (newcells=1)
  if((*prop)->wetdry) {
    (*prop)->conserveMomentum = 0;
    (*prop)->thetaM = 1;
    (*prop)->newcells = 1;
  }

  if((*prop)->nonlinear==2) {
    (*prop)->laxWendroff = MPI_GetValue(DATAFILE,"laxWendroff","ReadProperties",myproc);
    if((*prop)->laxWendroff!=0)
      (*prop)->laxWendroff_Vertical = MPI_GetValue(DATAFILE,"laxWendroff_Vertical","ReadProperties",myproc);
    else
      (*prop)->laxWendroff_Vertical = 0;
  } else {
    (*prop)->laxWendroff = 0;
    (*prop)->laxWendroff_Vertical = 0;
  }

  (*prop)->hprecond = MPI_GetValue(DATAFILE,"hprecond","ReadProperties",myproc);

  // addition for interpolation methods
  switch((int)MPI_GetValue(DATAFILE,"interp","ReadProperties",myproc)) {
    case 0: //Perot
      (*prop)->interp = PEROT;
      break;
    case 1: //Quad
      (*prop)->interp = QUAD;
      break;
    default:
      printf("ERROR: Specification of interpolation type is incorrect!\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
      break;
  }

  // additional data for pretty plot methods
  (*prop)->prettyplot = MPI_GetValue(DATAFILE,"prettyplot","ReadProperties",myproc);

  // addition for linearized free surface where dzz=dz
  (*prop)->linearFS = (int)MPI_GetValue(DATAFILE,"linearFS","ReadProperties",myproc);

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
 * triangles have not been corrected.  This affects the Coriolis term as currently 
 * implemented.
 *
 */
REAL InterpToFace(int j, int k, REAL **phi, REAL **u, gridT *grid) {
  int nc1, nc2;
  REAL def1, def2, Dj;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  Dj = grid->dg[j];
  def1 = sqrt(pow(grid->xv[nc1]-grid->xe[j],2)+
      pow(grid->yv[nc1]-grid->ye[j],2));
  def2 = Dj-def1;

  if(def1==0 || def2==0) {
    return UpWind(u[j][k],phi[nc1][k],phi[nc2][k]);
  }
  else {
    return (phi[nc1][k]*def2+phi[nc2][k]*def1)/(def1+def2);
  }
}

/*
 * Function: UFaceFlux
 * Usage: UFaceFlux(j,k,phi,phys->u,grid,prop->dt,prop->nonlinear);
 * ---------------------------------------------------------------------
 * Interpolation to obtain the flux of the scalar field phi (of type REAL **) on 
 * face j, k;  method==2: Central-differencing, method==4: Lax-Wendroff.
 *
 */
// note that we should not compute def1 and def2 in this function as they don't change
// these should be computed in grid.c and just looked up for efficiency
static REAL UFaceFlux(int j, int k, REAL **phi, REAL **u, gridT *grid, REAL dt, int method) {
  int nc1, nc2;
  REAL def1, def2, Dj, C=0;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  if(nc1==-1) nc1=nc2;
  if(nc2==-1) nc2=nc1;
  Dj = grid->dg[j];
  def1 = sqrt(pow(grid->xv[nc1]-grid->xe[j],2)+
      pow(grid->yv[nc1]-grid->ye[j],2));
  def2 = Dj-def1;

  if(method==4) C = u[j][k]*dt/Dj;

  if(def1==0 || def2==0) {
    // this happens on a boundary cell
    return UpWind(u[j][k],phi[nc1][k],phi[nc2][k]);
  }
  else {
    // on an interior cell (orthogonal cell) we can just distance interpolation the 
    // value to the face)
    // basically Eqn 51
    return (phi[nc1][k]*def2+phi[nc2][k]*def1)/(def1+def2)
      -C/2*(phi[nc1][k]-phi[nc2][k]);
  }
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
      p=RHO0*prop->grav*z;
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
      p=RHO0*prop->grav*z;
      phys->boundary_rho[jptr-grid->edgedist[2]][k]=
        StateEquation(prop,phys->boundary_s[jptr-grid->edgedist[2]][k],
            phys->boundary_T[jptr-grid->edgedist[2]][k],p);
      z+=0.5*grid->dzz[ib][k];
    }
  }
}

/*
 * Function: SetFluxHeight
 * Usage: SetFluxHeight(grid,phys,prop);
 * -------------------------------------
 * Set the value of the flux height dzf at time step n for use
 * in continuity and scalar transport.
 *
 */
void SetFluxHeight(gridT *grid, physT *phys, propT *prop) {
  //static void SetFluxHeight(gridT *grid, physT *phys, propT *prop) {
  int i, j, k, nc1, nc2;
  REAL dz_bottom, dzsmall=grid->dzsmall;

  //  // assuming upwinding
  //  for(j=0;j<grid->Ne;j++) {
  //    nc1 = grid->grad[2*j];
  //    nc2 = grid->grad[2*j+1];
  //    if(nc1==-1) nc1=nc2;
  //    if(nc2==-1) nc2=nc1;
  //
  //    for(k=0;k<grid->etop[j];k++)
  
  //      grid->dzf[j][k]=0;
  //    for(k=grid->etop[j];k<grid->Nke[j];k++) 
  //      grid->dzf[j][k]=UpWind(phys->u[j][k],grid->dzz[nc1][k],grid->dzz[nc2][k]);
  //
  //    k=grid->Nke[j]-1;
  //    /* This works with Wet/dry but not with cylinder case...*/
  //    if(grid->etop[j]==k)
  //      grid->dzf[j][k]=Max(0,UpWind(phys->u[j][k],phys->h[nc1],phys->h[nc2])+Min(grid->dv[nc1],grid->dv[nc2]));
  //    else 
  //      grid->dzf[j][k]=Min(grid->dzz[nc1][k],grid->dzz[nc2][k]);
  //    /*
  //       if(grid->Nk[nc1]!=grid->Nk[nc2])
  //       dz_bottom = Min(grid->dzz[nc1][k],grid->dzz[nc2][k]);
  //       else
  //       dz_bottom = 0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
  //       if(k==grid->etop[j])
  //       grid->dzf[j][k]=Max(0,dz_bottom + UpWind(phys->u[j][k],phys->h[nc1],phys->h[nc2]));
  //       else
  //       grid->dzf[j][k]=dz_bottom;
  //     */
  //
  //    for(k=grid->etop[j];k<grid->Nke[j];k++) 
  //      if(grid->dzf[j][k]<=DRYCELLHEIGHT)
  //        grid->dzf[j][k]=0;
  //  }
  // initialize dzf
  for(j=0;j<grid->Ne;j++) {
    for(k=0; k<grid->Nkc[j];k++)
      grid->dzf[j][k]=0;
  }

  // assuming central differencing
  if(grid->smoothbot)
    for(i=0;i<grid->Nc;i++) {
      //if(grid->dzz[i][grid->Nk[i]-1]>grid->dzbot[i])
      //printf("i=%d dzz=%e, dzbot=%e dzsmall=%e\n",i,grid->dzz[i][grid->Nk[i]-1],grid->dzbot[i],dzsmall);
      grid->dzz[i][grid->Nk[i]-1]=Max(grid->dzbot[i],grid->smoothbot*grid->dz[grid->Nk[i]-1]); 

    }


  for(j=0;j<grid->Ne;j++) {
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    for(k=0;k<grid->etop[j];k++)
      grid->dzf[j][k]=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
//      grid->dzf[j][k]=UFaceFlux(j,k, grid->dzz, phys->u, grid,prop->dt,prop->nonlinear);
      grid->dzf[j][k]=UpWind(phys->u[j][k],grid->dzz[nc1][k],grid->dzz[nc2][k]);

    k=grid->Nke[j]-1;
    /* This works with Wet/dry but not with cylinder case...*/
    if(grid->etop[j]==k) {
//      grid->dzf[j][k]=Max(0,HFaceFlux(j,k,phys->h,phys->u,grid,prop->dt,prop->nonlinear)
//          +Min(grid->dv[nc1],grid->dv[nc2]));
      grid->dzf[j][k]=Max(0,UpWind(phys->u[j][k],phys->h[nc1],phys->h[nc2])+Min(grid->dv[nc1],grid->dv[nc2]));
    }
    else 
      grid->dzf[j][k]=Min(grid->dzz[nc1][k],grid->dzz[nc2][k]);
    /*
       if(grid->Nk[nc1]!=grid->Nk[nc2])
       dz_bottom = Min(grid->dzz[nc1][k],grid->dzz[nc2][k]);
       else
       dz_bottom = 0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
       if(k==grid->etop[j])
       grid->dzf[j][k]=Max(0,dz_bottom + UpWind(phys->u[j][k],phys->h[nc1],phys->h[nc2]));
       else
       grid->dzf[j][k]=dz_bottom;
       */

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      if(grid->dzf[j][k]<=DRYCELLHEIGHT)
        grid->dzf[j][k]=0;
  }


  //set minimum dzz
  if(grid->smoothbot)
    for(i=0;i<grid->Nc;i++)
      grid->dzz[i][grid->Nk[i]-1]=Max(grid->dzz[i][grid->Nk[i]-1],dzsmall*grid->dz[grid->Nk[i]-1]);
}




/*
 * Function: ComputeUC
 * Usage: ComputeUC(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  
 *
 */
//inline static void ComputeUC(physT *phys, gridT *grid, int myproc) {
inline void ComputeUC(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc, interpolation interp) {

  switch(interp) {
    case QUAD:
      // using Wang et al 2011 methods
      ComputeUCRT(ui, vi, phys,grid, myproc);
      break;
    case PEROT:
      ComputeUCPerot(phys->u,ui,vi,grid);
      break;
  }

}

/*
 * Function: ComputeUCRT
 * Usage: ComputeUCRT(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  This function estimates the velocity vector with
 * methods outlined in Wang et al, 2011
 *
 */
inline static void ComputeUCRT(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc) {

  int k, n, ne, nf, iptr;
  REAL sum;

  // first we need to reconstruct the nodal velocities using the RT0 basis functions
  //  if(myproc==0) printf("ComputeNodalVelocity\n");
  ComputeNodalVelocity(phys, grid, nRT2, myproc);

  // now we can get the tangential velocities from these results
  //  if(myproc==0) printf("ComputeTangentialVelocity\n");
  ComputeTangentialVelocity(phys, grid, nRT2, tRT2, myproc);

  // for each computational cell (non-stage defined)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get cell pointer transfering from boundary coordinates 
    // to grid coordinates
    n=grid->cellp[iptr];

    // initialize over all depths
    for(k=0;k<grid->Nk[n];k++) {
      phys->uc[n][k]=0;
      phys->vc[n][k]=0;
    }

    // over the entire depth (cell depth)
    // possible mistake here- make sure to the get the right value for the start
    // of the loop (etop?)
    //    if(myproc==0) printf("ComputeQuadraticInterp\n");
    for(k=grid->ctop[n];k<grid->Nk[n];k++) {
      // now we can compute the quadratic interpolated velocity from these results
      ComputeQuadraticInterp(grid->xv[n], grid->yv[n], n, k, ui, 
          vi, phys, grid, nRT2, tRT2, myproc);
    }
  } 
}

/*
 * Function: ComputeQuadraticInterp
 * Usage: ComputeQuadraticInterp(U, V, un, grid, ninterp, tinterp, myproc)
 * -------------------------------------------
 * Compute the quadratic interpolation of the velocity based on
 * a choice for tinterp.
 * Two options presently exist based on choice of tinterp:
 *  1. tRT1
 *  2. tRT2
 * which are outlined in Wang et al, 2011.
 *
 */
static void  ComputeQuadraticInterp(REAL x, REAL y, int ic, int ik, REAL **uc, 
    REAL **vc, physT *phys, gridT *grid, interpolation ninterp, 
    interpolation tinterp, int myproc) {
  // this formulation considers a specific cell ic and a specific level ik
  // for computation of uc and vc

  // we need to have several pieces for this to work: 
  // we need to be able to quickly compute the area subsets in the triangle
  // A12, A13, A23; we need nodal velocity vectors (nRT2), and tangential
  // velocity vectors
  REAL points[2][NUMEDGECOLUMNS];
  REAL xt[NUMEDGECOLUMNS], yt[NUMEDGECOLUMNS];
  REAL SubArea[NUMEDGECOLUMNS];
  REAL nu[NUMEDGECOLUMNS], nv[NUMEDGECOLUMNS];
  REAL eu[NUMEDGECOLUMNS], ev[NUMEDGECOLUMNS];
  int np[NUMEDGECOLUMNS], ne[NUMEDGECOLUMNS], nf, ie, ip;
  const REAL TotalArea = grid->Ac[ic];

  // loop over each of the faces 
  // possible error if NUMEDGECOLUMNS != NFACES!!
  for( nf=0; nf < NUMEDGECOLUMNS; nf ++) {
    // get the index for the vertex of the cell
    np[nf] = grid->cells[NUMEDGECOLUMNS*ic + nf];
    // get the index for the cell edge
    ne[nf] = grid->face[NFACES*ic + nf];

    // get the coordinates for the points for the vertex
    // to compute the area of the subtriangle
    points[0][nf] = grid->xp[np[nf]];
    points[1][nf] = grid->yp[np[nf]];
  }
  // now we have the points for each node: points[x/y][nf=0,1,2]
  // and also the edges via grid->face[NFACES*ic + nf=0,1,2] for
  // each cell ic

  // now all that remains is to store the velocities and then
  // get A12, A13, A23 and we can 
  // use the formula after getting the vector values for the 
  // velocities at each of the points

  for( nf=0; nf < NUMEDGECOLUMNS; nf++) {
    // get the points for each vertex of the triangle
    xt[0] = points[0][nf];
    yt[0] = points[1][nf];
    // wrap the next one if it goes out of bounds
    xt[1] = points[0][(nf+1)%NUMEDGECOLUMNS]; 
    yt[1] = points[1][(nf+1)%NUMEDGECOLUMNS];
    // the interpolated center
    xt[2] = x;
    yt[2] = y;

    // get the area from the points where 
    // nf=0 => A12, nf=1 => A23, nf=2 => A13 from Wang et al 2011
    // note that this is normalized
    SubArea[nf] = GetArea(xt, yt, NUMEDGECOLUMNS)/TotalArea;
    //          printf("nf = %d Area=%e\n", nf, SubArea[nf]);
  }

  // get the nodal and edge velocities projected unto global x,y coords
  for (nf=0; nf < NUMEDGECOLUMNS; nf++) {
    // node and edge indexes
    ip = np[nf]; ie = ne[nf];
    // get the nodal velocity components 
    nu[nf] = phys->nRT2u[ip][ik];
    nv[nf] = phys->nRT2v[ip][ik];
    // get the tangential velocity components
    if(tinterp == tRT2) {
      eu[nf] = phys->u[ie][ik]*grid->n1[ie] + phys->tRT2[ie][ik]*grid->n2[ie];
      ev[nf] = phys->u[ie][ik]*grid->n2[ie] - phys->tRT2[ie][ik]*grid->n1[ie];
    }
    else if(tinterp == tRT1) {
      // need to have the specific cell neighbor here!  Current implementation 
      // is not correct.  This seems like the easiest way to make sure 
      // that this works
      eu[nf] = phys->u[ie][ik]*grid->n1[ie] + phys->tRT1[ie][ik]*grid->n2[ie];
      ev[nf] = phys->u[ie][ik]*grid->n2[ie] - phys->tRT1[ie][ik]*grid->n1[ie];
    }
  }

  // now perform interpolation from the results via Wang et al 2011 eq 12
  uc[ic][ik] = 
    (2*SubArea[1]-1)*SubArea[1]*nu[0] 
    + (2*SubArea[2]-1)*SubArea[2]*nu[1]
    + (2*SubArea[0]-1)*SubArea[0]*nu[2]
    + 4*SubArea[2]*SubArea[1]*eu[0]
    + 4*SubArea[0]*SubArea[2]*eu[1]
    + 4*SubArea[0]*SubArea[1]*eu[2];
  vc[ic][ik] = 
    (2*SubArea[1]-1)*SubArea[1]*nv[0] 
    + (2*SubArea[2]-1)*SubArea[2]*nv[1]
    + (2*SubArea[0]-1)*SubArea[0]*nv[2]
    + 4*SubArea[2]*SubArea[1]*ev[0]
    + 4*SubArea[0]*SubArea[2]*ev[1]
    + 4*SubArea[0]*SubArea[1]*ev[2];

}

/*
 * Function: ComputeTangentialVelocity
 * Usage: ComputeTangentialVelocity(phys, grid, ninterp, tinterp, myproc)
 * -------------------------------------------
 * Compute the tangential velocity from nodal velocities (interp determines type)
 * Two options presently exist based on choice of tinterp:
 *  1. tRT1
 *  2. tRT2
 * which are outlined in Wang et al, 2011.
 *
 */
static void  ComputeTangentialVelocity(physT *phys, gridT *grid, 
    interpolation ninterp, interpolation tinterp, int myproc) {
  int ie, in, ic, tn, ink;
  int nodes[2], cells[2];
  REAL tempu, tempv, tempA; 
  int tempnode, tempcell;


  if(tinterp == tRT2) {
    // for each edge, we must compute the value of its tangential velocity 
    for(ie = 0; ie < grid->Ne; ie++) {

      // get the nodes on either side of the edge
      nodes[0] = grid->edges[NUMEDGECOLUMNS*ie];
      nodes[1] = grid->edges[NUMEDGECOLUMNS*ie+1];
      //     printf("nodes = %d %d\n", nodes[0], nodes[1]);

      // for each layer at the edge
      for (ink = 0; ink < grid->Nkc[ie]; ink++) {
        // initialize temporary values
        tempu = tempv = tempA = 0;

        // for each node
        for (in = 0; in < 2; in++) {
          // get particular value for node 
          tempnode = nodes[in];

          // ensure that the node exists at the level 
          // before continuing
          if(ink < grid->Nkp[tempnode]) {
            // accumulate the velocity area-weighted values and 
            // cell area
            tempA += grid->Actotal[tempnode][ink];
            tempu += grid->Actotal[tempnode][ink]*phys->nRT2u[tempnode][ink];
            tempv += grid->Actotal[tempnode][ink]*phys->nRT2v[tempnode][ink];
          }
        }
        // area-average existing nodal velocities
        if(tempA == 0) {
          phys->tRT2[ie][ink] = 0;
        }
        else {
          tempu /= tempA;
          tempv /= tempA;
          // get the correct component in the tangential direction to store it
          phys->tRT2[ie][ink] = grid->n2[ie]*tempu - grid->n1[ie]*tempv;
        }
      }

    }

  }
}

/*
 * Function: ComputeNodalVelocity
 * Usage: ComputeNodalVelocity(phys, grid, interp, myproc)
 * -------------------------------------------
 * Compute the nodal velocity using RT0 basis functions.  
 * Two options presently exist based on choice of interp:
 *  1. nRT1
 *  2. nRT2
 * which are outlined in Wang et al, 2011.
 *
 */
static void ComputeNodalVelocity(physT *phys, gridT *grid, interpolation interp, int myproc) {
  //  int in, ink, e1, e2, cell, cp1, cp2;
  int in, ink, inpc, ie, intemp, cell, cp1, cp2,
      e1, e2, n1, n2, onode;
  REAL tempu, tempv, Atemp, tempAu, tempAv;


  /* compute the nodal velocity for RT1 elements */
  // for each node compute its nodal value
  for(in = 0; in < grid->Np; in++) {/*{{{*/
    // there are Nkp vertical layers for each node, so much compute over each 
    // of these

    for(ink = 0; ink < grid->Nkp[in]; ink++) {
      // there will be numpcneighs values for each node at a
      //particular layer so we must compute each separately

      if(interp == nRT2) { Atemp = tempAu = tempAv = 0; }
      for(inpc = 0; inpc < grid->numpcneighs[in]; inpc++) {

        // check to make sure that the cell under consideration exists at the given z-level
        if (ink < grid->Nk[grid->pcneighs[in][inpc]]) { 
          // if so get neighbors to the node
          e1 = grid->peneighs[in][2*inpc];
          e2 = grid->peneighs[in][2*inpc+1];

          // compute the RT0 reconstructed nodal value for the edges
          ComputeRT0Velocity(&tempu, &tempv, 
              grid->n1[e1], grid->n2[e1], grid->n1[e2], grid->n2[e2], phys->u[e1][ink], phys->u[e2][ink]);

          // store the computed values
          phys->nRT1u[in][ink][inpc] = tempu;
          phys->nRT1v[in][ink][inpc] = tempv;
          if(interp == nRT2) {
            // accumulate area and weighted values
            Atemp += grid->Ac[grid->pcneighs[in][inpc]];
            tempAu += grid->Ac[grid->pcneighs[in][inpc]]*tempu;
            tempAv += grid->Ac[grid->pcneighs[in][inpc]]*tempv;
          }

        } 
        else { // cell neighbor doesn't exist at this level
          phys->nRT1u[in][ink][inpc] = 0;
          phys->nRT1v[in][ink][inpc] = 0;
        }
      }
      // compute nRT2 from nRT1 values if desired  (area-weighted average of all cells)
      if(interp == nRT2) {
        if(Atemp == 0) {
          printf("Error as Atemp is 0 in nodal calc!! at in,ink,Nkp=%d,%d,%d\n",in, ink, grid->Nkp[in]);
          // print all nodal neighbors
          printf("cell neighbors = ");
          for(inpc = 0; inpc < grid->numpcneighs[in]; inpc++) {
            printf(" %d(%d)", grid->pcneighs[in][inpc], grid->Nk[grid->pcneighs[in][inpc]]);

          }
          printf("\n");

        }
        /* now we can compute the area-averaged values with nRT1 elements */
        phys->nRT2u[in][ink] = tempAu/Atemp;
        phys->nRT2v[in][ink] = tempAv/Atemp;
      }

    }

  }

}

/*
 * Function: ComputeRT0Velocity
 * Usage: ComputeRT0Velocity(REAL* tempu, REAL* tempv, int e1, int e2, physT* phys);
 * -------------------------------------------
 * Compute the nodal velocity using RT0 basis functions (Appendix B, Wang et al 2011)
 *
 */
inline static void ComputeRT0Velocity(REAL *tempu, REAL *tempv, REAL e1n1, REAL e1n2, 
    REAL e2n1, REAL e2n2, REAL Uj1, REAL Uj2) 
{
  const REAL det = e1n1*e2n2 - e1n2*e2n1;
  // compute using the basis functions directly with inverted 2x2 matrix
  *tempu = (e2n2*Uj1 - e1n2*Uj2)/det;
  *tempv = (e1n1*Uj2 - e2n1*Uj1)/det;
}

/*
 * Function: HFaceFlux
 * Usage: HFaceFlux(j,k,phi,phys->u,grid,prop->dt,prop->nonlinear);
 * ---------------------------------------------------------------------
 * Interpolation to obtain the flux of the scalar field phi (of type REAL **) on 
 * face j, k;  method==2: Central-differencing, method==4: Lax-Wendroff.
 *
 */
// note that we should not compute def1 and def2 in this function as they don't change
// these should be computed in grid.c and just looked up for efficiency
static REAL HFaceFlux(int j, int k, REAL *phi, REAL **u, gridT *grid, REAL dt, int method) {
  int nc1, nc2;
  REAL def1, def2, Dj, C=0;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  if(nc1==-1) nc1=nc2;
  if(nc2==-1) nc2=nc1;
  Dj = grid->dg[j];
  def1 = sqrt(pow(grid->xv[nc1]-grid->xe[j],2)+
      pow(grid->yv[nc1]-grid->ye[j],2));
  def2 = Dj-def1;

  if(method==4) C = u[j][k]*dt/Dj;

  if(def1==0 || def2==0) {
    // this happens on a boundary cell
    return UpWind(u[j][k],phi[nc1],phi[nc2]);
  }
  else {
    // on an interior cell (orthogonal cell) we can just distance interpolation the 
    // value to the face)
    // basically Eqn 51
    return (phi[nc1]*def2+phi[nc2]*def1)/(def1+def2)
      -C/2*(phi[nc1]-phi[nc2]);
  }
}

/*
 *
 */
static void GetMomentumFaceValues(REAL **uface, REAL **ui, REAL **boundary_ui, REAL **U, gridT *grid, physT *phys, propT *prop,
				  MPI_Comm comm, int myproc, int nonlinear) {
  int i, iptr, j, jptr, nc, nf, ne, nc1, nc2, k, kmin;
  REAL tempu;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];
    i = grid->grad[2*j];

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      if(U[j][k]>0)
	uface[j][k]=boundary_ui[jptr-grid->edgedist[2]][k];
      else
	uface[j][k]=ui[i][k];
    }
  }

  // type 4 boundary conditions used for no slip 
  // but typically we think of no-slip boundary conditions as 
  // not having flux across them, so this is more general
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      uface[j][k]=boundary_ui[jptr-grid->edgedist[2]][k];
  }

  if(prop->nonlinear==5) //use tvd for advection of momemtum
    HorizontalFaceScalars(grid,phys,prop,ui,boundary_ui,prop->TVDmomentum,comm,myproc);

  // over each of the "computational" cells
  // Compute the u-component fluxes at the faces
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // figure out which column adjacent to edge is limiter for the
    // top of the cell (cells much share face to have flux)
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];
    
    for(k=0;k<kmin;k++)
      uface[j][k]=0;

    // compute mass flow rate by interpolating the velocity to the face flux
    // note that we can use UFaceFlux to just compute the value directly if 
    // we have quadratic upwinding available to us
    // basically Eqn 47
    for(k=kmin;k<grid->Nke[j];k++) {
      // compute upwinding data
      if(U[j][k]>0)
	nc=nc2;
      else
	nc=nc1;

      switch(prop->nonlinear) {
      case 1:
	uface[j][k]=ui[nc][k];
            break;
          case 2:
            uface[j][k]=UFaceFlux(j,k,ui,U,grid,prop->dt,nonlinear);
            break;
          case 4:
            uface[j][k]=UFaceFlux(j,k,ui,U,grid,prop->dt,nonlinear);
            break;
          case 5:
            if(U[j][k]>0)
              tempu=phys->SfHp[j][k];
            else
              tempu=phys->SfHm[j][k];
            uface[j][k]=tempu;
            break;
          default:
            uface[j][k]=ui[nc][k];
            break;
        }
      }
    }

    // Faces on type 3 cells are always updated with first-order upwind

    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      for(nf=0;nf<NFACES;nf++) {
        if((ne=grid->neigh[i*NFACES+nf])!=-1) {
          j=grid->face[i*NFACES+nf];

          nc1 = grid->grad[2*j];
          nc2 = grid->grad[2*j+1];

          if(grid->ctop[nc1]>grid->ctop[nc2])
            kmin = grid->ctop[nc1];
          else
            kmin = grid->ctop[nc2];

          for(k=0;k<kmin;k++)
            uface[j][k]=0;

          for(k=kmin;k<grid->Nke[j];k++) {
            if(U[j][k]>0)
              nc=nc2;
            else
              nc=nc1;
            uface[j][k]=ui[nc][k];
          }
        }
      }
    }
}












