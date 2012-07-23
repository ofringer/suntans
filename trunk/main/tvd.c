/*
 * File: tvd.c
 * Author: Zhonghua Zhang and Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains the TVD scalar-computing functions.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "tvd.h"
#include "util.h"
#include <math.h>

// To prevent denominators from going to zero.
#define EPS 1e-12

// Local function
static REAL Psi(REAL r, int TVD);

/*
 * Function: HorizontalFaceScalars
 * Usage: HorizontalFaceScalars(grid, phys, boundary_scal);
 * ---------------------------------------------------------------------------
 * Calculate the horizontal face scalars with upwind/TVD schemes.  
 * SfHp[Ne][Nk] & SfHm[Ne][Nk] are used to store the scalar facial values.
 * where S--scalar, f--face, H--horizontal, p--plus, m--minus;
 *       Ne--the number of horizontal edges, Nk--the number of vertical layers. 
 *
 */
void HorizontalFaceScalars(gridT *grid, physT *phys, propT *prop, REAL **boundary_scal, int TVD, 
			   MPI_Comm comm, int myproc) 
{
  int i, k, nf, iptr;
  int normal, nc1, nc2, ne;

  REAL df, dg, *sp, Ac, dt=prop->dt;
  REAL *Cp, *Cm, *rp, *rm, **gradSx, **gradSy;   

  // For check!
  int iu, ku, nfu, gradflag, faceflag, Cflag, Rflag;

  // Pointer Variables for TVD scheme
  Cp = phys->Cp;
  Cm = phys->Cm;
  rp = phys->rp;
  rm = phys->rm;
  gradSx = phys->gradSx;
  gradSy = phys->gradSy;


  for(i=0; i<grid->Nc; i++) {
    Ac = grid->Ac[i];
   
    // Initialize the gradSx and gradSy
    for(k=0;k<grid->Nk[i];k++){
      gradSx[i][k] = 0;
      gradSy[i][k] = 0;
    }

    // Loop through all faces of the current cell
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) {
	nc2=nc1;
	if(boundary_scal && grid->mark[ne]>1)
	  sp=phys->stmp2[nc1];
	else
	  sp=phys->stmp[nc1];
      }else 
        sp=phys->stmp[nc2];

      for(k=0;k<grid->Nke[ne];k++) {
	gradSx[i][k]+=1/Ac*0.5*(phys->stmp[nc1][k]+sp[k])*grid->n1[ne]*normal*df; 
	gradSy[i][k]+=1/Ac*0.5*(phys->stmp[nc1][k]+sp[k])*grid->n2[ne]*normal*df; 
      }
      for(k=grid->Nke[ne];k<grid->Nk[i];k++) {
	gradSx[i][k]+=1/Ac*phys->stmp[i][k]*grid->n1[ne]*normal*df; 
	gradSy[i][k]+=1/Ac*phys->stmp[i][k]*grid->n2[ne]*normal*df; 	
      }
    } 
  } 
  ISendRecvCellData3D(gradSx,grid,myproc,comm);  
  ISendRecvCellData3D(gradSy,grid,myproc,comm);  
 
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) {
	nc2=nc1;
	if(boundary_scal && grid->mark[ne]>1)
	  sp=phys->stmp2[nc1];
	else
	  sp=phys->stmp[nc1];
      } else 
	sp=phys->stmp[nc2];


      for(k=0;k<grid->Nke[ne];k++) {
	Cp[k]= 0.5*(phys->utmp2[ne][k]+fabs(phys->utmp2[ne][k]))*dt/dg;
	Cm[k]= 0.5*(phys->utmp2[ne][k]-fabs(phys->utmp2[ne][k]))*dt/dg;
      }

      for(k=0;k<grid->Nke[ne];k++) {	
	rp[k]= 2*(gradSx[nc2][k]*grid->n1[ne]*dg + gradSy[nc2][k]*grid->n2[ne]*dg + EPS )/
	  (phys->stmp[nc1][k]-sp[k]+EPS)-1;
	rm[k]= 2*(gradSx[nc1][k]*grid->n1[ne]*dg + gradSy[nc1][k]*grid->n2[ne]*dg + EPS )/
	  (phys->stmp[nc1][k]-sp[k]+EPS)-1;
      }

      for(k=0;k<grid->Nke[ne];k++) {
	phys->SfHp[ne][k] = sp[k]+0.5*Psi(rp[k],TVD)*(1-Cp[k])*(phys->stmp[nc1][k]-sp[k]);
	phys->SfHm[ne][k] = phys->stmp[nc1][k]-0.5*Psi(rm[k],TVD)*(1+Cm[k])*(phys->stmp[nc1][k]-sp[k]);
      }

      for(k=grid->Nke[ne];k<grid->Nk[nc1];k++) 
	phys->SfHm[ne][k] = phys->stmp[nc1][k];

      for(k=grid->Nke[ne];k<grid->Nk[nc2];k++) 
	phys->SfHp[ne][k] = sp[k];
    } 
  }
}

/*
 * Function: Psi
 * Usage: Psi(r, TVD);
 * ---------------------------------------------------------------------------
 * Calculate the filter function Psi for TVD schemes 
 * TVD=1  first-order upwind,    TVD=2  Lax-Wendroff
 * TVD=3  Superbee,              TVD=4  Van Leer  
 */
static REAL Psi(REAL r, int TVD){
  switch(TVD) {
  case 1:
    return 0;
    break;
  case 2:
    return 1;
    break;
  case 3:
    return ( Max(0, Max(Min(2*(r),1), Min(r,2))));
    break;
  case 4:
    if(r>0)
      return Min(r,1);
    else
      return 0;
    break;
  default:
    return 0;
    break;
  }
}

/*
 * Function: GetApAm
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *	           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the fluxes for vertical advection using the TVD schemes.
 *
 */
void GetApAm(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	     REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD) {
  int k;

  // Implicit vertical advection terms
  for(k=0;k<Nk+1;k++) {
    wp[k] = 0.5*(w[i][k]+fabs(w[i][k]));
    wm[k] = 0.5*(w[i][k]-fabs(w[i][k]));
  }

  // Courant number: C=w*dt/dz
  for(k=1;k<Nk;k++) {
    Cp[k] = 2 * wp[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
    Cm[k] = 2 * wm[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
  }
  k=Nk;
  Cp[k] = wp[k]*dt/dzz[i][k-1];
  Cm[k] = wm[k]*dt/dzz[i][k-1];


  // Compute the upwind gradient ratios r
  for(k=ktop+2;k<Nk-1;k++) {
    rp[k-ktop] = (scal[i][k]-scal[i][k+1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);
    rm[k-ktop] = (scal[i][k-2]-scal[i][k-1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);
  }

  rp[1]= (scal[i][ktop+1]-scal[i][ktop+2]+EPS) / (scal[i][ktop]-scal[i][ktop+1]+EPS);
  rm[1]= EPS / (scal[i][ktop]-scal[i][ktop+1]+EPS);

  k=Nk-1;
  rp[k-ktop]=EPS / (scal[i][k-1]-scal[i][k]+EPS);
  rm[k-ktop]=(scal[i][k-2]-scal[i][k-1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);

  k=Nk;
  rp[k-ktop]=1;
  rm[k-ktop]=(scal[i][k-2]-scal[i][k-1]+EPS) / EPS;

  for(k=ktop+1;k<Nk+1;k++) {
    am[k]= 0.5*wp[k]*Psi(rp[k-ktop], TVD)*(1-Cp[k]) 
             + wm[k]*(1-0.5*Psi(rm[k-ktop], TVD)*(1+Cm[k]));
    ap[k]= wp[k]*(1-0.5*Psi(rp[k-ktop],TVD)*(1-Cp[k])) 
             + 0.5*wm[k]*Psi(rm[k-ktop],TVD)*(1+Cm[k]);
  }
}
/*
 * Function: HorizontalFaceU
 * Usage: HorizontalFaceScalars(uc, grid, phys, boundary_scal);
 * ---------------------------------------------------------------------------
 * Calculate the horizontal face values for uc and vc using TVD schemes.  
 * SfHp[Ne][Nk] & SfHm[Ne][Nk] are used to store the facial values.
 */
void HorizontalFaceU(REAL **uc, gridT *grid, physT *phys, propT *prop, int TVD, 
			   MPI_Comm comm, int myproc) 
{
  int i, k, nf, iptr;
  int normal, nc1, nc2, ne;

  REAL df, dg, *sp, Ac, dt=prop->dt;
  REAL *Cp, *Cm, *rp, *rm, **gradSx, **gradSy, **stmp;   

  // For check!
  int iu, ku, nfu, gradflag, faceflag, Cflag, Rflag;

  // Pointer Variables for TVD scheme
  Cp = phys->Cp;
  Cm = phys->Cm;
  rp = phys->rp;
  rm = phys->rm;
  gradSx = phys->gradSx;
  gradSy = phys->gradSy;
  stmp = uc;


  for(i=0; i<grid->Nc; i++) {
    Ac = grid->Ac[i];
   
    // Initialize the gradSx and gradSy
    for(k=0;k<grid->Nk[i];k++){
      gradSx[i][k] = 0;
      gradSy[i][k] = 0;
    }

    // Loop through all faces of the current cell
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) {
	nc2=nc1;
	//if(grid->mark[ne]>1) //(boundary_scal)
	//  sp=phys->stmp2[nc1];
	//else
	sp=stmp[nc1];
      }else 
        sp=stmp[nc2];


      for(k=0;k<grid->Nke[ne];k++) {
	gradSx[i][k]+=1/Ac*0.5*(stmp[nc1][k]+sp[k])*grid->n1[ne]*normal*df; 
	gradSy[i][k]+=1/Ac*0.5*(stmp[nc1][k]+sp[k])*grid->n2[ne]*normal*df;

      }

      for(k=grid->Nke[ne];k<grid->Nk[i];k++) {
	gradSx[i][k]+=1/Ac*stmp[i][k]*grid->n1[ne]*normal*df; 
	gradSy[i][k]+=1/Ac*stmp[i][k]*grid->n2[ne]*normal*df; 	
      }
    } 
  } 
  ISendRecvCellData3D(gradSx,grid,myproc,comm);  
  ISendRecvCellData3D(gradSy,grid,myproc,comm);  
 
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) {
	nc2=nc1;
	//if(boundary_scal)
	//  sp=phys->stmp2[nc1];
	//else
	  sp=stmp[nc1];
      } else 
	sp=stmp[nc2];


      for(k=0;k<grid->Nke[ne];k++) {
	Cp[k]= 0.5*(phys->utmp2[ne][k]+fabs(phys->utmp2[ne][k]))*dt/dg;
	Cm[k]= 0.5*(phys->utmp2[ne][k]-fabs(phys->utmp2[ne][k]))*dt/dg;
      }

      for(k=0;k<grid->Nke[ne];k++) {	
	rp[k]= 2*(gradSx[nc2][k]*grid->n1[ne]*dg + gradSy[nc2][k]*grid->n2[ne]*dg + EPS )/
	  (stmp[nc1][k]-sp[k]+EPS)-1;
	rm[k]= 2*(gradSx[nc1][k]*grid->n1[ne]*dg + gradSy[nc1][k]*grid->n2[ne]*dg + EPS )/
	  (stmp[nc1][k]-sp[k]+EPS)-1;
      }

      for(k=0;k<grid->Nke[ne];k++) {
	phys->SfHp[ne][k] = sp[k]+0.5*Psi(rp[k],TVD)*(1-Cp[k])*(stmp[nc1][k]-sp[k]);
	phys->SfHm[ne][k] = stmp[nc1][k]-0.5*Psi(rm[k],TVD)*(1+Cm[k])*(stmp[nc1][k]-sp[k]);
      }

      for(k=grid->Nke[ne];k<grid->Nk[nc1];k++) 
	phys->SfHm[ne][k] = stmp[nc1][k];

      for(k=grid->Nke[ne];k<grid->Nk[nc2];k++) 
	phys->SfHp[ne][k] = sp[k];
    } 
  }
}
