/*
 * File: turbulence.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains the Mellor-Yamada level 2.5 turbulence model.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "math.h"
#include "phys.h"
#include "grid.h"
#include "phys.h"
#include "sendrecv.h"
#include "util.h"
#include "turbulence.h"
#include "boundaries.h"
#include "scalars.h"

// Local function
static void StabilityFunctions(REAL *Sm, REAL *Sh, REAL Gh, REAL A1, REAL A2, REAL B1, REAL B2, REAL C1);

/*
 * Function: my25
 * Usage: my25(grid,phys,prop,wnew,phys->qT,phys->lT,phys->Cn_q,phys->Cn_l,phys->nu_tv,phys->kappa_tv);
 * -----------------------------------------------------------------------------------------------
 * Computes the eddy viscosity and scalar diffusivity via the MY25 closure.  Advection of the turbulent
 * quantities q^2 and q^2l is included with the use of UpdateScalars
 *
 */
void my25(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **q, REAL **l, REAL **Cn_q, REAL **Cn_l, 
	  REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc) {
  int i, ib, j, iptr, jptr, k, nf, nc1, nc2, ne;
  REAL thetaQ=1, CdAvgT, CdAvgB, *dudz, *dvdz, *drdz, z, *N, *Gh, tauAvgT;
  REAL A1, A2, B1, B2, C1, E1, E2, E3, Sq, Sm, Sh;

  N = dudz = phys->a;
  dvdz = phys->b;
  drdz = phys->c;
  Gh = phys->d;

  // Specification of constants
  A1 = 0.92;
  A2 = 0.74;
  B1 = 16.6;
  B2 = 10.1;
  C1 = 0.08;
  E1 = 1.8;
  E2 = 1.33;
  E3 = 0.25;
  Sq = 0.2;
  
  // First solve for q^2 and store its old value in stmp3
  for(i=0;i<grid->Nc;i++) {

    // dudz, dvdz, and drdz store gradients at k-1/2
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      dudz[k]=2.0*(phys->uc[i][k-1]-phys->uc[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
      dvdz[k]=2.0*(phys->vc[i][k-1]-phys->vc[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
      drdz[k]=2.0*(phys->rho[i][k-1]-phys->rho[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
    }
    dudz[grid->ctop[i]]=dudz[grid->ctop[i]+1];
    dvdz[grid->ctop[i]]=dvdz[grid->ctop[i]+1];
    drdz[grid->ctop[i]]=drdz[grid->ctop[i]+1];
    dudz[grid->Nk[i]]=dudz[grid->Nk[i]-1];
    dvdz[grid->Nk[i]]=dvdz[grid->Nk[i]-1];
    drdz[grid->Nk[i]]=drdz[grid->Nk[i]-1];
    
    // uold will store src1 for q^2, which is the 2q/B1 l term
    // wtmp will store src2 for q^2, which is the 2 (Ps+Pb) term
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->uold[i][k]=2.0*q[i][k]/B1/(l[i][k]+SMALL);
      phys->wtmp[i][k]=2.0*fabs((prop->nu+nuT[i][k])*(pow(0.5*(dudz[k]+dudz[k+1]),2)+pow(0.5*(dvdz[k]+dvdz[k+1]),2))+
      				prop->grav*(prop->kappa_s+kappaT[i][k])*0.5*(drdz[k]+drdz[k+1]));
    }

    // kappaT will store the diffusion coefficient for q^2
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      kappaT[i][k]=q[i][k]*l[i][k]*Sq;

    // q will store q^2 and stmp3 will store the old value of q
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->stmp3[i][k]=q[i][k];
      q[i][k]*=q[i][k];
    }

    // htmp will store the value at the top boundary for q^2
    // hold will store it at the bottom boundary
    // The drag coefficient is the average of the coefficients on the cell faces
    CdAvgT=0;
    CdAvgB=0;
    tauAvgT=0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      CdAvgT+=phys->CdT[ne]/3;
      CdAvgB+=phys->CdB[ne]/3;
      tauAvgT+=fabs(phys->tau_T[ne])/3;
    }
    phys->htmp[i]=pow(B1,2.0/3.0)*(CdAvgT*(pow(phys->uc[i][grid->ctop[i]],2)+pow(phys->vc[i][grid->ctop[i]],2))+
				   tauAvgT);
    phys->hold[i]=pow(B1,2.0/3.0)*CdAvgB*(pow(phys->uc[i][grid->Nk[i]-1],2)+pow(phys->vc[i][grid->Nk[i]-1],2));
  }
  // Specify turbulence at boundaries for use in updatescalars.  Assume that all incoming turbulence is zero and let outgoing
  // turbulence flow outward.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    ib = grid->grad[2*j];
    
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
      phys->boundary_tmp[jptr-grid->edgedist[2]][k]=q[ib][k];
  }    
  UpdateScalars(grid,phys,prop,wnew,q,phys->boundary_tmp,phys->Cn_q,0,0,kappaT,thetaQ,phys->uold,phys->wtmp,
		phys->htmp,phys->hold,1,1,comm,myproc,0,prop->TVDturb);

  // q now contains q^2
  for(i=0;i<grid->Nc;i++) {

    // uold will store src1 for q^2 l, which is the q/B1 l*(1+E2(l/kz)^2+E3(l/k(H-z))^2) term
    // wtmp will store src2 for q^2 l, which is the l E1 (Ps+Pb) term
    z = phys->h[i];
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      z-=grid->dzz[i][k]/2;
      phys->uold[i][k]*=0.5*(1+E2*pow(l[i][k]/KAPPA_VK/(z-phys->h[i]),2)+E3*pow(l[i][k]/KAPPA_VK/(grid->dv[i]+z),2));
      phys->wtmp[i][k]*=0.5*l[i][k]*E1;
      z-=grid->dzz[i][k]/2;
    }

    // kappaT already stores q l Sq from before
    // l will store q^2 l 
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      l[i][k]*=pow(phys->stmp3[i][k],2);
    }

    // htmp will store the value at the top boundary for q^2 l
    // hold will store it at the bottom boundary (both are 0)
    phys->htmp[i]=0;
    phys->hold[i]=0;
  }

  // Specify turbulence at boundaries for use in updatescalars.  Assume that all incoming turbulence is zero and let outgoing
  // turbulence flow outward.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    ib = grid->grad[2*j];
    
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
      phys->boundary_tmp[jptr-grid->edgedist[2]][k]=l[ib][k];
  }
  UpdateScalars(grid,phys,prop,wnew,l,phys->boundary_tmp,phys->Cn_l,0,0,kappaT,thetaQ,phys->uold,phys->wtmp,
		phys->htmp,phys->hold,1,1,comm,myproc,0,prop->TVDturb);

  // Set l to a background value if it gets too small.
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i=grid->cellp[iptr];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(l[i][k]<LBACKGROUND) l[i][k]=LBACKGROUND;
    }
  }

  // Send/Recv q and l data to neighboring processors
  ISendRecvCellData3D(q,grid,myproc,comm);
  ISendRecvCellData3D(l,grid,myproc,comm);

  // l stores q^2 l
  // q stores q^2
  // Extract q and l from their stored quantities
  // and then set the values of nuT and kappaT
  for(i=0;i<grid->Nc;i++) {

    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      drdz[k]=-2.0*prop->grav*(phys->rho[i][k-1]-phys->rho[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
      if(drdz[k]<0) drdz[k]=0;
      N[k]=sqrt(drdz[k]);
    }
    if(grid->ctop[i]<grid->Nk[i]-1)
      N[grid->ctop[i]]=N[grid->ctop[i]+1];
    else
      N[grid->ctop[i]]=0;

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(l[i][k]<0) l[i][k]=0;
      if(q[i][k]<0) q[i][k]=0;
      l[i][k]=l[i][k]/(q[i][k]+SMALL);
      q[i][k]=sqrt(q[i][k]);
      l[i][k]=Min(0.53*q[i][k]/(N[k]+SMALL),l[i][k]);
      Gh[k]=-pow(N[k]*l[i][k]/(q[i][k]+SMALL),2);
      StabilityFunctions(&Sm,&Sh,Gh[k],A1,A2,B1,B2,C1);

      nuT[i][k]=Sm*q[i][k]*l[i][k];
      kappaT[i][k]=Sh*q[i][k]*l[i][k];
    }
    for(k=0;k<grid->ctop[i];k++)
      nuT[i][k]=kappaT[i][k]=l[i][k]=q[i][k]=0;
  }
}

/*
 * Function: Stability Functions
 * Usage:  StabilityFunctions(&Sm,&Sh,b[k],A1,A2,B1,B2,C1);
 * --------------------------------------------------------
 * Computes the Stability functions of Blumberg et al. (1992) and 
 * places them into Sm and Sh.
 *
 */
static void StabilityFunctions(REAL *Sm, REAL *Sh, REAL Gh, REAL A1, REAL A2, REAL B1, REAL B2, REAL C1) {
  *Sm = (pow(B1,-1.0/3.0)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1)))/
    ((1-3*A2*Gh*(6*A1+B2))*(1-9*A1*A2*Gh));
  *Sh = A2*(1-6*A1/B1)/(1-3*A2*Gh*(6*A1+B2));
}

