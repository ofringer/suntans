/*
 * File: turbulence.c
 * Description:  Contains the Mellor-Yamad level 2.5 turbulence model.
 *
 * $Id: turbulence.c,v 1.4 2004-09-27 01:04:31 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2004/09/22 06:29:33  fringer
 * Latest version of my25 function with several changes.  The previous
 * version was checked in too soon and was not yet working properly.
 *
 * Revision 1.2  2004/09/15 01:10:47  fringer
 * Added MY25 but it seems to have problems predicting the log layer
 * correctly for a channel flow.  For some reason the velocity is
 * underpredicted by the code because the eddy viscosity is overpredicted
 * when compared to its theoretical value.
 *
 * Revision 1.1  2004/09/13 04:14:36  fringer
 * Contains functions that compute the eddy viscosity and scalar
 * diffusivity.
 *
 *
 */
#include "math.h"
#include "phys.h"
#include "grid.h"
#include "phys.h"
#include "util.h"
#include "turbulence.h"

// Local function
static void StabilityFunctions(REAL *Sm, REAL *Sh, REAL Gh, REAL A1, REAL A2, REAL B1, REAL B2, REAL C1);

/*
 * Function: my25
 * Usage: my25(grid,phys,prop,phys->qT,phys->lT,phys->Cn_q,phys->Cn_l,phys->nu_tv,phys->kappa_tv);
 * -----------------------------------------------------------------------------------------------
 * Computes the eddy viscosity and scalar diffusivity via the MY25 closure.  Advection of the turbulent
 * quantities q^2 and q^2l is included with the use of UpdateScalars
 *
 */
void my25(gridT *grid, physT *phys, propT *prop, REAL **q, REAL **l, REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc) {
  int i, j, iptr, k, nf, nc1, nc2;
  REAL thetaQ=0.75, CdAvgT, CdAvgB, *dudz, *dvdz, *drdz, z, *N, *Gh;
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
      drdz[k]=2.0*prop->beta*(phys->s[i][k-1]-phys->s[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
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
      				GRAV*(prop->kappa_s+kappaT[i][k])*0.5*(drdz[k]+drdz[k+1]));
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
    for(nf=0;nf<NFACES;nf++) {
      CdAvgT+=phys->CdT[grid->face[i*NFACES+nf]]/3;
      CdAvgB+=phys->CdB[grid->face[i*NFACES+nf]]/3;
    }
    phys->htmp[i]=pow(B1,2/3)*CdAvgT*pow(phys->uc[i][grid->ctop[i]],2)+pow(phys->vc[i][grid->ctop[i]],2);
    phys->hold[i]=pow(B1,2/3)*CdAvgB*pow(phys->uc[i][grid->Nk[i]-1],2)+pow(phys->vc[i][grid->Nk[i]-1],2);
  }
  UpdateScalars(grid,phys,prop,q,NULL,phys->Cn_q,0,0,kappaT,thetaQ,phys->uold,phys->wtmp,phys->htmp,phys->hold,1,1);

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
  if(prop->n<3)
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];
      
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	l[i][k]+=prop->dt*pow(q[i][k],3/2)/B1;
      }
    }
  else
    UpdateScalars(grid,phys,prop,l,NULL,phys->Cn_l,0,0,kappaT,thetaQ,phys->uold,phys->wtmp,phys->htmp,phys->hold,1,1);

  // Send/Recv q and l data to neighboring processors
  ISendRecvCellData3D(q,grid,myproc,comm);
  ISendRecvCellData3D(l,grid,myproc,comm);

  // l stores q^2 l
  // q stores q^2
  // Extract q and l from their stored quantities
  // and then set the values of nuT and kappaT
  for(i=0;i<grid->Nc;i++) {

    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      drdz[k]=-2.0*GRAV*prop->beta*(phys->s[i][k-1]-phys->s[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);
      if(drdz[k]<0) drdz[k]=0;
      N[k]=sqrt(drdz[k]);
    }
    N[grid->Nk[i]]=N[grid->Nk[i]-1];

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
  *Sm = (pow(B1,-1/3)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1)))/
    ((1-3*A2*Gh*(6*A1+B2))*(1-9*A1*A2*Gh));
  *Sh = A2*(1-6*A1/B1)/(1-3*A2*Gh*(6*A1+B2));
}

