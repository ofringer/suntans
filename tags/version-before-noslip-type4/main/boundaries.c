/*
 * File: boundaries.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains functions to impose the boundary conditions on u.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "boundaries.h"

// Local functions
static void GetBoundaryVelocity(REAL *ub, int *forced, REAL x, REAL y, REAL t, REAL h, REAL d, REAL omega, REAL amp);
static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag);

/*
 * Function: GetBoundaryVelocity
 * Usage: GetBoundaryVelocity(&ub0,&forced,grid->xv[ib],grid->yv[ib],
 *                            prop->rtime,phys->h[ib],grid->dv[ib],prop->omega,prop->amp);
 * ---------------------------------------------------------------------------------------
 * Set the boundary velocity based on the location.  If this is a partially-clamped-free bc, then
 * set forced to 1, otherwise, if this is a free open bc, then set forced to 0.
 *
 */
static void GetBoundaryVelocity(REAL *ub, int *forced, REAL x, REAL y, REAL t, REAL h, REAL d, REAL omega, REAL amp) {

  // For Monterey
  if(x<1000) {
    *ub=0.002445*cos(omega*t)+.00182*cos(2*PI/(24*3600)*t+.656);
    //*ub=-2.9377*sqrt(GRAV/d)*(cos(omega*t)+0.00182/.002445*cos(2*PI/(24*3600)*t+.656));
    *forced=1;
  } else 
    *forced=0;

  // New 3MS boundary forcing (cleaner) 9/27/2004 OBF/DAF
  *forced=1; // always
  /*
  if(y>1500) // Northern boundary
    if(cos(omega*t)>0) {
      *ub = -h*sqrt(GRAV/d);
    } else {
      *ub = amp*fabs(cos(omega*t));
    }
  else
    if(cos(omega*t)<0) {
      *ub = -h*sqrt(GRAV/d);
    } else {
      *ub = amp*fabs(cos(omega*t));
    }
  // end 3MS new forcing
  */
  *forced = 1;

}

/*
 * Function: OpenBoundaryFluxes
 * Usage: OpenBoundaryFluxes(q,ubnew,ubn,grid,phys,prop);
 * ----------------------------------------------------
 * This will update the boundary flux at the edgedist[2] to edgedist[3] edges.
 * 
 * Note that phys->uold,vold contain the velocity at time step n-1 and 
 * phys->uc,vc contain it at time step n.
 *
 * The radiative open boundary condition does not work yet!!!  For this reason c[k] is
 * set to 0
 *
 */
void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop) {
  int j, jptr, ib, k, forced;
  REAL *uboundary = phys->a, **u = phys->uc, **v = phys->vc, **uold = phys->uold, **vold = phys->vold;
  REAL z, c0, c1, C0, C1, dt=prop->dt, u0, u0new, uc0, vc0, uc0old, vc0old, ub0;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    GetBoundaryVelocity(&ub0,&forced,grid->xv[ib],grid->yv[ib],prop->rtime,phys->h[ib],grid->dv[ib],prop->omega,prop->amp);

    c0 = sqrt(GRAV*grid->dv[ib]);
    c1 = 1.94;//0.5533; <- for openbc NH/pi
    
    // First compute u0, uc0, vc0, uc0old, vcd0ld;
    u0=uc0=vc0=uc0old=vc0old=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      u0+=ub[j][k]*grid->dzz[ib][k];
      uc0+=u[ib][k]*grid->dzz[ib][k];
      vc0+=v[ib][k]*grid->dzz[ib][k];
      uc0old+=uold[ib][k]*grid->dzz[ib][k];
      vc0old+=vold[ib][k]*grid->dzz[ib][k];
    }
    u0/=(grid->dv[ib]+phys->h[ib]);
    uc0/=(grid->dv[ib]+phys->h[ib]);
    vc0/=(grid->dv[ib]+phys->h[ib]);
    uc0old/=(grid->dv[ib]+phys->h[ib]);
    vc0old/=(grid->dv[ib]+phys->h[ib]);

    // Set the Courant numbers
    C0=2.0*c0*dt/grid->dg[j];
    C1=2.0*c1*dt/grid->dg[j];

    // Now update u0 with c0 to get u0new
    u0new=((1-C0/2-forced*dt/2/prop->timescale)*u0+0.5*C0*(3.0*(grid->n1[j]*uc0+grid->n2[j]*vc0)-
							   (grid->n1[j]*uc0old+grid->n2[j]*vc0old))
	   +forced*ub0*prop->dt/prop->timescale)/(1+C0/2+dt/2/prop->timescale);

    // Now update u1 with c1 but internal waves are only radiated and not forced
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      ub[j][k]=u0new+((1-C1/2)*(ub[j][k]-u0)+0.5*C1*(3.0*(grid->n1[j]*(u[ib][k]-uc0)+grid->n2[j]*(v[ib][k]-vc0))-
						     (grid->n1[j]*(uold[ib][k]-uc0old)+grid->n2[j]*(vold[ib][k]-vc0old))))/(1+C1/2);
  }
}

/*
 * Function: BoundaryScalars
 * Usage: BoundaryScalars(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the scalars at the open boundaries.
 * 
 */
void BoundaryScalars(gridT *grid, physT *phys, propT *prop) {
  int jptr, j, ib, k;
  REAL z;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      z=0;
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
	z-=grid->dzz[ib][k];
	phys->boundary_s[jptr-grid->edgedist[2]][k]=phys->s[ib][k];
	if(z>-10)
	  phys->boundary_T[jptr-grid->edgedist[2]][k]=0;
	else
	  phys->boundary_T[jptr-grid->edgedist[2]][k]=1;
	z-=grid->dzz[ib][k];
      }
  }
}

/*
 * Function: BoundaryVelocities
 * Usage: BoundaryVelocities(grid,phys,prop,myproc);
 * -------------------------------------------------
 * This will set the values of u,v,w, and h at the boundaries.
 * 
 */
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc) {
  int jptr, j, ib, k, boundary_index;
  REAL z, amp=prop->amp, rtime=prop->rtime, omega=prop->omega, boundary_flag;

  /* For three-mile slough */
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    boundary_index = jptr-grid->edgedist[2];
    j=grid->edgep[jptr];
    ib=grid->grad[2*j];

    if(cos(omega*rtime)>=0) { // Ebb (northward)
      if(grid->yv[ib]>1500) 
	boundary_flag=open;
      else
	boundary_flag=specified;
    } else {
      if(grid->yv[ib]>1500) 
	boundary_flag=specified;
      else
	boundary_flag=open;
    }

    phys->boundary_flag[boundary_index]=boundary_flag;

    SetUVWH(grid,phys,prop,ib,j,boundary_index,boundary_flag);
  }
}

static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag) {
  int k;

  if(boundary_flag==open) {
    phys->boundary_h[boundary_index]=phys->h[ib];
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_u[boundary_index][k]=phys->uc[ib][k];
      phys->boundary_v[boundary_index][k]=phys->vc[ib][k];
      phys->boundary_w[boundary_index][k]=0.5*(phys->w[ib][k]+phys->w[ib][k+1]);
    }
  } else {
    phys->boundary_h[boundary_index]=prop->amp*fabs(cos(prop->omega*prop->rtime));
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_u[boundary_index][k]=phys->u[j][k]*grid->n1[j];
      phys->boundary_v[boundary_index][k]=phys->u[j][k]*grid->n2[j];
      phys->boundary_w[boundary_index][k]=0.5*(phys->w[ib][k]+phys->w[ib][k+1]);
    }
  }
}
	
/*
 * Function: WindStress
 * Usage: WindStress(grid,phys,prop,myproc);
 * -----------------------------------------
 * Set the wind stress as well as the bottom stress.
 * tau_B is not currently in use (4/1/05).
 *
 */
void WindStress(gridT *grid, physT *phys, propT *prop, int myproc) {
  int j, jptr;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    phys->tau_T[j]=grid->n2[j]*prop->tau_T;
    phys->tau_B[j]=0;
  }
}
