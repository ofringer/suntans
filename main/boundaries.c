/*
 * File: boundaries.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 03/11/04
 * --------------------------------
 * This file contains functions to impose the boundary conditions on u.
 *
 * $Id: boundaries.c,v 1.8 2004-09-27 01:28:15 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.7  2004/07/27 20:31:13  fringer
 * Added SetBoundaryScalars function, which allows the specification of
 * the salinity or temperature on the boundaries.
 *
 * Revision 1.6  2004/06/23 06:20:36  fringer
 * This is the form of boundaries.c used to force the Monterey Bay run
 * over the Spring-Neap cycle.
 *
 * Revision 1.5  2004/06/17 02:53:58  fringer
 * Added river plume forcing code.
 *
 * Revision 1.4  2004/06/15 18:24:23  fringer
 * Cleaned up and removed extraneous code used for testing.
 *
 * Revision 1.3  2004/06/13 07:08:57  fringer
 * Changes after testing of open boundaries.  No physical breakthroughs...
 *
 * Revision 1.2  2004/05/29 20:25:02  fringer
 * Revision before converting to CVS.
 *
 * Revision 1.1  2004/03/12 06:21:56  fringer
 * Initial revision
 *
 */
#include "boundaries.h"

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

    // For Monterey
    if(grid->xv[ib]<1000) {
      ub0=0.002445*cos(prop->omega*prop->rtime)+.00182*cos(2*PI/(24*3600)*prop->rtime+.656);
      //ub0=-2.9377*sqrt(GRAV/grid->dv[ib])*(cos(prop->omega*prop->rtime)+0.00182/.002445*cos(2*PI/(24*3600)*prop->rtime+.656));
      forced=1;
    } else 
      forced=0;
    
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
    C0=0*2.0*c0*dt/grid->dg[j];
    C1=2.0*c1*dt/grid->dg[j];

    // Now update u0 with c0 to get u0new
    u0new=ub0;((1-C0/2-forced*dt/2/prop->timescale)*u0+0.5*C0*(3.0*(grid->n1[j]*uc0+grid->n2[j]*vc0)-
						   (grid->n1[j]*uc0old+grid->n2[j]*vc0old))
	   -forced*ub0/prop->timescale)/(1+C0/2+dt/2/prop->timescale);

    // Now update u1 with c1 but internal waves are only radiated and not forced
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      ub[j][k]=u0new+((1-C1/2)*(ub[j][k]-u0)+0.5*C1*(3.0*(grid->n1[j]*(u[ib][k]-uc0)+grid->n2[j]*(v[ib][k]-vc0))-
						     (grid->n1[j]*(uold[ib][k]-uc0old)+grid->n2[j]*(vold[ib][k]-vc0old))))
	/(1+C1/2);
  }
}

/*
 * Function: SetBoundaryScalars
 * Usage: SetBoundaryScalars(boundary_s,grid,phys,prop,"salt");
 * ------------------------------------------------------------
 * This will set the values of the scalars at the open boundaries.
 * 
 */
void SetBoundaryScalars(REAL **boundary_scal, gridT *grid, physT *phys, propT *prop, char *type) {
  int jptr, j, ib, k;

  if(!strcmp(type,"salt")) {
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      for(k=grid->ctop[ib];k<grid->Nk[ib];k++)
	boundary_scal[jptr-grid->edgedist[2]][k]=-0.4;
    }
  } else if(!strcmp(type,"temperature")) {
    // Set the temperature here!
  } else if(WARNING)
    printf("Warning! Invalid boundary scalar specification in SetBoundaryScalars.\n");  
}
	
      
