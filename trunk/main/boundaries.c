/*
 * File: boundaries.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 03/11/04
 * --------------------------------
 * This file contains functions to impose the boundary conditions on u.
 *
 * $Id: boundaries.c,v 1.1 2004-03-12 06:21:56 fringer Exp $
 * $Log: not supported by cvs2svn $
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
 */
void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop) {

  int forced, ib, j, jptr, k, nf, nf0, ne, nc1, nc2, numneighs, neigh;
  REAL z, uf, vf, cmax, *uboundary = phys->c, *c = phys->b, *phin = phys->a, xc, yc, width,
    **u = phys->uc, **v = phys->vc, **uold = phys->uold, **vold = phys->vold;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];
    
    // For Three-mile slough
    if(grid->yv[ib]>2500) {
      if(cos(prop->omega*prop->rtime)>0) {
	forced=0;
	for(k=grid->etop[j];k<grid->Nke[j];k++) 
	  uboundary[k] = -phys->h[ib]*sqrt(GRAV/(grid->dv[ib]));
      } else {
	forced=1;
	for(k=grid->etop[j];k<grid->Nke[j];k++) 
	  uboundary[k] = -prop->amp*cos(prop->omega*prop->rtime);
      }
    } else {
      if(cos(prop->omega*prop->rtime)<0) {
	forced=0;
	for(k=grid->etop[j];k<grid->Nke[j];k++) 
	  uboundary[k] = -phys->h[ib]*sqrt(GRAV/(grid->dv[ib]));
      } else {
	forced=1;
	for(k=grid->etop[j];k<grid->Nke[j];k++) 
	  uboundary[k] = prop->amp*cos(prop->omega*prop->rtime);
      }
    }
    forced=1;

    // For Huntington Beach
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      uboundary[k] = prop->amp*cos(prop->omega*prop->rtime);

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      c[k] = 0;
      phin[k] = 0;
    }
    numneighs=0;

    for(nf0=0;nf0<NFACES;nf0++) {
      neigh = grid->neigh[ib*NFACES+nf0];

      if(neigh!=-1) {

	for(k=grid->etop[j];k<grid->Nke[j];k++) 
	  c[k]+=(grid->n1[j]*(u[neigh][k]-uold[neigh][k])+
		 grid->n2[j]*(v[neigh][k]-vold[neigh][k]))/prop->dt;

	for(nf=0;nf<NFACES;nf++) {
	  ne=grid->face[neigh*NFACES+nf];
	  nc1 = grid->grad[2*ne];
	  nc2 = grid->grad[2*ne+1];
	  if(nc1==-1) nc1=nc2;
	  if(nc2==-1) nc2=nc1;
	  
	  for(k=grid->etop[j];k<grid->Nke[j];k++) {
	    //	    uf = u[nc1][k]*grid->def[neigh*NFACES+grid->gradf[2*ne]]/grid->dg[ne]+
	    //	      u[nc2][k]*(1-grid->def[neigh*NFACES+grid->gradf[2*ne]]/grid->dg[ne]);
	    //	    vf = v[nc1][k]*grid->def[neigh*NFACES+grid->gradf[2*ne]]/grid->dg[ne]+
	    //	      v[nc2][k]*(1-grid->def[neigh*NFACES+grid->gradf[2*ne]]/grid->dg
	    uf = 0.5*(u[nc1][k]+u[nc2][k]);
	    vf = 0.5*(v[nc1][k]+v[nc2][k]);
	    phin[k]+=(grid->n1[j]*uf+grid->n2[j]*vf)*
	      (grid->n1[ne]*grid->n1[j]+grid->n2[ne]*grid->n2[j])*
	      grid->normal[neigh*NFACES+nf]*grid->df[ne]/grid->Ac[neigh];
	    //	    phin[k] += (vf+0*grid->n1[j]*uf+0*grid->n2[j]*vf)*
	    //	      (grid->n1[ne]*grid->n1[j]+grid->n2[ne]*grid->n2[j])*
	    //	      grid->normal[neigh*NFACES+nf]*grid->df[ne]/grid->Ac[neigh];
	    //	    phintemp=(grid->n1[j]*uf+grid->n2[j]*vf);
	    //	    if(fabs(phintemp)>fabs(phin[k]))
	    //	      phin[k]=phintemp;
	  }
	}
	numneighs++;
      }
    }

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phin[k]/=numneighs;
      c[k]/=numneighs;
    }

    cmax=sqrt(GRAV*(grid->dv[ib]+phys->h[ib]));
    for(k=grid->etop[j];k<grid->Nke[j];k++) {

      if(phin[k]==0)
	c[k]=0;
      else 
	c[k]/=phin[k];

      if(c[k]<0)
      	c[k]=0;
      if(c[k]>cmax)
      	c[k]=cmax;

      c[k]=cmax;
    }
    //    printf("c[%d]=%f,%f,phin[%d]=%f,(%f,%f)\n",0,c[0],cmax,0,phin[0],grid->xv[ib],grid->yv[ib]);
    //    printf("numneighs = %d\n",numneighs);

    if(forced) {
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	ub[j][k]=(1-prop->dt*(1-prop->theta)/prop->timescale)*ub[j][k]+prop->dt*uboundary[k]/prop->timescale
      	  +0*2.0*prop->dt*c[k]*(phys->uc[ib][k]*grid->n1[j]+phys->vc[ib][k]*grid->n2[j]-ub[j][k])/
      	  grid->dg[j];
      for(k=grid->etop[j];k<grid->Nke[j];k++)
      	ub[j][k]/=(1+prop->theta*prop->dt/prop->timescale);
    } else {
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	ub[j][k]=(ub[j][k]+2.0*prop->dt*c[k]*(phys->uc[ib][k]*grid->n1[j]+phys->vc[ib][k]*grid->n2[j])/grid->dg[j])
	  /(1+2*prop->dt*c[k]/grid->dg[j]);
      //	ub[j][k]+=2.0*prop->dt*c[k]*(phys->uc[ib][k]*grid->n1[j]+phys->vc[ib][k]*grid->n2[j]-ub[j][k])/
      //	  grid->dg[j];
     }

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      ub[j][k]=prop->amp*cos(prop->omega*prop->rtime);
  }
}


