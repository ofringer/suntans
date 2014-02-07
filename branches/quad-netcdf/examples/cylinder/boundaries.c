/*
 * Boundaries test file.
 *
 */
#include "boundaries.h"
#include "sediments.h"

static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag);

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

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      ub[j][k]=0;

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      ub[j][k]=phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->n1[j]+
	phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->n2[j];
  }
}

/*
 * Function: BoundaryScalars
 * Usage: BoundaryScalars(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the scalars at the open boundaries.
 * 
 */
void BoundaryScalars(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int jptr, j, ib, k;
  REAL z, Tb;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      if(grid->yv[ib]<0.2) 
	Tb = 1.0;
      else
	Tb = 0.0;
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
	phys->boundary_T[jptr-grid->edgedist[2]][k]=Tb;
	phys->boundary_s[jptr-grid->edgedist[2]][k]=phys->s[ib][k];
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
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int jptr, j, ib, k;
  REAL z, utheta, r, u, v, x0, y0;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];


    if(grid->xv[ib]>1)
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	// Right boundary
	phys->boundary_u[jptr-grid->edgedist[2]][k]=1;
	phys->boundary_v[jptr-grid->edgedist[2]][k]=phys->vc[ib][k];
	phys->boundary_w[jptr-grid->edgedist[2]][k]=phys->w[ib][k];
      }
    else
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	// Left boundary
	phys->boundary_u[jptr-grid->edgedist[2]][k]=1;
	phys->boundary_v[jptr-grid->edgedist[2]][k]=0.0;
	phys->boundary_w[jptr-grid->edgedist[2]][k]=0;
      }
  }

  // Add circulation with utheta nonzero to simulate cylinder with lift
  r = 0.1;
  utheta = 0.0;
  x0 = 0.3;
  y0 = 0.15;
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];
    u = utheta*(grid->ye[j]-y0)/r;
    v = -utheta*(grid->xe[j]-x0)/r;

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->boundary_u[jptr-grid->edgedist[2]][k]=u;
      phys->boundary_v[jptr-grid->edgedist[2]][k]=v;
      phys->boundary_w[jptr-grid->edgedist[2]][k]=0;
    }
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
 * Set the wind stress.
 *
 */
void WindStress(gridT *grid, physT *phys, propT *prop, metT *met, int myproc) {
  int j, jptr;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    phys->tau_T[j]=grid->n2[j]*prop->tau_T;
    phys->tau_B[j]=0;
  }
}

void InitBoundaryData(propT *prop, gridT *grid, int myproc){}
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc){}
void BoundarySediment(gridT *grid, physT *phys, propT *prop) {}
