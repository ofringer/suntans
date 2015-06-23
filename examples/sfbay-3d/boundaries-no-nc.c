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

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
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
  int i, ib, iptr, k, j, jptr;

  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++)  {
    i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	phys->s[i][k]=32;
  }

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++)  {
    j = grid->edgep[jptr];
    ib = grid->grad[2*j];
  
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
      phys->boundary_s[jptr-grid->edgedist[2]][k]=0;
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
  int i, ib, iptr, k, j, jptr;
  REAL A_sac, A_san, Q_sac, Q_san, ub;

  Q_sac = 300;
  Q_san = 300;
  A_sac = 0;
  A_san = 0;

  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++)  {
    i = grid->cellp[iptr];

    phys->h[i]=-5+prop->amp*sin(prop->omega*prop->rtime);
  }

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++)  {
    j = grid->edgep[jptr];
    ib = grid->grad[2*j];
  
    if(grid->yv[ib]>4.21e6) 
      A_sac+=(grid->dv[ib]+phys->h[ib])*grid->df[j];
    else
      A_san+=(grid->dv[ib]+phys->h[ib])*grid->df[j];
  }
  
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++)  {
    j = grid->edgep[jptr];
    ib = grid->grad[2*j];

    if(grid->yv[ib]>4.21e6) 
      ub=Q_sac/A_sac;
    else
      ub=Q_san/A_san;

    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_u[jptr-grid->edgedist[2]][k]=-ub;
      phys->boundary_v[jptr-grid->edgedist[2]][k]=0;
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
  REAL omegaWind = prop->omega/2;  // Fake diurnal wind (2*M2)

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    phys->tau_T[j]=grid->n2[j]*prop->tau_T*cos(omegaWind*prop->rtime);
    phys->tau_B[j]=0;
  }
}

void InitBoundaryData(propT *prop, gridT *grid, int myproc, MPI_Comm comm){}
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc, MPI_Comm comm){}

/*
 * Function: BoundarySediment
 * Usage: BoundarySediment(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the suspended sediment concentration
 * at the open boundaries.
 * 
 */
void BoundarySediment(gridT *grid, physT *phys, propT *prop) {
  int jptr, j, ib, k,nosize,i,iptr;
  REAL z;

  // At the upstream boundary
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j=grid->edgep[jptr];
    ib=grid->grad[2*j];
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
        sediments->boundary_sediC[nosize][jptr-grid->edgedist[2]][k]=200;
      }
    }
  }

  // At the ocean boundary
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      for(k=0;k<grid->ctop[i];k++) {
        sediments->SediC[nosize][i][k]=0;
      } 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        sediments->SediC[nosize][i][k]=0;
      }
    }
  }
}
