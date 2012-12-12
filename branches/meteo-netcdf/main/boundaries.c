
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

#ifdef USENETCDF
static void ReadBndNCcoord(int ncid, propT *prop, gridT *grid, int myproc);
static void MatchBndPoints(propT *prop, gridT *grid, int myproc);
void ReadBdyNC(propT *prop, gridT *grid, int myproc);
void UpdateBdyNC(propT *prop, gridT *grid, int myproc);
static size_t returndimlenBC(int ncid, char *dimname);

#endif 

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
  REAL **uc = phys->uc, **vc = phys->vc, **ucold = phys->uold, **vcold = phys->vold;
  REAL z, c0, c1, C0, C1, dt=prop->dt, u0, u0new, uc0, vc0, uc0old, vc0old, ub0;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    for(k=grid->etop[j];k<grid->Nke[j];k++) {    
      ub[j][k] = phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->n1[j] 
    	+ phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->n2[j]; 
    }
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
  int jptr, j, ib, k, jind;
  int iptr, i, ii;
  int nf,ne,neigh;
  REAL z;
  
  //Type-2 zero gradient (Neumann) boundary condition
  /* 
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
	phys->boundary_T[jptr-grid->edgedist[2]][k]=phys->T[ib][k];
	phys->boundary_s[jptr-grid->edgedist[2]][k]=phys->s[ib][k];
      }
  }
  */

  // Type-2
  ii=-1;
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    jind = jptr-grid->edgedist[2];
    j = grid->edgep[jptr];
    ib=grid->grad[2*j];
    ii+=1;

    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_T[jind][k]=bound->boundary_T[bound->ind2[ii]][k];
      phys->boundary_s[jind][k]=bound->boundary_S[bound->ind2[ii]][k];
    }
  }

  //Type-3 with Dirichlet tracer boundary (set edge value to cell value in boundary structure)
  /* ???? (NEEDS TESTING)
  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    ii+=1;
    // Find the edge index of the boundary cell
    for(nf=0;nf<NFACES;nf++){ 
	if(neigh=grid->neigh[i*NFACES+nf]==-1) {
	     ne = grid->edgep[grid->face[i*NFACES+nf]];
	     phys->boundary_T[ne] = bound->T[bound->ind3[ii]][k];
	}
    }
  } 
  */

} // End funciton

/*
 * Function: BoundaryVelocities
 * Usage: BoundaryVelocities(grid,phys,prop);
 * ------------------------------------------
 * This will set the values of u,v,w, and h at the boundaries.
 * 
 */
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc) {
  int i, ii, j, jj, jind, iptr, jptr, n, k;
  REAL u,v,w,h;

  REAL rampfac = 1-exp(-prop->rtime/prop->thetaramptime);//Tidal rampup factor 

   // Update the netcdf boundary data
#ifdef USENETCDF
   if(prop->netcdfBdy==1) 
       UpdateBdyNC(prop,grid,myproc);
#endif

  // Type-2
  ii=-1;
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    jind = jptr-grid->edgedist[2];
    j = grid->edgep[jptr];
    ii+=1;

    phys->boundary_h[jind]=0.0; // Not used??
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
     //printf("ii=%d, bound->ind2[ii]=%d\n",ii,bound->ind2[ii]);
     phys->boundary_u[jind][k]=bound->boundary_u[bound->ind2[ii]][k]*rampfac;
     phys->boundary_v[jind][k]=bound->boundary_v[bound->ind2[ii]][k]*rampfac;
     phys->boundary_w[jind][k]=bound->boundary_w[bound->ind2[ii]][k]*rampfac;
    }
  }
 
  // Type-3
  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    //jind = iptr-grid->celldist[1]+grid->edgedist[3]-grid->edgedist[2];
    i = grid->cellp[iptr];
    ii+=1;

    phys->h[i]=bound->h[bound->ind3[ii]]*rampfac;
    //printf("phys->h[%d] = %f (%f)\n",i,h,phys->h[i]);
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //phys->uc[i][k]=u*(1-exp(-prop->rtime/prop->thetaramptime));
      //phys->vc[i][k]=v*(1-exp(-prop->rtime/prop->thetaramptime));
      //phys->w[i][k]=0;
      phys->uc[i][k]=bound->uc[bound->ind3[ii]][k]*rampfac;
      phys->vc[i][k]=bound->vc[bound->ind3[ii]][k]*rampfac;
      phys->vc[i][k]=bound->wc[bound->ind3[ii]][k]*rampfac;
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
 * !!!!!!!!!!! This should be moved to met.c !!!!!!!!!!!!!!!
 */
void WindStress(gridT *grid, physT *phys, propT *prop, metT *met, int myproc) {
  int j, jptr;
  int Nc=grid->Nc; 
  int i, nf, ne, nc1, nc2, neigh;
  REAL def1, def2;

   if(prop->metmodel>=2){// Interpoalte the spatially variable wind stress onto the edges
       for(i=0;i<Nc;i++){
	  for(nf=0;nf<NFACES;nf++){ 
            if((neigh=grid->neigh[i*NFACES+nf])!=-1) {
	     ne = grid->face[i*NFACES+nf];
	     nc1 = grid->grad[2*ne];
	     nc2 = grid->grad[2*ne+1];
	     
	     def1 = grid->def[nc1*NFACES+grid->gradf[2*ne]];
	     def2 = grid->def[nc2*NFACES+grid->gradf[2*ne+1]];

	     phys->tau_T[ne] = (met->tau_x[nc1]*def1/grid->dg[ne] + met->tau_x[nc2]*def2/grid->dg[ne])*grid->n1[ne] + 
		(met->tau_y[nc1]*def1/grid->dg[ne] + met->tau_y[nc2]*def2/grid->dg[ne])*grid->n2[ne];  
	     phys->tau_T[ne] /= RHO0; 
//	     printf("%3.6e, %3.6e, %3.6e, %3.6f, %3.6f, %3.6f,%3.6f\n",
//	       phys->tau_T[ne], met->tau_x[nc1],met->tau_x[nc2],grid->dg[ne],def1,def2,grid->n1[ne]);
	     //printf("%3.6e, %3.6e, %3.6e\n",met->tau_y[nc1],met->tau_y[nc2],phys->tau_T[ne]);
            }
	  }
       }
   }else{// Set stress to constant

      for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
	j = grid->edgep[jptr];

	phys->tau_T[j]=grid->n2[j]*prop->tau_T;
	phys->tau_B[j]=0;
      }
    }
}

/****************************************************************
 * Beginning of the boundary netCDF-IO functions
 ****************************************************************/
#ifdef USENETCDF
 /*
  * Function: InitBoundaryData()
  * -----------------------------
  * Initialise the boundary condition data for the model
  *
  * This is called from phys.c
  */
 void InitBoundaryData(propT *prop, gridT *grid, int myproc){
   
    // Step 1) Allocate the structure array data
	// Moved to phys.c

    // Step 2) Read in the coordinate info
    if(VERBOSE>2 && myproc==0) printf("Reading netcdf boundary data...\n");
    ReadBndNCcoord(prop->netcdfBdyFileID, prop, grid, myproc);

    // Step 3) Match each boundary point with its local grid point
    MatchBndPoints(prop, grid, myproc);

    // Step 4) Read in the forward and backward time steps into the boundary arrays
    ReadBdyNC(prop, grid, myproc);

   
    
 }//end function

 /*
  * Function: UpdateBdyNC()
  * -----------------------------
  * Update the boundary netcdf data and temporally interpolate onto the model time step
  *
  */     
 void UpdateBdyNC(propT *prop, gridT *grid, int myproc){
     int j, k, t0; 
     REAL dt, r1, r2, mu;
   
     t0 = getTimeRec(prop->nctime,bound->time,bound->Nt);
    
     /* Only interpolate the data onto the grid if need to*/
     if (bound->t0!=t0){
	if(VERBOSE>3 && myproc==0) printf("Updating netcdf boundary data at nc timestep: %d\n",t0);
	/* Read in the data two time steps*/
	ReadBdyNC(prop, grid, myproc);
	bound->t0=t0;
	bound->t1=t0+1;
      }

   /*Linear temporal interpolation coefficients*/
    //dt = bound->time[bound->t1]-bound->time[bound->t0];
    //r2 = (prop->nctime - bound->time[bound->t0])/dt;
    //r1 = 1.0-r2;

    /*Cosine temporal interpolation coefficients*/
    dt = bound->time[bound->t1]-bound->time[bound->t0];
    mu = (prop->nctime - bound->time[bound->t0])/dt;
    r2 = (1.0 - cos(PI*mu))/2.0;
    r1 = 1.0-r2;

    if(bound->hasType2>0){
	for (j=0;j<bound->Ntype2;j++){
	  for (k=0;k<bound->Nk;k++){
	    bound->boundary_u[j][k] = bound->boundary_u_b[j][k]*r1 + bound->boundary_u_f[j][k]*r2;
	    bound->boundary_v[j][k] = bound->boundary_v_b[j][k]*r1 + bound->boundary_v_f[j][k]*r2;
	    bound->boundary_w[j][k] = bound->boundary_w_b[j][k]*r1 + bound->boundary_w_f[j][k]*r2;
	    bound->boundary_T[j][k] = bound->boundary_T_b[j][k]*r1 + bound->boundary_T_f[j][k]*r2;
	    bound->boundary_S[j][k] = bound->boundary_S_b[j][k]*r1 + bound->boundary_S_f[j][k]*r2;
	  }
	}
    }

    if(bound->hasType3>0){
	for (j=0;j<bound->Ntype3;j++){
	  for (k=0;k<bound->Nk;k++){
	    bound->uc[j][k] = bound->uc_b[j][k]*r1 + bound->uc_f[j][k]*r2;
	    bound->vc[j][k] = bound->vc_b[j][k]*r1 + bound->vc_f[j][k]*r2;
	    bound->wc[j][k] = bound->wc_b[j][k]*r1 + bound->wc_f[j][k]*r2;
	    bound->T[j][k] = bound->T_b[j][k]*r1 + bound->T_f[j][k]*r2;
	    bound->S[j][k] = bound->S_b[j][k]*r1 + bound->S_f[j][k]*r2;
	  }
	  bound->h[j] = bound->h_b[j]*r1 + bound->h_f[j]*r2;
	}
    }
 }//End function

  /*
  * Function: ReadBdyNC()
  * -----------------------------
  * Reads in boundary netcdf data into the forward and back time steps 
  *
  */     
 void ReadBdyNC(propT *prop, gridT *grid, int myproc){
    int retval, j,k;
    int t0;
    int varid;
    char *vname;
    size_t start[3];
    size_t start2[2];
    size_t count[]={1,1,1};
    size_t count2[]={1,1};
    int ncid = prop->netcdfBdyFileID;  
    int Nk = bound->Nk;
    int Ntype3 = bound->Ntype3;
    int Ntype2 = bound->Ntype2;

    // Read the first and second time steps
    t0 = getTimeRec(prop->nctime,bound->time,(int)bound->Nt); //this is in met.c
    bound->t0=t0;
    //if(myproc==0) printf("t0 = %d [Nt = %d]\n",t0,bound->Nt);    
    if(bound->hasType2){

	vname = "boundary_u";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype2;j++){
	  for (k=0;k<Nk;k++){
	    //if(myproc==0) printf("t0[%d],k[%d] of %d,j[%d]\n",t0,k,bound->Nk,j);
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_u_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_u_f[j][k]))) 
		ERR(retval); 
	  }
	}

	vname = "boundary_v";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype2;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_v_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_v_f[j][k]))) 
		ERR(retval); 
	  }
	}
	
	vname = "boundary_w";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype2;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_w_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_w_f[j][k]))) 
		ERR(retval); 
	  }
	}

	vname = "boundary_T";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype2;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_T_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_T_f[j][k]))) 
		ERR(retval); 
	  }
	}

	vname = "boundary_S";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype2;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_S_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->boundary_S_f[j][k]))) 
		ERR(retval); 
	  }
	}
    }

    if(bound->hasType3){
	vname = "uc";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype3;j++){
	  for (k=0;k<Nk;k++){
	    //if(myproc==0) printf("t0[%d],k[%d] of %d,j[%d]\n",t0,k,bound->Nk,j);
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->uc_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->uc_f[j][k]))) 
		ERR(retval); 
	  }
	}

	vname = "vc";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype3;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->vc_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->vc_f[j][k]))) 
		ERR(retval); 
	  }
	}
	
	vname = "wc";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype3;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->wc_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->wc_f[j][k]))) 
		ERR(retval); 
	  }
	}

	vname = "T";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype3;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->T_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->T_f[j][k]))) 
		ERR(retval); 
	  }
	}

	vname = "S";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype3;j++){
	  for (k=0;k<Nk;k++){
	    start[0]=t0;
	    start[1]=k;
	    start[2]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->S_b[j][k]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start, count, &bound->S_f[j][k]))) 
		ERR(retval); 
	  }
	}

	vname = "h";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	for (j=0;j<Ntype3;j++){
	    start2[0]=t0;
	    start2[1]=j;
	    if ((retval = nc_get_vara_double(ncid, varid, start2, count2, &bound->h_b[j]))) 
		ERR(retval); 
	    start[0]=t0+1;
	    if ((retval = nc_get_vara_double(ncid, varid, start2, count2, &bound->h_f[j]))) 
		ERR(retval); 
	}
     }// End read type-3
 }//End function

/*
 * Function: MatchBndPoints()
 * -------------------------
 * Checks that the boundary arrays match the grid sizes
 */
 static void MatchBndPoints(propT *prop, gridT *grid, int myproc){
 int iptr, jptr, jj, ii, ne, nc1, nc2, nc;

    if(myproc==0) printf("SUNTANS grid # type 2 points = %d\n",grid->edgedist[3] - grid->edgedist[2]);
    if(myproc==0) printf("SUNTANS grid # type 3 points = %d\n",grid->celldist[2] - grid->celldist[1]);
    if(myproc==0) printf("Boundary NetCDF file grid # type 2 points = %d\n",(int)bound->Ntype2);
    if(myproc==0) printf("Boundary NetCDF file grid # type 3 points = %d\n",(int)bound->Ntype3);

     //Type-2
     ii=-1;
     for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
	 //printf("Type 2 : Processor = %d, jptr = %d, edgep[jptr]\n",myproc,jptr,grid->edgep[jptr]);
	 // Need to find the edge index from the unpartitioned grid
	 /*ne = grid->edgep[jpr]; // Edge index for partitioned grid
	 nc1 = grid->grad[2*ne];
         nc2 = grid->grad[2*ne+1];
	 if(nc1==-1) nc1=nc2; // nc1 is the cell index for the partitioned grid
	 
	 nc = grid->mnptr[nc1]; // Cell index for the unpartitioned grid
        // Find the edge index of the boundary cell
	for(nf=0;nf<NFACES;nf++){ 
	    if(neigh=grid->neigh[i*NFACES+nf]==-1) {
		 ne = grid->edgep[grid->face[i*NFACES+nf]];
	}
	*/
	 ii+=1;
	 // Match suntans edge cell with the type-2 boundary file point
	 for(jj=0;jj<bound->Ntype2;jj++){
	     if(grid->eptr[grid->edgep[jptr]]==bound->edgep[jj])
		bound->ind2[ii]=jj;
	 }
	 printf("Type 2 : Processor = %d, jptr = %d, edgep[jptr]=%d, bound->ind2[ii]=%d\n",myproc,jptr,grid->edgep[jptr],bound->ind2[ii]);
     }	 
     // Type-3
     ii=-1;
     for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
         ii+=1;
	 // Match suntans grid cell with the type-3 boundary file point
	 for(jj=0;jj<bound->Ntype3;jj++){
	     if(grid->mnptr[grid->cellp[iptr]]==bound->cellp[jj])
		bound->ind3[ii]=jj;
	 }
	 //printf("Type 3 : Processor = %d, jptr = %d, cellp[jptr]=%d, bound->ind3[ii]=%d\n",myproc,jptr,grid->cellp[jptr],bound->ind3[ii]);
     }
     // Check that ind2 and ind3 do not contain any -1 (non-matching points)

     // Check that the number of vertical grid points match

 } // End function

/*
 * Function: ReadBndNCcoord()
 * --------------------------
 * Reads the coordinate information from the netcdf file into the boundary structure
 *
 */
static void ReadBndNCcoord(int ncid, propT *prop, gridT *grid, int myproc){

    int retval, j;
    int varid;
    char *vname;
    //int ncid = prop->netcdfBdyFileID; 

    vname = "time";
    if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,bound->time))) 
      ERR(retval); 
    if(VERBOSE>2 && myproc==0) printf("done.\n");

    vname = "z";
    if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,bound->z))) 
      ERR(retval); 
    if(VERBOSE>2 && myproc==0) printf("done.\n");

    if(bound->hasType3>0){

	vname = "xv";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->xv))) 
	  ERR(retval); 
	if(VERBOSE>2 && myproc==0) printf("done.\n");

	vname = "yv";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->yv))) 
	      ERR(retval); 
	if(VERBOSE>2 && myproc==0) printf("done.\n");

	vname = "cellp";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->cellp))) 
	  ERR(retval); 
	if(VERBOSE>2 && myproc==0) printf("done.\n");
    }//end if

    if(bound->hasType2>0){

	vname = "xe";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->xe))) 
	  ERR(retval); 

	vname = "ye";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->ye))) 
	  ERR(retval); 

	vname = "edgep";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->edgep))) 
	  ERR(retval); 
    }//end if
 }//End function

/*
 * Function: AllocateBoundaryData()
 * -------------------------------
 * Allocate boundary structure arrays.
 */
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc){

     int Ntype2, Ntype3, Nt, Nk;
     int j, k, i;
     int n3, n2; // Number of type 2 and 3 on a processor's grid

    //Allocate memory
    if(VERBOSE>1 && myproc==0) printf("Allocating boundary data structure...\n");
      *bound = (boundT *)SunMalloc(sizeof(boundT),"AllocateBoundaryData");

    // Read in the dimensions
    if(myproc==0) printf("Reading boundary netcdf dimensions...\n");
    if(myproc==0) printf("Reading dimension: Ntype3...\n");
    (*bound)->Ntype3 = returndimlenBC(prop->netcdfBdyFileID,"Ntype3");
        if(myproc==0) printf("Reading dimension: Ntype2...\n");
    (*bound)->Ntype2 = returndimlenBC(prop->netcdfBdyFileID,"Ntype2");
        if(myproc==0) printf("Reading dimension: Nt...\n");
    (*bound)->Nt = returndimlenBC(prop->netcdfBdyFileID,"Nt");
    if(myproc==0) printf("Reading dimension: Nk...\n");
    (*bound)->Nk = returndimlenBC(prop->netcdfBdyFileID,"Nk");
    Ntype3 = (*bound)->Ntype3;
    Ntype2 = (*bound)->Ntype2;
    Nt = (*bound)->Nt;
    Nk = (*bound)->Nk;

    // Check if boundary types are in the file
    if ((*bound)->Ntype3==0){
    	(*bound)->hasType3=0;
    }else{
    	(*bound)->hasType3=1;
    }

    if ((*bound)->Ntype2==0){
    	(*bound)->hasType2=0;
    }else{
    	(*bound)->hasType2=1;
    }
    
    // Print the array sizes
    if(VERBOSE>1 && myproc==0){
	printf("Ntype 3 = %d\n",Ntype3);
	printf("Ntype 2 = %d\n",Ntype2);
	printf("Nt = %d\n",Nt);
	printf("Nk = %d\n",Nk);
	printf("hasType2 = %d\n",(*bound)->hasType2);
	printf("hasType3 = %d\n",(*bound)->hasType3);
    }
    // Allocate the type2 arrays
    if((*bound)->hasType2==1){
	(*bound)->edgep = (int *)SunMalloc(Ntype2*sizeof(int),"AllocateBoundaryData");
	(*bound)->xe = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->ye = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->boundary_u = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_v = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_w = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_T = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_S = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_u_f = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_v_f = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_w_f = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_T_f = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_S_f = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_u_b = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_v_b = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_w_b = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_T_b = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_S_b = (REAL **)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");

	for(j=0;j<Ntype2;j++){
	    (*bound)->boundary_u[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_v[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_w[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_T[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_S[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_u_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_v_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_w_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_T_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_S_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_u_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_v_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_w_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_T_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_S_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	}
	// Allocate the pointer arrays for each processor's grid
	n2 = grid->edgedist[3]-grid->edgedist[2];
	(*bound)->ind2 = (int *)SunMalloc(n2*sizeof(REAL),"AllocateBoundaryData");
    }//endif

    // Allocate the type3 arrays   
    if((*bound)->hasType3==1){
	(*bound)->cellp = (int *)SunMalloc(Ntype3*sizeof(int),"AllocateBoundaryData");
	(*bound)->xv = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->yv = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->h = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->h_f = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->h_b = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->uc = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->vc = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->wc = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->T = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->S = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->uc_f = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->vc_f = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->wc_f = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->T_f = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->S_f = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->uc_b = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->vc_b = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->wc_b = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->T_b = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->S_b = (REAL **)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");

	for(j=0;j<Ntype3;j++){
	    (*bound)->uc[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->vc[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->wc[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->T[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->S[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->uc_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->vc_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->wc_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->T_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->S_f[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->uc_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->vc_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->wc_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->T_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->S_b[j] = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	}
	// Allocate the pointer arrays for each processor's grid
	n3 = grid->celldist[2]-grid->celldist[1];
	(*bound)->ind3 = (int *)SunMalloc(n3*sizeof(REAL),"AllocateBoundaryData");

    }//endif

    (*bound)->time = (REAL *)SunMalloc(Nt*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->z = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");

    // Zero the boundary arrays

    for(j=0;j<(*bound)->Nt;j++){
	(*bound)->time[j]=0.0;
    }
    for(k=0;k<(*bound)->Nk;k++){
	(*bound)->z[k]=0.0;
    }

    if((*bound)->hasType2==1){
	for(j=0;j<Ntype2;j++){
	    (*bound)->edgep[j]=0;
	    (*bound)->xe[j]=0.0;
	    (*bound)->ye[j]=0.0;
	    for(k=0;k<Nk;k++){
		(*bound)->boundary_u[j][k]=0.0;
		(*bound)->boundary_v[j][k]=0.0;
		(*bound)->boundary_w[j][k]=0.0;
		(*bound)->boundary_T[j][k]=0.0;
		(*bound)->boundary_S[j][k]=0.0;
		(*bound)->boundary_u_f[j][k]=0.0;
		(*bound)->boundary_v_f[j][k]=0.0;
		(*bound)->boundary_w_f[j][k]=0.0;
		(*bound)->boundary_T_f[j][k]=0.0;
		(*bound)->boundary_S_f[j][k]=0.0;
		(*bound)->boundary_u_b[j][k]=0.0;
		(*bound)->boundary_v_b[j][k]=0.0;
		(*bound)->boundary_w_b[j][k]=0.0;
		(*bound)->boundary_T_b[j][k]=0.0;
		(*bound)->boundary_S_b[j][k]=0.0;
	    }
	}
	for(i=0;i<n2;i++){
	    (*bound)->ind2[i]=-1;
	}
    }

    if((*bound)->hasType3==1){
	for(j=0;j<Ntype3;j++){
	    (*bound)->cellp[j]=0;
	    (*bound)->xv[j]=0.0;
	    (*bound)->yv[j]=0.0;

	    (*bound)->h[j]=0.0;
	    (*bound)->h_f[j]=0.0;
            (*bound)->h_b[j]=0.0;
	    for(k=0;k<Nk;k++){
		(*bound)->uc[j][k]=0.0;
		(*bound)->vc[j][k]=0.0;
		(*bound)->wc[j][k]=0.0;
		(*bound)->T[j][k]=0.0;
		(*bound)->S[j][k]=0.0;
		(*bound)->uc_f[j][k]=0.0;
		(*bound)->vc_f[j][k]=0.0;
		(*bound)->wc_f[j][k]=0.0;
		(*bound)->T_f[j][k]=0.0;
		(*bound)->S_f[j][k]=0.0;
		(*bound)->uc_b[j][k]=0.0;
		(*bound)->vc_b[j][k]=0.0;
		(*bound)->wc_b[j][k]=0.0;
		(*bound)->T_b[j][k]=0.0;
		(*bound)->S_b[j][k]=0.0;
	    }
	}
	for(i=0;i<n3;i++){
	    (*bound)->ind3[i]=-1;
	}
    }
 }//End function


 /*
 * Function: returndimlenBC()
 * --------------------------
 * Returns the length of a dimension 
 * Returns a zero if the dimension is not found and does not raise an error
 */
static size_t returndimlenBC(int ncid, char *dimname){
 int retval;
 int dimid;
 size_t dimlen;
 
 if ((retval =nc_inq_dimid(ncid,dimname,&dimid)))
    return 0;
 
 if ((retval = nc_inq_dimlen(ncid,dimid, &dimlen)))
    ERR(retval);
 return dimlen;
} // End function
#endif
