
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
#include "mynetcdf.h"

// Local functions
static void GetBoundaryVelocity(REAL *ub, int *forced, REAL x, REAL y, REAL t, REAL h, REAL d, REAL omega, REAL amp);
static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag);

static void MatchBndPoints(propT *prop, gridT *grid, int myproc);
static void FluxtoUV(propT *prop, gridT *grid, int myproc,MPI_Comm comm);
static void SegmentArea(propT *prop, gridT *grid, int myproc, MPI_Comm comm);
int isGhostEdge(int j, gridT *grid, int myproc);
int getTimeRecBnd(REAL nctime, REAL *time, int nt);

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
void BoundaryScalars(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
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
      phys->boundary_T[jind][k]=bound->boundary_T[k][bound->ind2[ii]];
      phys->boundary_s[jind][k]=bound->boundary_S[k][bound->ind2[ii]];
    }
  }

  //Type-3 
  // ???? (NEEDS TESTING)
  // Set boundary cell to value in the boundary structure

  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    ii+=1;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
//    for(k=0;k<grid->Nk[i];k++){//Go from the very top
	 phys->T[i][k] = bound->T[k][bound->ind3[ii]];
	 phys->s[i][k] = bound->S[k][bound->ind3[ii]];
    }
   }
  // Need to communicate the cell data for type 3 boundaries
  ISendRecvCellData3D(phys->T,grid,myproc,comm);
  ISendRecvCellData3D(phys->s,grid,myproc,comm);

   // Set the edge array to the value in the boundary array
  /* ii=-1;
 int myproc,int myproc,   for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++) {
    jind = jptr-grid->edgedist[2];
    j = grid->edgep[jptr];
    ib=grid->grad[2*j];
    ii+=1;// This is the same as jind...
    //printf("ii: %d, jptr: %d, jind: %d, j: %d,ind3edge[ii]: %d\n",ii, jptr,jind,j,bound->ind3edge[ii]);
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_T[jind][k]=bound->T[bound->ind3edge[ii]][k];
      phys->boundary_s[jind][k]=bound->S[bound->ind3edge[ii]][k];
    }
  }
  /*
    // Find the edge index of the boundary cell
    /*
    for(nf=0;nf<NFACES;nf++){ 
	if(neigh=grid->neigh[i*NFACES+nf]==-1) {
	     ne = grid->edgep[grid->face[i*NFACES+nf]];
	     printf("ne: %d, grid->edgedist[3]: %d\n",ne,grid->edgedist[3]);
	     for(k=grid->ctop[i];k<grid->Nk[i];k++) {
		 phys->boundary_T[grid->edgedist[3]-ne][k] = bound->T[bound->ind3[ii]][k];
		 phys->boundary_s[grid->edgedist[3]-ne][k] = bound->S[bound->ind3[ii]][k];
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
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int i, ii, j, jj, jind, iptr, jptr, n, k;
  REAL u,v,w,h;

  REAL rampfac = 1-exp(-prop->rtime/prop->thetaramptime);//Tidal rampup factor 

   // Test
   REAL amp = 0.25;
   REAL omega = 7.27e-5;

   // Update the netcdf boundary data
   if(prop->netcdfBdy==1) 
       UpdateBdyNC(prop,grid,myproc,comm);

  // Type-2
  ii=-1;
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    jind = jptr-grid->edgedist[2];
    j = grid->edgep[jptr];
    ii+=1;

    phys->boundary_h[jind]=0.0; // Not used??
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
     phys->boundary_u[jind][k]=bound->boundary_u[k][bound->ind2[ii]]*rampfac;
     phys->boundary_v[jind][k]=bound->boundary_v[k][bound->ind2[ii]]*rampfac;
     phys->boundary_w[jind][k]=bound->boundary_w[k][bound->ind2[ii]]*rampfac;
    }
  }
 
  // Type-3
/*
  // Flather condition. This updates the boundary velocity values
  REAL dh, Umod, Ubc, gfac; 
  int ne, nf, neigh;
  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    ii+=1;
    gfac = sqrt(9.81/grid->dv[i]);
    dh = phys->h[i] - bound->h[bound->ind3[ii]]*rampfac;
    for(nf=0;nf<NFACES;nf++){ 
      if((neigh=grid->neigh[i*NFACES+nf])==-1) {
        ne = grid->face[i*NFACES+nf];
	for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	    Ubc = bound->uc[k][bound->ind3[ii]]*rampfac*grid->n1[ne] + bound->vc[k][bound->ind3[ii]]*rampfac*grid->n2[ne]; 
	    Umod = Ubc - gfac*dh;
	    bound->uc[k][bound->ind3[ii]]=Umod*grid->n1[ne]; 
	    bound->vc[k][bound->ind3[ii]]=Umod*grid->n2[ne]; 
	}
      }
    }
  }//End Flather condition. Everything should be the same after here...
*/
  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    ii+=1;

    phys->h[i]=bound->h[bound->ind3[ii]]*rampfac;
    //    phys->h[i]=bound->h[bound->ind3[ii]]*sin(omega*prop->rtime)*rampfac;// Test harmonic netcdf type bc
    //    phys->h[i]=amp*sin(omega*prop->rtime)*rampfac; // Test analytical BC
  }

  // Try recalculating the cell tops here
  ISendRecvCellData2D(phys->h,grid,myproc,comm);
  UpdateDZ(grid,phys,prop,0);

  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    ii+=1;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->uc[i][k]=bound->uc[k][bound->ind3[ii]]*rampfac;
      phys->vc[i][k]=bound->vc[k][bound->ind3[ii]]*rampfac;
      phys->wc[i][k]=bound->wc[k][bound->ind3[ii]]*rampfac;
    }
  }
  // Need to communicate the cell data for type 3 boundaries
  ISendRecvCellData3D(phys->uc,grid,myproc,comm);
  ISendRecvCellData3D(phys->vc,grid,myproc,comm);
  ISendRecvCellData3D(phys->wc,grid,myproc,comm);

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

    REAL rampfac = 1-exp(-prop->rtime/prop->thetaramptime);// rampup factor 
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
	     phys->tau_T[ne] *= rampfac;
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
    if(VERBOSE>2 && myproc==0) printf("Reading netcdf boundary coordinate data...\n");
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
 void UpdateBdyNC(propT *prop, gridT *grid, int myproc, MPI_Comm comm){
     int j, k, t0, t1, t2; 
     REAL dt, r1, r2, mu;
   
     t1 = getTimeRecBnd(prop->nctime,bound->time,bound->Nt);
     t0=t1-1;
     t2=t1+1;

     /* Only update the boundary data if need to*/
     if (bound->t1!=t1){
	if(VERBOSE>3 && myproc==0) printf("Updating netcdf boundary data at nc timestep: %d\n",t1);
	/* Read in the data two time steps*/
	bound->t1=t1;
	bound->t2=t2;
	bound->t0=t0;
        //printf("myproc: %d, bound->t0: %d, nctime: %f, rtime: %f \n",myproc,bound->t0, prop->nctime, prop->rtime);
	ReadBdyNC(prop, grid, myproc);
      }

   /*Linear temporal interpolation coefficients*/
    //dt = bound->time[bound->t1]-bound->time[bound->t0];
    //r2 = (prop->nctime - bound->time[bound->t0])/dt;
    //r1 = 1.0-r2;

    /*Cosine temporal interpolation coefficients*/
    //dt = bound->time[bound->t1]-bound->time[bound->t0];
    //mu = (prop->nctime - bound->time[bound->t0])/dt;
    //r2 = (1.0 - cos(PI*mu))/2.0;
    //r1 = 1.0-r2;

//    if (myproc==0) printf("t1: %f, %t2: %f, tnow:%f, ,mu %f, r2: %f, r1: %f\n",bound->time[bound->t0],bound->time[bound->t1], prop->nctime, mu, r2, r1);

    if(bound->hasType2>0){
	for (j=0;j<bound->Ntype2;j++){
	  for (k=0;k<bound->Nk;k++){
	    //Quadratic temporal interpolation
	    bound->boundary_u[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->boundary_u_t[0][k][j],bound->boundary_u_t[1][k][j],bound->boundary_u_t[2][k][j] );
	    bound->boundary_v[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->boundary_v_t[0][k][j],bound->boundary_v_t[1][k][j],bound->boundary_v_t[2][k][j] );
	    bound->boundary_w[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->boundary_w_t[0][k][j],bound->boundary_w_t[1][k][j],bound->boundary_w_t[2][k][j] );
	    bound->boundary_T[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->boundary_T_t[0][k][j],bound->boundary_T_t[1][k][j],bound->boundary_T_t[2][k][j] );
	    bound->boundary_S[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->boundary_S_t[0][k][j],bound->boundary_S_t[1][k][j],bound->boundary_S_t[2][k][j] );
	  }
	}
    }

    if(bound->hasType3>0){
	for (j=0;j<bound->Ntype3;j++){
	  for (k=0;k<bound->Nk;k++){
	    // Quadratic temporal interpolation
	    bound->uc[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->uc_t[0][k][j],bound->uc_t[1][k][j],bound->uc_t[2][k][j]);
	    bound->vc[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->vc_t[0][k][j],bound->vc_t[1][k][j],bound->vc_t[2][k][j]);
	    bound->wc[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->wc_t[0][k][j],bound->wc_t[1][k][j],bound->wc_t[2][k][j]);
	    bound->T[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->T_t[0][k][j],bound->T_t[1][k][j],bound->T_t[2][k][j]);
	    bound->S[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->S_t[0][k][j],bound->S_t[1][k][j],bound->S_t[2][k][j]);
	  }
	  bound->h[j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->h_t[0][j],bound->h_t[1][j],bound->h_t[2][j]);
	}
    }
    // Interpolate Q and find the velocities based on the (dynamic) segment area
    if(bound->hasSeg>0){
	for (j=0;j<bound->Nseg;j++){
	  bound->boundary_Q[j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->boundary_Q_t[0][j],bound->boundary_Q_t[1][j],bound->boundary_Q_t[2][j]);
	}
	// This function does the actual conversion
	FluxtoUV(prop,grid,myproc,comm);
    }

 }//End function

 /*
  * Function: FluxtoUV()
  * -----------------------------
  * Converts boundary flux information into u and v information 
  *
  */     
static void FluxtoUV(propT *prop, gridT *grid, int myproc, MPI_Comm comm){
    int ii, j, k, jptr, jind;
    int ss,n;
    REAL dz;

    // Step 1) Find the total area of each segment
    SegmentArea(prop,grid,myproc,comm);
   //for (n=0;n<bound->Nseg;n++){
   //	 printf("Processor: %d,  segment #: %d, segment ID: %d, segment area: %10.6f [m2]\n",myproc,n,bound->segp[n],bound->segarea[n]); 
   // }
    //Step 2) Loop through again and calculate the velocity for each type-2 cell 
    for(j=0;j<bound->Ntype2;j++){
    	jind = bound->localedgep[j];

	//Only calculate if the segment flag >0 
	if(bound->segedgep[j]>0 && jind!=-1 && !isGhostEdge(jind,grid,myproc) ){
	    //Find the segment index
	    ss = -1;
	    for (n=0;n<bound->Nseg;n++){
		if( bound->segp[n] == bound->segedgep[j]) 
		    ss=n;
	    }

	    //printf("Processor: %d, j: %d, jind: %d, mark: %d, segment #: %d, segment ID: %d, segment area: %10.6f [m2]\n",myproc,j,jind,grid->mark[jind],ss,bound->segp[ss],bound->segarea[ss]); 
	    //isGhostEdge(jind,grid,myproc);

	    for(k=grid->etop[jind];k<grid->Nke[jind];k++) {
		bound->boundary_u[k][j] = grid->n1[jind] * bound->boundary_Q[ss] / bound->segarea[ss]; 
		bound->boundary_v[k][j] = grid->n2[jind] * bound->boundary_Q[ss] / bound->segarea[ss]; 
	    }
	}//End if
    }
}//End function FluxtoUV

/*
 * Function: SegmentArea()
 * ----------------------------
 * Calculates the area of each boundary segment
 * Note that this needs to be done across processors...
 */
static void SegmentArea(propT *prop, gridT *grid, int myproc, MPI_Comm comm){
    int ii, j, k, jptr, jind;
    int ss,n;
    REAL dz;
    
    //Zero the area first
    for (n=0;n<bound->Nseg;n++){
	bound->localsegarea[n]=0.0;
	//bound->segarea[n]=0.0;
    }
    for(j=0;j<bound->Ntype2;j++){
        jind = bound->localedgep[j]; // the local edge pointer will = -1 if the boundary cell is not on the processor's grid
	//Only calculate if the segment flag >0 
	if(bound->segedgep[j]>0 && jind != -1 && !isGhostEdge(jind,grid,myproc) ){
	    //Find the segment index
	    ss = -1;
	    for (n=0;n<bound->Nseg;n++){
		if(bound->segp[n] ==bound->segedgep[j]) 
		    ss=n;
	    }

	    //Loop through all vertical cells and sum the area
	    for(k=grid->etop[jind];k<grid->Nke[jind];k++) {
	    	//dz = grid->dzf[jind][k];
		dz = grid->dzz[grid->grad[2*jind]][k];//Updwind cell height
		bound->localsegarea[ss]+= grid->df[jind] * dz ;
	//	 printf("Processor: %d, k: %d, j: %d,  segment #: %d, segment ID: %d, df = %f, dzf = %f, segment area: %10.6f [m2]\n",myproc,k, j, ss,bound->segp[ss],grid->df[j],dz,bound->segarea[ss]); 
	    }
	}//end if
    }

    //Sum the area across processors
    MPI_Reduce(bound->localsegarea, bound->segarea,(int)bound->Nseg, MPI_DOUBLE, MPI_SUM, 0,comm);
    MPI_Bcast(bound->segarea,(int)bound->Nseg, MPI_DOUBLE,0,comm);

}//End function SegmentArea()

/*
* Function: isGhostEdge()
* ----------------------------
* Find whether or not an edge is part of a ghost cell on a particular processor
* I think the quickest test is to find whether one of the other edges has mark==6.
*/

int isGhostEdge(int j, gridT *grid, int myproc){
    int ib, isGhost, nf, ne;
    
    isGhost=0;
    // Cell index
    ib = grid->grad[2*j];
    // Loop through the cell edges
    //printf("Processor: %d, j: %d, ",myproc,j);
    for(nf=0;nf<NFACES;nf++){ 
	// Edge index
	ne = grid->face[ib*NFACES+nf];
	//printf("ne: %d, mark: %d, ",ne,grid->mark[ne]);
	if (grid->mark[ne]==6)
	    isGhost=1;
    }
    //printf("\n");

    return isGhost;

}//End function isGhostEdge


/*
 * Function: MatchBndPoints()
 * -------------------------
 * Checks that the boundary arrays match the grid sizes
 */
 static void MatchBndPoints(propT *prop, gridT *grid, int myproc){
 int iptr, jptr, jj, ii, ne, nc1, nc2, nc, j, ib;

    if(myproc==0) printf("Boundary NetCDF file grid # type 2 points = %d\n",(int)bound->Ntype2);
    if(myproc==0) printf("Boundary NetCDF file grid # type 3 points = %d\n",(int)bound->Ntype3);

     //Type-2
     ii=-1;
     for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
	 ii+=1;
	 // Match suntans edge cell with the type-2 boundary file point
	 for(jj=0;jj<bound->Ntype2;jj++){
	     if(grid->eptr[grid->edgep[jptr]]==bound->edgep[jj]){
		bound->ind2[ii]=jj;
		bound->localedgep[jj]=grid->edgep[jptr]; 
	        //printf("grid->eptr:%d, bound->edgep[jj]: %d, jj: %d, jptr: %d, grid->edgep[jptr]: %d, localedgep[jj]: %d\n",grid->eptr[grid->edgep[jptr]],bound->edgep[jj],jj,jptr,grid->edgep[jptr],bound->localedgep[jj]); 
	     }
	 }	 
     }
     // Type-3
     ii=-1;
     for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
         ii+=1;
	 // Match suntans grid cell with the type-3 boundary file point
	 for(jj=0;jj<bound->Ntype3;jj++){
	     if(grid->mnptr[grid->cellp[iptr]]==bound->cellp[jj]){
		bound->ind3[ii]=jj;
	     }
	 }
	 //printf("Type 3 : Processor = %d, jptr = %d, cellp[jptr]=%d, bound->ind3[ii]=%d\n",myproc,jptr,grid->cellp[jptr],bound->ind3[ii]);
     }
     //Type-3 edges
     ii=-1;
     for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++) {
	 ii+=1;
	 j = grid->edgep[jptr];
	 ib=grid->mnptr[grid->grad[2*j]];
	 // Match suntans edge cell with the type-3 boundary file point
	 for(jj=0;jj<bound->Ntype3;jj++){
	     if(ib==bound->cellp[jj]){
		bound->ind3edge[ii]=jj;
	     }
	 }	 
     }
     // Check that ind2 and ind3 do not contain any -1 (non-matching points)

     // Check that the number of vertical grid points match
    if(bound->Nk != grid->Nkmax){
	printf("Error! Number of layers in open boundary file (%d) not equal to Nkmax (%d).\n",bound->Nk,grid->Nkmax); 
	MPI_Finalize();
        exit(EXIT_FAILURE);
    }
 } // End function


/*
 * Function: AllocateBoundaryData()
 * -------------------------------
 * Allocate boundary structure arrays.
 */
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc){

     int Ntype2, Ntype3, Nseg, Nt, Nk;
     int j, k, i, n;
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
    if(myproc==0) printf("Reading dimension: Nseg...\n");
    (*bound)->Nseg = returndimlenBC(prop->netcdfBdyFileID,"Nseg");
    if(myproc==0) printf("Reading dimension: Nt...\n");
    (*bound)->Nt = returndimlenBC(prop->netcdfBdyFileID,"Nt");
    if(myproc==0) printf("Reading dimension: Nk...\n");
    (*bound)->Nk = returndimlenBC(prop->netcdfBdyFileID,"Nk");
    Ntype3 = (*bound)->Ntype3;
    Nseg = (*bound)->Nseg;
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
    
    if ((*bound)->Nseg==0){
    	(*bound)->hasSeg=0;
    }else{
    	(*bound)->hasSeg=1;
    }
    // Print the array sizes
    if(VERBOSE>1 && myproc==0){
	printf("Ntype 3 = %d\n",Ntype3);
	printf("Ntype 2 = %d\n",Ntype2);
	printf("Nseg= %d\n",Nseg);
	printf("Nt = %d\n",Nt);
	printf("Nk = %d\n",Nk);
	printf("hasType2 = %d\n",(*bound)->hasType2);
	printf("hasType3 = %d\n",(*bound)->hasType3);
	printf("hasSeg = %d\n",(*bound)->hasSeg);
    }
    // Allocate the type2 arrays
    if((*bound)->hasType2==1){
	(*bound)->edgep = (int *)SunMalloc(Ntype2*sizeof(int),"AllocateBoundaryData");
	(*bound)->localedgep = (int *)SunMalloc(Ntype2*sizeof(int),"AllocateBoundaryData");
	(*bound)->xe = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->ye = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->boundary_u = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_v = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_w = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_T = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_S = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->boundary_u_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_v_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_w_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_T_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_S_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");

	for(k=0;k<Nk;k++){
	    (*bound)->boundary_u[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_v[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_w[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_T[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_S[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	}
	for(n=0;n<NT;n++){
	    (*bound)->boundary_u_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_v_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_w_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_T_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->boundary_S_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    for(k=0;k<Nk;k++){
		 (*bound)->boundary_u_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
		 (*bound)->boundary_v_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
		 (*bound)->boundary_w_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
		 (*bound)->boundary_T_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
		 (*bound)->boundary_S_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
	    }
	}
	// Allocate the pointer arrays for each processor's grid
	n2 = grid->edgedist[3]-grid->edgedist[2];
	(*bound)->ind2 = (int *)SunMalloc(n2*sizeof(REAL),"AllocateBoundaryData");
    }//endif

    // Allocate the segment arrays (subset of type2)
    if ((*bound)->hasSeg==1){
	(*bound)->segp = (int *)SunMalloc(Nseg*sizeof(int),"AllocateBoundaryData");
	(*bound)->segedgep = (int *)SunMalloc(Ntype2*sizeof(int),"AllocateBoundaryData");

	(*bound)->segarea = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->localsegarea = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_Q  = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->boundary_Q_t  = (REAL **)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	for(n=0;n<NT;n++){
	    (*bound)->boundary_Q_t[n] = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
	}
    }
    // Allocate the type3 arrays   
    if((*bound)->hasType3==1){
	(*bound)->cellp = (int *)SunMalloc(Ntype3*sizeof(int),"AllocateBoundaryData");
	(*bound)->xv = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->yv = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->h = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->h_t = (REAL **)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->uc = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->vc = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->wc = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->T = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->S = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");

	(*bound)->uc_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->vc_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->wc_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->T_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->S_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");

	for(k=0;k<Nk;k++){
	    (*bound)->uc[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->vc[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->wc[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->T[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->S[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	}
	for(n=0;n<NT;n++){
	    (*bound)->uc_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->vc_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->wc_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->T_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    (*bound)->S_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
	    for(k=0;k<Nk;k++){
		(*bound)->uc_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
		(*bound)->vc_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
		(*bound)->wc_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
		(*bound)->T_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
		(*bound)->S_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	    }
	    (*bound)->h_t[n] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
	}
	// Allocate the pointer arrays for each processor's grid
	n3 = grid->celldist[2]-grid->celldist[1];
	(*bound)->ind3 = (int *)SunMalloc(n3*sizeof(REAL),"AllocateBoundaryData");
	(*bound)->ind3edge = (int *)SunMalloc(n3*sizeof(REAL),"AllocateBoundaryData");

    }//endif

    (*bound)->time = (REAL *)SunMalloc(Nt*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->z = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");

    // Zero the boundary arrays
    if(myproc==0) printf("Zeroing boundary arrays...\n");

    for(j=0;j<(*bound)->Nt;j++){
	(*bound)->time[j]=0.0;
    }
    (*bound)->t0=-1;
    (*bound)->t1=-1;
    (*bound)->t2=-1;
    for(k=0;k<(*bound)->Nk;k++){
	(*bound)->z[k]=0.0;
    }

    if((*bound)->hasType2==1){
	for(j=0;j<Ntype2;j++){
	    (*bound)->edgep[j]=0;
	    (*bound)->localedgep[j]=-1;
	    (*bound)->xe[j]=0.0;
	    (*bound)->ye[j]=0.0;
	    for(k=0;k<Nk;k++){
		(*bound)->boundary_u[k][j]=0.0;
		(*bound)->boundary_v[k][j]=0.0;
		(*bound)->boundary_w[k][j]=0.0;
		(*bound)->boundary_T[k][j]=0.0;
		(*bound)->boundary_S[k][j]=0.0;
		for(n=0;n>NT;n++){
		    (*bound)->boundary_u_t[n][k][j]=0.0;
		    (*bound)->boundary_v_t[n][k][j]=0.0;
		    (*bound)->boundary_w_t[n][k][j]=0.0;
		    (*bound)->boundary_T_t[n][k][j]=0.0;
		    (*bound)->boundary_S_t[n][k][j]=0.0;
		}
	    }
	}
	for(i=0;i<n2;i++){
	    (*bound)->ind2[i]=-1;
	}
    }

    if((*bound)->hasSeg==1){
	for(j=0;j<Ntype2;j++){
	    (*bound)->segedgep[j]=0;
	}
	for(j=0;j<Nseg;j++){
	    (*bound)->segp[j]=0;
	    (*bound)->segarea[j]=0.0;
	    (*bound)->localsegarea[j]=0.0;
	    (*bound)->boundary_Q[j]=0.0;
	    for(n=0;n<NT;n++){
		(*bound)->boundary_Q_t[n][j]=0.0;
	    }
	}
    }
    if(myproc==0) printf("Finished Zeroing Type2 boundary arrays...\n");

    if((*bound)->hasType3==1){
	for(j=0;j<Ntype3;j++){
	    (*bound)->cellp[j]=0;
	    (*bound)->xv[j]=0.0;
	    (*bound)->yv[j]=0.0;

	    (*bound)->h[j]=0.0;
	    for(n=0;n<NT;n++){
		 (*bound)->h_t[n][j]=0.0;
	    }
	    for(k=0;k<Nk;k++){
		(*bound)->uc[k][j]=0.0;
		(*bound)->vc[k][j]=0.0;
		(*bound)->wc[k][j]=0.0;
		(*bound)->T[k][j]=0.0;
		(*bound)->S[k][j]=0.0;
		for(n=0;n<NT;n++){
		    (*bound)->uc_t[n][k][j]=0.0;
		    (*bound)->vc_t[n][k][j]=0.0;
		    (*bound)->wc_t[n][k][j]=0.0;
		    (*bound)->T_t[n][k][j]=0.0;
		    (*bound)->S_t[n][k][j]=0.0;
		}
	    }
	}
	for(i=0;i<n3;i++){
	    (*bound)->ind3[i]=-1;
	    (*bound)->ind3edge[i]=-1;
	}
    }
    if(myproc==0) printf("Finished Zeroing Type 3 boundary arrays...\n");
 }//End function


/*
* Function: GetTimeRecBnd()
* ------------------
* Retuns the index of the first preceding time step in the vector time
*/
int getTimeRecBnd(REAL nctime, REAL *time, int nt){
    int j;

    for(j=0;j<nt;j++){
       if (time[j]>=nctime)
	 return j-1;
    }
    return nt;
}

