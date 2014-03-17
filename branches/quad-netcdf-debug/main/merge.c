/*
 * File: merge.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Functions for merging data onto one processor for writing.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "merge.h"
#include "memory.h"

/*
 * Private Functions
 */
static void MergeGridVariables(gridT *grid, int numprocs, int myproc, MPI_Comm comm);
static void InitializeMergeEdges(gridT *grid, int numprocs, int myproc, MPI_Comm comm);

/*
 * Function: InitializeMerging
 * Usage: InitializeMerging(grid,numprocs,myproc,comm);
 * ----------------------------------------------------
 * Allocate space needed for merging of data in global arrays defined in merge.h:
 *
 *  mergedGrid: Contains Nk which stores the number of vertical layers in each cell on the merged grid.
 *  **localTempMergeArray: 2D array of size Nc_max X Nkmax, where Nc_max is the maximum number of computational 
 *    cells in 2d, i.e. cells of type 0,1,2,or 3, on any of the processors.
 *  *merged2DArray: 2D array of size mergedGrid->Nc, where mergedGrid->Nc contains the total number
 *    of computational cells on the merged grid.
 *  **merged3DArray: 3D array of size mergedGrid->Nc X Nkmax, where mergedGrid->Nc contains the total number
 *    of computational cells on the merged grid.
 *  **mnptr_all: Contains the pointers mapping local grid indices to the merged grid indices.
 *  *Nc_all: Contains the number of computational cells on each processor.
 *  *send3DSize: Size of 3D array to send,i.e. sum(Nk[i=0:Nc_computational]).
 *  Nc_max: Maximum of Nc_all.
 *
 */
void InitializeMerging(gridT *grid, int mergeedges, int numprocs, int myproc, MPI_Comm comm) {
  int i, iptr, p, Nc_computational, Ne_computational, size3D;
  int *mnptr_temp, *Nk_temp, *eptr_temp;
  MPI_Status status;

  size3D = 0;
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i=grid->cellp[iptr];
    
    size3D+=grid->Nk[i];
  }
  send3DSize=(int *)SunMalloc(numprocs*sizeof(int),"InitializeMerging");
  
  if(myproc!=0) {
    Nc_computational = grid->celldist[2]-grid->celldist[0];
    MPI_Send(&(Nc_computational),1,MPI_INT,0,1,comm);     
    MPI_Send(&(size3D),1,MPI_INT,0,1,comm);     
  } else {
    mergedGrid=(gridT *)SunMalloc(sizeof(gridT),"InitializeMerging");

    Nc_all=(int *)SunMalloc(numprocs*sizeof(int),"InitializeMerging");
    Nc_all[0]=grid->celldist[2]-grid->celldist[0];
    
    Nc_max=Nc_all[0];
    mergedGrid->Nc=Nc_all[0];
    send3DSize[0]=size3D;

    for(p=1;p<numprocs;p++) {
      //Cells
      MPI_Recv(&(Nc_all[p]),1,MPI_INT,p,1,comm,&status);         
      MPI_Recv(&(send3DSize[p]),1,MPI_INT,p,1,comm,&status);         

      if(Nc_all[p]>Nc_max)
	Nc_max=Nc_all[p];

      mergedGrid->Nc+=Nc_all[p];

      // These variables are constant across processors
      mergedGrid->maxfaces=grid->maxfaces;
      mergedGrid->Np=grid->Np;
    }
  }
  MPI_Bcast(&(Nc_max),1,MPI_INT,0,comm);
  MPI_Bcast(send3DSize,numprocs,MPI_INT,0,comm);

  mnptr_temp=(int *)SunMalloc(Nc_max*sizeof(int),"InitializeMerging");

  //Allocate temp variables for merging grid variables
  Nk_temp=(int *)SunMalloc(Nc_max*sizeof(int),"InitializeMerging");

  
  if(myproc==0) {
    mnptr_all=(int **)SunMalloc(numprocs*sizeof(int *),"InitializeMerging");
    for(p=0;p<numprocs;p++) {
      mnptr_all[p]=(int *)SunMalloc(Nc_all[p]*sizeof(int),"InitializeMerging");
      for(i=0;i<Nc_all[p];i++)
	mnptr_all[p][i]=0;

    }

    mergedGrid->Nk=(int *)SunMalloc(mergedGrid->Nc*sizeof(int),"MergeGridVariables");
  }

  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];

      mnptr_temp[iptr-grid->celldist[0]]=grid->mnptr[i];
      Nk_temp[iptr-grid->celldist[0]]=grid->Nk[i];
    }
    MPI_Send(mnptr_temp,grid->celldist[2]-grid->celldist[0],MPI_INT,0,1,comm);       
    MPI_Send(Nk_temp,grid->celldist[2]-grid->celldist[0],MPI_INT,0,1,comm);       
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];

      mnptr_all[0][iptr-grid->celldist[0]]=grid->mnptr[i];
      mergedGrid->Nk[grid->mnptr[i]]=grid->Nk[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(mnptr_all[p],Nc_all[p],MPI_INT,p,1,comm,&status);         
      MPI_Recv(Nk_temp,Nc_all[p],MPI_INT,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++){
	mergedGrid->Nk[mnptr_all[p][i]]=Nk_temp[i];
      }
    }
  }

  // Allocate space for temporary arrays used in merging
  // All processors need a temporary array used to send data
  localTempMergeArray=(REAL *)SunMalloc(Nc_max*grid->Nkmax*sizeof(REAL),"InitializeMerging");

  // Only processor 0 needs the temporary array storing the entire grid
  if(myproc==0) {
    merged2DArray=(REAL *)SunMalloc(mergedGrid->Nc*sizeof(REAL),"InitializeMerging");
    merged3DArray=(REAL **)SunMalloc(mergedGrid->Nc*sizeof(REAL *),"InitializeMerging");
    for(i=0;i<mergedGrid->Nc;i++)
      merged3DArray[i]=(REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializeMerging");
    //This is necessary for netcdf write
    //Allocate extra vertical layer to allow for w  
    //merged3DVector=(REAL *)SunMalloc(mergedGrid->Nc*(grid->Nkmax+1)*sizeof(REAL),"InitializeMerging");

  }


  //Merge the edges and grid if output netcdf only
  // Now do the edge arrays
  if (mergeedges>0){
      InitializeMergeEdges(grid,numprocs,myproc,comm);
      // Merge all of the grid variables onto one processor
      MergeGridVariables(grid, numprocs, myproc, comm);
  }else if(myproc==0) {
      merged3DVector=(REAL *)SunMalloc(mergedGrid->Nc*(grid->Nkmax+1)*sizeof(REAL),"InitializeMerging");
  }

  SunFree(mnptr_temp,Nc_max*sizeof(int),"InitializeMerging");
  SunFree(Nk_temp,Nc_max*sizeof(int),"InitializeMerging");

}

static void InitializeMergeEdges(gridT *grid, int numprocs, int myproc, MPI_Comm comm){
  int j, jptr, p, Ne_computational, size3D;
  int *eptr_temp, *Nke_temp;
  MPI_Status status;

  size3D = 0;
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[EDGEMAX];jptr++) {
    j=grid->edgep[jptr];
    
    size3D+=grid->Nke[j];
  }
  send3DESize=(int *)SunMalloc(numprocs*sizeof(int),"InitializeMerging");
 
  if(myproc!=0) {
    Ne_computational = grid->edgedist[EDGEMAX]-grid->edgedist[0];
    MPI_Send(&(Ne_computational),1,MPI_INT,0,1,comm);     
    MPI_Send(&(size3D),1,MPI_INT,0,1,comm);     
  } else {

    Ne_all=(int *)SunMalloc(numprocs*sizeof(int),"InitializeMerging");
    Ne_all[0]=grid->edgedist[EDGEMAX]-grid->edgedist[0];


    Ne_max=Ne_all[0];
    mergedGrid->Ne=Ne_all[0];
    send3DESize[0]=size3D;
    
    for(p=1;p<numprocs;p++) {
      //Edges
      MPI_Recv(&(Ne_all[p]),1,MPI_INT,p,1,comm,&status);         
      MPI_Recv(&(send3DESize[p]),1,MPI_INT,p,1,comm,&status);         

      if(Ne_all[p]>Ne_max)
        Ne_max=Ne_all[p];

      mergedGrid->Ne+=Ne_all[p];

    }
  }
  MPI_Bcast(&(Ne_max),1,MPI_INT,0,comm);
  MPI_Bcast(send3DESize,numprocs,MPI_INT,0,comm);
  eptr_temp=(int *)SunMalloc(Ne_max*sizeof(int),"InitializeMerging");
//  Nke_temp=(int *)SunMalloc(Ne_max*sizeof(int),"InitializeMerging");

  if(myproc==0) {
    eptr_all=(int **)SunMalloc(numprocs*sizeof(int *),"InitializeMerging");
    for(p=0;p<numprocs;p++) {

      eptr_all[p]=(int *)SunMalloc(Ne_all[p]*sizeof(int),"InitializeMerging");
      for(j=0;j<Ne_all[p];j++)
	eptr_all[p][j]=0;

    }
//    mergedGrid->Nke=(int *)SunMalloc(mergedGrid->Ne*sizeof(int),"MergeGridVariables");
  }

  // Merge the eptr and Nke arrays
  if(myproc!=0) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[EDGEMAX];jptr++) {
      j=grid->edgep[jptr];

      eptr_temp[jptr-grid->edgedist[0]]=grid->eptr[j];
//     Nke_temp[jptr-grid->edgedist[0]]=grid->Nke[j];
    }
    MPI_Send(eptr_temp,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_INT,0,1,comm);       
//    MPI_Send(Nke_temp,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_INT,0,1,comm);       
  } else {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[EDGEMAX];jptr++) {
      j=grid->edgep[jptr];

      eptr_all[0][jptr-grid->edgedist[0]]=grid->eptr[j];
//      mergedGrid->Nke[grid->eptr[j]]=grid->Nke[j];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(eptr_all[p],Ne_all[p],MPI_INT,p,1,comm,&status);         
//      MPI_Recv(Nke_temp,Ne_all[p],MPI_INT,p,1,comm,&status);         
//      for(j=0;j<Ne_all[p];j++)
//	mergedGrid->Nke[eptr_all[p][j]]=Nke_temp[j];
    }
  }
//// Find the maximum value in mnptr for testing
  //if(myproc==0){
  //    mergedGrid->Ne=0;
  //    // processor 0
  //    for(j=0;j<grid->Nc;j++){
  //       if(grid->mnptr[j]>mergedGrid->Ne)
  //           mergedGrid->Ne = grid->mnptr[j];
  //    }
  //    // other processors
  //    for(p=1;p<numprocs;p++){
  //       for(j=0;j<Nc_all[p];j++){
  //          if(mnptr_all[p][j]>>mergedGrid->Ne)
  //          	mergedGrid->Ne = mnptr_all[p][j]+1;
  //       }
  //    }
  //}
  // Now go back through and get the maximum value in eptr_all
  if(myproc==0){
      mergedGrid->Ne=0;
      // processor 0
      for(j=0;j<grid->Ne;j++){
         if(grid->eptr[j]>=mergedGrid->Ne)
             mergedGrid->Ne = grid->eptr[j]+1;
      }
      // other processors
      for(p=1;p<numprocs;p++){
         for(j=0;j<Ne_all[p];j++){
            if(eptr_all[p][j]>=mergedGrid->Ne)
            	mergedGrid->Ne = eptr_all[p][j]+1;
         }
      }
  }
  localTempEMergeArray=(REAL *)SunMalloc(Ne_max*grid->Nkmax*sizeof(REAL),"InitializeMerging");
  // Only processor 0 needs the temporary array storing the entire grid
  if(myproc==0) {
    merged3DEArray=(REAL **)SunMalloc(mergedGrid->Ne*sizeof(REAL *),"InitializeMerging");
    for(j=0;j<mergedGrid->Ne;j++)
      merged3DEArray[j]=(REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializeMerging");
    //This is necessary for netcdf write
    merged3DVector=(REAL *)SunMalloc(mergedGrid->Ne*(grid->Nkmax)*sizeof(REAL),"InitializeMerging");

  }

  SunFree(eptr_temp,Ne_max*sizeof(int),"InitializeMerging");
  //SunFree(Nke_temp,Ne_max*sizeof(int),"InitializeMerging");


}// End function

void MergeCellCentered2DArray(REAL *localArray, gridT *grid, int numprocs, int myproc, MPI_Comm comm) {
  int i, p, iptr;
  MPI_Status status;

  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];

      localTempMergeArray[iptr-grid->celldist[0]]=localArray[i];
    }
    MPI_Send(&(localTempMergeArray[0]),grid->celldist[2]-grid->celldist[0],MPI_DOUBLE,0,1,comm);       
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];

      merged2DArray[grid->mnptr[i]]=localArray[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(&(localTempMergeArray[0]),Nc_all[p],MPI_DOUBLE,p,1,comm,&status);         

      for(i=0;i<Nc_all[p];i++) 
	merged2DArray[mnptr_all[p][i]]=localTempMergeArray[i];
    }
  }
}

void MergeCellCentered3DArray(REAL **localArray, gridT *grid, int numprocs, int myproc, MPI_Comm comm) {
  int i, k, m, p, iptr;
  MPI_Status status;

  if(myproc!=0) {
    m=0;
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];

      
      for(k=0;k<grid->Nk[i];k++)
      	localTempMergeArray[m++]=localArray[i][k];
    }
    MPI_Send(&(localTempMergeArray[0]),send3DSize[myproc],MPI_DOUBLE,0,1,comm);       
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];

      for(k=0;k<grid->Nk[i];k++) 
	merged3DArray[grid->mnptr[i]][k]=localArray[i][k];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(&(localTempMergeArray[0]),send3DSize[p],MPI_DOUBLE,p,1,comm,&status);         

      m=0;
      for(i=0;i<Nc_all[p];i++) {
	for(k=0;k<mergedGrid->Nk[mnptr_all[p][i]];k++)
	  merged3DArray[mnptr_all[p][i]][k]=localTempMergeArray[m++];
      }
    }
  }
}

void MergeEdgeCentered3DArray(REAL **localArray, gridT *grid, int numprocs, int myproc, MPI_Comm comm) {
  int j, k, m, p, jptr;
  MPI_Status status;

  if(myproc!=0) {
    m=0;
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[EDGEMAX];jptr++) {
      j=grid->edgep[jptr];
      for(k=0;k<grid->Nke[j];k++)
      	localTempEMergeArray[m++]=localArray[j][k];
    }
    MPI_Send(&(localTempEMergeArray[0]),send3DESize[myproc],MPI_DOUBLE,0,1,comm);       
  } else {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[EDGEMAX];jptr++) {
      j=grid->edgep[jptr];

      for(k=0;k<grid->Nke[j];k++) 
	//merged3DEArray[eptr_all[0][j]][k]=localArray[j][k];
	merged3DEArray[grid->eptr[j]][k]=localArray[j][k];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(&(localTempEMergeArray[0]),send3DESize[p],MPI_DOUBLE,p,1,comm,&status);      
      m=0;
      for(j=0;j<Ne_all[p];j++) {
	for(k=0;k<mergedGrid->Nke[eptr_all[p][j]];k++)
	  merged3DEArray[eptr_all[p][j]][k]=localTempEMergeArray[m++];
      }
    }
  }
}
/*
 *MergeGridVariables()()
 *----------------------------
 * Merge all of the grid variables into the mergedGrid structure.
 * 
 */
static void MergeGridVariables(gridT *grid, int numprocs, int myproc, MPI_Comm comm){

  int iptr, i, p, nf, nff, jptr, j;
  int *Ne_int, *Nc_int, *NcNs_int, *Ne2_int;
  REAL *Nc_real, *Ne_real, NcNs_real;

  MPI_Status status;

  // Allocate space for the temp arrays on each processor
  // Some of these variables are used more than once
  Nc_int=(int *)SunMalloc(Nc_max*sizeof(int),"MergeGridVariables");
  Ne_int=(int *)SunMalloc(Ne_max*sizeof(int),"MergeGridVariables");
  NcNs_int=(int *)SunMalloc(grid->maxfaces*Nc_max*sizeof(int),"MergeGridVariables");
  Ne2_int=(int *)SunMalloc(2*Ne_max*sizeof(int),"MergeGridVariables");

  Nc_real=(REAL *)SunMalloc(Nc_max*sizeof(REAL),"MergeGridVariables");
  Ne_real=(REAL *)SunMalloc(Ne_max*sizeof(REAL),"MergeGridVariables");

  if(myproc==0){
    // Allocate arrays
    mergedGrid->Nkmax = grid->Nkmax;
    mergedGrid->maxfaces = grid->maxfaces;

    // These are all needed to write to netcdf
    mergedGrid->nfaces = (int *)SunMalloc(mergedGrid->Nc*sizeof(REAL),"MergeGridVariables");	
    mergedGrid->xv = (REAL *)SunMalloc(mergedGrid->Nc*sizeof(REAL),"MergeGridVariables");
    mergedGrid->yv = (REAL *)SunMalloc(mergedGrid->Nc*sizeof(REAL),"MergeGridVariables");
    mergedGrid->dv = (REAL *)SunMalloc(mergedGrid->Nc*sizeof(REAL),"MergeGridVariables");
    mergedGrid->Ac = (REAL *)SunMalloc(mergedGrid->Nc*sizeof(REAL),"MergeGridVariables");
  
    mergedGrid->neigh = (int *)SunMalloc(mergedGrid->maxfaces*mergedGrid->Nc*sizeof(int),"MergeGridVariables");
    mergedGrid->face = (int *)SunMalloc(mergedGrid->maxfaces*mergedGrid->Nc*sizeof(int),"MergeGridVariables");
    mergedGrid->normal = (int *)SunMalloc(mergedGrid->maxfaces*mergedGrid->Nc*sizeof(int),"MergeGridVariables");
    mergedGrid->def = (REAL *)SunMalloc(mergedGrid->maxfaces*mergedGrid->Nc*sizeof(REAL),"MergeGridVariables");
    mergedGrid->cells = (int *)SunMalloc(mergedGrid->maxfaces*mergedGrid->Nc*sizeof(REAL),"MergeGridVariables");

    mergedGrid->df = (REAL *)SunMalloc(mergedGrid->Ne*sizeof(REAL),"MergeGridVariables");
    mergedGrid->dg = (REAL *)SunMalloc(mergedGrid->Ne*sizeof(REAL),"MergeGridVariables");
    mergedGrid->n1 = (REAL *)SunMalloc(mergedGrid->Ne*sizeof(REAL),"MergeGridVariables");
    mergedGrid->n2 = (REAL *)SunMalloc(mergedGrid->Ne*sizeof(REAL),"MergeGridVariables");
    mergedGrid->xe = (REAL *)SunMalloc(mergedGrid->Ne*sizeof(REAL),"MergeGridVariables");
    mergedGrid->ye = (REAL *)SunMalloc(mergedGrid->Ne*sizeof(REAL),"MergeGridVariables");
  
    mergedGrid->grad = (int *)SunMalloc(2*mergedGrid->Ne*sizeof(int),"MergeGridVariables");
    mergedGrid->gradf = (int *)SunMalloc(2*mergedGrid->Ne*sizeof(int),"MergeGridVariables");
    mergedGrid->edges = (int *)SunMalloc(2*mergedGrid->Ne*sizeof(int),"MergeGridVariables");
    mergedGrid->mark = (int *)SunMalloc(mergedGrid->Ne*sizeof(int),"MergeGridVariables");
    mergedGrid->Nke=(int *)SunMalloc(mergedGrid->Ne*sizeof(int),"MergeGridVariables");
  }
  
  // Now go through and merge each variable individually. This cannot
  // really go into a function as the mergedGrid array is only on processor 0
  // and cannot be passed to functions.

  /* Integer variables size: Nc */
  //Nk is already written

  //nfaces
  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      Nc_int[iptr-grid->celldist[0]]=grid->nfaces[i];
    }
    MPI_Send(Nc_int,grid->celldist[2]-grid->celldist[0],MPI_INT,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      mergedGrid->nfaces[grid->mnptr[i]]=grid->nfaces[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Nc_int,Nc_all[p],MPI_INT,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++)
	mergedGrid->nfaces[mnptr_all[p][i]]=Nc_int[i];
    }
  }

  /*Real variables size: Nc */
  //xv
  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      Nc_real[iptr-grid->celldist[0]]=grid->xv[i];
    }
    MPI_Send(Nc_real,grid->celldist[2]-grid->celldist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      mergedGrid->xv[grid->mnptr[i]]=grid->xv[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Nc_real,Nc_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++)
	mergedGrid->xv[mnptr_all[p][i]]=Nc_real[i];
    }
  }

  //yv
  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      Nc_real[iptr-grid->celldist[0]]=grid->yv[i];
    }
    MPI_Send(Nc_real,grid->celldist[2]-grid->celldist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      mergedGrid->yv[grid->mnptr[i]]=grid->yv[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Nc_real,Nc_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++)
	mergedGrid->yv[mnptr_all[p][i]]=Nc_real[i];
    }
  }

  //Ac
   if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      Nc_real[iptr-grid->celldist[0]]=grid->Ac[i];
    }
    MPI_Send(Nc_real,grid->celldist[2]-grid->celldist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      mergedGrid->Ac[grid->mnptr[i]]=grid->Ac[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Nc_real,Nc_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++)
	mergedGrid->Ac[mnptr_all[p][i]]=Nc_real[i];
    }
  }  

  //dv
   if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      Nc_real[iptr-grid->celldist[0]]=grid->dv[i];
    }
    MPI_Send(Nc_real,grid->celldist[2]-grid->celldist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      mergedGrid->dv[grid->mnptr[i]]=grid->dv[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Nc_real,Nc_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++)
	mergedGrid->dv[mnptr_all[p][i]]=Nc_real[i];
    }
  }

  /* Integer variables size: [Nc, maxfaces]*/

  //cells - node_connectivity
  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      for(nf=0;nf<grid->nfaces[i];nf++){
	  NcNs_int[(iptr-grid->celldist[0])*grid->maxfaces+nf]=grid->cells[i*grid->maxfaces+nf];
      }
    }
    MPI_Send(NcNs_int,(grid->celldist[2]-grid->celldist[0])*grid->maxfaces,MPI_INT,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      for(nf=0;nf<grid->nfaces[i];nf++){
	  //mergedGrid->cells[mnptr_all[0][i]*grid->maxfaces+nf]=grid->cells[i*grid->maxfaces+nf];
	  mergedGrid->cells[grid->mnptr[i]*grid->maxfaces+nf]=grid->cells[i*grid->maxfaces+nf];
      }
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(NcNs_int,Nc_all[p]*grid->maxfaces,MPI_INT,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++){
         nff = mergedGrid->nfaces[mnptr_all[p][i]];
	 for(nf=0;nf<nff;nf++){
	    mergedGrid->cells[mnptr_all[p][i]*grid->maxfaces+nf]=NcNs_int[i*grid->maxfaces+nf];
	    //mergedGrid->cells[mnptr_all[p][i]*grid->maxfaces+nf]=nf;
	 }
      }
    }
  }

  // face - face_edge_connectivity
  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      for(nf=0;nf<grid->nfaces[i];nf++){
	  //NcNs_int[(iptr-grid->celldist[0])*grid->maxfaces+nf]=grid->face[i*grid->maxfaces+nf];
	  NcNs_int[(iptr-grid->celldist[0])*grid->maxfaces+nf]=grid->eptr[grid->face[i*grid->maxfaces+nf]];
      }
    }
    MPI_Send(NcNs_int,(grid->celldist[2]-grid->celldist[0])*grid->maxfaces,MPI_INT,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      for(nf=0;nf<grid->nfaces[i];nf++){
	  mergedGrid->face[grid->mnptr[i]*grid->maxfaces+nf]=grid->eptr[grid->face[i*grid->maxfaces+nf]];
      }
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(NcNs_int,Nc_all[p]*grid->maxfaces,MPI_INT,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++){
         nff = mergedGrid->nfaces[mnptr_all[p][i]];
	 for(nf=0;nf<nff;nf++){
	    mergedGrid->face[mnptr_all[p][i]*grid->maxfaces+nf]=NcNs_int[i*grid->maxfaces+nf];
	 }
      }
    }
  }

  // normal - outward normal for each edge
  if(myproc!=0) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      for(nf=0;nf<grid->nfaces[i];nf++){
	  NcNs_int[(iptr-grid->celldist[0])*grid->maxfaces+nf]=grid->normal[i*grid->maxfaces+nf];
      }
    }
    MPI_Send(NcNs_int,(grid->celldist[2]-grid->celldist[0])*grid->maxfaces,MPI_INT,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      for(nf=0;nf<grid->nfaces[i];nf++){
	  mergedGrid->normal[grid->mnptr[i]*grid->maxfaces+nf]=grid->normal[i*grid->maxfaces+nf];
      }
    }
    for(p=1;p<numprocs;p++) {
      MPI_Recv(NcNs_int,Nc_all[p]*grid->maxfaces,MPI_INT,p,1,comm,&status);         
      for(i=0;i<Nc_all[p];i++){
         nff = mergedGrid->nfaces[mnptr_all[p][i]];
	 for(nf=0;nf<nff;nf++){
	    mergedGrid->normal[mnptr_all[p][i]*grid->maxfaces+nf]=NcNs_int[i*grid->maxfaces+nf];
	 }
      }
    }
  }

  /* Integer variables size: [Ne, 2]*/
  
  //grad - edge_face_connectivity
  if(myproc!=0) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[EDGEMAX];jptr++) {
      j=grid->edgep[jptr];
      for(nf=0;nf<2;nf++){
	  Ne2_int[(jptr-grid->edgedist[0])*2+nf]=grid->grad[j*2+nf];
      }
    }
    MPI_Send(Ne2_int,(grid->edgedist[EDGEMAX]-grid->edgedist[0])*2,MPI_INT,0,1,comm);       
  } else {//myproc==0
    	
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[EDGEMAX];jptr++) {
      j=grid->edgep[jptr];
      for(nf=0;nf<2;nf++){
	  //mergedGrid->grad[grid->eptr[i]*2+nf]=grid->mnptr[grid->grad[i*2+nf]];
	  if(grid->grad[j*2+nf]==-1){
	      mergedGrid->grad[grid->eptr[j]*2+nf]=-1;
	  }else{
	      mergedGrid->grad[grid->eptr[j]*2+nf]=grid->mnptr[grid->grad[j*2+nf]];
	  }
      }
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne2_int,Ne_all[p]*2,MPI_INT,p,1,comm,&status);         
      for(j=0;j<Ne_all[p];j++){
	 for(nf=0;nf<2;nf++){
	    //mergedGrid->grad[eptr_all[p][i]*2+nf]=mnptr_all[p][Ne2_int[i*2+nf]];
	    if(Ne2_int[j*2+nf]==-1){
	      mergedGrid->grad[eptr_all[p][j]*2+nf]=-1;
	    }else{
	      i = Ne2_int[j*2+nf];
	      // Type-5 marked edges will cause the index to overflow the mnptr array
	      if(i<Nc_all[p]) 
		  mergedGrid->grad[eptr_all[p][j]*2+nf]=mnptr_all[p][Ne2_int[j*2+nf]];
	    }
	 }
      }
    }
  }

  //edges - edge_node_connectivity
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      for(nf=0;nf<2;nf++){
	  Ne2_int[(iptr-grid->edgedist[0])*2+nf]=grid->edges[i*NUMEDGECOLUMNS+nf];
      }
    }
    MPI_Send(Ne2_int,(grid->edgedist[EDGEMAX]-grid->edgedist[0])*2,MPI_INT,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      for(nf=0;nf<2;nf++){
	  mergedGrid->edges[grid->eptr[i]*2+nf]=grid->edges[i*NUMEDGECOLUMNS+nf];
      }
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne2_int,Ne_all[p]*2,MPI_INT,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++){
	 for(nf=0;nf<2;nf++){
	    mergedGrid->edges[eptr_all[p][i]*2+nf]=Ne2_int[i*2+nf];
	 }
      }
    }
  }


  /*Real variables size: Ne */
  //xe
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_real[iptr-grid->edgedist[0]]=grid->xe[i];
    }
    MPI_Send(Ne_real,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->xe[grid->eptr[i]]=grid->xe[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_real,Ne_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->xe[eptr_all[p][i]]=Ne_real[i];
    }
  }

  //ye
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_real[iptr-grid->edgedist[0]]=grid->ye[i];
    }
    MPI_Send(Ne_real,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->ye[grid->eptr[i]]=grid->ye[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_real,Ne_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->ye[eptr_all[p][i]]=Ne_real[i];
    }
  }

  //df
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_real[iptr-grid->edgedist[0]]=grid->df[i];
    }
    MPI_Send(Ne_real,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->df[grid->eptr[i]]=grid->df[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_real,Ne_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->df[eptr_all[p][i]]=Ne_real[i];
    }
  }

  //dg
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_real[iptr-grid->edgedist[0]]=grid->dg[i];
    }
    MPI_Send(Ne_real,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->dg[grid->eptr[i]]=grid->dg[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_real,Ne_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->dg[eptr_all[p][i]]=Ne_real[i];
    }
  }

  //n1
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_real[iptr-grid->edgedist[0]]=grid->n1[i];
    }
    MPI_Send(Ne_real,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->n1[grid->eptr[i]]=grid->n1[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_real,Ne_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->n1[eptr_all[p][i]]=Ne_real[i];
    }
  }

  //n2
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_real[iptr-grid->edgedist[0]]=grid->n2[i];
    }
    MPI_Send(Ne_real,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_DOUBLE,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->n2[grid->eptr[i]]=grid->n2[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_real,Ne_all[p],MPI_DOUBLE,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->n2[eptr_all[p][i]]=Ne_real[i];
    }
  }
  /* Integer variables - size Ne */
  //Nke
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_int[iptr-grid->edgedist[0]]=grid->Nke[i];
    }
    MPI_Send(Ne_int,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_INT,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->Nke[grid->eptr[i]]=grid->Nke[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_int,Ne_all[p],MPI_INT,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->Nke[eptr_all[p][i]]=Ne_int[i];
    }
  }


  //mark
  if(myproc!=0) {
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      Ne_int[iptr-grid->edgedist[0]]=grid->mark[i];
    }
    MPI_Send(Ne_int,grid->edgedist[EDGEMAX]-grid->edgedist[0],MPI_INT,0,1,comm);       
  } else {//myproc==0
    for(iptr=grid->edgedist[0];iptr<grid->edgedist[EDGEMAX];iptr++) {
      i=grid->edgep[iptr];
      mergedGrid->mark[grid->eptr[i]]=grid->mark[i];
    }

    for(p=1;p<numprocs;p++) {
      MPI_Recv(Ne_int,Ne_all[p],MPI_INT,p,1,comm,&status);         
      for(i=0;i<Ne_all[p];i++)
	mergedGrid->mark[eptr_all[p][i]]=Ne_int[i];
    }
  }

  // Free up the swap arrays
  SunFree(Nc_int,Nc_max*sizeof(int),"MergeGridVariables");
  SunFree(Ne_int,Ne_max*sizeof(int),"MergeGridVariables");
  SunFree(NcNs_int,grid->maxfaces*Nc_max*sizeof(int),"MergeGridVariables");
  SunFree(Ne2_int,2*Ne_max*sizeof(int),"MergeGridVariables");
  SunFree(Nc_real,Nc_max*sizeof(REAL),"MergeGridVariables");
  SunFree(Ne_real,Ne_max*sizeof(REAL),"MergeGridVariables");

}

/*
 * Function: FreeMergingArrays
 * Usage: FreeMergingArrays(grid,myproc);
 * --------------------------------------
 * Free space associated with array merging.
 *
 */
void FreeMergingArrays(gridT *grid, int myproc) {
  int i;
  
  // Space was allocated for temporary arrays used in merging
  // All processors needed a temporary array used to send data
  SunFree(localTempMergeArray,Nc_max*grid->Nkmax*sizeof(REAL *),"FreeMergingArrays");  

  // Only processor 0 needed the temporary array storing the entire grid
  if(myproc==0) {
    for(i=0;i<grid->Nkmax;i++)
      SunFree(merged3DArray[i],grid->Nkmax*sizeof(REAL),"FreeMergingArrays");
    SunFree(merged3DArray,mergedGrid->Nc*sizeof(REAL *),"FreeMergingArrays");
    SunFree(merged2DArray,mergedGrid->Nc*sizeof(REAL),"FreeMergingArrays");

    SunFree(mergedGrid,sizeof(gridT *),"FreeMergingArrays");
  }
}
