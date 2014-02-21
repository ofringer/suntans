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
void InitializeMerging(gridT *grid, int numprocs, int myproc, MPI_Comm comm) {
  int i, iptr, p, Nc_computational, size3D;
  int *mnptr_temp, *Nk_temp;
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
      MPI_Recv(&(Nc_all[p]),1,MPI_INT,p,1,comm,&status);         
      MPI_Recv(&(send3DSize[p]),1,MPI_INT,p,1,comm,&status);         

      if(Nc_all[p]>Nc_max)
	Nc_max=Nc_all[p];

      mergedGrid->Nc+=Nc_all[p];
    }
  }
  MPI_Bcast(&(Nc_max),1,MPI_INT,0,comm);
  MPI_Bcast(send3DSize,numprocs,MPI_INT,0,comm);

  mnptr_temp=(int *)SunMalloc(Nc_max*sizeof(int),"InitializeMerging");
  Nk_temp=(int *)SunMalloc(Nc_max*sizeof(int),"InitializeMerging");

  if(myproc==0) {
    mnptr_all=(int **)SunMalloc(numprocs*sizeof(int *),"InitializeMerging");
    for(p=0;p<numprocs;p++) {
      mnptr_all[p]=(int *)SunMalloc(Nc_all[p]*sizeof(int),"InitializeMerging");
      for(i=0;i<Nc_all[p];i++)
	mnptr_all[p][i]=0;
    }
    mergedGrid->Nk=(int *)SunMalloc(mergedGrid->Nc*sizeof(int),"InitializeMerging");
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
      for(i=0;i<Nc_all[p];i++)
	mergedGrid->Nk[mnptr_all[p][i]]=Nk_temp[i];
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
  }

  SunFree(mnptr_temp,Nc_max*sizeof(int),"InitializeMerging");
  SunFree(Nk_temp,Nc_max*sizeof(int),"InitializeMerging");
}

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
