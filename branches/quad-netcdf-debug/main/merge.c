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
    //Ne_computational = grid->edgedist[1]-grid->edgedist[0];
    MPI_Send(&(Nc_computational),1,MPI_INT,0,1,comm);     
    //MPI_Send(&(Ne_computational),1,MPI_INT,0,1,comm);     
    MPI_Send(&(size3D),1,MPI_INT,0,1,comm);     
  } else {
    mergedGrid=(gridT *)SunMalloc(sizeof(gridT),"InitializeMerging");

    Nc_all=(int *)SunMalloc(numprocs*sizeof(int),"InitializeMerging");
    Nc_all[0]=grid->celldist[2]-grid->celldist[0];
    
    //Ne_all=(int *)SunMalloc(numprocs*sizeof(int),"InitializeMerging");
    //Ne_all[0]=grid->edgedist[1]-grid->edgedist[0];

    Nc_max=Nc_all[0];
    mergedGrid->Nc=Nc_all[0];
    send3DSize[0]=size3D;

    //Ne_max=Ne_all[0];
    //mergedGrid->Ne=Ne_all[0];
    mergedGrid->Ne=0;
    
    for(p=1;p<numprocs;p++) {
      //Cells
      MPI_Recv(&(Nc_all[p]),1,MPI_INT,p,1,comm,&status);         
      MPI_Recv(&(send3DSize[p]),1,MPI_INT,p,1,comm,&status);         

      if(Nc_all[p]>Nc_max)
	Nc_max=Nc_all[p];

      mergedGrid->Nc+=Nc_all[p];

      //Edges
      //MPI_Recv(&(Ne_all[p]),1,MPI_INT,p,1,comm,&status);         

      //if(Ne_all[p]>Ne_max)
      //  Ne_max=Ne_all[p];

      //mergedGrid->Ne+=Ne_all[p];

      // These variables are constant across processors
      mergedGrid->maxfaces=grid->maxfaces;
      mergedGrid->Np=grid->Np;
    }
  }
  MPI_Bcast(&(Nc_max),1,MPI_INT,0,comm);
  MPI_Bcast(send3DSize,numprocs,MPI_INT,0,comm);

  //MPI_Bcast(&(Ne_max),1,MPI_INT,0,comm);

  mnptr_temp=(int *)SunMalloc(Nc_max*sizeof(int),"InitializeMerging");

  //eptr_temp=(int *)SunMalloc(Ne_max*sizeof(int),"InitializeMerging");

  //Allocate temp variables for merging grid variables
  Nk_temp=(int *)SunMalloc(Nc_max*sizeof(int),"InitializeMerging");

  
  if(myproc==0) {
    mnptr_all=(int **)SunMalloc(numprocs*sizeof(int *),"InitializeMerging");
    //eptr_all=(int **)SunMalloc(numprocs*sizeof(int *),"InitializeMerging");
    for(p=0;p<numprocs;p++) {
      mnptr_all[p]=(int *)SunMalloc(Nc_all[p]*sizeof(int),"InitializeMerging");
      for(i=0;i<Nc_all[p];i++)
	mnptr_all[p][i]=0;

//      eptr_all[p]=(int *)SunMalloc(Ne_all[p]*sizeof(int),"InitializeMerging");
//      for(i=0;i<Ne_all[p];i++)
//	eptr_all[p][i]=0;

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
      for(i=0;i<Nc_all[p];i++)
	mergedGrid->Nk[mnptr_all[p][i]]=Nk_temp[i];
    }
  }

  // Merge all of the grid variables onto one processor
  MergeGridVariables(grid, numprocs, myproc, comm);

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
    merged3DVector=(REAL *)SunMalloc(mergedGrid->Nc*(grid->Nkmax+1)*sizeof(REAL),"InitializeMerging");

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
  
    mergedGrid->Nke = (int *)SunMalloc(mergedGrid->Ne*sizeof(int),"MergeGridVariables");
    mergedGrid->grad = (int *)SunMalloc(2*mergedGrid->Ne*sizeof(int),"MergeGridVariables");
    mergedGrid->gradf = (int *)SunMalloc(2*mergedGrid->Ne*sizeof(int),"MergeGridVariables");
    mergedGrid->edges = (int *)SunMalloc(mergedGrid->Ne*NUMEDGECOLUMNS*sizeof(int),"MergeGridVariables");
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
	 }
      }
    }
  }



  // Free up the swap arrays
  SunFree(Nc_int,Nc_max*sizeof(int),"MergeGridVariables");
  SunFree(NcNs_int,grid->maxfaces*Nc_max*sizeof(int),"MergeGridVariables");
  SunFree(Nc_real,Nc_max*sizeof(REAL),"MergeGridVariables");

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
