/*
 * File: physio.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Functions for reading/writing physical data.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "physio.h"
#include "merge.h"
#include "sendrecv.h"
#include "mynetcdf.h"

/************************************************************************/
/*                                                                      */
/*                    All functions are public.                         */
/*                                                                      */
/************************************************************************/

/*
 * Function: Write2DData
 * Usage: Write2DData(phys->h,prop->mergeArrays,fid,"Error outputting h data!",grid,numprocs,myproc,comm);
 * -------------------------------------------------------------------------------------------------------
 * Write a 2D array in the pointer *array to the file pointer fid, usually for free-surface data.  
 *
 */
void Write2DData(REAL *array, int merge, FILE *fid, char *error_message, 
		 gridT *grid, int numprocs, int myproc, MPI_Comm comm) {
  int arraySize, writeProc, nwritten;
  REAL *array2DPointer;

  if(merge) {
    MergeCellCentered2DArray(array,grid,numprocs,myproc,comm);

    if(myproc==0) {
      arraySize=mergedGrid->Nc;
      array2DPointer=merged2DArray;
    }
    writeProc=0;
  } else {
    arraySize=grid->Nc;
    array2DPointer=array;
    writeProc=myproc;
  }

  if(myproc==writeProc) {
    nwritten=fwrite(array2DPointer,sizeof(REAL),arraySize,fid);
    if(nwritten!=arraySize) {
      printf("%s",error_message);
      exit(EXIT_WRITING);
    }
  }
  fflush(fid);
}

/*
 * Function: Write3DData
 * Usage: Write3DData(phys->s,prop->mergeArrays,fid,"Error outputting salinity data!",grid,numprocs,myproc,comm);
 * --------------------------------------------------------------------------------------------------------------
 * Write a 3D array in the pointer *array to the file pointer fid.  This array can be any 3D array with the
 * same size as phys->s[Nc][Nkmax].
 *
 */
void Write3DData(REAL **array, REAL *temp_array, int merge, FILE *fid, char *error_message, 
		 gridT *grid, int numprocs, int myproc, MPI_Comm comm) {
  int i, k, nwritten;

  if(merge) {
    MergeCellCentered3DArray(array,grid,numprocs,myproc,comm);

    if(myproc==0) {
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<mergedGrid->Nc;i++) {
	  if(k<mergedGrid->Nk[i])
	    merged2DArray[i]=merged3DArray[i][k];
	  else
	    merged2DArray[i]=EMPTY;
	}
	nwritten=fwrite(merged2DArray,sizeof(REAL),mergedGrid->Nc,fid);
	if(nwritten!=mergedGrid->Nc) {
	  printf("%s",error_message);
	  exit(EXIT_WRITING);
	}
      }
    }
  } else {
    for(k=0;k<grid->Nkmax;k++) {
      for(i=0;i<grid->Nc;i++) {
	if(k<grid->Nk[i])
	  temp_array[i]=array[i][k];
	else
	  temp_array[i]=EMPTY;
      }
      nwritten=fwrite(temp_array,sizeof(REAL),grid->Nc,fid);
      if(nwritten!=grid->Nc) {
	printf("%s",error_message);
	exit(EXIT_WRITING);
      }
    }
  }
  fflush(fid);
}

/* 
 * Function: OpenFiles
 * Usage: OpenFiles(prop,myproc);
 * ------------------------------
 * Open all of the files used for i/o to store the file pointers.
 *
 */
void OpenFiles(propT *prop, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  if(prop->readSalinity && prop->readinitialnc == 0) {
    MPI_GetFile(filename,DATAFILE,"InitSalinityFile","OpenFiles",myproc);
    prop->InitSalinityFID = MPI_FOpen(filename,"r","OpenFiles",myproc);
  }
  if(prop->readTemperature && prop->readinitialnc == 0) {
    MPI_GetFile(filename,DATAFILE,"InitTemperatureFile","OpenFiles",myproc);
    prop->InitTemperatureFID = MPI_FOpen(filename,"r","OpenFiles",myproc);
  }
  if(prop->readinitialnc>0){
    MPI_GetFile(filename,DATAFILE,"initialNCfile","OpenFiles",myproc);
    prop->initialNCfileID = MPI_NCOpen(filename,NC_NOWRITE,"OpenFiles",myproc);
  }
  if(prop->netcdfBdy>0){
    MPI_GetFile(filename,DATAFILE,"netcdfBdyFile","OpenFiles",myproc);
    //if(myproc==0){
	prop->netcdfBdyFileID = MPI_NCOpen(filename,NC_NOWRITE,"OpenFiles",myproc);
    //}else{
    //    prop->netcdfBdyFileID = -1;
    //}
  }
  if(prop->metmodel>0){
    MPI_GetFile(filename,DATAFILE,"metfile","OpenFiles",myproc);
    prop->metncid = MPI_NCOpen(filename,NC_NOWRITE,"OpenFiles",myproc);
  }

  if(prop->calcaverage && prop->mergeArrays==0){
    MPI_GetFile(filename,DATAFILE,"averageNetcdfFile","OpenFiles",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    prop->averageNetcdfFileID = MPI_NCOpen(str,NC_NETCDF4,"OpenFiles",myproc);
  }
  
  if(prop->outputNetcdf==0) {
    MPI_GetFile(filename,DATAFILE,"FreeSurfaceFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->FreeSurfaceFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"HorizontalVelocityFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->HorizontalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"VerticalVelocityFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->VerticalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"SalinityFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->SalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"BGSalinityFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->BGSalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"TemperatureFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->TemperatureFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"PressureFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->PressureFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"EddyViscosityFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->EddyViscosityFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    MPI_GetFile(filename,DATAFILE,"ScalarDiffusivityFile","OpenFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    prop->ScalarDiffusivityFID = MPI_FOpen(str,"w","OpenFiles",myproc);
    
    // No longer writing to verticalgridfile
    
  }else {
    if(prop->mergeArrays==0){
	MPI_GetFile(filename,DATAFILE,"outputNetcdfFile","OpenFiles",myproc);
	sprintf(str,"%s.nc.%d",filename,myproc);
	prop->outputNetcdfFileID = MPI_NCOpen(str,NC_NETCDF4,"OpenFiles",myproc);
    }
  }

  if(RESTART) {
    MPI_GetFile(filename,DATAFILE,"StartFile","OpenFiles",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    prop->StartFID = MPI_FOpen(str,"r","OpenFiles",myproc);
  }

  if(myproc==0) {
    MPI_GetFile(filename,DATAFILE,"ConserveFile","OpenFiles",myproc);
    sprintf(str,"%s",filename);
    prop->ConserveFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  }
}

/*
 * Function: OutputPhysicalVariables
 * Usage: OutputPhysicalVariables(grid,phys,prop,myproc,numprocs,blowup,comm);
 * ---------------------------------------------------------------------------
 * Output the data every ntout steps as specified in suntans.dat
 * If this is the last time step or if the run is blowing up (blowup==1),
 * then output the data to the restart file specified by the file pointer
 * prop->StoreFID.
 *
 * Note that ASCII output is no longer implemented.
 *
 */
void OutputPhysicalVariables(gridT *grid, physT *phys, propT *prop,int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int i, j, jptr, k, nwritten, arraySize, writeProc;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  REAL *tmp = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"OutputData"), 
    *array2DPointer, **array3DPointer;

  if(!(prop->n%prop->ntconserve) && !blowup) {
    ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);
    if(myproc==0)
      fprintf(prop->ConserveFID,"%e %e %e %e %e %e %e %e\n",prop->rtime,phys->mass,phys->volume,
          phys->Ep-phys->Ep0,phys->Eflux1,phys->Eflux2,phys->Eflux3,phys->Eflux4);
  }

  if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) {

    if(myproc==0 && VERBOSE>1) {
      if(!blowup) {
        printf("Outputting data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      } else {
        printf("Outputting blowup data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      }
    }
    Write2DData(phys->h,prop->mergeArrays,prop->FreeSurfaceFID,"Error outputting free-surface data!\n",
    		grid,numprocs,myproc,comm);

    // compute quadratic interpolated estimates for velocity using pretty plot and don't redo work
    if(prop->prettyplot==1 && prop->interp != QUAD) {
      ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
      ComputeUC(phys->uc, phys->vc, phys,grid, myproc, QUAD);
    }

    // Output u, v, w.  Interpolate w from faces to obtain it at the cell-centers first
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
	phys->stmp2[i][k]=0.5*(phys->w[i][k]+phys->w[i][k+1]);

    Write3DData(phys->uc,phys->htmp,prop->mergeArrays,prop->HorizontalVelocityFID,
		"Error outputting uc-data!\n",grid,numprocs,myproc,comm);
    Write3DData(phys->vc,phys->htmp,prop->mergeArrays,prop->HorizontalVelocityFID,
		"Error outputting vc-data!\n",grid,numprocs,myproc,comm);
    Write3DData(phys->stmp2,phys->htmp,prop->mergeArrays,prop->HorizontalVelocityFID,
		"Error outputting wc-data!\n",grid,numprocs,myproc,comm);

    // Output face-centered vertical velocity data
    Write3DData(phys->w,phys->htmp,prop->mergeArrays,prop->VerticalVelocityFID,
		"Error outputting vertical velocity data!\n",grid,numprocs,myproc,comm);

    // Background salinity field
    if(prop->n==1+prop->nstart)
      Write3DData(phys->s0,phys->htmp,prop->mergeArrays,prop->BGSalinityFID,
		  "Error outputting background salinity data!\n",grid,numprocs,myproc,comm);

    // Salinity field
    Write3DData(phys->s,phys->htmp,prop->mergeArrays,prop->SalinityFID,
		"Error outputting salinity data!\n",grid,numprocs,myproc,comm);

    // Temperature field
    Write3DData(phys->T,phys->htmp,prop->mergeArrays,prop->TemperatureFID,
		"Error outputting temperature data!\n",grid,numprocs,myproc,comm);

    // Nonhydrostatic pressure
    Write3DData(phys->q,phys->htmp,prop->mergeArrays,prop->PressureFID,
		"Error outputting nonhydrostatic pressure data!\n",grid,numprocs,myproc,comm);

    if(prop->turbmodel) {
      // Eddy-viscosity
      Write3DData(phys->nu_tv,phys->htmp,prop->mergeArrays,prop->EddyViscosityFID,
		  "Error outputting eddy-viscosity data!\n",grid,numprocs,myproc,comm);
      Write3DData(phys->kappa_tv,phys->htmp,prop->mergeArrays,prop->ScalarDiffusivityFID,
		  "Error outputting scalar-diffusivity data!\n",grid,numprocs,myproc,comm);
    }
    
    // No longer outputting vertical grid data
  }

  if(prop->n==1)
    fclose(prop->BGSalinityFID);

  if(prop->n==prop->nsteps+prop->nstart) {
    fclose(prop->FreeSurfaceFID);
    fclose(prop->HorizontalVelocityFID);
    fclose(prop->VerticalVelocityFID);
    fclose(prop->SalinityFID);
    // No longer writing to vertical grid file
    if(myproc==0) fclose(prop->ConserveFID);
  }

  // probably should change to make a distinction between blowup and restarts
  if(!(prop->n%prop->ntoutStore) || blowup) {
    if(VERBOSE>1 && myproc==0) 
      printf("Outputting restart data at step %d\n",prop->n);

    MPI_GetFile(filename,DATAFILE,"StoreFile","OutputData",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    prop->StoreFID = MPI_FOpen(str,"w","OpenFiles",myproc);

    nwritten=fwrite(&(prop->n),sizeof(int),1,prop->StoreFID);

    fwrite(phys->h,sizeof(REAL),grid->Nc,prop->StoreFID);
    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->Cn_U2[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_W2[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    if(prop->turbmodel>=1) {
      for(i=0;i<grid->Nc;i++) 
        fwrite(phys->Cn_q[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
      for(i=0;i<grid->Nc;i++) 
        fwrite(phys->Cn_l[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

      for(i=0;i<grid->Nc;i++) 
        fwrite(phys->qT[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
      for(i=0;i<grid->Nc;i++) 
        fwrite(phys->lT[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    }
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->nu_tv[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->kappa_tv[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->qc[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);

    fclose(prop->StoreFID);
  }

  SunFree(tmp,grid->Ne*sizeof(REAL),"OutputData");
}

/*
 * Function: ReadPhysicalVariables
 * Usage: ReadPhysicalVariables(grid,phys,prop,myproc,comm);
 * ---------------------------------------------------------
 * This function reads in physical variables for a restart run
 * from the restart file defined by prop->StartFID.
 *
 */
void ReadPhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {

  int i, j;

  if(VERBOSE>1 && myproc==0) printf("Reading from rstore...\n");
  //fixdzz
  UpdateDZ(grid,phys,prop,-1); 

  if(fread(&(prop->nstart),sizeof(int),1,prop->StartFID) != 1)
    printf("Error reading prop->nstart\n");

  if(fread(phys->h,sizeof(REAL),grid->Nc,prop->StartFID) != grid->Nc)
    printf("Error reading phys->h\n");
  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j])
      printf("Error reading phys->Cn_U[j]\n");
  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->Cn_U2[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j]) //AB3
      printf("Error reading phys->Cn_U2[j]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_W[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_W2[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_W[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_R[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_T[i]\n");

  if(prop->turbmodel>=1) {
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->Cn_q[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->Cn_q[i]\n");
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->Cn_l[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->Cn_l[i]\n");

    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->qT[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->qT[i]\n");
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->lT[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->lT[i]\n");
  }
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->nu_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->nu_tv[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->kappa_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->kappa_tv[i]\n");

  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j])
      printf("Error reading phys->u[j]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StartFID) != grid->Nk[i]+1)
      printf("Error reading phys->w[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->q[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->qc[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->qc[i]\n");

  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->s[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->T[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->s0[i]\n");
  fclose(prop->StartFID);

  UpdateDZ(grid,phys,prop, 0);

  // cell centered velocity computed so that this does not 
  // need to be reconsidered 
  ComputeUC(phys->uc, phys->vc, phys,grid, myproc, prop->interp);

  ISendRecvCellData3D(phys->uc,grid,myproc,comm);
  ISendRecvCellData3D(phys->vc,grid,myproc,comm);

  // Set the density from s and T using the equation of state
  SetDensity(grid,phys,prop);
}
