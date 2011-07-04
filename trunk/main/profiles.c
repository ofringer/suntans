/*
 * File: profiles.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Outputs data in vertical profiles at locations specified in an input file. 
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include<stdio.h>
#include<string.h>
#include "util.h"
#include "memory.h"
#include "grid.h"
#include "phys.h"
#include "suntans.h"
#include "profiles.h"

/*
 * All functions are private except for InterpData which is called from phys.c
 *
 */
void InitializeOutputIndices(gridT *grid, MPI_Comm comm, int numprocs, int myproc);
static int IsOnGrid(REAL x, REAL y, REAL *xp, REAL *yp, 
		    int *cells, int *celldist, int *cellp, int Np, int Nc);
static int InPolygon(REAL x, REAL y, REAL *xg, REAL *yg, int N);
static void OpenDataFiles(int myproc);
static void CloseDataFiles(void);
static void AllInitialWriteToFiles(gridT *grid, propT *prop, MPI_Comm comm, int numprocs, int myproc);
static void WriteAllProfileData(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int numprocs, int myproc);
static void Write2DData(REAL *data, REAL *tmp, REAL *tmp2, FILE *fid, MPI_Comm comm, int numprocs, int myproc);
static void Write3DData(REAL **data, int *Nk, int Nkmax, REAL *tmp, REAL *tmp2, FILE *fid, MPI_Comm comm, int numprocs, int myproc);
static void FreeProfileArrays(gridT *grid, int numprocs);
static int Merge2dIntData(int *local, int N, int *global, int *allN, MPI_Comm comm, int numprocs, int myproc);
static int Merge2dRealData(REAL *local, int N, REAL *global, int *allN, MPI_Comm comm, int numprocs, int myproc);
static void ResortIntData(int *in, int *out, int *order, int N, int n);
static void ResortRealData(REAL *in, REAL *out, int *order, int N, int n);
static int GetProfileVariables(void);
static int ContainsCharacter(char *string, char c);

/*
 * Function: InterpData
 * Usage: InterpData(grid,phys,prop,comm,numprocs,myproc);
 * -------------------------------------------------------
 * This function first determines the indices of the nearest neighbors to the given
 * input points on each processor, and then uses this to output the data at the
 * specified intervals in suntans.dat.  The variables in suntans.dat are given by:
 * 
 * ProfileVariables         default, all, none, or one or all of husbTqnk
 *    h = free-surface output to file FreeSurfaceFile.prof
 *    u = u,v,w output to file HorizontalVelocityFile.prof
 *    s = salinity output to file SalinityFile.prof
 *    b = background salinity output to file BGSalinityFile.prof
 *    T = temperature output to file TemperatureFile.prof
 *    q = nonhydrostatic pressure output to file PressureFile.prof
 *    n = eddy-viscosity output to file EddyViscosityFile.prof
 *    k = scalar-diffusivity output to file ScalarDiffusivityFile.prof
 *    all = husbTqnk
 *    none = no profile output (same as omitting ProfileVariables)
 *    default = husb
 * DataLocations            name of file containing the x-y coordinates of the locations
 *                          at which output profiles are desired.
 * ProfileDataFile          name of file containing information about the profile data.
 *                          See AllInitialWriteToFiles() below.
 * ntoutProfs               frequency at which profile data is written.
 * NkmaxProfs               number of vertical levels to output (starting from k=0).
 *                          Nkmax=0 implies Nkmax as specified in suntans.dat.
 * numInterpPoints          number of nearest neighbors to output in the vicinity of the
 *                          points specified in the DataLocations file.  If numInterpPoints=1
 *                          then the output will effectively be a nearest-neighbor output.
 *                          Note that the user must interpolate the data after it is output; this
 *                          code does not perform the interpolation.
 *
 * Omitting the DataLocations variable from suntans.dat ignores all code in this file.
 *
 * Here is an example of entries in suntans.dat:
 *
 * ProfileVariables	default          # Output u, s, s0, and h
 * DataLocations	dataxy.dat       # dataxy.dat contains column x-y data
 * ProfileDataFile	profdata.dat     # Information about profiles is in profdata.dat
 * ntoutProfs		10               # Output profile data every 10 time steps
 * NkmaxProfs		10               # Only output the top 10 z-levels
 * numInterpPoints	3                # Output data at the three nearest neighbors to each input point.
 *
 */
void InterpData(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int numprocs, int myproc) {
  char filename[BUFFERLENGTH], str[BUFFERLENGTH];
  FILE *ofid;

  if(prop->n==prop->nstart+1) {
    existProfs=GetProfileVariables();
    if(existProfs) {
      
      // Compute nearest-neighbor indices and determine the processor on which each
      // point is located.
      InitializeOutputIndices(grid,comm,numprocs,myproc);

      // Open the data files only if specified in the variable ProfileVariables in suntans.dat.
      OpenDataFiles(myproc);

      // Output information about number of points, Nkmax, ntoutProfs, and other information
      // relevant to the profile data.
      AllInitialWriteToFiles(grid,prop,comm,numprocs,myproc);
    }
  }

  // Output the profile data at the frequency specified by ntoutProfs
  if(existProfs && !(prop->n%ntoutProfs)) 
    WriteAllProfileData(grid,phys,prop,comm,numprocs,myproc);

  // Free up arrays associated with the profile data and close the
  // data files.
  if(existProfs && prop->n==prop->nstart+prop->nsteps) {
    FreeProfileArrays(grid,numprocs);
    if(myproc==0) CloseDataFiles();
  }
}    

/*
 * Function: InitializeOutputIndices
 * Usage: InitializeOutputIndices(grid,comm,numprocs,myproc);
 * ----------------------------------------------------------
 * This function loads the points from the data file onto the current processor and 
 * determines if those points lie on the current processor.  Nearest neighbor
 * indices are then computed for those points and stored for future use.
 *
 */
void InitializeOutputIndices(gridT *grid, MPI_Comm comm, int numprocs, int myproc) {
  int i, ni, Ndata, Np, *cells, nf, total2dtemp, total3dtemp, proc, tempSum;
  REAL x, y, *xp, *yp;
  char filename[BUFFERLENGTH], str[BUFFERLENGTH];
  FILE *ifid;

  // Get ntoutProfs, which is the frequency that output profiles is desired
  ntoutProfs = MPI_GetValue(DATAFILE,"ntoutProfs","InitializeOutputIndices",myproc);

  // Get NkmaxProfs, which allows the specification of the number of output levels starting
  // from k=0.  A value of 0 implies Nkmax.
  NkmaxProfs = MPI_GetValue(DATAFILE,"NkmaxProfs","InitializeOutputIndices",myproc);
  if(NkmaxProfs==0 || NkmaxProfs>grid->Nkmax)
    NkmaxProfs=grid->Nkmax;

  // Number of nearest neighbors required for the interpolation
  numInterpPoints = MPI_GetValue(DATAFILE,"numInterpPoints","InitializeOutputIndices",myproc);

  // Read in the cells data for this proc since it is not available anymore
  // as part of the grid struct.
  MPI_GetFile(str,DATAFILE,"cells","InitializeOutputIndices",myproc);
  sprintf(filename,"%s.%d",str,myproc);

  cells = (int *)SunMalloc(3*grid->Nc*sizeof(int),"InitializeOutputIndices");

  ifid = fopen(filename,"r");
  for(i=0;i<grid->Nc;i++) {
    getfield(ifid,str);
    getfield(ifid,str);
    for(nf=0;nf<NFACES;nf++) 
      cells[i*NFACES+nf]=(int)getfield(ifid,str);
    getfield(ifid,str);
    getfield(ifid,str);
    getfield(ifid,str);
  }
  fclose(ifid);

  // Points are also not avialable anymore in the grid struct.
  MPI_GetFile(filename,DATAFILE,"points","InitializeOutputIndices",myproc);
  Np = MPI_GetSize(filename,"InitializeOutputIndices",myproc);

  xp = (REAL *)SunMalloc(Np*sizeof(REAL),"InitializeOutputIndices");
  yp = (REAL *)SunMalloc(Np*sizeof(REAL),"InitializeOutputIndices");

  ifid=fopen(filename,"r");
  for(i=0;i<Np;i++) {
    xp[i]=getfield(ifid,str);
    yp[i]=getfield(ifid,str);
    getfield(ifid,str);
  }
  fclose(ifid);

  // Load in the x-y coordinates of the profile locations.
  MPI_GetFile(filename,DATAFILE,"DataLocations","InitializeOutputIndices",myproc);
  numTotalDataPoints = MPI_GetSize(filename,"InitializeOutputIndices",myproc);
  
  // dataIndices contains the index of the current data point in its original location
  // in the DataLocations file.  This is in order to output the data in the same order
  // in which it was loaded.  dataXY contains the x-y coordinates loaded in from the file.
  dataIndices = (int *)SunMalloc(numTotalDataPoints*sizeof(int),"InitializeOutputIndices");
  dataXY = (REAL *)SunMalloc(2*numTotalDataPoints*sizeof(REAL),"InitializeOutputIndices");
  
  ifid = fopen(filename,"r");
  numLocalDataPoints=0;
  for(i=0;i<numTotalDataPoints;i++) {
    x = getfield(ifid,str);
    y = getfield(ifid,str);
    if(IsOnGrid(x,y,xp,yp,cells,grid->celldist,grid->cellp,Np,grid->Nc)) {
      dataIndices[numLocalDataPoints]=i;
      dataXY[2*numLocalDataPoints]=x;
      dataXY[2*numLocalDataPoints+1]=y;
      numLocalDataPoints++;
    }
  }
  if(VERBOSE>3) printf("Profile output: found %d points on processor %d\n",numLocalDataPoints,myproc);

  // Determine the total number of points that were found on all processors, since some points
  // may lie outside of the entire domain.
  MPI_Reduce(&numLocalDataPoints,&(tempSum),1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&tempSum,1,MPI_INT,0,comm);
  if(myproc==0 && VERBOSE>3 && numTotalDataPoints!=tempSum) 
    printf("Found %d interpolation point(s) not within computational domain.\n",numTotalDataPoints-tempSum);
  numTotalDataPoints = tempSum;

  // Now determine the indices for the interpolation using a nearest-neighbor search.
  interpIndices = (int *)SunMalloc(numInterpPoints*numLocalDataPoints*sizeof(int),"InitializeOutputIndices");

  for(i=0;i<numLocalDataPoints;i++) 
    FindNearest(&(interpIndices[i*numInterpPoints]),grid->xv,grid->yv,grid->Nc,
    		numInterpPoints,dataXY[2*i],dataXY[2*i+1]);

  // Processor 0 needs to know about everyones sizes.  total2d stores the
  // number of points that are output in two dimensions on each processor, while total3d stores
  // the same for three dimensions.
  total2d = (int *)SunMalloc(numprocs*sizeof(int),"InitializeOutputIndices");
  total3d = (int *)SunMalloc(numprocs*sizeof(int),"InitializeOutputIndices");
  total2dtemp = numLocalDataPoints*numInterpPoints;
  total3dtemp = NkmaxProfs*numLocalDataPoints*numInterpPoints;

  // All processors need to know how much data they store
  MPI_Gather(&total2dtemp,1,MPI_INT,total2d,1,MPI_INT,0,comm); 
  MPI_Bcast(total2d,numprocs,MPI_INT,0,comm);

  MPI_Gather(&total3dtemp,1,MPI_INT,total3d,1,MPI_INT,0,comm); 
  MPI_Bcast(total3d,numprocs,MPI_INT,0,comm);

  // all2d is the total number of points to be output on all processors.
  if(myproc==0) {
    all2d = 0;
    for(proc=0;proc<numprocs;proc++) 
      all2d+=total2d[proc];
  }
  MPI_Bcast(&(all2d),1,MPI_INT,0,comm);
  
  // These are used to merge and resort data distributed over the processors
  // for writing to file.
  merge_tmp = (REAL *)SunMalloc(3*all2d*NkmaxProfs*sizeof(REAL),"WriteAllProfileData");
  merge_tmp2 = (REAL *)SunMalloc(3*all2d*NkmaxProfs*sizeof(REAL),"WriteAllProfileData");

  // Only free up space associated with the grid since this is no longer needed.
  SunFree(xp,Np*sizeof(REAL),"InitializeOutputIndices");
  SunFree(yp,Np*sizeof(REAL),"InitializeOutputIndices");
  SunFree(cells,3*grid->Nc*sizeof(int),"InitializeOutputIndices");
}

/*
 * Function: IsOnGrid
 * Usage: if(IsOnGrid(x,y,xp,yp,cells,grid->celldist,grid->cellp,Np,grid->Nc)) {} ...
 * ----------------------------------------------------------------------------------
 * Determines if the point (x,y) is on the grid associated with the current processor.
 * If so, it returns 1, otherwise 0.
 *
 */
static int IsOnGrid(REAL x, REAL y, REAL *xp, REAL *yp, int *cells, int *celldist, int *cellp, int Np, int Nc) {
  int i, iptr, nf;
  REAL xg[NFACES], yg[NFACES];

  // Need to make sure only check the computational cells, not the interprocessor boundary
  // cells.
  for(iptr=celldist[0];iptr<celldist[2];iptr++) {
    i = cellp[iptr];

    for(nf=0;nf<NFACES;nf++) {
      xg[nf]=xp[cells[i*NFACES+nf]];
      yg[nf]=yp[cells[i*NFACES+nf]];
    }
    if(InPolygon(x,y,xg,yg,NFACES))
      return 1;
  }
  return 0;
}

/*
 * Function: InPolygon
 * Usage: if(InPolygon(x,y,xg,yg,N))
 * ---------------------------------
 * Returns true if the point (x,y) is contained within the polygon defined by the
 * N points defined by xg,yg.
 *
 * Copyright (c) 1970-2003, Wm. Randolph Franklin
 * http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
 *
 */
static int InPolygon(REAL x, REAL y, REAL *xg, REAL *yg, int N) {

  int i, j=0, in=0;

  for (i=0; i<N; i++) {
    j++; if (j==N) j=0;
    if (yg[i]<y && yg[j]>=y || yg[j]<y && yg[i]>=y) 
      if (xg[i]+(y-yg[i])/(yg[j]-yg[i])*(xg[j]-xg[i])<x) 
        in=!in;
  }

  return in;
}

/*
 * Function: OpenDataFiles
 * Usage: OpenDataFiles(myproc);
 * -----------------------------
 * Open the data files on processor 0 only if they are specified in the 
 * variable ProfileVariables in suntans.dat.
 *
 */
static void OpenDataFiles(int myproc) {
  int status;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  if(strlen(ProfileVariables)>0) {
    if(ContainsCharacter(ProfileVariables,'h')) {
      MPI_GetFile(filename,DATAFILE,"FreeSurfaceFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	FreeSurfaceProfFID = fopen(str,"w");
      }
    }
    
    if(ContainsCharacter(ProfileVariables,'u')) {
      MPI_GetFile(filename,DATAFILE,"HorizontalVelocityFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	HorizontalVelocityProfFID = fopen(str,"w");
      }
    }
    
    if(ContainsCharacter(ProfileVariables,'s')) {
      MPI_GetFile(filename,DATAFILE,"SalinityFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	SalinityProfFID = fopen(str,"w");
      }
    }

    if(ContainsCharacter(ProfileVariables,'b')) {
      MPI_GetFile(filename,DATAFILE,"BGSalinityFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	BGSalinityProfFID = fopen(str,"w");
      }
    }
    
    if(ContainsCharacter(ProfileVariables,'T')) {
      MPI_GetFile(filename,DATAFILE,"TemperatureFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	TemperatureProfFID = fopen(str,"w");
      }
    }
    
    if(ContainsCharacter(ProfileVariables,'q')) {
      MPI_GetFile(filename,DATAFILE,"PressureFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	PressureProfFID = fopen(str,"w");
      }
    }
    
    if(ContainsCharacter(ProfileVariables,'n')) {
      MPI_GetFile(filename,DATAFILE,"EddyViscosityFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	EddyViscosityProfFID = fopen(str,"w");
      }
    }
    
    if(ContainsCharacter(ProfileVariables,'k')) {
      MPI_GetFile(filename,DATAFILE,"ScalarDiffusivityFile","OpenDataFiles",myproc);
      if(myproc==0) {
	sprintf(str,"%s.prof",filename);
	ScalarDiffusivityProfFID = fopen(str,"w");
      }
    }
    
    MPI_GetFile(filename,DATAFILE,"ProfileDataFile","OpenDataFiles",myproc);
    if(myproc==0) ProfileDataFID = fopen(filename,"w");
  }
}

/*
 * Function: CloseDataFiles
 * Usage: if(myproc==0) CloseDataFiles();
 * --------------------------------------
 * Close the data files only if they were opened in OpenDataFiles, making
 * sure to do so on processor 0 only.
 * 
 */
static void CloseDataFiles(void) {
  if(ContainsCharacter(ProfileVariables,'h'))
    fclose(FreeSurfaceProfFID);
  if(ContainsCharacter(ProfileVariables,'u'))
    fclose(HorizontalVelocityProfFID);
  if(ContainsCharacter(ProfileVariables,'s'))
    fclose(SalinityProfFID);
  if(ContainsCharacter(ProfileVariables,'b'))
    fclose(BGSalinityProfFID);
  if(ContainsCharacter(ProfileVariables,'T'))
    fclose(TemperatureProfFID);
  if(ContainsCharacter(ProfileVariables,'q'))
    fclose(PressureProfFID);
  if(ContainsCharacter(ProfileVariables,'n'))
    fclose(EddyViscosityProfFID);
  if(ContainsCharacter(ProfileVariables,'k'))
    fclose(ScalarDiffusivityProfFID);
}

/*
 * Function: AllInitialWriteToFiles
 * Usage: AllInitialWriteToFiles(grid,prop,comm,numprocs,myproc);
 * --------------------------------------------------------------
 * Write data pertaining to the profiles in binary form to the file specified by the
 * variable ProfileDataFile in suntans.dat.  The format of the output is given by:
 * 
 * (4 byte int)numTotalDataPoints: Number of data points found on all processors.  Note that
 *     that this could be different from the number specified since some may lie outside the domain.
 * (4 byte int)numInterpPoints: Number of nearest neighbors to each point used for interpolation.
 * (4 byte int)NkmaxProfs: Number of vertical levels output in the profiles.
 * (4 byte int)nsteps: Total number of time steps in the simulation.
 * (4 byte int)ntoutProfs: Frequency of profile output.  This implies a total of nsteps/ntoutProfs are output.
 * (8 byte double)dt: Time step size
 * (8 byte double array X NkmaxProfs)dz: Contains the vertical grid spacings.
 * (4 byte int array X numTotalDataPoints)allIndices: Contains the indices of each point that determines its
 *     original location in the data file.  This is mostly for debugging since the output data is resorted
 *     so that it is in the same order as it appeared in the data file.
 * (4 byte int array X 2*numTotalDataPoints)dataXY: Contains the original data points at (or near) which profiles
 *     are output.
 * (8 byte double array X numTotalDataPoints*numInterpPoints)xv: Array containing the x-locations of the nearest
 *     neighbors to the dataXY points.  If numInterpPoints=3, then the 3 closest neighbors to the point
 *     (dataXY[2*i],dataXY[2*i+1]) are (xv[3*i],yv[3*i]), (xv[3*i+1],yv[3*i+1]), (xv[3*i+2],yv[3*i+2]).
 * (8 byte double array X numTotalDataPoints*numInterpPoints)yv: Array containing the y-locations of the nearest
 *     neighbors to the dataXY points (see xv above).
 * 
 */
static void AllInitialWriteToFiles(gridT *grid, propT *prop, MPI_Comm comm, int numprocs, int myproc) {
  int i, ni, *allN, *tmp_int, *tmp_int2, size, proc, totalDataPoints;
  REAL *tmp_real, *tmp_real2;
  tmp_real = (REAL *)SunMalloc(2*all2d*sizeof(REAL),"AllInitialWriteToFiles");
  tmp_real2 = (REAL *)SunMalloc(2*all2d*sizeof(REAL),"AllInitialWriteToFiles");
  tmp_int = (int *)SunMalloc(2*all2d*sizeof(int),"AllInitialWriteToFiles");
  tmp_int2 = (int *)SunMalloc(2*all2d*sizeof(int),"AllInitialWriteToFiles");
  allN = (int *)SunMalloc(numprocs*sizeof(int),"AllInitialWriteToFiles");
  allIndices = (int *)SunMalloc(numTotalDataPoints*sizeof(int),"AllInitialWriteToFiles");

  if(myproc==0) {
    fwrite(&numTotalDataPoints,sizeof(int),1,ProfileDataFID);
    fwrite(&numInterpPoints,sizeof(int),1,ProfileDataFID);
    fwrite(&(NkmaxProfs),sizeof(int),1,ProfileDataFID);
    fwrite(&(prop->nsteps),sizeof(int),1,ProfileDataFID);
    fwrite(&ntoutProfs,sizeof(int),1,ProfileDataFID);
    fwrite(&prop->dt,sizeof(REAL),1,ProfileDataFID);
    fwrite(grid->dz,sizeof(REAL),NkmaxProfs,ProfileDataFID);
  }
  
  // Indices of data points indicating their locations in the original data file.  This
  // is used to reorder the data once it is sent back to processor 0.
  allN[myproc]=numLocalDataPoints;
  MPI_Gather(&(allN[myproc]),1,MPI_INT,allN,1,MPI_INT,0,comm);   
  size=Merge2dIntData(dataIndices,numLocalDataPoints,allIndices,allN,comm,numprocs,myproc);
  if(myproc==0) fwrite(allIndices,sizeof(int),size,ProfileDataFID);

  // Original data points specified in the file ProfileData in suntans.dat.
  allN[myproc]=2*numLocalDataPoints;
  MPI_Gather(&(allN[myproc]),1,MPI_INT,allN,1,MPI_INT,0,comm);   
  size=Merge2dRealData(dataXY,2*numLocalDataPoints,tmp_real,allN,comm,numprocs,myproc);
  if(myproc==0) {
    ResortRealData(tmp_real,tmp_real2,allIndices,numTotalDataPoints,2);
    fwrite(tmp_real2,sizeof(REAL),size,ProfileDataFID);
  }

  // x-coordinates of nearest-neighbors.
  allN[myproc]=numInterpPoints*numLocalDataPoints;
  MPI_Gather(&(allN[myproc]),1,MPI_INT,allN,1,MPI_INT,0,comm);   
  for(i=0;i<numLocalDataPoints;i++)
    for(ni=0;ni<numInterpPoints;ni++) 
      tmp_real2[i*numInterpPoints+ni] = grid->xv[interpIndices[i*numInterpPoints+ni]];
  size=Merge2dRealData(tmp_real2,numLocalDataPoints*numInterpPoints,tmp_real,allN,comm,numprocs,myproc);
  if(myproc==0) {
    ResortRealData(tmp_real,tmp_real2,allIndices,numTotalDataPoints,numInterpPoints);
    fwrite(tmp_real2,sizeof(REAL),size,ProfileDataFID);
  }

  // y-coordinates of nearest-neighbors.
  for(i=0;i<numLocalDataPoints;i++)
    for(ni=0;ni<numInterpPoints;ni++) 
      tmp_real2[i*numInterpPoints+ni] = grid->yv[interpIndices[i*numInterpPoints+ni]];
  size=Merge2dRealData(tmp_real2,numLocalDataPoints*numInterpPoints,tmp_real,allN,comm,numprocs,myproc);
  if(myproc==0) {
    ResortRealData(tmp_real,tmp_real2,allIndices,numTotalDataPoints,numInterpPoints);
    fwrite(tmp_real2,sizeof(REAL),size,ProfileDataFID);
  }

  if(myproc==0)
    fclose(ProfileDataFID);

  SunFree(tmp_real,2*all2d*sizeof(REAL),"AllInitialWriteToFiles");
  SunFree(tmp_real2,2*all2d*sizeof(REAL),"AllInitialWriteToFiles");
  SunFree(tmp_int,2*all2d*sizeof(int),"AllInitialWriteToFiles");
  SunFree(tmp_int2,2*all2d*sizeof(int),"AllInitialWriteToFiles");
  SunFree(allN,numprocs*sizeof(int),"AllInitialWriteToFiles");
}

/*
 * Function: WriteAllProfileData
 * Usage: WriteAllProfileData(grid,phys,prop,comm,numprocs,myproc);
 * ----------------------------------------------------------------
 * If the variable is specified in the ProfileData variable in suntans.dat, all processors
 * send their data for interpolation to processor 0, which resorts the data so that it is
 * in the same order as it was upon loading from the data file.  Then, it is written to
 * a raw binary file.  Write2DData writes free surface data to a file, while Write3DData writes
 * all other 3d arrays to a file.
 *
 */
static void WriteAllProfileData(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int numprocs, int myproc) {
  int i, ni, k, n, iloc, size;

  if(myproc==0 && VERBOSE>2)
    printf("Outputting profile data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);

  if(ContainsCharacter(ProfileVariables,'h')) 
    Write2DData(phys->h,merge_tmp,merge_tmp2,FreeSurfaceProfFID,comm,numprocs,myproc);

  if(ContainsCharacter(ProfileVariables,'u')) {
    Write3DData(phys->uc,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,HorizontalVelocityProfFID,comm,numprocs,myproc);
    Write3DData(phys->vc,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,HorizontalVelocityProfFID,comm,numprocs,myproc);
    Write3DData(phys->w,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,HorizontalVelocityProfFID,comm,numprocs,myproc);
  }

  if(ContainsCharacter(ProfileVariables,'s')) 
    Write3DData(phys->s,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,SalinityProfFID,comm,numprocs,myproc);

  if(ContainsCharacter(ProfileVariables,'b') && prop->n==prop->nstart+ntoutProfs) 
    Write3DData(phys->s0,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,BGSalinityProfFID,comm,numprocs,myproc);

  if(ContainsCharacter(ProfileVariables,'T')) 
    Write3DData(phys->T,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,TemperatureProfFID,comm,numprocs,myproc);

  if(ContainsCharacter(ProfileVariables,'q')) 
    Write3DData(phys->q,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,PressureProfFID,comm,numprocs,myproc);

  if(ContainsCharacter(ProfileVariables,'n')) 
    Write3DData(phys->nu_tv,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,EddyViscosityProfFID,comm,numprocs,myproc);

  if(ContainsCharacter(ProfileVariables,'k')) 
    Write3DData(phys->kappa_tv,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,ScalarDiffusivityProfFID,comm,numprocs,myproc);
}

/*
 * Function: Write2DData
 * Usage: Write2DData(phys->h,merge_tmp,merge_tmp2,FreeSurfaceProfFID,comm,numprocs,myproc);
 * -----------------------------------------------------------------------------------------
 * Gather the two-dimensional data into a temporary array and then send it all to processor 0 for
 * resorting, then write it to the specified file descriptor.  Data is written as one contiguous
 * block of 8-byte doubles containing numTotalDataPoints*numInterpPoints elements.
 * 
 */
static void Write2DData(REAL *data, REAL *tmp, REAL *tmp2, FILE *fid, MPI_Comm comm, int numprocs, int myproc) {
  int i, iloc, size, ni, k, n=0;

  for(i=0;i<numLocalDataPoints;i++) 
    for(ni=0;ni<numInterpPoints;ni++) {
      iloc = interpIndices[numInterpPoints*i+ni];
      tmp[n++] = data[iloc];
    }
  size=Merge2dRealData(tmp,total2d[myproc],tmp2,total2d,comm,numprocs,myproc);
  if(myproc==0) {
    ResortRealData(tmp2,tmp,allIndices,numTotalDataPoints,numInterpPoints);
    fwrite(tmp,sizeof(REAL),size,fid);
    fflush(fid);
  }
}

/*
 * Function: Write3DData
 * Usage: Write3DData(phys->kappa_tv,grid->Nk,NkmaxProfs,merge_tmp,merge_tmp2,ScalarDiffusivityProfFID,comm,numprocs,myproc);
 * --------------------------------------------------------------------------------------------------------------------------
 * Gather the three-dimensional data into a temporary array and then send it all to processor 0 for
 * resorting, then write it to the specified file descriptor.  Data is written as one contiguous
 * block of 8-byte doubles containing Nkmax*numTotalDataPoints*numInterpPoints elements.
 * 
 */
static void Write3DData(REAL **data, int *Nk, int Nkmax, REAL *tmp, REAL *tmp2, FILE *fid, MPI_Comm comm, int numprocs, int myproc) {
  int i, iloc, size, ni, k, n=0;

  for(i=0;i<numLocalDataPoints;i++) 
    for(ni=0;ni<numInterpPoints;ni++) {
      iloc = interpIndices[numInterpPoints*i+ni];
      if(Nkmax<Nk[iloc]) {
	for(k=0;k<Nkmax;k++) 
	  tmp[n++]=data[iloc][k];
      } else {
	for(k=0;k<Nk[iloc];k++) 
	  tmp[n++] = data[iloc][k];
	for(k=Nk[iloc];k<Nkmax;k++) 
	  tmp[n++] = EMPTY;
      }
    }
  size=Merge2dRealData(tmp,total3d[myproc],tmp2,total3d,comm,numprocs,myproc);
  if(myproc==0) {
    ResortRealData(tmp2,tmp,allIndices,numTotalDataPoints,Nkmax*numInterpPoints);
    fwrite(tmp,sizeof(REAL),size,fid);
    fflush(fid);
  }
}

/*
 * Function: Merge2dIntData
 * Usage: size=Merge2dIntData(localdata,nlocal,globaldata,sizes,comm,numprocs,myproc);
 * -----------------------------------------------------------------------------------
 * Merges integer data on all processors, each of length nlocal, into one array of length sum(sizes), 
 * where sizes is an array containing the length of the arrays on each processor.
 *
 */
static int Merge2dIntData(int *local, int N, int *global, int *allN, MPI_Comm comm, int numprocs, int myproc) {
  int i, offset, proc;

  if(myproc!=0) 
    MPI_Send((void *)local,N,MPI_INT,0,1,comm);     
  else {
    for(i=0;i<N;i++) 
      global[i]=local[i];

    offset=N;
    for(proc=1;proc<numprocs;proc++) {
      MPI_Recv((void *)(&(global[offset])),allN[proc],MPI_INT,proc,1,comm,MPI_STATUS_IGNORE);
      offset+=allN[proc];
    }
  }
  return offset;
}

/*
 * Function: Merge2dRealData
 * Usage: size=Merge2dRealData(localdata,nlocal,globaldata,sizes,comm,numprocs,myproc);
 * -----------------------------------------------------------------------------------
 * Merges floating-point data on all processors, each of length nlocal, into one array of length sum(sizes), 
 * where sizes is an array containing the length of the arrays on each processor.
 *
 */
static int Merge2dRealData(REAL *local, int N, REAL *global, int *allN, MPI_Comm comm, int numprocs, int myproc) {
  int i, offset, proc;

  if(myproc!=0) {
    MPI_Send((void *)local,N,MPI_DOUBLE,0,1,comm);     
  } else {
    for(i=0;i<N;i++) 
      global[i]=local[i];

    offset=N;
    for(proc=1;proc<numprocs;proc++) {
      MPI_Recv((void *)(&(global[offset])),allN[proc],MPI_DOUBLE,proc,1,comm,MPI_STATUS_IGNORE);
      offset+=allN[proc];
    }
  }
  return offset;
}

/*
 * Function: FreeProfileArrays
 * Usage: FreeProfileArrays(grid,numprocs);
 * ----------------------------------------
 * Free up space associated with global arrays.
 *
 */
static void FreeProfileArrays(gridT *grid, int numprocs) {
  SunFree(dataIndices,numTotalDataPoints*sizeof(int),"FreeProfileArrays");
  SunFree(dataXY,2*numTotalDataPoints*sizeof(REAL),"FreeProfileArrays");
  SunFree(interpIndices,numInterpPoints*numLocalDataPoints*sizeof(int),"FreeProfileArrays");
  SunFree(total2d,numprocs*sizeof(int),"FreeProfileArrays");
  SunFree(total3d,numprocs*sizeof(int),"FreeProfileArrays");
  SunFree(merge_tmp,3*all2d*NkmaxProfs*sizeof(REAL),"FreeProfileArrays");
  SunFree(merge_tmp2,3*all2d*NkmaxProfs*sizeof(REAL),"FreeProfileArrays");
}

/*
 * Function: ResortIntData
 * Usage: ResortIntData(in,out,order,N,n);
 * ---------------------------------------
 * Resorts a two-dimensional integer array of size N by n that is
 * stored in a one-dimensional array using the order[] array.
 *
 */
static void ResortIntData(int *in, int *out, int *order, int N, int n) {
  int i, j, ind;

  for(i=0;i<N;i++) {
    ind = order[i];

    for(j=0;j<n;j++) 
      out[ind*n+j]=in[i*n+j];
  }
}

/*
 * Function: ResortRealData
 * Usage: ResortRealData(in,out,order,N,n);
 * ---------------------------------------
 * Resorts a two-dimensional floating-point array of size N by n that is
 * stored in a one-dimensional array using the order[] array.
 *
 */
static void ResortRealData(REAL *in, REAL *out, int *order, int N, int n) {
  int i, j, ind;

  for(i=0;i<N;i++) {
    ind = order[i];

    for(j=0;j<n;j++) 
      out[ind*n+j]=in[i*n+j];
  }
}

/*
 * Function: GetProfileVariables
 * Usage: GetProfileVariables();
 * -----------------------------
 * Read the ProfileVariables array from suntans.dat and determine which profile
 * variables are to be output during the simulation.  For details see the description
 * for the function InterpData.
 *
 */
static int GetProfileVariables(void) {
  int i, length, status;

  GetString(ProfileVariables,DATAFILE,"ProfileVariables",&status);
  if(!status)
    return 0;

  if(!strcmp(ProfileVariables,"none")) {
    sprintf(ProfileVariables,"\0");
    return 0;
  } else if(!strcmp(ProfileVariables,"all")) {
    sprintf(ProfileVariables,"%s",ALLPROFILEVARIABLES);
  } else if(!strcmp(ProfileVariables,"default")) {
    sprintf(ProfileVariables,"%s",DEFAULTPROFILEVARIABLES);
  } else {
    length = strlen(ProfileVariables);
    for(i=0;i<length;i++)
      if(!ContainsCharacter(ALLPROFILEVARIABLES,ProfileVariables[i])) {
	printf("Error in format of ProfileVariables in suntans.dat: Invalid string \"%s\"\n",ProfileVariables);
	printf("ProfileVariables must be either \"all\", \"default\", \"none\", or one or all of \"%s\".\n",ALLPROFILEVARIABLES);
	MPI_Finalize();
	exit(EXIT_FAILURE);
      }
  }
  return 1;
}

/*
 * Function: ContainsCharacter
 * Usage: if(ContainsCharacter("hello",'h')) ...
 * ---------------------------------------------
 * Returns 1 if the given string contains the given character, and 0 otherwise.
 *
 */
static int ContainsCharacter(char *string, char c) {
  int i, length = strlen(string);
  for(i=0;i<length;i++)
    if(string[i]==c)
      return 1;
  return 0;
}
  
