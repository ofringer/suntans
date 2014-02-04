/*
 * File: tides.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains functions that read/write data for specification of tidal
 * components at boundaries of type 2.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "mympi.h"
#include "grid.h"
#include "tides.h"
#include "memory.h"

static void InitTidalArrays(int N);
static void ReadTidalArrays(FILE *ifid, char *istr, int N);

/*
 * Function: SetTideComponents
 * Usage: SetTideComponents(grid,myproc);
 * --------------------------------------
 * This function looks for the tidal data in the file specified by
 * TideInput in suntans.dat and reads in the appropriate tidal data.
 * If that file does not exist, then the locations of the x-y points
 * of the boundaries are output to the file specified by TideOutput
 * in suntans.dat.
 *
 * Note that velocity mags must be in m/s and phases are in rad/s!!
 *
 * Format of the binary TideInput file:
 * ------------------------------------
 * numtides (1 X int) number of tidal constituents.
 * numboundaryedges (1 X int) number of boundary edges.
 * omegas (numtides X REAL) frequencies of tidal components.
 * u_amp[0] (numtides X REAL) amplitude of easting velocity at location 0.
 * u_phase[0] (numtides X REAL) phase of easting velocity at location 0.
 * v_amp[0] (numtides X REAL) amplitude of northing velocity at location 0.
 * v_phase[0] (numtides X REAL) phase of northing velocity at location 0.
 * h_amp[0] (numtides X REAL) amplitude of free surface at location 0.
 * h_phase[0] (numtides X REAL) phase of free surface at location 0.
 * u_amp[1] (numtides X REAL) amplitude of easting velocity at location 1.
 * u_phase[1] (numtides X REAL) phase of easting velocity at location 1.
 * v_amp[1] (numtides X REAL) amplitude of northing velocity at location 1.
 * v_phase[1] (numtides X REAL) phase of northing velocity at location 1.
 * h_amp[1] (numtides X REAL) amplitude of free surface at location 1.
 * h_phase[1] (numtides X REAL) phase of free surface at location 1.
 * ...
 * u_amp[numboundaryedges-1] (numtides X REAL) amplitude of easting velocity at last location.
 * u_phase[numboundaryedges-1] (numtides X REAL) phase of easting velocity at last location.
 * v_amp[numboundaryedges-1] (numtides X REAL) amplitude of northing velocity at last location.
 * v_phase[numboundaryedges-1] (numtides X REAL) phase of northing velocity at last location.
 * h_amp[numboundaryedges-1] (numtides X REAL) amplitude of free surface at last location.
 * h_phase[numboundaryedges-1] (numtides X REAL) phase of free surface at last location 
 *
 */
void SetTideComponents(gridT *grid, int myproc) {
  int i, j, iptr, jptr, numboundaryedges;
  char istr[BUFFERLENGTH], ostr[BUFFERLENGTH], filename[BUFFERLENGTH];
  FILE *ifid, *ofid;

  MPI_GetFile(filename,DATAFILE,"TideInput","SetTideComponents",myproc);
  sprintf(istr,"%s.%d",filename,myproc);
  MPI_GetFile(filename,DATAFILE,"TideOutput","SetTideComponents",myproc);
  sprintf(ostr,"%s.%d",filename,myproc);
  
  if(VERBOSE>2) printf("Set Tidecomponents on proc %d\n",myproc);
  if((ifid=fopen(istr,"r"))==NULL) {
    printf("Error opening %s!\n",istr);
    printf("Writing x-y boundary locations to %s instead.\n",ostr);

    ofid=fopen(ostr,"w");
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      
      fprintf(ofid,"%f %f\n",grid->xe[j],grid->ye[j]);
    }

    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i=grid->cellp[iptr];
      
      fprintf(ofid,"%f %f\n",grid->xv[i],grid->yv[i]);
    }
    fclose(ofid);

    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else {
    if(VERBOSE>2 && myproc==0) printf("Reading tidal data.\n");

    // Read in number of tidal components from file
    fread(&numtides,sizeof(int),1,ifid);
    // Read in number of edges in file
    fread(&numboundaryedges,sizeof(int),1,ifid);

    if(numboundaryedges!=grid->edgedist[3]-grid->edgedist[2]+grid->celldist[2]-grid->celldist[1]) {
      printf("Error reading %s.  Number of edges does not match current run!\n",istr);

      MPI_Finalize();
      exit(EXIT_FAILURE);
    } else {
      InitTidalArrays(numboundaryedges);
      ReadTidalArrays(ifid,istr,numboundaryedges);
    }
  }
}

static void InitTidalArrays(int N) {
  int j;

  u_amp = (REAL **)SunMalloc(N*sizeof(REAL *),"InitTidalArrays");
  v_amp = (REAL **)SunMalloc(N*sizeof(REAL *),"InitTidalArrays");
  h_amp = (REAL **)SunMalloc(N*sizeof(REAL *),"InitTidalArrays");
  u_phase = (REAL **)SunMalloc(N*sizeof(REAL *),"InitTidalArrays");
  v_phase = (REAL **)SunMalloc(N*sizeof(REAL *),"InitTidalArrays");
  h_phase = (REAL **)SunMalloc(N*sizeof(REAL *),"InitTidalArrays");

  omegas = (REAL *)SunMalloc(numtides*sizeof(REAL),"InitTidalArrays");
  for(j=0;j<N;j++) {
    u_amp[j] = (REAL *)SunMalloc(numtides*sizeof(REAL),"InitTidalArrays");
    u_phase[j] = (REAL *)SunMalloc(numtides*sizeof(REAL),"InitTidalArrays");
    v_amp[j] = (REAL *)SunMalloc(numtides*sizeof(REAL),"InitTidalArrays");
    v_phase[j] = (REAL *)SunMalloc(numtides*sizeof(REAL),"InitTidalArrays");
    h_amp[j] = (REAL *)SunMalloc(numtides*sizeof(REAL),"InitTidalArrays");
    h_phase[j] = (REAL *)SunMalloc(numtides*sizeof(REAL),"InitTidalArrays");
  }
}

static void ReadTidalArrays(FILE *ifid, char *istr, int N) {
  int j, t, flag=0;
  
  if(fread(omegas,sizeof(REAL),numtides,ifid)!=numtides) flag=1;

  for(j=0;j<N && !flag;j++) {
    if(fread(u_amp[j],sizeof(REAL),numtides,ifid)!=numtides ||
       fread(u_phase[j],sizeof(REAL),numtides,ifid)!=numtides ||
       fread(v_amp[j],sizeof(REAL),numtides,ifid)!=numtides ||
       fread(v_phase[j],sizeof(REAL),numtides,ifid)!=numtides ||
       fread(h_amp[j],sizeof(REAL),numtides,ifid)!=numtides ||
       fread(h_phase[j],sizeof(REAL),numtides,ifid)!=numtides)
      flag=1;
  }
  if(flag)
    printf("Error reading tidal data.  Only %d of %d edges in %s.\n",j,N,istr);
}
    
  

  
    
