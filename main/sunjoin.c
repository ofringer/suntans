/*
 * File: sunjoin.c
 * Institution: Stanford University
 * --------------------------------
 * Suntans netcdf joining post-processing program
 *
 */
#include "sunjoin.h"
#include "suntans.h"
#include "mympi.h"
#include "grid.h"
#include "phys.h"
#include "report.h"
#include "mynetcdf.h"

 
void JoinNetcdf(propT *prop, gridT *grid, physT *phys, metT *met);
static void alloc_mnptr(int ncproc, int Nclocal, ptrT *ptr, int myproc);
static void alloc_eptr(int ncproc, int Nelocal, ptrT *ptr, int myproc);
static void alloc_ptr(int Nc[NUMPROCS], int Ne[NUMPROCS], gridT *grid, ptrT **ptr);
void read_cell_2D(int ncproc, int Nc, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc);
void read_cell_3D(int ncproc, int Nc, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL **tmp3d, int myproc);
void read_cell_3Dw(int ncproc, int Nc, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL **tmp3d, int myproc);
void read_edge_3D(int ncproc, int Ne, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL **tmp3d, int myproc);
void read_cell_int_2D(int ncproc, int Nrows, int Ndim, char *vname, ptrT *ptr, gridT *grid, int *tmp2d, int myproc);
void read_edge_int_2D(int ncproc, int Nrows, int Ndim, char *vname, ptrT *ptr, gridT *grid, int *tmp2d, int myproc);
void read_edge_1D(int ncproc, int Nrows,  char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc);
void read_cell_1D(int ncproc, int Nrows,  char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc);
void read_cell_int_1D(int ncproc, int nrows, char *vname, ptrT *ptr, gridT *grid, int *tmp2d, int myproc);
void read_cell_2Dall(int ncproc, int Nc, int Ndim, char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc);

main(int argc, char *argv[])
{
  int myproc, numprocs, j;
  MPI_Comm comm;
  gridT *grid;
  physT *phys;
  propT *prop;
  metT *met;
  REAL t0,t1;

  t0=time(NULL);
  printf("###############################################\n");
  printf(" Initiating SUNTANS NetCDF joining program...\n");
  printf("###############################################\n");

  StartMpi(&argc,&argv,&comm,&myproc,&numprocs);
  // Same steps as suntans
  ParseFlags(argc,argv,myproc);
  GetGrid(&grid,myproc,numprocs,comm);
  ReadProperties(&prop,myproc);
  InitializeVerticalGrid(&grid,myproc);
  AllocatePhysicalVariables(grid,&phys,prop);
  AllocateMet(prop,grid,&met,myproc);
  // Joining function
  JoinNetcdf(prop,grid,phys,met);
  
  t1=time(NULL);
  printf("###############################################\n");
  printf(" Finished joining output files.\n");
  printf("\nTotal elapsed time: %.2f s\n",t1-t0);
  printf("###############################################\n");
  EndMpi(&comm);
} // End of main

/*
 * Function: JoinNetcdf()
 * ----------------------
 * Main function for joining netcdf files together
 *
 */
void JoinNetcdf(propT *prop, gridT *grid, physT *phys, metT *met){

    int myproc, i, j, n;
    int ntime, varid, retval;
    int timectr, filenum;
    int ncid[NUMPROCS], Nc[NUMPROCS], Ne[NUMPROCS];
    char str[BUFFERLENGTH], filename[BUFFERLENGTH], outfile[BUFFERLENGTH]; 
    char *basefile;
    size_t start1[1], count1[1];
    int blowup=0;

    // Step 1: Open all of the netcdf files
    for(myproc=0;myproc<NUMPROCS;myproc++){
	MPI_GetFile(filename,DATAFILE,"outputNetcdfFile","OpenFiles",myproc);
	sprintf(str,"%s.%d",filename,myproc);
	printf("Opening netcdf file: %s...\n",str);
	ncid[myproc] = MPI_NCOpen(str,NC_NOWRITE,"OpenFiles",myproc);

	// Read the length of the dimensions for each file
	Nc[myproc] = returndimlen(ncid[myproc],"Nc");
	Ne[myproc] = returndimlen(ncid[myproc],"Ne");
    }

    // Step 2: Read the number of time steps from the first file
    ntime = returndimlen(ncid[0],"time"); 
    printf("Total time steps = %d.\n",ntime);

    // Step 3: Update the mnptr and eptr arrays 
    alloc_ptr(Nc,Ne,grid,&ptr);
    for(myproc=0;myproc<NUMPROCS;myproc++){ 
        printf("Reading partitioned array pointers (proc: %d)\n",myproc);
    	alloc_mnptr(ncid[myproc],Nc[myproc],ptr,myproc);
    	alloc_eptr(ncid[myproc],Ne[myproc],ptr,myproc);
    }

    // Step 4: Reallocate the grid variables
    for(myproc=0;myproc<NUMPROCS;myproc++){ 
	read_cell_int_2D(ncid[myproc], Nc[myproc], NFACES, "cells", ptr, grid, grid->cells, myproc);
	read_cell_int_2D(ncid[myproc], Nc[myproc], NFACES, "face", ptr, grid, grid->face, myproc);
	read_cell_int_2D(ncid[myproc], Nc[myproc], NFACES, "neigh", ptr, grid, grid->neigh, myproc);
	read_cell_int_2D(ncid[myproc], Nc[myproc], NFACES, "normal", ptr, grid, grid->normal, myproc);
	read_cell_int_1D(ncid[myproc], Nc[myproc], "Nk",ptr, grid, grid->Nk, myproc);

	read_edge_int_2D(ncid[myproc], Ne[myproc], 2, "edges", ptr, grid, grid->edges, myproc);
	read_edge_int_2D(ncid[myproc], Ne[myproc], 2, "grad", ptr, grid, grid->grad, myproc);
	read_edge_int_2D(ncid[myproc], Ne[myproc], 1, "Nke", ptr, grid, grid->Nke, myproc);

	read_cell_1D(ncid[myproc], Nc[myproc], "xv", ptr, grid, grid->xv, myproc);
	read_cell_1D(ncid[myproc], Nc[myproc], "yv", ptr, grid, grid->yv, myproc);
	read_cell_1D(ncid[myproc], Nc[myproc], "Ac", ptr, grid, grid->Ac, myproc);
	read_cell_1D(ncid[myproc], Nc[myproc], "dv", ptr, grid, grid->dv, myproc);

	read_edge_1D(ncid[myproc], Ne[myproc], "xe", ptr, grid, grid->xe, myproc);
	read_edge_1D(ncid[myproc], Ne[myproc], "ye", ptr, grid, grid->ye, myproc);
	read_edge_1D(ncid[myproc], Ne[myproc], "n1", ptr, grid, grid->n1, myproc);
	read_edge_1D(ncid[myproc], Ne[myproc], "n2", ptr, grid, grid->n2, myproc);
	read_edge_1D(ncid[myproc], Ne[myproc], "df", ptr, grid, grid->df, myproc);
	read_edge_1D(ncid[myproc], Ne[myproc], "dg", ptr, grid, grid->dg, myproc);

	read_cell_2Dall(ncid[myproc], Nc[myproc], NFACES, "def", ptr, grid, grid->def, myproc);
    }

    // Step 5: Initialise the first combined file
    filenum=1;
    basefile=strndup(filename,strlen(filename)-3);
    //printf("Filename=%s\n",basefile);

    // This generates the file
    sprintf(outfile,"%s_%03d.nc",basefile,filenum);
    printf("\tCreating filename=%s\n",outfile);
    prop->outputNetcdfFileID = MPI_NCOpen(outfile,NC_NETCDF4|NC_CLASSIC_MODEL,"OpenFiles",myproc);
    InitialiseOutputNCugrid(prop, grid, phys, met, myproc);

    // Step 6: Loop through the time steps
    count1[0]=1;
    timectr=0;
    for(n=0;n<ntime;n++){
        timectr += 1;
    	// Read the time variable
	start1[0]=n;
	if ((retval = nc_inq_varid(ncid[0], "time", &varid)))
	    ERR(retval);
	if ((retval = nc_get_vara_double(ncid[0], varid, start1, count1, &prop->nctime )))
	    ERR(retval);

	 printf("Joining step: %d of %d, nctime: %f...\n",n,ntime,prop->nctime);

	 // Update the file
	 if(timectr==STEPSPERFILE){
	    timectr=0;	 
	    filenum+=1;
	    //Close the old file
	    MPI_NCClose(prop->outputNetcdfFileID);

	    //Create the new file	    
	    sprintf(outfile,"%s_%03d.nc",basefile,filenum);
	    printf("\tCreating filename=%s\n",outfile);
	    prop->outputNetcdfFileID = MPI_NCOpen(outfile,NC_NETCDF4|NC_CLASSIC_MODEL,"OpenFiles",myproc);
	    InitialiseOutputNCugrid(prop, grid, phys, met, myproc);
	 }



	 // Read the other variables for each processor...
	for(myproc=0;myproc<NUMPROCS;myproc++){ 
	   read_cell_2D(ncid[myproc], Nc[myproc], n, "eta", ptr, grid, phys->h, myproc);
	   read_cell_3D(ncid[myproc], Nc[myproc], n, "uc", ptr, grid, phys->uc, myproc);
	   read_cell_3D(ncid[myproc], Nc[myproc], n, "vc", ptr, grid, phys->vc, myproc);
	   read_cell_3Dw(ncid[myproc], Nc[myproc], n, "w", ptr, grid, phys->w, myproc);
	   read_cell_3D(ncid[myproc], Nc[myproc], n, "nu_v", ptr, grid, phys->nu_tv, myproc);
	   read_cell_3D(ncid[myproc], Nc[myproc], n, "dzz", ptr, grid, grid->dzz, myproc);
	   if(prop->beta>0)
	       read_cell_3D(ncid[myproc], Nc[myproc], n, "salt", ptr, grid, phys->s, myproc);
	   if(prop->gamma>0)
	       read_cell_3D(ncid[myproc], Nc[myproc], n, "temp", ptr, grid, phys->T, myproc);
	   if( (prop->gamma>0) || (prop->beta>0) )
	       read_cell_3D(ncid[myproc], Nc[myproc], n, "rho", ptr, grid, phys->rho, myproc);
	   read_edge_3D(ncid[myproc], Ne[myproc], n, "dzf", ptr, grid, grid->dzf, myproc);
	   read_edge_3D(ncid[myproc], Ne[myproc], n, "U", ptr, grid, phys->u, myproc);
	   if(prop->metmodel>0){
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Uwind", ptr, grid, met->Uwind, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Vwind", ptr, grid, met->Vwind, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Tair", ptr, grid, met->Tair, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Pair", ptr, grid, met->Pair, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "rain", ptr, grid, met->rain, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "RH", ptr, grid, met->RH, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "cloud", ptr, grid, met->cloud, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Hs", ptr, grid, met->Hs, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Hl", ptr, grid, met->Hl, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Hlw", ptr, grid, met->Hlw, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "Hsw", ptr, grid, met->Hsw, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "tau_x", ptr, grid, met->tau_x, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "tau_y", ptr, grid, met->tau_y, myproc);
	       read_cell_2D(ncid[myproc], Nc[myproc], n, "EP", ptr, grid, met->EP, myproc);
	   }
	}

	 // Write the data to the file
	  WriteOuputNC(prop, grid, phys, met, blowup, myproc);
 

    }

}// End JoinNetcdf

/*
 * Function: read_edge_3D()
 * -------------------------
 * Read a 3D edge based double array onto the phys variable
 */
void read_edge_3D(int ncproc, int Ne, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL **tmp3d, int myproc){
    int j, jloc, ind, k, retval, varid;
    size_t start[]={tstep,0,0};
    size_t count[]={1,grid->Nkmax,Ne};
    //REAL tmpvar[1][grid->Nkmax][Ne];

    //printf("Reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncproc, varid, start, count, &ptr->ncscratch[0] )))
	ERR(retval);

    for(j=0;j<Ne;j++){
    	jloc=ptr->eptr[myproc][j];
	for(k=0;k<grid->Nke[jloc];k++){
	    ind = k*Ne+j;
	    tmp3d[jloc][k]=ptr->ncscratch[ind];
	    //tmp3d[jloc][k]=tmpvar[0][k][j];
	}
    }

}//End function

/*
 * Function: read_cell_3D()
 * -------------------------
 * Read a 3D cell based double array onto the phys variable
 */
void read_cell_3D(int ncproc, int Nc, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL **tmp3d, int myproc){
    int i, iloc, k, ind, retval, varid;
    size_t start[]={tstep,0,0};
    size_t count[]={1,grid->Nkmax,Nc};
    //REAL tmpvar[1][grid->Nkmax][Nc];

    //printf("Reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncproc, varid, start, count, &ptr->ncscratch[0] )))
    //if ((retval = nc_get_vara_double(ncproc, varid, start, count, &tmpvar[0][0][0] )))
	ERR(retval);

    for(i=0;i<Nc;i++){
    	iloc=ptr->mnptr[myproc][i];
	for(k=0;k<grid->Nkmax;k++){
	    ind = k*Nc+i;
	    if(k < grid->Nk[iloc])
		tmp3d[iloc][k]=ptr->ncscratch[ind];
	    //tmp3d[iloc][k]=tmpvar[0][k][i];
	}
    }

}//End function

/*
 * Function: read_cell_3Dw()
 * -------------------------
 * Read a 3D cell based double array onto the phys variable
 */
void read_cell_3Dw(int ncproc, int Nc, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL **tmp3d, int myproc){
    int i, iloc, ind, k, retval, varid;
    size_t start[]={tstep,0,0};
    size_t count[]={1,grid->Nkmax+1,Nc};
    //REAL tmpvar[1][grid->Nkmax+1][Nc];

    //printf("Reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncproc, varid, start, count, &ptr->ncscratch[0] )))
	ERR(retval);

    for(i=0;i<Nc;i++){
    	iloc=ptr->mnptr[myproc][i];
	for(k=0;k<grid->Nk[iloc]+1;k++){
	    ind = k*Nc+i;
	    tmp3d[iloc][k]=ptr->ncscratch[ind];
	    //tmp3d[iloc][k]=tmpvar[0][k][i];
	}
    }

}//End function


/*
 * Function: read_cell_2D()
 * -------------------------
 * Read a 2D cell based double array onto the phys variable
 */
void read_cell_2D(int ncproc, int Nc, int tstep, char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc){
    int i, iloc, retval, varid;
    size_t start[]={tstep,0};
    size_t count[]={1,Nc};
    //REAL tmpvar[1][Nc];

    //printf("Reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncproc, varid, start, count, &ptr->ncscratch[0] )))
	ERR(retval);

    for(i=0;i<Nc;i++){
    	iloc=ptr->mnptr[myproc][i];
	tmp2d[iloc]=ptr->ncscratch[i];
    }

}//End function

/*
 * Function: read_cell_1D()
 * -------------------------
 * Read a 2D cell based double array onto the phys variable
 */
void read_cell_1D(int ncproc, int Nrows,  char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc){
    int i, iloc, k, ind, retval, varid;
    //REAL tmpvar[1][Nc];
    size_t start[]={0};
    size_t count[]={Nrows};

    //printf("Reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncproc, varid, start, count, &ptr->ncscratch[0] )))
	ERR(retval);

    for(i=0;i<Nrows;i++){
    	iloc=ptr->mnptr[myproc][i];
	tmp2d[iloc]=ptr->ncscratch[i];
    }

}//End function

/*
 * Function: read_edge_1D()
 * -------------------------
 * Read a 2D cell based double array onto the phys variable
 */
void read_edge_1D(int ncproc, int Nrows,  char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc){
    int i, iloc, k, ind, retval, varid;
    //REAL tmpvar[1][Nc];
    size_t start[]={0};
    size_t count[]={Nrows};

    //printf("Reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncproc, varid, start, count, &ptr->ncscratch[0] )))
	ERR(retval);

    for(i=0;i<Nrows;i++){
    	iloc=ptr->eptr[myproc][i];
	tmp2d[iloc]=ptr->ncscratch[i];
    }

}//End function

/*
 * Function: read_cell_2Dall()
 * -------------------------
 * Read a 2D cell based double array onto the phys variable
 */
void read_cell_2Dall(int ncproc, int Nc, int Ndim, char *vname, ptrT *ptr, gridT *grid, REAL *tmp2d, int myproc){
    int i, iloc, k, retval, varid;
    size_t start[]={0,0};
    size_t count[]={Nc,Ndim};
    //REAL tmpvar[1][Nc];

    //printf("Reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncproc, varid, start, count, &ptr->ncscratch[0] )))
	ERR(retval);

    for(i=0;i<Nc;i++){
    	iloc=ptr->mnptr[myproc][i];
	for(k=0;k<Ndim;k++){
	     tmp2d[iloc*Ndim+k]=ptr->ncscratch[i*Ndim+k];
	 }
    }

}//End function

/*
 * function: read_cell_int_2d()
 * -------------------------
 * read a 2d cell based double array onto the phys variable
 */
void read_cell_int_2D(int ncproc, int nrows, int ndim, char *vname, ptrT *ptr, gridT *grid, int *tmp2d, int myproc){
    int i, iloc, k, ind, retval, varid;
    //real tmpvar[1][nc];
    size_t start[]={0,0};
    size_t count[]={nrows,ndim};

    //printf("reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	err(retval);
    if ((retval = nc_get_vara_int(ncproc, varid, start, count, &ptr->nctmpint[0] )))
	err(retval);

    for(i=0;i<nrows;i++){
    	iloc=ptr->mnptr[myproc][i];
	for(k=0;k<ndim;k++){
	     tmp2d[iloc*ndim+k]=ptr->nctmpint[i*ndim+k];
	 }
    }

}//end function

/*
 * function: read_cell_int_1d()
 * -------------------------
 * read a 2d cell based double array onto the phys variable
 */
void read_cell_int_1D(int ncproc, int nrows, char *vname, ptrT *ptr, gridT *grid, int *tmp2d, int myproc){
    int i, iloc, k, ind, retval, varid;
    //real tmpvar[1][nc];
    size_t start[]={0};
    size_t count[]={nrows};

    //printf("reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	err(retval);
    if ((retval = nc_get_vara_int(ncproc, varid, start, count, &ptr->nctmpint[0] )))
	err(retval);

    for(i=0;i<nrows;i++){
    	iloc=ptr->mnptr[myproc][i];
	tmp2d[iloc]=ptr->nctmpint[i];
    }

}//end function


/*
 * function: read_edge_int_2d()
 * -------------------------
 * read a 2d cell based double array onto the phys variable
 */
void read_edge_int_2D(int ncproc, int nrows, int ndim, char *vname, ptrT *ptr, gridT *grid, int *tmp2d, int myproc){
    int i, iloc, k, ind, retval, varid;
    //real tmpvar[1][nc];
    size_t start[]={0,0};
    size_t count[]={nrows,ndim};

    //printf("reading variable: %s on processor: %d...\n",vname,myproc);
    if ((retval = nc_inq_varid(ncproc, vname, &varid)))
	err(retval);
    if ((retval = nc_get_vara_int(ncproc, varid, start, count, &ptr->nctmpint[0] )))
	err(retval);

    for(i=0;i<nrows;i++){
    	iloc=ptr->eptr[myproc][i];
	for(k=0;k<ndim;k++){
	     tmp2d[iloc*ndim+k]=ptr->nctmpint[i*ndim+k];
	 }
    }

}//end function



/*
 * Function: alloc_ptr()
 * -------------------------
 *
 */
static void alloc_ptr(int Nc[NUMPROCS], int Ne[NUMPROCS], gridT *grid, ptrT **ptr){
    int myproc, i, j;
    *ptr = (ptrT *)SunMalloc(sizeof(ptrT),"JoinNetcdf");
    (*ptr)->mnptr = (int **)SunMalloc(NUMPROCS*sizeof(int *),"JoinNetcdf");
    (*ptr)->eptr = (int **)SunMalloc(NUMPROCS*sizeof(int *),"JoinNetcdf");


    for(myproc=0;myproc<NUMPROCS;myproc++){ 
        printf("myproc: %d, Nc: %d, Ne: %d\n",myproc,Nc[myproc],Ne[myproc]);
	(*ptr)->mnptr[myproc] = (int *)SunMalloc(Nc[myproc]*sizeof(int),"JoinNetcdf");
	(*ptr)->eptr[myproc] = (int *)SunMalloc(Ne[myproc]*sizeof(int),"JoinNetcdf");
    }
    
    for(myproc=0;myproc<NUMPROCS;myproc++){ 
    	for(i=0;i<Nc[myproc];i++) (*ptr)->mnptr[myproc][i] = 0.0;
    	for(j=0;j<Ne[myproc];j++) (*ptr)->eptr[myproc][j] = 0.0;
    }


    // Allocate memory to the netcdf read in arrays
    (*ptr)->ncscratch = (REAL *)SunMalloc(grid->Nkmax*grid->Ne*sizeof(REAL),"JoinNetcdf");
    (*ptr)->nctmpint = (int *)SunMalloc(3*grid->Ne*sizeof(int),"JoinNetcdf");



}//End function

 /*
 * Function: alloc_mnptr()
 * -------------------------
 *
 */
static void alloc_mnptr(int ncproc, int Nclocal, ptrT *ptr, int myproc){
    int i, varid,retval;
    size_t start[]={0};
    size_t count[]={Nclocal};
    //int tmp[Nclocal];
    //printf("myproc: %d, Nc: %d, Nc(file): %d\n",myproc,Nclocal,returndimlen(ncproc,"Nc"));
    // Read the data into the array
    if ((retval = nc_inq_varid(ncproc, "mnptr", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_int(ncproc, varid,start, count, &ptr->nctmpint[0] )))
	ERR(retval);

    for(i=0;i<Nclocal;i++){
	ptr->mnptr[myproc][i]=ptr->nctmpint[i];
	//printf("myproc: %d, i: %i, mnptr: %i\n",myproc,i,tmp[i]);
    }
} // end alloc_mnptr

/*
 * Function: alloc_eptr()
 * -------------------------
 *
 */
static void alloc_eptr(int ncproc, int Nelocal, ptrT *ptr, int myproc){
    int j, varid,retval;
    size_t start[]={0};
    size_t count[]={Nelocal};
 
    //int tmp[Nelocal];

    // Read the data into the array
    if ((retval = nc_inq_varid(ncproc, "eptr", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_int(ncproc, varid, start, count, &ptr->nctmpint[0] )))
	ERR(retval);

    for(j=0;j<Nelocal;j++){
    	ptr->eptr[myproc][j]=ptr->nctmpint[j];
    }
}
