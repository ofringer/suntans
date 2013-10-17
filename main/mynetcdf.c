/*
* NetCDF IO functions
* -------------------
* All functions either generic or specific involved in netcdf io should go here
* 
*/ 

#include "mynetcdf.h"

/***********************************************
* Private function
***********************************************/

const void* FillValue(int empty);
static void ravel(REAL **tmparray, REAL *tmpvec, gridT *grid);
static void ravelW(REAL **tmparray, REAL *tmpvec,gridT *grid);
static void ravelEdge(REAL **tmparray, REAL *tmpvec, gridT *grid);
static void nc_addattr(int ncid, int varid, char *attname, char *attvalue);
static void nc_addattr_int(int ncid, int varid, char *attname, int *attvalue);
static void nc_addattr_real(int ncid, int varid, char *attname, REAL *attvalue);

void nc_read_3D(int ncid, char *vname, size_t *start, size_t *count, REAL ***tmparray);
void nc_read_2D(int ncid, char *vname, size_t *start, size_t *count, REAL **tmparray, int myproc);
void nc_write_double(int ncid, char *vname, REAL *tmparray, int myproc);
void nc_write_int(int ncid, char *vname, int *tmparray, int myproc);
void nc_write_intvar(int ncid, char *vname, gridT *grid, int *tmparray, int myproc);
void nc_write_doublevar(int ncid, char *vname, gridT *grid, REAL *tmparray, int myproc);

/*########################################################
*
* General functions
*
*#########################################################*/

/*
* Function: nc_read_3D()
* ----------------------
* Reads a 3D array from a netcdf file and returns the output in an array (not a vector)
*
* Warning: there are no dimension checks performed here so be careful.
* The size of the array dimension should equal 'count' ie [n, k ,j]
*/

void nc_read_3D(int ncid, char *vname, size_t start[3], size_t count[3], REAL ***tmparray){

    int j, k, n, ii;
    int varid, retval;
    //REAL tmpvec[ (int)count[0] * (int)count[1] * (int)count[2] ];
    REAL outdata[(int)count[0]][(int)count[1]][(int)count[2]];

    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
//    if ((retval = nc_get_vara_double(ncid, varid, start, count, &tmparray[0][0][0]))) 
//	ERR(retval);    
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &outdata[0][0][0]))) 
	ERR(retval); 

    // Loop through and insert the vector values into an array
    for(n=0;n<(int)count[0];n++){
	for(k=0;k<(int)count[1];k++){
	    for(j=0;j<(int)count[2];j++){
		//Linear index
		//ii = n*(int)count[1]*(int)count[2]+k*(int)count[2]+j;
		//tmparray[n][k][j]=tmpvec[ii];
		tmparray[n][k][j]=outdata[n][k][j];
	    }
	}
    }

}// End function


/*
* Function: nc_read_2D()
* ----------------------
* Reads a 2D array from a netcdf file and returns the output in an array (not a vector)
*
* Warning: there are no dimension checks performed here so be careful.
* The size of the array dimension should equal 'count' ie [n,,j]
*/

void nc_read_2D(int ncid, char *vname, size_t start[2], size_t count[2], REAL **tmparray, int myproc){

    int j, n, ii;
    int varid, retval;
    //REAL tmpvec[ (int)count[0] * (int)count[1] ];
    REAL outdata[(int)count[0]][(int)count[1]];

  //  printf("nc_read_2d -- Proc: %d, vname: %s, start[%d][%d], count[%d][%d]",myproc,vname,(int)start[0],(int)start[1],(int)count[0],(int)count[1]);
  
    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    //if ((retval = nc_get_vara_double(ncid, varid, start, count, &tmparray[0][0]))) 
    //  ERR(retval); 
    
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &outdata[0][0]))) 
	ERR(retval); 
    // Loop through and insert the vector values into an array
    for(n=0;n<(int)count[0];n++){
	for(j=0;j<(int)count[1];j++){
	    //Linear index
	    //ii = n*(int)count[1]+j;
	    //tmparray[n][j]=tmpvec[ii];
	    tmparray[n][j]=outdata[n][j];
	    //printf("myproc: %d, start[0]: %d, n: %d of %d, j: %d of %d, outdata[n][j]: %f, tmparray[n][j]: %f\n",myproc, (int)start[0], n,(int)count[0],j,(int)count[1],outdata[n][j], tmparray[n][j]);
	}
    }
    //printf(" Done\n");
}// End function

/*
 * Function: nc_write_double()
 * --------------------------
 *
 * Wrapper function for writing a double variable
 *
 */
void nc_write_double(int ncid, char *vname, REAL *tmparray, int myproc){
    int varid, retval;

    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_double(ncid, varid, tmparray)))
      ERR(retval);

} //end function

/*
 * Function: nc_write_int()
 * --------------------------
 *
 * Wrapper function for writing an integer variable
 *
 */
void nc_write_int(int ncid, char *vname, int *tmparray, int myproc){
    int varid, retval;

    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_int(ncid, varid, tmparray)))
      ERR(retval);

} //end function

/*
 * Function: nc_write_intvar()
 * --------------------------
 *
 * Wrapper function for writing a variable length integer variable
 *
 */
void nc_write_intvar(int ncid, char *vname, gridT *grid, int *tmparray, int myproc){
    int varid, retval,j,nf;
    int *tmpvar;
    tmpvar = (int *)SunMalloc((grid->Nc*grid->maxfaces)*sizeof(int),"nc_write_intvar");

    for(j=0;j<grid->Nc;j++){
	for(nf=0;nf<grid->maxfaces;nf++){
	    if(nf < grid->nfaces[j])
	    	tmpvar[j*grid->maxfaces+nf]=tmparray[j*grid->maxfaces+nf];
	    else
	    	tmpvar[j*grid->maxfaces+nf]=(int)EMPTY;
	}
    }
    
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_int(ncid, varid, tmpvar)))
      ERR(retval);

} //end function

/*
 * Function: nc_write_doublevar()
 * --------------------------
 *
 * Wrapper function for writing a variable length double variable
 *
 */
void nc_write_doublevar(int ncid, char *vname, gridT *grid, REAL *tmparray, int myproc){
    int varid, retval,j,nf;
    REAL *tmpvar;
    tmpvar = (REAL *)SunMalloc(grid->Nc*grid->maxfaces*sizeof(REAL),"nc_write_intvar");

    for(j=0;j<grid->Nc;j++){
	for(nf=0;nf<grid->maxfaces;nf++){
	    if(nf < grid->nfaces[j])
	    	tmpvar[j*grid->maxfaces+nf]=tmparray[j*grid->maxfaces+nf];
	    else
	    	tmpvar[j*grid->maxfaces+nf]=(REAL)EMPTY;
	}
    }
    
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_double(ncid, varid, tmpvar)))
      ERR(retval);

} //end function



/*
* Function: getTimeRec()
* -----------------------------
*  Retuns the index of the first preceding time step in the vector time
*/
int getTimeRec(REAL nctime, REAL *time, int nt){
   int j;
   
   for(j=0;j<nt;j++){
      if (time[j]>=nctime)
	return j-1;
   }
   return nt;
}

/*###############################################################
*
* SUNTANS output file functions
*
#################################################################*/

/*
* Function: WriteOutputNC()
* -----------------------------
* Main function for writing SUNTANS output to netcdf file/s
* 
*/
void WriteOuputNC(propT *prop, gridT *grid, physT *phys, metT *met, int blowup, int myproc){
   int ncid = prop->outputNetcdfFileID;
   int varid, retval, k;
   // Start and count vectors for one, two and three dimensional arrays
   const size_t startone[] = {prop->nctimectr};
   const size_t countone[] = {1};
   const size_t starttwo[] = {prop->nctimectr,0};
   const size_t counttwo[] = {1,grid->Nc};
   size_t startthree[] = {prop->nctimectr,0,0};
   size_t countthree[] = {1,grid->Nkmax,grid->Nc};
   const size_t countthreew[] = {1,grid->Nkmax+1,grid->Nc};
   const REAL time[] = {prop->nctime};

   nc_set_log_level(3); // This helps with debugging errors
   
   //REAL *tmpvar, *tmpvarE;
   // Need to write the 3-D arrays as vectors
   //tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"WriteOutputNC");
   //tmpvarE = (REAL *)SunMalloc(grid->Ne*grid->Nkmax*sizeof(REAL),"WriteOutputNC");
   
   if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) {

    if(myproc==0 && VERBOSE>1){ 
      if(!blowup) 
        printf("Outputting data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      else
        printf("Outputting blowup data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
    }
    
    /* Write the time data*/
    if ((retval = nc_inq_varid(ncid, "time", &varid)))
	ERR(retval);
    if ((retval = nc_put_vara_double(ncid, varid, startone, countone, time )))
	ERR(retval);
    
    /* Write to the physical variables*/
    if ((retval = nc_inq_varid(ncid, "eta", &varid)))
	ERR(retval);
    if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, phys->h )))
 	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "uc", &varid)))
	ERR(retval);
    ravel(phys->uc, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "vc", &varid)))
	ERR(retval);
    ravel(phys->vc, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	ERR(retval);
      
    // write w at cell top and bottom
    if ((retval = nc_inq_varid(ncid, "w", &varid)))
	ERR(retval);
    ravelW(phys->w, phys->tmpvarW, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthreew, phys->tmpvarW )))
	ERR(retval);

    if ((retval = nc_inq_varid(ncid, "nu_v", &varid)))
	ERR(retval);
    ravel(phys->nu_tv, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	ERR(retval);
    
    // Tracers
     if(prop->beta>0){
       if ((retval = nc_inq_varid(ncid, "salt", &varid)))
	  ERR(retval);
      ravel(phys->s, phys->tmpvar, grid);
      if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	  ERR(retval);
     }
     
     if(prop->gamma>0){
	if ((retval = nc_inq_varid(ncid, "temp", &varid)))
	  ERR(retval);
	ravel(phys->T, phys->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	  ERR(retval);
     }
      
     if( (prop->gamma>0) || (prop->beta>0) ){ 
	if ((retval = nc_inq_varid(ncid, "rho", &varid)))
	  ERR(retval);
	ravel(phys->rho, phys->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	  ERR(retval);
     }

     if(prop->calcage>0){ 
	if ((retval = nc_inq_varid(ncid, "agec", &varid)))
	  ERR(retval);
	ravel(phys->agec, phys->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	  ERR(retval);

	if ((retval = nc_inq_varid(ncid, "agealpha", &varid)))
	  ERR(retval);
	ravel(phys->agealpha, phys->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	  ERR(retval);
     }

     // Vertical grid spacing
     if ((retval = nc_inq_varid(ncid, "dzz", &varid)))
       ERR(retval);
     ravel(grid->dzz, phys->tmpvar, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
       ERR(retval);

     countthree[2] = grid->Ne;
     if ((retval = nc_inq_varid(ncid, "dzf", &varid)))
       ERR(retval);
     ravelEdge(grid->dzf, phys->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvarE )))
        ERR(retval);

     // Edge normal velocity
     if ((retval = nc_inq_varid(ncid, "U", &varid)))
	ERR(retval);
     ravelEdge(phys->u, phys->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvarE )))
        ERR(retval);

     // Wind variables
     if(prop->metmodel>0){
       if ((retval = nc_inq_varid(ncid, "Uwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Uwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Vwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Vwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Tair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Tair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Pair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Pair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "rain", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->rain )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "RH", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->RH )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "cloud", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->cloud )))
	 ERR(retval);
       
       // Heat flux variables
       if ((retval = nc_inq_varid(ncid, "Hs", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hs )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hl", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hl )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hlw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hlw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hsw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hsw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_x", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->tau_x )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_y", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->tau_y )))
	 ERR(retval);
       
       if(prop->beta > 0.0){
	  if ((retval = nc_inq_varid(ncid, "EP", &varid)))
	     ERR(retval);
	  if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->EP )))
	     ERR(retval);
       }
     }
     
    /* Update the time counter*/
    prop->nctimectr += 1;  
   }
   
   // Free the temporary vector
   //SunFree(tmpvar,grid->Nc*grid->Nkmax*sizeof(REAL),"WriteOuputNC");
   //SunFree(tmpvarE,grid->Ne*grid->Nkmax*sizeof(REAL),"WriteOuputNC");
  
} // End of function



/*
* Function: InitialiseOutputNCugrid()
* ------------------------------------
*
* Initialises the output netcdf file/s
* Files conform to the UGRID-0.9 CF conventions
* 
* One file per processor
* The pointer to each file is stored in prop->outputNetcdfFileID
* 
*/
void InitialiseOutputNCugrid(propT *prop, gridT *grid, physT *phys, metT *met, int myproc){
   int ncid = prop->outputNetcdfFileID;
   int retval, k;
   int varid;
   int dimid_Nc, dimid_Ne , dimid_Np, dimid_time, dimid_numsides, dimid_Two, dimid_Nkw, dimid_Nk; 
   int dimidone[1];
   int dimidtwo[2];
   int dimidthree[3];
   int nofill=0;
   const size_t starttwo[] = {0,0};
   const size_t counttwo[] = {grid->Nkmax,grid->Nc};
   REAL *z_r;
   REAL *z_w;
   const int DEFLATE=1;
   const int DEFLATELEVEL=2;
   const REAL FILLVALUE = (REAL)EMPTY;


   //REAL *tmpvar;
   // Need to write the 3-D arrays as vectors
   //tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"InitialiseOutputNC");  
   
   /* Initialise the depth arrays */
   z_r = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"InitialiseOutputNCugrid");
   z_w = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");
   
   /**************
    *
    * Start writing...
    *
    **************/
   if(VERBOSE>1 && myproc==0) printf("Initialising output netcdf files...");
   
   // Set the netcdf time ctr to 0
   prop->nctimectr=0;
   
   /* Define the global attributes - this should be expanded to include model input parameters*/
   nc_addattr(ncid, NC_GLOBAL,"title","SUNTANS NetCDF output file");
   
   /********************************************************************** 
    *
    * Define the dimensions
    *
    **********************************************************************/
   if ((retval = nc_def_dim(ncid, "Nc", grid->Nc, &dimid_Nc)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Np", grid->Np, &dimid_Np)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Ne", grid->Ne, &dimid_Ne)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nk", grid->Nkmax, &dimid_Nk)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nkw", grid->Nkmax+1, &dimid_Nkw)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "numsides", grid->maxfaces, &dimid_numsides)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Two", 2, &dimid_Two)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid_time)))
	ERR(retval);
   
    /********************************************************************** 
    *
    * Define the grid topology variables and attributes
    *
    **********************************************************************/

    //suntans_mesh
    if ((retval = nc_def_var(ncid,"suntans_mesh",NC_INT,0,0,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","mesh_topology");
    nc_addattr(ncid, varid,"long_name","Topology data of 2D unstructured mesh");
    nc_addattr(ncid, varid,"topology_dimension","2");
    nc_addattr(ncid, varid,"node_coordinates","xp yp");
    nc_addattr(ncid, varid,"face_node_connectivity","cells");
    nc_addattr(ncid, varid,"edge_node_connectivity","edges");
    nc_addattr(ncid, varid,"face_coordinates","xv yv");
    nc_addattr(ncid, varid,"edge_coordinates","xe ye");
    nc_addattr(ncid, varid,"face_edge_connectivity","face");
    nc_addattr(ncid, varid,"edge_face_connectivity","grad");

    // cells
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"cells",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its corner nodes");
    //if ((retval = nc_put_var_int(ncid,varid, grid->cells)))
    //  ERR(retval);

    //face
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"face",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_edge_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its edges");
    //if ((retval = nc_put_var_int(ncid,varid, grid->face)))
    //  ERR(retval);

    //nfaces
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"nfaces",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Number of cell faces");


    //edges
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"edges",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two nodes it connects");
    //if ((retval = nc_put_var_int(ncid,varid, grid->edges)))
    //  ERR(retval);

    //neigh
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"neigh",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its neighbouring faces");
    //if ((retval = nc_put_var_int(ncid,varid, grid->neigh)))
    //  ERR(retval);

    //grad
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"grad",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two faces it connects ");
    //if ((retval = nc_put_var_int(ncid,varid, grid->grad)))
    //  ERR(retval);

    //mnptr
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"mnptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps face indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->mnptr)))
    //  ERR(retval);

    //eptr
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"eptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps edge indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->eptr)))
    //  ERR(retval);

   /********************************************************************** 
    *
    * Define the grid coordinate variables and attributes 
    *
    **********************************************************************/

    //xv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"xv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xv)))
    //   ERR(retval);
      
    //yv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"yv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yv)))
    //   ERR(retval);
       
    //xp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"xp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xp)))
    //   ERR(retval);
        
    //yp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"yp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yp)))
    //   ERR(retval);
         
    //xe
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"xe",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xe)))
    //   ERR(retval);
          
    //ye
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"ye",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->ye)))
    //   ERR(retval);
        
    /********************************************************************** 
    *
    * Define the grid metric variables 
    *
    **********************************************************************/

    //normal
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"normal",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Dot product of unique normal with outward normal of each edge");
    //if ((retval = nc_put_var_int(ncid,varid, grid->normal)))
    //  ERR(retval);

    //n1
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n1",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","x-component of the edge normal");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n1)))
    //  ERR(retval);

    //n2
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n2",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","y-component of the edge normal");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n2)))
    //  ERR(retval);

    //df
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"df",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","edge length");
    nc_addattr(ncid, varid,"units","m");
    //if ((retval = nc_put_var_double(ncid,varid, grid->df)))
    //  ERR(retval);

    //dg
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dg",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","distance between faces on either side of edge");
    nc_addattr(ncid, varid,"units","m");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dg)))
    //  ERR(retval);

    //def
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"def",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Distance between faces and edges");
    nc_addattr(ncid, varid,"units","m");
    //if ((retval = nc_put_var_double(ncid,varid, grid->def)))
    //  ERR(retval);

    //Ac
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"Ac",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Horizontal area of 2D mesh");
    nc_addattr(ncid, varid,"units","m2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->Ac)))
    //  ERR(retval);

    /********************************************************************** 
    *
    * Define the vertical grid variables and attributes 
    *
    **********************************************************************/
   //dz
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"dz",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","z layer spacing");
   nc_addattr(ncid, varid,"units","m");
   //if ((retval = nc_put_var_double(ncid,varid, grid->dz)))
   //  ERR(retval);

   // Calculate and write the vertical coordinate z levels 
   z_w[0]=0.0;
   for(k=0;k<grid->Nkmax;k++){
      z_w[k+1] = z_w[k] + grid->dz[k];
      if(k==0){
	 z_r[k] = grid->dz[k]*0.5;
      }else{
	 z_r[k] = z_r[k-1]+grid->dz[k];
      }
   }
   //z_r
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"z_r",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer mid points");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_r)))
   //  ERR(retval);

   //z_w
   dimidone[0] = dimid_Nkw;   
   if ((retval = nc_def_var(ncid,"z_w",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer edges");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_w)))
   //  ERR(retval);

   //Nk
   dimidone[0] = dimid_Nc;   
   if ((retval = nc_def_var(ncid,"Nk",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at face"); 
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nk)))
   //  ERR(retval);

   //Nke
   dimidone[0] = dimid_Ne;   
   if ((retval = nc_def_var(ncid,"Nke",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at edge");
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nke)))
   //  ERR(retval);

    //dv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"stanford_name","sea_floor_depth_below_geoid");
    nc_addattr(ncid, varid,"long_name","seafloor depth");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dv)))
    //  ERR(retval);

    //dzz
    dimidthree[0] = dimid_time;
    dimidthree[1] = dimid_Nk;
    dimidthree[2] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dzz",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","z layer spacing at faces");
    nc_addattr(ncid, varid,"units","m");

    //dzf
    dimidthree[0] = dimid_time;
    dimidthree[1] = dimid_Nk;
    dimidthree[2] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dzf",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","z layer spacing at edges");
    nc_addattr(ncid, varid,"units","m");

   
    
    //time
    dimidone[0] = dimid_time;
    if ((retval = nc_def_var(ncid,"time",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","time");
    nc_addattr(ncid, varid,"units","seconds since 1990-01-01 00:00:00");  

   /********************************************************************** 
    * 
    * Define the physical variables and attributes 
    *
    **********************************************************************/
   
   dimidtwo[0] = dimid_time;
   dimidtwo[1] = dimid_Nc;
   
   dimidthree[0] = dimid_time;
   dimidthree[1] = dimid_Nk;
   dimidthree[2] = dimid_Nc;
   
   // eta
   if ((retval = nc_def_var(ncid,"eta",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Sea surface elevation");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   //u
   if ((retval = nc_def_var(ncid,"uc",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Eastward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //v
   if ((retval = nc_def_var(ncid,"vc",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval);   
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Northward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //w
   dimidthree[1] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"w",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Vertical water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_w yv xv");  

   dimidthree[1] = dimid_Nk;
   
   //nu_v
   if ((retval = nc_def_var(ncid,"nu_v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Vertical eddy viscosity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //salinity
   if(prop->beta>0){
     if ((retval = nc_def_var(ncid,"salt",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Salinity");
    nc_addattr(ncid, varid,"units","ppt");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //temperature
   if(prop->gamma>0){
     if ((retval = nc_def_var(ncid,"temp",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval); 
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Water temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //rho
   if( (prop->gamma>0) || (prop->beta>0) ){
     if ((retval = nc_def_var(ncid,"rho",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Water density");
    nc_addattr(ncid, varid,"units","kg m-3");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //age
   if(prop->calcage>0){
     if ((retval = nc_def_var(ncid,"agec",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Age concentration");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
    
    if ((retval = nc_def_var(ncid,"agealpha",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Age alpha parameter");
    nc_addattr(ncid, varid,"units","seconds");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }

   //U
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"U",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Edge normal velocity");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

   // Meteorological variables (2-D) //
   
   if(prop->metmodel>0){
    // Uwind
    if ((retval = nc_def_var(ncid,"Uwind",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Vwind
    if ((retval = nc_def_var(ncid,"Vwind",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Tair
    if ((retval = nc_def_var(ncid,"Tair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval);
    nc_addattr(ncid, varid,"long_name","Air temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Pair
    if ((retval = nc_def_var(ncid,"Pair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Air pressure");
    nc_addattr(ncid, varid,"units","millibar");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // rain
    if ((retval = nc_def_var(ncid,"rain",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Rain fall rate");
    nc_addattr(ncid, varid,"units","kg m2 s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //RH
    if ((retval = nc_def_var(ncid,"RH",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Relative humidity");
    nc_addattr(ncid, varid,"units","Percent (%)");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //cloud
    if ((retval = nc_def_var(ncid,"cloud",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Cloud cover fraction");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
    // Surface flux variables //
    // Hs
    if ((retval = nc_def_var(ncid,"Hs",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Sensible heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hl
    if ((retval = nc_def_var(ncid,"Hl",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Latent heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hlw
    if ((retval = nc_def_var(ncid,"Hlw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Net longwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");      
    nc_addattr(ncid, varid,"positive","down");   

    // Hsw
    if ((retval = nc_def_var(ncid,"Hsw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Net shortwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // tau_x
    if ((retval = nc_def_var(ncid,"tau_x",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    // tau_y
    if ((retval = nc_def_var(ncid,"tau_y",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Northward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    //EP
    if(prop->beta > 0.0){
	if ((retval = nc_def_var(ncid,"EP",NC_DOUBLE,2,dimidtwo,&varid)))
	    ERR(retval); 
	nc_addattr(ncid, varid,"long_name","Evaporation minus precipiaton");
	nc_addattr(ncid, varid,"units","m s-1");
	nc_addattr(ncid, varid,"mesh","suntans_mesh");
	nc_addattr(ncid, varid,"location","face");
	nc_addattr(ncid, varid,"coordinates","time yv xv");   
    }

   }

   //End file definition mode
   if ((retval = nc_enddef(ncid)))
	ERR(retval);

   
   /**********************************************************
   *
   * Write data (needs to be done out of definition mode for classic model)
   *
   ****************************************************************/
   nc_write_intvar(ncid,"cells",grid,grid->cells,myproc);
   nc_write_intvar(ncid,"face",grid,grid->face,myproc);
   nc_write_int(ncid,"nfaces",grid->nfaces,myproc);
   nc_write_int(ncid,"edges",grid->edges,myproc);
   nc_write_intvar(ncid,"neigh",grid,grid->neigh,myproc);
   nc_write_int(ncid,"grad",grid->grad,myproc);
   //nc_write_int(ncid,"mnptr",grid->mnptr,myproc);
   //nc_write_int(ncid,"eptr",grid->eptr,myproc);
   
   nc_write_double(ncid,"xv",grid->xv,myproc);
   nc_write_double(ncid,"yv",grid->yv,myproc);
   nc_write_double(ncid,"xe",grid->xe,myproc);
   nc_write_double(ncid,"ye",grid->ye,myproc);
   nc_write_double(ncid,"xp",grid->xp,myproc);
   nc_write_double(ncid,"yp",grid->yp,myproc);

   nc_write_intvar(ncid,"normal",grid,grid->normal,myproc);
   nc_write_double(ncid,"n1",grid->n1,myproc);
   nc_write_double(ncid,"n2",grid->n2,myproc);
   nc_write_double(ncid,"df",grid->df,myproc);
   nc_write_double(ncid,"dg",grid->dg,myproc);
   nc_write_doublevar(ncid,"def",grid,grid->def,myproc);
   nc_write_double(ncid,"Ac",grid->Ac,myproc);

   nc_write_double(ncid,"dz",grid->dz,myproc);
   nc_write_double(ncid,"z_r",z_r,myproc);
   nc_write_double(ncid,"z_w",z_w,myproc);
   nc_write_int(ncid,"Nk",grid->Nk,myproc);
   nc_write_int(ncid,"Nke",grid->Nke,myproc);
   nc_write_double(ncid,"dv",grid->dv,myproc);

   // Free the temporary vectors
   //SunFree(z_r,grid->Nkmax*sizeof(REAL),"InitialiseOutputNCugrid");
   //SunFree(z_w,(grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");
   //SunFree(tmpvar,grid->Nc*grid->Nkmax,"InitialiseOutputNC");
   if(VERBOSE>1 && myproc==0) printf("Done.\n");

}// End function


/*
* Function: InitialiseAverageNCugrid()
* ------------------------------------
*
* Initialises the average netcdf file/s
* Files conform to the UGRID-0.9 CF conventions
* 
* One file per processor
* The pointer to each file is stored in prop->outputNetcdfFileID
* 
*/
void InitialiseAverageNCugrid(propT *prop, gridT *grid, averageT *average, int myproc){
   int ncid = prop->averageNetcdfFileID;
   int retval, k;
   int varid;
   int dimid_Nc, dimid_Ne , dimid_Np, dimid_time, dimid_numsides, dimid_Two, dimid_Nkw, dimid_Nk; 
   int dimidone[1];
   int dimidtwo[2];
   int dimidthree[3];
   int nofill=0;
   const size_t starttwo[] = {0,0};
   const size_t counttwo[] = {grid->Nkmax,grid->Nc};
   REAL *z_r;
   REAL *z_w;
   const int DEFLATE=1;
   const int DEFLATELEVEL=2;
   const REAL FILLVALUE = (REAL)EMPTY;


   //REAL *tmpvar;
   // Need to write the 3-D arrays as vectors
   //tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"InitialiseOutputNC");  
   
   /* Initialise the depth arrays */
   z_r = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"InitialiseOutputNCugrid");
   z_w = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");
   
   /**************
    *
    * Start writing...
    *
    **************/
   if(VERBOSE>1 && myproc==0) printf("Initialising output netcdf files...\n");
   
   // Set the netcdf time ctr to 0
   prop->avgtimectr=0;
   prop->avgctr=0;
   
   /* Define the global attributes - this should be expanded to include model input parameters*/
   nc_addattr(ncid, NC_GLOBAL,"title","SUNTANS NetCDF time-averaged file");

   nc_addattr_int(ncid,NC_GLOBAL,"ntaverage",&prop->ntaverage);
   nc_addattr_real(ncid,NC_GLOBAL,"dt",&prop->dt);
   
   /********************************************************************** 
    *
    * Define the dimensions
    *
    **********************************************************************/
   if ((retval = nc_def_dim(ncid, "Nc", grid->Nc, &dimid_Nc)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Np", grid->Np, &dimid_Np)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Ne", grid->Ne, &dimid_Ne)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nk", grid->Nkmax, &dimid_Nk)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nkw", grid->Nkmax+1, &dimid_Nkw)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "numsides", grid->maxfaces, &dimid_numsides)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Two", 2, &dimid_Two)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid_time)))
	ERR(retval);
   
    /********************************************************************** 
    *
    * Define the grid topology variables and attributes
    *
    **********************************************************************/

    //suntans_mesh
    if ((retval = nc_def_var(ncid,"suntans_mesh",NC_INT,0,0,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","mesh_topology");
    nc_addattr(ncid, varid,"long_name","Topology data of 2D unstructured mesh");
    nc_addattr(ncid, varid,"topology_dimension","2");
    nc_addattr(ncid, varid,"node_coordinates","xp yp");
    nc_addattr(ncid, varid,"face_node_connectivity","cells");
    nc_addattr(ncid, varid,"edge_node_connectivity","edges");
    nc_addattr(ncid, varid,"face_coordinates","xv yv");
    nc_addattr(ncid, varid,"edge_coordinates","xe ye");
    nc_addattr(ncid, varid,"face_edge_connectivity","face");
    nc_addattr(ncid, varid,"edge_face_connectivity","grad");

    // cells
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"cells",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its corner nodes");
    //if ((retval = nc_put_var_int(ncid,varid, grid->cells)))
    //  ERR(retval);

    //face
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"face",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_edge_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its edges");
    //if ((retval = nc_put_var_int(ncid,varid, grid->face)))
    //  ERR(retval);

    //edges
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"edges",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two nodes it connects");
    //if ((retval = nc_put_var_int(ncid,varid, grid->edges)))
    //  ERR(retval);

    //neigh
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"neigh",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its neighbouring faces");
    //if ((retval = nc_put_var_int(ncid,varid, grid->neigh)))
    //  ERR(retval);

    //grad
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"grad",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two faces it connects ");
    //if ((retval = nc_put_var_int(ncid,varid, grid->grad)))
    //  ERR(retval);

    //mnptr
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"mnptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps face indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->mnptr)))
    //  ERR(retval);

    //eptr
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"eptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps edge indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->eptr)))
    //  ERR(retval);

   /********************************************************************** 
    *
    * Define the grid coordinate variables and attributes 
    *
    **********************************************************************/

    //xv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"xv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xv)))
    //   ERR(retval);
      
    //yv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"yv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yv)))
    //   ERR(retval);
       
    //xp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"xp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xp)))
    //   ERR(retval);
        
    //yp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"yp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yp)))
    //   ERR(retval);
         
    //xe
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"xe",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xe)))
    //   ERR(retval);
          
    //ye
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"ye",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->ye)))
    //   ERR(retval);
        
    /********************************************************************** 
    *
    * Define the grid metric variables 
    *
    **********************************************************************/

    //normal
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"normal",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Dot product of unique normal with outward normal of each edge");
    //if ((retval = nc_put_var_int(ncid,varid, grid->normal)))
    //  ERR(retval);

    //n1
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n1",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","x-component of the edge normal");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n1)))
    //  ERR(retval);

    //n2
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n2",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","y-component of the edge normal");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n2)))
    //  ERR(retval);

    //df
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"df",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","edge length");
    nc_addattr(ncid, varid,"units","m");
    //if ((retval = nc_put_var_double(ncid,varid, grid->df)))
    //  ERR(retval);

    //dg
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dg",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","distance between faces on either side of edge");
    nc_addattr(ncid, varid,"units","m");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dg)))
    //  ERR(retval);

    //def
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"def",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Distance between faces and edges");
    nc_addattr(ncid, varid,"units","m");
    //if ((retval = nc_put_var_double(ncid,varid, grid->def)))
    //  ERR(retval);

    //Ac
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"Ac",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Horizontal area of 2D mesh");
    nc_addattr(ncid, varid,"units","m2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->Ac)))
    //  ERR(retval);

    /********************************************************************** 
    *
    * Define the vertical grid variables and attributes 
    *
    **********************************************************************/
   //dz
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"dz",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","z layer spacing");
   nc_addattr(ncid, varid,"units","m");
   //if ((retval = nc_put_var_double(ncid,varid, grid->dz)))
   //  ERR(retval);

   // Calculate and write the vertical coordinate z levels 
   z_w[0]=0.0;
   for(k=0;k<grid->Nkmax;k++){
      z_w[k+1] = z_w[k] + grid->dz[k];
      if(k==0){
	 z_r[k] = grid->dz[k]*0.5;
      }else{
	 z_r[k] = z_r[k-1]+grid->dz[k];
      }
   }
   //z_r
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"z_r",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer mid points");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_r)))
   //  ERR(retval);

   //z_w
   dimidone[0] = dimid_Nkw;   
   if ((retval = nc_def_var(ncid,"z_w",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer edges");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_w)))
   //  ERR(retval);

   //Nk
   dimidone[0] = dimid_Nc;   
   if ((retval = nc_def_var(ncid,"Nk",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at face"); 
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nk)))
   //  ERR(retval);

   //Nke
   dimidone[0] = dimid_Ne;   
   if ((retval = nc_def_var(ncid,"Nke",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at edge");
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nke)))
   //  ERR(retval);

    //dv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"stanford_name","sea_floor_depth_below_geoid");
    nc_addattr(ncid, varid,"long_name","seafloor depth");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dv)))
    //  ERR(retval);

    //time
    dimidone[0] = dimid_time;
    if ((retval = nc_def_var(ncid,"time",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","time");
    nc_addattr(ncid, varid,"units","seconds since 1990-01-01 00:00:00");  

  //  //average time
  //  dimidone[0] = 0;
  //  if ((retval = nc_def_var(ncid,"average_time",NC_DOUBLE,1,dimidone,&varid)))
  //     ERR(retval);
  //  nc_addattr(ncid, varid,"long_name","Averaging time interval");
  //  nc_addattr(ncid, varid,"units","seconds");  


   /********************************************************************** 
    * 
    * Define the physical variables and attributes 
    *
    **********************************************************************/
   
   dimidtwo[0] = dimid_time;
   dimidtwo[1] = dimid_Nc;
   
   dimidthree[0] = dimid_time;
   dimidthree[1] = dimid_Nk;
   dimidthree[2] = dimid_Nc;
   
   // eta
   if ((retval = nc_def_var(ncid,"eta",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Sea surface elevation");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   //u
   if ((retval = nc_def_var(ncid,"uc",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Eastward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //v
   if ((retval = nc_def_var(ncid,"vc",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval);   
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Northward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //w
   dimidthree[1] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"w",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Vertical water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_w yv xv");  

   dimidthree[1] = dimid_Nk;
   
   //nu_v
   if ((retval = nc_def_var(ncid,"nu_v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Vertical eddy viscosity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //salinity
   if(prop->beta>0){
     if ((retval = nc_def_var(ncid,"salt",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Salinity");
    nc_addattr(ncid, varid,"units","ppt");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

    // Depth-integrated salinity
    if ((retval = nc_def_var(ncid,"s_dz",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged depth-integrated salinity");
    nc_addattr(ncid, varid,"units","psu m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   }
   
   //temperature
   if(prop->gamma>0){
     if ((retval = nc_def_var(ncid,"temp",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval); 
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Water temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

     // Depth-integrated temperature
    if ((retval = nc_def_var(ncid,"T_dz",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged depth-integrated temperature");
    nc_addattr(ncid, varid,"units","degrees C m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   }
   
   //rho
   if( (prop->gamma>0) || (prop->beta>0) ){
     if ((retval = nc_def_var(ncid,"rho",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Water density");
    nc_addattr(ncid, varid,"units","kg m-3");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //age
   if(prop->calcage>0){
     if ((retval = nc_def_var(ncid,"agemean",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged mean age");
    nc_addattr(ncid, varid,"units","seconds");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
    
   }

   //U_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"U_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge flux rate");
   nc_addattr(ncid, varid,"units","m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

    //s_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"s_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge salt flux rate");
   nc_addattr(ncid, varid,"units","psu m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

    //T_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"T_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge temperature flux rate");
   nc_addattr(ncid, varid,"units","degreesC m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  


   // Meteorological variables (2-D) //
   
   if(prop->metmodel>0){
    // Uwind
    if ((retval = nc_def_var(ncid,"Uwind",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Vwind
    if ((retval = nc_def_var(ncid,"Vwind",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Tair
    if ((retval = nc_def_var(ncid,"Tair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Air temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Pair
    if ((retval = nc_def_var(ncid,"Pair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Air pressure");
    nc_addattr(ncid, varid,"units","millibar");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // rain
    if ((retval = nc_def_var(ncid,"rain",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Rain fall rate");
    nc_addattr(ncid, varid,"units","kg m2 s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //RH
    if ((retval = nc_def_var(ncid,"RH",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Relative humidity");
    nc_addattr(ncid, varid,"units","Percent (%)");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //cloud
    if ((retval = nc_def_var(ncid,"cloud",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Cloud cover fraction");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
    // Surface flux variables //
    // Hs
    if ((retval = nc_def_var(ncid,"Hs",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Sensible heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hl
    if ((retval = nc_def_var(ncid,"Hl",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Latent heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hlw
    if ((retval = nc_def_var(ncid,"Hlw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Net longwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");      
    nc_addattr(ncid, varid,"positive","down");   

    // Hsw
    if ((retval = nc_def_var(ncid,"Hsw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Net shortwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // tau_x
    if ((retval = nc_def_var(ncid,"tau_x",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    // tau_y
    if ((retval = nc_def_var(ncid,"tau_y",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Northward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    //EP
    if(prop->beta > 0.0){
	if ((retval = nc_def_var(ncid,"EP",NC_DOUBLE,2,dimidtwo,&varid)))
	    ERR(retval); 
	nc_addattr(ncid, varid,"long_name","Time-averaged Evaporation minus precipiaton");
	nc_addattr(ncid, varid,"units","m s-1");
	nc_addattr(ncid, varid,"mesh","suntans_mesh");
	nc_addattr(ncid, varid,"location","face");
	nc_addattr(ncid, varid,"coordinates","time yv xv");   
    }

   }

   //End file definition mode
   if ((retval = nc_enddef(ncid)))
	ERR(retval);

   
   /**********************************************************
   *
   * Write data (needs to be done out of definition mode for classic model)
   *
   ****************************************************************/
   nc_write_int(ncid,"cells",grid->cells,myproc);
   nc_write_int(ncid,"face",grid->face,myproc);
   nc_write_int(ncid,"edges",grid->edges,myproc);
   nc_write_int(ncid,"neigh",grid->neigh,myproc);
   nc_write_int(ncid,"grad",grid->grad,myproc);
   nc_write_int(ncid,"mnptr",grid->mnptr,myproc);
   nc_write_int(ncid,"eptr",grid->eptr,myproc);
   
   nc_write_double(ncid,"xv",grid->xv,myproc);
   nc_write_double(ncid,"yv",grid->yv,myproc);
   nc_write_double(ncid,"xe",grid->xe,myproc);
   nc_write_double(ncid,"ye",grid->ye,myproc);
   nc_write_double(ncid,"xp",grid->xp,myproc);
   nc_write_double(ncid,"yp",grid->yp,myproc);

   nc_write_int(ncid,"normal",grid->normal,myproc);
   nc_write_double(ncid,"n1",grid->n1,myproc);
   nc_write_double(ncid,"n2",grid->n2,myproc);
   nc_write_double(ncid,"df",grid->df,myproc);
   nc_write_double(ncid,"dg",grid->dg,myproc);
   nc_write_double(ncid,"def",grid->def,myproc);
   nc_write_double(ncid,"Ac",grid->Ac,myproc);

   nc_write_double(ncid,"dz",grid->dz,myproc);
   nc_write_double(ncid,"z_r",z_r,myproc);
   nc_write_double(ncid,"z_w",z_w,myproc);
   nc_write_int(ncid,"Nk",grid->Nk,myproc);
   nc_write_int(ncid,"Nke",grid->Nke,myproc);
   nc_write_double(ncid,"dv",grid->dv,myproc);

  // nc_write_double(ncid,"average_time",(float)prop->ntaverage*prop->dt,myproc);

}// End function

/*
* Function: WriteAverageNC()
* -----------------------------
* Main function for writing SUNTANS output to netcdf file/s
* 
*/
void WriteAverageNC(propT *prop, gridT *grid, averageT *average, physT *phys, metT *met, int blowup, MPI_Comm comm, int myproc){
   int ncid = prop->averageNetcdfFileID;
   int varid, retval, k;
   // Start and count vectors for one, two and three dimensional arrays
   const size_t startone[] = {prop->avgtimectr};
   const size_t countone[] = {1};
   const size_t starttwo[] = {prop->avgtimectr,0};
   const size_t counttwo[] = {1,grid->Nc};
   size_t startthree[] = {prop->avgtimectr,0,0};
   size_t countthree[] = {1,grid->Nkmax,grid->Nc};
   const size_t countthreew[] = {1,grid->Nkmax+1,grid->Nc};
   const REAL time[] = {prop->nctime};
   int ntaverage=prop->ntaverage;

   nc_set_log_level(3); // This helps with debugging errors
   
   //REAL *tmpvar, *tmpvarE;
   // Need to write the 3-D arrays as vectors
   //tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"WriteOutputNC");
   //tmpvarE = (REAL *)SunMalloc(grid->Ne*grid->Nkmax*sizeof(REAL),"WriteOutputNC");
   
//   if(!(prop->n%prop->ntaverage) || prop->n==1+prop->nstart || blowup) {
//
    prop->avgctr+=1;
   // Output the first time step but don't compute the average 
   //if(!(prop->n%ntaverage) || prop->n==1+prop->nstart) {
//    if(prop->avgctr==ntaverage || prop->n==1+prop->nstart) {
   if(!(prop->n%ntaverage)) {
     //printf("prop->n/prop->ntaverage=%d\n",prop->n/prop->ntaverage);
     
    //Compute the averages 
    //printf("prop->avgctr=%d\n",prop->avgctr);
     //if(!(prop->n%ntaverage)) 
    ComputeAverageVariables(grid,average,phys,met,prop->avgctr,prop);

    //Communicate the values
    SendRecvAverages(prop,grid,average,comm,myproc); 

    //Reset the counter
    prop->avgctr=0;

    if(myproc==0 && VERBOSE>1){ 
      if(!blowup) 
        printf("Outputting average data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      else
        printf("Outputting blowup averagedata to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
    }
    
    /* Write the time data*/
    if ((retval = nc_inq_varid(ncid, "time", &varid)))
	ERR(retval);
    if ((retval = nc_put_vara_double(ncid, varid, startone, countone, time )))
	ERR(retval);
    
    /* Write to the physical variables*/
    if ((retval = nc_inq_varid(ncid, "eta", &varid)))
	ERR(retval);
    if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->h )))
 	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "uc", &varid)))
	ERR(retval);
    ravel(average->uc, average->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "vc", &varid)))
	ERR(retval);
    ravel(average->vc, average->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	ERR(retval);
      
    // write w at cell top and bottom
    if ((retval = nc_inq_varid(ncid, "w", &varid)))
	ERR(retval);
    ravelW(average->w, average->tmpvarW, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthreew, average->tmpvarW )))
	ERR(retval);

    if ((retval = nc_inq_varid(ncid, "nu_v", &varid)))
	ERR(retval);
    ravel(average->nu_v, average->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	ERR(retval);
    
    // Tracers
     if(prop->beta>0){
       if ((retval = nc_inq_varid(ncid, "salt", &varid)))
	  ERR(retval);
      ravel(average->s, average->tmpvar, grid);
      if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);

	if ((retval = nc_inq_varid(ncid, "s_dz", &varid)))
	    ERR(retval);
	if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->s_dz )))
	    ERR(retval);
     }
     
     if(prop->gamma>0){
	if ((retval = nc_inq_varid(ncid, "temp", &varid)))
	  ERR(retval);
	ravel(average->T, average->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);

	if ((retval = nc_inq_varid(ncid, "T_dz", &varid)))
	    ERR(retval);
	if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->T_dz )))
	    ERR(retval);
     }
      
     if( (prop->gamma>0) || (prop->beta>0) ){ 
	if ((retval = nc_inq_varid(ncid, "rho", &varid)))
	  ERR(retval);
	ravel(average->rho, average->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);
     }

     if(prop->calcage>0){ 
	if ((retval = nc_inq_varid(ncid, "agemean", &varid)))
	  ERR(retval);
	ravel(average->agemean, average->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);

    }

     // Edge fluxes
     countthree[2] = grid->Ne;
     if ((retval = nc_inq_varid(ncid, "U_F", &varid)))
	ERR(retval);
     ravelEdge(average->U_F, average->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvarE )))
        ERR(retval);

     if ((retval = nc_inq_varid(ncid, "s_F", &varid)))
	ERR(retval);
     ravelEdge(average->s_F, average->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvarE )))
        ERR(retval);

     if ((retval = nc_inq_varid(ncid, "T_F", &varid)))
	ERR(retval);
     ravelEdge(average->T_F, average->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvarE )))
        ERR(retval);

     // Wind variables
     if(prop->metmodel>0){
       if ((retval = nc_inq_varid(ncid, "Uwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Uwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Vwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Vwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Tair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Tair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Pair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Pair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "rain", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->rain )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "RH", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->RH )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "cloud", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->cloud )))
	 ERR(retval);
       
       // Heat flux variables
       if ((retval = nc_inq_varid(ncid, "Hs", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hs )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hl", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hl )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hlw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hlw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hsw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hsw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_x", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->tau_x )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_y", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->tau_y )))
	 ERR(retval);
       
       if(prop->beta > 0.0){
	  if ((retval = nc_inq_varid(ncid, "EP", &varid)))
	     ERR(retval);
	  if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->EP )))
	     ERR(retval);
       }
     }
     
    // Zero the arrays after they have been written(don't do it for the initial step)
    if(prop->avgctr>1)
	ZeroAverageVariables(grid,average,prop);

    /* Update the time counter*/
    prop->avgtimectr += 1;  
   }
  
} // End of function



/*
 * Function: nc_addattr()
 * -------------------------
 *
 * Wrapper function to add a text attribute into a netcdf file
 *
 */
static void nc_addattr(int ncid, int varid, char *attname, char *attvalue){
    int retval;
    if ((retval = nc_put_att_text(ncid, varid, attname ,strlen(attvalue),attvalue)))
	ERR(retval);
	
}//End function

/*
 * Function: nc_addattr_int()
 * -------------------------
 *
 * Wrapper function to add a int attribute into a netcdf file
 *
 */
static void nc_addattr_int(int ncid, int varid, char *attname, int *attvalue){
    int retval;
    if ((retval = nc_put_att_int(ncid, varid, attname ,NC_INT,1,attvalue)))
	ERR(retval);
	
}

/*
 * Function: nc_addattr_real()
 * -------------------------
 *
 * Wrapper function to add a double attribute into a netcdf file
 *
 */
static void nc_addattr_real(int ncid, int varid, char *attname, REAL *attvalue){
    int retval;
    if ((retval = nc_put_att_double(ncid, varid, attname ,NC_DOUBLE,1,attvalue)))
	ERR(retval);
	
}
/* 
* Function: ravel()
* -----------------
* Unravel a 2-D SUNTANS array [Nc, Nk] into a vector 
* This is necessary for writing a 2-D array to netcdf as the missing cells need to be filled
*
*/
static void ravel(REAL **tmparray, REAL *tmpvec,gridT *grid){
  int j,k;
  int nk=grid->Nkmax, nc=grid->Nc;
  
  for(j=0;j<nc;j++){
    for(k=0;k<nk;k++){
      if(k<grid->Nk[j]){
        tmpvec[k*nc+j] = tmparray[j][k];
      }else{
	tmpvec[k*nc+j] = (REAL)EMPTY;
      }
    }
  }
}//End of function


/* 
* Function: ravelW()
* -----------------
* Unravel a 2-D SUNTANS array [Nc, Nk+1] into a vector 
* This is necessary for writing a 2-D array to netcdf as the missing cells need to be filled
*
*/
static void ravelW(REAL **tmparray, REAL *tmpvec,gridT *grid){
  int j,k;
  int nk=grid->Nkmax+1, nc=grid->Nc;
  
  for(j=0;j<nc;j++){
    for(k=0;k<nk;k++){
      if(k<grid->Nk[j]+1){
        tmpvec[k*nc+j] = tmparray[j][k];
      }else{
	tmpvec[k*nc+j] = (REAL)EMPTY;
      }
    }
  }
}//End of function

/* 
* Function: ravelEdge()
* -----------------
* Unravel a 2-D SUNTANS array [Ne, Nk] into a vector 
* This is necessary for writing a 2-D array to netcdf as the missing edges need to be filled
*
*/
static void ravelEdge(REAL **tmparray, REAL *tmpvec,gridT *grid){
  int i,k;
  int nk=grid->Nkmax, ne=grid->Ne;
  
  for(i=0;i<ne;i++){
    for(k=0;k<nk;k++){
      if(k<grid->Nke[i]){
        tmpvec[k*ne+i] = tmparray[i][k];
      }else{
	tmpvec[k*ne+i] = (REAL)EMPTY;
      }
    }
  }
}//End of function
const void* FillValue(int empty){
  /* Converts the EMPTY value expression type to match the type expected by nc_def_var_fill*/
  empty = (REAL)empty;
}
 
/*###############################################################
*
* Meteorological input NetCDF functions
*
#################################################################*/

/*
* Function: ReadMetNC()
* ---------------------
* Main function for reading in the meteorological data from the netcdf file
*
*/

void ReadMetNC(propT *prop, gridT *grid, metinT *metin,int myproc){
    int retval, j,k;
    int t0;
    int varid;
    char *vname;
    size_t start[2];
    size_t count[]={1,1};
    int ncid = prop->metncid;

    if(metin->t0==-1){
	metin->t1 = getTimeRec(prop->nctime,metin->time,(int)metin->nt);
	metin->t0 = metin->t1-1;
	metin->t2 = metin->t1+1;
    }
    t0 = metin->t0;
    
    //printf("Model time(0) = %f, time index = %d of %d\n",prop->nctime,t0,metin->nt);
    start[0] = t0;
    start[1] = 0;
    count[0] = NTmet;

    vname = "Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NUwind;
    nc_read_2D(ncid,vname,start,count, metin->Uwind, myproc);

    vname = "Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NVwind;
    nc_read_2D(ncid,vname,start,count, metin->Vwind, myproc);

    vname = "Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NTair;
    nc_read_2D(ncid,vname,start,count, metin->Tair, myproc); 

    vname = "Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NPair;
    nc_read_2D(ncid,vname,start,count, metin->Pair, myproc);

    vname = "rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->Nrain;
    nc_read_2D(ncid,vname,start,count, metin->rain, myproc);

    vname = "RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NRH;
    nc_read_2D(ncid,vname,start,count, metin->RH, myproc);

    vname = "cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->Ncloud;
    nc_read_2D(ncid,vname,start,count, metin->cloud, myproc);
} //End function

/*
* Function: ReadMetNCcoord()
* --------------------------
* Read the coordinate information from the netcdf file 
*
*/
void ReadMetNCcoord(propT *prop, gridT *grid, metinT *metin, int myproc){
  
    /* Read the data from the meteorological netcdf file into the metin structure */
    int retval, j;
    int varid;
    char *vname;
    int ncid = prop->metncid;

    /* Get the horizontal coordintates*/
    vname = "x_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Uwind))) 
      ERR(retval); 
    vname = "y_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Uwind))) 
      ERR(retval); 
    vname = "x_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Vwind))) 
      ERR(retval); 
    vname = "y_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Vwind))) 
      ERR(retval); 
    vname = "x_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Tair))) 
      ERR(retval); 
    vname = "y_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Tair))) 
      ERR(retval); 
    vname = "x_Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Pair))) 
      ERR(retval); 
    vname = "y_Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Pair))) 
      ERR(retval); 
    vname = "x_rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_rain))) 
      ERR(retval); 
    vname = "y_rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_rain))) 
      ERR(retval); 
    vname = "x_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_RH))) 
      ERR(retval); 
    vname = "y_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_RH))) 
      ERR(retval); 
    vname = "x_cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_cloud))) 
      ERR(retval); 
    vname = "y_cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_cloud))) 
      ERR(retval); 
    
    /* Vertical coordinates */
    vname = "z_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_Uwind))) 
      ERR(retval); 
    vname = "z_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_Vwind))) 
      ERR(retval); 
    vname = "z_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_Tair))) 
      ERR(retval); 
    vname = "z_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_RH))) 
      ERR(retval); 
    
    /* Time */
    vname = "Time";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->time))) 
      ERR(retval); 
    
    if(VERBOSE>2 && myproc==0) printf("Finished Reading met netcdf coordinates.\n");

}// End function

/*
* Function: returndimlen()
* ------------------------
* Returns the length of a netcdf dimension 
*/

size_t returndimlen(int ncid, char *dimname){
 int retval;
 int dimid;
 size_t dimlen;
 
 if ((retval =nc_inq_dimid(ncid,dimname,&dimid)))
    ERR(retval);
 
 if ((retval = nc_inq_dimlen(ncid,dimid, &dimlen)))
    ERR(retval);
 return dimlen;
} //End function

/*###############################################################
*
* Boundary contition input NetCDF functions
*
#################################################################*/

/*
 * Function: ReadBdyNC()
 * -----------------------------
 * Reads in boundary netcdf data into the forward and back time steps 
 *
 */     
void ReadBdyNC(propT *prop, gridT *grid, int myproc){
    int retval, j, k, n;
    int t0, t1;
    int varid;
    char *vname;
    size_t start[]={0,0,0};
    size_t start2[]={0,0};
    size_t count[]={0,0,0};
    size_t count2[]={0,0};
    int ncid = prop->netcdfBdyFileID;  
    size_t Nk = bound->Nk;
    size_t Ntype3 = bound->Ntype3;
    size_t Ntype2 = bound->Ntype2;
    size_t Nseg = bound->Nseg;

    //Find the time index of the middle time step (t1) 
    if(bound->t0==-1){
       bound->t1 = getTimeRecBnd(prop->nctime,bound->time,(int)bound->Nt); //this is in met.c
       bound->t0=bound->t1-1;
       bound->t2=bound->t1+1;
       printf("myproc: %d, bound->t0: %d, nctime: %f\n",myproc,bound->t0, prop->nctime);
    }
    t0 = bound->t0;

    count[0]=NT;
    count[1]=Nk;

    count2[0]=NT;

    start[0]=t0;
    start[1]=0;
    start[2]=0;

    start2[0]=t0;
    start2[1]=0;

    //if(myproc==0) printf("t0 = %d [Nt = %d]\n",t0,bound->Nt);    
    if(bound->hasType2){

	count[2]=Ntype2;

	vname = "boundary_u";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->boundary_u_t );

	vname = "boundary_v";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->boundary_v_t );

	vname = "boundary_w";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->boundary_w_t );

	vname = "boundary_T";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->boundary_T_t );

	vname = "boundary_S";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->boundary_S_t );

    }

    if(bound->hasType3){
	count[0]=NT;
	count[1]=Nk;
	count[2]=Ntype3;
	count2[1]=Ntype3;

	vname = "uc";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
        nc_read_3D(ncid, vname, start, count, bound->uc_t );

	vname = "vc";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->vc_t );

	vname = "wc";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->wc_t );

	vname = "T";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->T_t );

	vname = "S";
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_3D(ncid, vname, start, count, bound->S_t);
	//printf("bound->S[1][0][0] = %f\n",bound->S_t[1][0][0]);

	vname = "h";//2D array
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_2D(ncid, vname, start2, count2, bound->h_t, myproc );

     }// End read type-3

     //Flux boundary data
     if(bound->hasType2 && bound->hasSeg){

	count2[1]=Nseg;
	vname = "boundary_Q";//2D array
	if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
	nc_read_2D(ncid, vname, start2, count2, bound->boundary_Q_t, myproc);

     }//End flux read
 }//End function

/*
 * Function: ReadBndNCcoord()
 * --------------------------
 * Reads the coordinate information from the netcdf file into the boundary structure
 *
 */
void ReadBndNCcoord(int ncid, propT *prop, gridT *grid, int myproc){

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

    if(bound->hasSeg>0){
	vname = "segedgep";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->segedgep))) 
	  ERR(retval); 

	vname = "segp";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->segp))) 
	  ERR(retval); 
    }

 }//End function

/*
* Function: returndimlenBC()
* --------------------------
* Returns the length of a dimension 
* Returns a zero if the dimension is not found and does not raise an error
*/
size_t returndimlenBC(int ncid, char *dimname){
 int retval;
 int dimid;
 size_t dimlen;
 
 if ((retval =nc_inq_dimid(ncid,dimname,&dimid)))
    return 0;
 
 if ((retval = nc_inq_dimlen(ncid,dimid, &dimlen)))
    ERR(retval);
 return dimlen;
} // End function

/*###############################################################
*
* Initial contition input NetCDF functions
*/

/*
 * Function: ReadInitialNCcoord()
 * -----------------------------
 * Reads the dimensions from the initial condition netcdf file
 *
 */
 void ReadInitialNCcoord(propT *prop, gridT *grid, int *Nci, int *Nki, int *T0, int myproc){

    int Nt;

   // Read the spatial dimension sizes
    *Nci = (int)returndimlenBC(prop->initialNCfileID,"Nc");
    *Nki = (int)returndimlenBC(prop->initialNCfileID,"Nk");

    // Check the dimension with the grid
    if(*Nki != grid->Nkmax){
	printf("Error! Number of layers in initial condition file (%d) not equal to Nkmax (%d).\n",*Nki,grid->Nkmax); 
	MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Find the index of the closest time point, T0
    Nt = (int)returndimlenBC(prop->initialNCfileID,"time");
    *T0 = getICtime(prop,Nt, myproc);
    if (*T0>=Nt) *T0 = Nt-1;
    //*T0=0;
    return;

 } // End function

/*
 * Function: getICtime()
 * -----------------------------
 * Return the closest time index from the initial condition file
 *
 */
int getICtime(propT *prop, int Nt, int myproc){

   int retval, varid;
   int ncid = prop->initialNCfileID;
   //REAL time[Nt]; 
   REAL *ictime;
   char *vname;

   ictime = (REAL *)SunMalloc(Nt*sizeof(REAL),"getICtime");

   vname = "time";
    if(VERBOSE>2 && myproc==0) printf("Reading initial condition %s...",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid, &ictime[0] ))) 
      ERR(retval); 
    if(VERBOSE>2 && myproc==0) printf("done.\n");

    return getTimeRecBnd(prop->nctime, ictime, (int)Nt);

} // End function

/*
 * Function: ReturnFreeSurfaceNC()
 * -------------------------------
 * Reads the free surface from the initial condition netcdf array
 *
 */
void ReturnFreeSurfaceNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int T0, int myproc){
   int i;
   size_t start[] = {T0, 0};
   size_t count[] = {1,Nci};
   //REAL htmp[Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading free-surface initial condition from netcdf file...\n");
   //printf("Initial condition file: T0 = %d, Nci = %d\n",T0,Nci);
   //nc_read_2D(prop->initialNCfileID, "eta", start, count, htmp , myproc);
    if ((retval = nc_inq_varid(ncid, "eta", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {

     phys->h[i]=htmp[grid->mnptr[i]];
     //phys->h[i]=0;
     if(phys->h[i]<-grid->dv[i] + DRYCELLHEIGHT) 
       phys->h[i]=-grid->dv[i] + DRYCELLHEIGHT;
  }
} // End function


/*
 * Function: ReturnSalinityNC()
 * -------------------------------
 * Reads the salinity from the initial condition netcdf array
 *
 */
void ReturnSalinityNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){
   int i,k,ind;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   //REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading salinity initial condition from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "salt", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //for(k=0;k<grid->Nk[i];k++) {
	 ind = k*Nci + grid->mnptr[i]; 
	 phys->s[i][k]=htmp[ind];
	 phys->s0[i][k]=htmp[ind];
      }
  }
} // End function


/*
 * Function: ReturnTemperatureNC()
 * -------------------------------
 * Reads the salinity from the initial condition netcdf array
 *
 */
void ReturnTemperatureNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){
   int i,k,ind;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   //REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading temperature initial condition from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "temp", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //for(k=0;k<grid->Nk[i];k++) {
	 ind = k*Nci + grid->mnptr[i]; 
	 phys->T[i][k]=htmp[ind];
	
      }
  }
} // End function

/*
 * Function: ReturnAgeNC()
 * -------------------------------
 * Reads the age variables (agec & agealpha) from the initial condition netcdf array
 *
 */
void ReturnAgeNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){
   int i,k,ind;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   //REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading agec initial condition from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "agec", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //for(k=0;k<grid->Nk[i];k++) {
	 ind = k*Nci + grid->mnptr[i]; 
	 phys->agec[i][k]=htmp[ind];
      }
  }

   if(VERBOSE>1 && myproc==0) printf("Reading agealpha initial condition from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "agealpha", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //for(k=0;k<grid->Nk[i];k++) {
	 ind = k*Nci + grid->mnptr[i]; 
	 phys->agealpha[i][k]=htmp[ind];
      }
  }

} // End function


/*
* Function: GetTimeRecBnd()
* ------------------
* Retuns the index of the first preceding time step in the vector time
*/
int getTimeRecBnd(REAL nctime, REAL *time, int nt){
    int j;

    for(j=0;j<nt;j++){
       if (time[j]>=nctime)
	 //return j-1;
	 return j;
    }
    return nt;
}

