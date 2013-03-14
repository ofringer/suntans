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
void ravel(REAL **tmparray, REAL *tmpvec, gridT *grid);


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

void nc_read_3D(int ncid, char *vname, size_t *start, size_t *count, REAL ***tmparray){

    int j, k, n, ii;
    int varid, retval;
    //REAL tmpvec[ (int)count[0] * (int)count[1] * (int)count[2] ];
    REAL outdata[(int)count[0]][(int)count[1]][(int)count[2]];

    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
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

void nc_read_2D(int ncid, char *vname, size_t *start, size_t *count, REAL **tmparray, int myproc){

    int j, n, ii;
    int varid, retval;
    //REAL tmpvec[ (int)count[0] * (int)count[1] ];
    REAL outdata[(int)count[0]][(int)count[1]];

    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
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
  
}// End function

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
   
   REAL *tmpvar;
   // Need to write the 3-D arrays as vectors
   tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"WriteOutputNC");
   
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
    
    if ((retval = nc_inq_varid(ncid, "u", &varid)))
	ERR(retval);
    ravel(phys->uc, tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, tmpvar )))
	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "v", &varid)))
	ERR(retval);
    ravel(phys->vc, tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, tmpvar )))
	ERR(retval);
      
    // write w at cell top and bottom
    if ((retval = nc_inq_varid(ncid, "w", &varid)))
	ERR(retval);
    ravel(phys->w, tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthreew, tmpvar )))
	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "nu_v", &varid)))
	ERR(retval);
    ravel(phys->nu_tv, tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, tmpvar )))
	ERR(retval);
    
    // Tracers
     if(prop->beta>0){
       if ((retval = nc_inq_varid(ncid, "salt", &varid)))
	  ERR(retval);
      ravel(phys->s, tmpvar, grid);
      if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, tmpvar )))
	  ERR(retval);
     }
     
     if(prop->gamma>0){
	if ((retval = nc_inq_varid(ncid, "temp", &varid)))
	  ERR(retval);
	ravel(phys->T, tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, tmpvar )))
	  ERR(retval);
     }
      
     if( (prop->gamma>0) || (prop->beta>0) ){ 
	if ((retval = nc_inq_varid(ncid, "rho", &varid)))
	  ERR(retval);
	ravel(phys->rho, tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, tmpvar )))
	  ERR(retval);
     }
     
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
   SunFree(tmpvar,grid->Nc*grid->Nkmax,"WriteOuputNC");
  
} // End of function


/*
* Function: InitialiseOutputNC()
* ------------------------------------
*
* Initialises the output netcdf file/s
* 
* One file per processor
* The pointer to each file is stored in prop->outputNetcdfFileID
* 
*/
void InitialiseOutputNC(propT *prop, gridT *grid, physT *phys, metT *met, int myproc){
   int ncid = prop->outputNetcdfFileID;
   int retval, i ,j;
   int varid;
   int dimid_Nc, dimid_Nk, dimid_Np, dimid_nt, dimid_numsides, dimid_Nkw; 
   int dimidone[1];
   int dimidtwo[2];
   int dimidthree[3];
   int nofill=0;
   const size_t starttwo[] = {0,0};
   const size_t counttwo[] = {grid->Nkmax,grid->Nc};
   REAL *z_r;
   REAL *z_w;
   const int DEFLATE=1;
   const int DEFLATELEVEL=1;
   
   /********************************************************************** 
    *
    * Define all of the attribute strings manually
    *		-- this is ugly but will do for now... 
    *
    **********************************************************************/
   
   /* Global attributes */ 
   static char globalatt[]="SUNTANS netcdf output file";
   
   /* Grid variable attributes */ 
   static char xv_long_name[]="x coordinate of grid cell centre point";
   static char xv_units[]="metres";
   static char yv_long_name[]="y coordinate of grid cell centre point";
   static char yv_units[]="metres";
   static char xp_long_name[]="x coordinate of grid nodes";
   static char xp_units[]="metres";
   static char yp_long_name[]="y coordinate of grid nodes";
   static char yp_units[]="metres";
   static char dv_long_name[]="depth at grid cell centre point";
   static char dv_units[]="metres";
   static char dz_long_name[]="Vertical grid z-layer spacing";
   static char dz_units[]="metres";
   static char z_r_long_name[]="z coordinate at grid cell mid-height";
   static char z_r_units[]="metres";
   static char z_w_long_name[]="z coordinate at grid cell top/bottom";
   static char z_w_units[]="metres";
   static char dzz_long_name[]="Cell-centred vertical grid spacing";
   static char dzz_units[]="metres";
   static char area_long_name[]="Grid cell area";
   static char area_units[]="metre2";  
   static char Nk_long_name[]="Number of vertical layers at grid cell centre";
   static char cells_long_name[]="Indices to the nodes of each grid cell";
   static char mnptr_long_name[]="Indices to the grid cells in the unpartitioned grid";
   
   /* Time variable attributes 
    *!!! need to read in the units from the input file!!! */ 
   static char time_long_name[]="Simulation time";
   static char time_units[]="seconds since 1990-01-01 00:00:00";
   
   /* Physical variable attributes */ 
   static char eta_long_name[]="Sea surface elevation";
   static char eta_units[]="metre";
   static char u_long_name[]="Eastward water velocity component";
   static char u_units[]="metre second-1";
   static char v_long_name[]="Northward water velocity component";
   static char v_units[]="metre second-1";
   static char w_long_name[]="Vertical water velocity component";
   static char w_units[]="metre second-1";
   static char nu_v_long_name[]="Vertical eddy viscosity";
   static char nu_v_units[]="metre2 second-1";
   static char salt_long_name[]="Salinity";
   static char salt_units[]="ppt";
   static char temp_long_name[]="Water temperature";
   static char temp_units[]="degrees C";
   static char rho_long_name[]="Water density";
   static char rho_units[]="kg m-3";
   
   /* Meteorological variable attribute */
   static char Uwind_long_name[]="Eastward wind velocity component";
   static char Uwind_units[]="metre second-1";
   static char Vwind_long_name[]="Northward wind velocity component";
   static char Vwind_units[]="metre second-1";
   static char Tair_long_name[]="Air temperature";
   static char Tair_units[]="degrees C";
   static char Pair_long_name[]="Air pressure";
   static char Pair_units[]="millibar";
   static char rain_long_name[]="Rain fall rate";
   static char rain_units[]="kg m-2 s-1";
   static char RH_long_name[]="Relative humidity";
   static char RH_units[]="percent (%)";
   static char cloud_long_name[]="Cloud cover fractions";
   static char cloud_units[]="fraction (0-1)";
   
   /* Model calculated air-sea flux parameters */
   static char Hs_long_name[]="Sensible heat flux";
   static char Hs_units[]="Watts metre-2";
   static char Hl_long_name[]="Latent heat flux";
   static char Hl_units[]="Watts metre-2";
   static char Hlw_long_name[]="Longwave radiation";
   static char Hlw_units[]="Watts metre-2";
   static char Hsw_long_name[]="Shortwave radiation";
   static char Hsw_units[]="Watts metre-2";
   static char tau_x_long_name[]="Eastward component surface wind stress";
   static char tau_x_units[]="Newton metre-2";
   static char tau_y_long_name[]="Northward component surface wind stress";
   static char tau_y_units[]="Newton metre-2";
   static char EP_long_name[]="Evaporation minus precipitation";
   static char EP_units[]="meter second-1";
   
   REAL *tmpvar;
   // Need to write the 3-D arrays as vectors
   tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"InitialiseOutputNC");  
   
   /* Initialise the depth arrays */
   z_r = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"InitialiseOutputNC");
   z_w = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNC");
   
   /**************
    *
    * Start writing...
    *
    **************/
   if(VERBOSE>1 && myproc==0) printf("Initialising output netcdf files...\n");
   
   // Set the netcdf time ctr to 0
   prop->nctimectr=0;
   
   /* Define the global attributes */
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "title",strlen(globalatt),globalatt)))
	ERR(retval);
   
   /********************************************************************** 
    *
    * Define the dimensions
    *
    **********************************************************************/
   if ((retval = nc_def_dim(ncid, "Nc", grid->Nc, &dimid_Nc)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Np", grid->Np, &dimid_Np)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nk", grid->Nkmax, &dimid_Nk)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nkw", grid->Nkmax+1, &dimid_Nkw)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "numsides", NUMEDGECOLUMNS, &dimid_numsides)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "nt", NC_UNLIMITED, &dimid_nt)))
	ERR(retval);
   
   /********************************************************************** 
    *
    * Define the grid variables and attributes 
    *
    **********************************************************************/
   
   //xv
   dimidone[0] = dimid_Nc;
   if ((retval = nc_def_var(ncid,"xv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   /*Write some attributes*/
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(xv_long_name),xv_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(xv_units),xv_units)))
	ERR(retval);

    /* Write the data */
   if ((retval = nc_put_var_double(ncid,varid, grid->xv)))
     ERR(retval);
   
   //yv
   dimidone[0] = dimid_Nc;
   if ((retval = nc_def_var(ncid,"yv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(yv_long_name),yv_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(yv_units),yv_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, grid->yv)))
     ERR(retval);
   
   //dv
   dimidone[0] = dimid_Nc;
   if ((retval = nc_def_var(ncid,"dv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(dv_long_name),dv_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(dv_units),dv_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, grid->dv)))
     ERR(retval);
   
   //area
   dimidone[0] = dimid_Nc;
   if ((retval = nc_def_var(ncid,"Ac",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(area_long_name),area_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(area_units),area_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, grid->Ac)))
     ERR(retval);
   //Nk
   dimidone[0] = dimid_Nc;
   if ((retval = nc_def_var(ncid,"Nk",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Nk_long_name),Nk_long_name)))
	ERR(retval);

   if ((retval = nc_put_var_int(ncid,varid, grid->Nk)))
     ERR(retval);
   
   //xp
   dimidone[0] = dimid_Np;   
   if ((retval = nc_def_var(ncid,"xp",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(xp_long_name),xp_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(xp_units),xp_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, grid->xp)))
     ERR(retval);
   
   //yp
   dimidone[0] = dimid_Np;
   if ((retval = nc_def_var(ncid,"yp",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(yp_long_name),yp_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(yp_units),yp_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, grid->yp)))
     ERR(retval);
   
   //dz
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"dz",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(dz_long_name),dz_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(dz_units),dz_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, grid->dz)))
     ERR(retval);
    
   //mnptr
   if(grid->mnptr){//Only write if mnptr exists
       dimidone[0] = dimid_Nc;
       if ((retval = nc_def_var(ncid,"mnptr",NC_INT,1,dimidone,&varid)))
	  ERR(retval);
   
       if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(mnptr_long_name),mnptr_long_name)))
	    ERR(retval);

       if ((retval = nc_put_var_int(ncid,varid, grid->mnptr)))
	 ERR(retval);
   }
   
   //cells
   dimidtwo[0] = dimid_Nc;
   dimidtwo[1] = dimid_numsides;
   if ((retval = nc_def_var(ncid,"cells",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(cells_long_name),cells_long_name)))
	ERR(retval);

   //**** Need to check that these are written in the right order - may need to change the order of the dimensions above ****/
   if ((retval = nc_put_var_int(ncid,varid, grid->cells)))
     ERR(retval);
   
    //time
   dimidone[0] = dimid_nt;
   if ((retval = nc_def_var(ncid,"time",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(time_long_name),time_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(time_units),time_units)))
	ERR(retval);
   
   /* Calculate and write the vertical coordinate data */
   
   z_w[0]=0.0;
   for(j=0;j<grid->Nkmax;j++){
      z_w[j+1] = z_w[j] + grid->dz[j];
      if(j==0){
	 z_r[j] = grid->dz[j]*0.5;
      }else{
	 z_r[j] = z_r[j-1]+grid->dz[j];
      }
   }
   
   //z_r
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"z_r",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(z_r_long_name),z_r_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(z_r_units),z_r_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, z_r)))
     ERR(retval);
   
   //z_w
   dimidone[0] = dimid_Nkw;   
   if ((retval = nc_def_var(ncid,"z_w",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(z_w_long_name),z_w_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(z_w_units),z_w_units)))
	ERR(retval);
   
   if ((retval = nc_put_var_double(ncid,varid, z_w)))
     ERR(retval);
   
    //dzz
   dimidtwo[0] = dimid_Nk;
   dimidtwo[1] = dimid_Nc;
   if ((retval = nc_def_var(ncid,"dzz",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(dzz_long_name),dzz_long_name)))
     ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(dzz_units),dzz_units)))
     ERR(retval);
   
   ravel(grid->dzz, tmpvar, grid);
   if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, tmpvar )))
	ERR(retval);

   /********************************************************************** 
    * 
    * Define the physical variables and attributes 
    *
    **********************************************************************/
   dimidtwo[0] = dimid_nt;
   dimidtwo[1] = dimid_Nc;
   
   dimidthree[0] = dimid_nt;
   dimidthree[1] = dimid_Nk;
   dimidthree[2] = dimid_Nc;
   
   // eta
   if ((retval = nc_def_var(ncid,"eta",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(eta_long_name),eta_long_name)))
	ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(eta_units),eta_units)))
	ERR(retval);
   
   //u
   if ((retval = nc_def_var(ncid,"u",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(u_long_name),u_long_name)))
     ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(u_units),u_units)))
     ERR(retval);
   
   //v
   if ((retval = nc_def_var(ncid,"v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval);   
   if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(v_long_name),v_long_name)))
     ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(v_units),v_units)))
     ERR(retval);
   
   //w
   dimidthree[1] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"w",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(w_long_name),w_long_name)))
     ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(w_units),w_units)))
     ERR(retval);
   
   dimidthree[1] = dimid_Nk;
   
   //nu_v
   if ((retval = nc_def_var(ncid,"nu_v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(nu_v_long_name),nu_v_long_name)))
     ERR(retval);
   if ((retval = nc_put_att_text(ncid, varid, "units",strlen(nu_v_units),nu_v_units)))
     ERR(retval);
   
   //salinity
   if(prop->beta>0){
     if ((retval = nc_def_var(ncid,"salt",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
     if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(salt_long_name),salt_long_name)))
       ERR(retval);
     if ((retval = nc_put_att_text(ncid, varid, "units",strlen(salt_units),salt_units)))
       ERR(retval);
   }
   
   //temperature
   if(prop->gamma>0){
     if ((retval = nc_def_var(ncid,"temp",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval); 
     if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
     if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(temp_long_name),temp_long_name)))
       ERR(retval);
     if ((retval = nc_put_att_text(ncid, varid, "units",strlen(temp_units),temp_units)))
       ERR(retval);     
   }
   
   //rho
   if( (prop->gamma>0) || (prop->beta>0) ){
     if ((retval = nc_def_var(ncid,"rho",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,FillValue(EMPTY)))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
     if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(rho_long_name),rho_long_name)))
       ERR(retval);
     if ((retval = nc_put_att_text(ncid, varid, "units",strlen(rho_units),rho_units)))
       ERR(retval);   
   }
   
   /* Meteorological variables (2-D)*/
   
   if(prop->metmodel>0){
    // Uwind
    if ((retval = nc_def_var(ncid,"Uwind",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Uwind_long_name),Uwind_long_name)))
      ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Uwind_units),Uwind_units)))
      ERR(retval);
    
    // Vwind
      if ((retval = nc_def_var(ncid,"Vwind",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
      if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Vwind_long_name),Vwind_long_name)))
	ERR(retval);
      if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Vwind_units),Vwind_units)))
	ERR(retval);
      
    // Tair
    if ((retval = nc_def_var(ncid,"Tair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Tair_long_name),Tair_long_name)))
	ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Tair_units),Tair_units)))
	ERR(retval);
    
    // Pair
    if ((retval = nc_def_var(ncid,"Pair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Pair_long_name),Pair_long_name)))
	ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Pair_units),Pair_units)))
	ERR(retval);
    
    // rain
    if ((retval = nc_def_var(ncid,"rain",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(rain_long_name),rain_long_name)))
	ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(rain_units),rain_units)))
	ERR(retval);
    
    //RH
    if ((retval = nc_def_var(ncid,"RH",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(RH_long_name),RH_long_name)))
	ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(RH_units),RH_units)))
	ERR(retval);
    
    //cloud
    if ((retval = nc_def_var(ncid,"cloud",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(cloud_long_name),cloud_long_name)))
	ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(cloud_units),cloud_units)))
	ERR(retval);
    
    /* Surface flux variables */
    // Hs
    if ((retval = nc_def_var(ncid,"Hs",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Hs_long_name),Hs_long_name)))
      ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Hs_units),Hs_units)))
      ERR(retval);
    
    // Hl
    if ((retval = nc_def_var(ncid,"Hl",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Hl_long_name),Hl_long_name)))
      ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Hl_units),Hl_units)))
      ERR(retval);
    
    // Hlw
    if ((retval = nc_def_var(ncid,"Hlw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Hlw_long_name),Hlw_long_name)))
      ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Hlw_units),Hlw_units)))
      ERR(retval);
    
    // Hsw
    if ((retval = nc_def_var(ncid,"Hsw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(Hsw_long_name),Hsw_long_name)))
      ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(Hsw_units),Hsw_units)))
      ERR(retval);
    
    // tau_x
    if ((retval = nc_def_var(ncid,"tau_x",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(tau_x_long_name),tau_x_long_name)))
      ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(tau_x_units),tau_x_units)))
      ERR(retval);
    
    // tau_y
    if ((retval = nc_def_var(ncid,"tau_y",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(tau_y_long_name),tau_y_long_name)))
      ERR(retval);
    if ((retval = nc_put_att_text(ncid, varid, "units",strlen(tau_y_units),tau_y_units)))
      ERR(retval);
     
    //EP
    if(prop->beta > 0.0){
	if ((retval = nc_def_var(ncid,"EP",NC_DOUBLE,2,dimidtwo,&varid)))
	    ERR(retval); 
	if ((retval = nc_put_att_text(ncid, varid, "long_name",strlen(EP_long_name),EP_long_name)))
	    ERR(retval);
	if ((retval = nc_put_att_text(ncid, varid, "units",strlen(EP_units),EP_units)))
	    ERR(retval);
    }
   }
   
   /****** 
    *End file definition mode
    ******/
   if ((retval = nc_enddef(ncid)))
	ERR(retval);

   // Free the temporary vector
   SunFree(tmpvar,grid->Nc*grid->Nkmax,"InitialiseOutputNC");


}// End of InitialiseOutputNC

/* 
* Function: ravel()
* -----------------
* Unravel a 2-D SUNTANS array [Nc, Nk] into a vector 
* This is necessary for writing a 2-D array to netcdf as the missing cells need to be filled
*
*/
void ravel(REAL **tmparray, REAL *tmpvec,gridT *grid){
  int j,k;
  int nk=grid->Nkmax, nc=grid->Nc;
  
  for(j=0;j<nc;j++){
    for(k=0;k<nk;k++){
      if(k<grid->Nk[j]){
        tmpvec[k*nc+j] = tmparray[j][k];
      }else{
	tmpvec[k*nc+j] = EMPTY;
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
    
    if(metin->t0==-1){
	metin->t1 = getTimeRec(prop->nctime,metin->time,metin->nt);
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
    nc_read_2D(prop->metncid,vname,start,count, metin->Uwind, myproc);

    vname = "Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NVwind;
    nc_read_2D(prop->metncid,vname,start,count, metin->Vwind, myproc);

    vname = "Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NTair;
    nc_read_2D(prop->metncid,vname,start,count, metin->Tair, myproc); 

    vname = "Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NPair;
    nc_read_2D(prop->metncid,vname,start,count, metin->Pair, myproc);

    vname = "rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->Nrain;
    nc_read_2D(prop->metncid,vname,start,count, metin->rain, myproc);

    vname = "RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NRH;
    nc_read_2D(prop->metncid,vname,start,count, metin->RH, myproc);

    vname = "cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->Ncloud;
    nc_read_2D(prop->metncid,vname,start,count, metin->cloud, myproc);
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

    /* Get the horizontal coordintates*/
    vname = "x_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->x_Uwind))) 
      ERR(retval); 
    vname = "y_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->y_Uwind))) 
      ERR(retval); 
    vname = "x_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->x_Vwind))) 
      ERR(retval); 
    vname = "y_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->y_Vwind))) 
      ERR(retval); 
    vname = "x_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->x_Tair))) 
      ERR(retval); 
    vname = "y_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->y_Tair))) 
      ERR(retval); 
    vname = "x_Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->x_Pair))) 
      ERR(retval); 
    vname = "y_Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->y_Pair))) 
      ERR(retval); 
    vname = "x_rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->x_rain))) 
      ERR(retval); 
    vname = "y_rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->y_rain))) 
      ERR(retval); 
    vname = "x_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->x_RH))) 
      ERR(retval); 
    vname = "y_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->y_RH))) 
      ERR(retval); 
    vname = "x_cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->x_cloud))) 
      ERR(retval); 
    vname = "y_cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->y_cloud))) 
      ERR(retval); 
    
    /* Vertical coordinates */
    vname = "z_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->z_Uwind))) 
      ERR(retval); 
    vname = "z_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->z_Vwind))) 
      ERR(retval); 
    vname = "z_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->z_Tair))) 
      ERR(retval); 
    vname = "z_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->z_RH))) 
      ERR(retval); 
    
    /* Time */
    vname = "Time";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(prop->metncid, varid,metin->time))) 
      ERR(retval); 
    
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
    size_t start[3];
    size_t start2[2];
    size_t count[3];
    size_t count2[2];
    int ncid = prop->netcdfBdyFileID;  
    int Nk = bound->Nk;
    int Ntype3 = bound->Ntype3;
    int Ntype2 = bound->Ntype2;
    int Nseg = bound->Nseg;

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
	nc_read_3D(ncid, vname, start, count, bound->S_t );

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
    Nt = (int)returndimlenBC(prop->initialNCfileID,"Nt");
    *T0 = getICtime(prop,Nt, myproc);
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
   REAL time[Nt]; 
   char *vname;

   vname = "time";
    if(VERBOSE>2 && myproc==0) printf("Reading initial condition %s...",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid, time))) 
      ERR(retval); 
    if(VERBOSE>2 && myproc==0) printf("done.\n");

    return getTimeRecBnd(prop->nctime, time, (int)Nt);

} // End function

/*
 * Function: ReturnFreeSurfaceNC()
 * -------------------------------
 * Reads the free surface from the initial condition netcdf array
 *
 */
void ReturnFreeSurfaceNC(propT *prop, physT *phys, gridT *grid, int Nci, int T0, int myproc){
   int i;
   size_t start[] = {T0, 0};
   size_t count[] = {1,Nci};
   REAL htmp[Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading free-surface initial condition from netcdf file...");
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
void ReturnSalinityNC(propT *prop, physT *phys, gridT *grid, int Nci, int Nki, int T0, int myproc){
   int i,k;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading salinity initial condition from netcdf file...");
    if ((retval = nc_inq_varid(ncid, "S", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0][0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	 phys->s[i][k]=htmp[k][grid->mnptr[i]];
	 phys->s0[i][k]=htmp[k][grid->mnptr[i]];
      }
  }
} // End function


/*
 * Function: ReturnTemperatureNC()
 * -------------------------------
 * Reads the salinity from the initial condition netcdf array
 *
 */
void ReturnTemperatureNC(propT *prop, physT *phys, gridT *grid, int Nci, int Nki, int T0, int myproc){
   int i,k;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading temperature initial condition from netcdf file...");
    if ((retval = nc_inq_varid(ncid, "T", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0][0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	 phys->T[i][k]=htmp[k][grid->mnptr[i]];
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
	 return j-1;
    }
    return nt;
}

