/*
 * Routines for parsing meteorological input data onto the SUNTANS grid
 * 
 */

#include "met.h"
#include "phys.h"

/* Private functions */
#ifdef USENETCDF
void ReadMetNCcoord(propT *prop, gridT *grid, metinT *metin,int myproc);
void ReadMetNC(propT *prop, gridT *grid, metinT *metin,int myproc);
int getTimeRec(REAL nctime, REAL *time, int nt);
const void* FillValue(int empty);
void ravel(REAL **tmparray, REAL *tmpvec, gridT *grid);
static size_t returndimlen(int ncid, char *dimname);
#endif
void calcInterpWeights(gridT *grid, propT *prop, REAL *xo, REAL *yo, int Ns, REAL **klambda,int myproc);
static REAL semivariogram(int varmodel, REAL nugget, REAL sill, REAL range, REAL D);
void weightInterpArray(REAL **D, REAL **klambda, int Nc, int Ns, int nt, REAL **Dout);
void weightInterpField(REAL *D, REAL **klambda, int Nc, int Ns, REAL *Dout);
void linsolve(REAL **A, REAL *b, int N);
static REAL specifichumidity(REAL RH, REAL Ta, REAL Pair);
static REAL qsat(REAL Tw, REAL Pair);
static REAL longwave(REAL Ta, REAL Tw, REAL C_cloud);
static void cor30a(REAL *y);
static REAL psiu_30(REAL zet);
static REAL psit_30(REAL zet);

/* Start of functions */


#ifdef USENETCDF
void WriteOuputNC(propT *prop, gridT *grid, physT *phys, metT *met, int blowup, int myproc){
  /*********************************************************************** 
   * 
   * Writes output to netcdf file/s
   * 
   ************************************************************************/
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
  
} // End of 

void InitialiseOutputNC(propT *prop, gridT *grid, physT *phys, metT *met, int myproc){
  /*********************************************************************** 
   * 
   * Initialises the output netcdf file/s
   * 
   * One file per processor
   * The pointer to each file is stored in prop->outputNetcdfFileID
   * 
   ************************************************************************/
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
   static char salt_long_name[]="Water salinity";
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

const void* FillValue(int empty){
  /* Converts the EMPTY value expression type to match the type expected by nc_def_var_fill*/
 empty = (REAL)empty;
}

void InitialiseMetFields(propT *prop, gridT *grid, metinT *metin, metT *met, int myproc){
  /* Wrapper function for initialising all meterological fields */
  int retval;
  int i,j;
  

 /*  Read in the coordinate data*/
 if(VERBOSE>3 && myproc==0) printf("Reading netcdf coordinate data...\n");
 ReadMetNCcoord(prop,grid,metin, myproc);
 
 /* Calculating the interpolation weights for each variable*/
 calcInterpWeights(grid,prop, metin->x_Uwind, metin->y_Uwind, metin->NUwind, metin->WUwind, myproc);
 calcInterpWeights(grid,prop, metin->x_Vwind, metin->y_Vwind, metin->NVwind, metin->WVwind, myproc);
 calcInterpWeights(grid,prop, metin->x_Tair, metin->y_Tair, metin->NTair, metin->WTair, myproc);
 calcInterpWeights(grid,prop, metin->x_Pair, metin->y_Pair, metin->NPair, metin->WPair, myproc);
 calcInterpWeights(grid,prop, metin->x_rain, metin->y_rain, metin->Nrain, metin->Wrain, myproc);
 calcInterpWeights(grid,prop, metin->x_RH, metin->y_RH, metin->NRH, metin->WRH, myproc);
 calcInterpWeights(grid,prop, metin->x_cloud, metin->y_cloud, metin->Ncloud, metin->Wcloud, myproc);
 
 if(VERBOSE>3 && myproc==0){
    printf("Uwind weights:\n");
    for(i=0;i<grid->Nc;i++){
	    printf("xv=%f, yv=%f, Weights: ",grid->xv[i],grid->yv[i]);
	    for(j=0;j<metin->NUwind;j++){
	      printf("%f, ",metin->WUwind[i][j]);
	    }
	    printf("\n");
    }
 }

 
 /*  Interpolate the heights of some variables */
 if(VERBOSE>1 && myproc==0) printf("Interpolating height coordinates onto grid...\n");
 weightInterpField(metin->z_Uwind, metin->WUwind, grid->Nc, metin->NUwind, met->z_Uwind);
 weightInterpField(metin->z_Vwind, metin->WVwind, grid->Nc, metin->NVwind, met->z_Vwind);
 weightInterpField(metin->z_Tair, metin->WTair, grid->Nc, metin->NTair, met->z_Tair);
 weightInterpField(metin->z_RH, metin->WRH, grid->Nc, metin->NRH, met->z_RH);
 

  
} // End of InitialiseMetFields

void ravel(REAL **tmparray, REAL *tmpvec,gridT *grid){
  /* 
   * Unravel a 2-D array into a vector 
   * This is necessary for writing a 2-D array to netcdf
   */
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
//   for(k=0;k<nk;k++){
//     for(j=0;j<nc;j++){
//       tmpvec[j*k+j] = tmparray[j][k]; 
//     }
//   }
//   
}

void ReadMetNC(propT *prop, gridT *grid, metinT *metin,int myproc){
    /* Read the met data from a netcddf file for the initial two time steps*/ 
    int retval, j,k;
    int t0;
    int varid;
    char *vname;
    size_t start[2];
    size_t count[]={1,1};
    ptrdiff_t stride[]={1,1};
    
    t0 = getTimeRec(prop->nctime,metin->time,metin->nt);
    
    //printf("Model time(0) = %f, time index = %d of %d\n",prop->nctime,t0,metin->nt);
    
    vname = "Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    for (j=0;j<metin->NUwind;j++){
      for (k=0;k<2;k++){
	start[0]=t0+k;
	start[1]=j;
	if ((retval = nc_get_vara_double(prop->metncid, varid, start, count, &metin->Uwind[k][j]))) 
	    ERR(retval); 
	if(VERBOSE>3 && myproc==0) printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->Uwind[k][j]); 
      }
    }
    
    vname = "Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    for (j=0;j<metin->NVwind;j++){
      for (k=0;k<2;k++){
	start[0]=t0+k;
	start[1]=j;
	if ((retval = nc_get_vara_double(prop->metncid, varid, start, count, &metin->Vwind[k][j]))) 
	    ERR(retval); 
	 if(VERBOSE>3 && myproc==0) printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->Vwind[k][j]); 
      }
    }
    
    vname = "Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    for (j=0;j<metin->NTair;j++){
      for (k=0;k<2;k++){
	start[0]=t0+k;
	start[1]=j;
	if ((retval = nc_get_vara_double(prop->metncid, varid, start, count, &metin->Tair[k][j]))) 
	    ERR(retval); 
	if(VERBOSE>3 && myproc==0) printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->Tair[k][j]); 
      }
    }
    
    vname = "Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    for (j=0;j<metin->NPair;j++){
      for (k=0;k<2;k++){
	start[0]=t0+k;
	start[1]=j;
	if ((retval = nc_get_vara_double(prop->metncid, varid, start, count, &metin->Pair[k][j]))) 
	    ERR(retval); 
	if(VERBOSE>3 && myproc==0) printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->Pair[k][j]); 
      }
    }
    
    vname = "rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    for (j=0;j<metin->Nrain;j++){
      for (k=0;k<2;k++){
	start[0]=t0+k;
	start[1]=j;
	if ((retval = nc_get_vara_double(prop->metncid, varid, start, count, &metin->rain[k][j]))) 
	    ERR(retval); 
	if(VERBOSE>3 && myproc==0) printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->rain[k][j]); 
      }
    }
    
    vname = "RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    for (j=0;j<metin->NRH;j++){
      for (k=0;k<2;k++){
	start[0]=t0+k;
	start[1]=j;
	if ((retval = nc_get_vara_double(prop->metncid, varid, start, count, &metin->RH[k][j]))) 
	    ERR(retval); 
	if(VERBOSE>3 && myproc==0) printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->RH[k][j]); 
      }
    }
    
    vname = "cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
	ERR(retval);
    for (j=0;j<metin->Ncloud;j++){
      for (k=0;k<2;k++){
	start[0]=t0+k;
	start[1]=j;
	if ((retval = nc_get_vara_double(prop->metncid, varid, start, count, &metin->cloud[k][j]))) 
	    ERR(retval); 
	if(VERBOSE>3 && myproc==0) printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->cloud[k][j]); 
      }
    }
//     count[0]=2;
//     count[1]=metin->NRH;
//     start[0]=t0;
//     start[1]=0;
//     if ((retval = nc_inq_varid(prop->metncid, vname, &varid)))
// 	ERR(retval);
//     if ((retval = nc_get_vars_double(prop->metncid, varid, start, count,stride, metin->RH))) 
// 	ERR(retval); 
//     for (j=0;j<metin->NRH;j++){
//        for (k=0;k<2;k++){
// 	  printf("%s[%d][%d] = %10.6f .\n",vname,j,k,metin->RH[k][j]); 
//        }
//     }
    
}

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
    
}


static size_t returndimlen(int ncid, char *dimname){
 /* Returns the length of a dimension */
 int retval;
 int dimid;
 size_t dimlen;
 
 if ((retval =nc_inq_dimid(ncid,dimname,&dimid)))
    ERR(retval);
 
 if ((retval = nc_inq_dimlen(ncid,dimid, &dimlen)))
    ERR(retval);
 return dimlen;
}
#endif

void updateMetData(propT *prop, gridT *grid, metinT *metin, metT *met, int myproc){
  
  /* Main function for updating the met structure and interpolating onto the model time step */
  
  int j, t0; 
  REAL dt, r1, r2;
   
  t0 = getTimeRec(prop->nctime,metin->time,metin->nt);
    
    /* Only interpolate the data onto the grid if need to*/
    if (metin->t0!=t0){
      if(VERBOSE>3 && myproc==0) printf("Updating netcdf variable at nc timestep: %d\n",t0);
      /* Read in the data two time steps*/
#ifdef USENETCDF	    
      ReadMetNC(prop, grid, metin, myproc);
#endif
      metin->t0=t0;
      metin->t1=t0+1;
      
      /* Interpolate the two time steps onto the grid*/
      weightInterpArray(metin->Uwind, metin->WUwind, grid->Nc, metin->NUwind, 2, met->Uwind_t);
      weightInterpArray(metin->Vwind, metin->WVwind, grid->Nc, metin->NVwind, 2, met->Vwind_t);
      weightInterpArray(metin->Tair, metin->WTair, grid->Nc, metin->NTair, 2, met->Tair_t);
      weightInterpArray(metin->Pair, metin->WPair, grid->Nc, metin->NPair, 2, met->Pair_t);
      weightInterpArray(metin->rain, metin->Wrain, grid->Nc, metin->Nrain, 2, met->rain_t);
      weightInterpArray(metin->RH, metin->WRH, grid->Nc, metin->NRH, 2, met->RH_t);
      weightInterpArray(metin->cloud, metin->Wcloud, grid->Nc, metin->Ncloud, 2, met->cloud_t);
    }
    
    /* Do a linear temporal interpolation */
    dt = metin->time[metin->t1]-metin->time[metin->t0];
    r2 = (prop->nctime - metin->time[metin->t0])/dt;
    r1 = 1.0-r2;
    
     //printf("tmod = %f, tlow = %f (r1=%f), thigh = %f (r2=%f)\n",prop->nctime, metin->time[metin->t0],r1,metin->time[metin->t1],r2);
    
    for (j=0;j<grid->Nc;j++){
      met->Uwind[j] = met->Uwind_t[0][j]*r1 + met->Uwind_t[1][j]*r2;
      met->Vwind[j] = met->Vwind_t[0][j]*r1 + met->Vwind_t[1][j]*r2;
      met->Tair[j] = met->Tair_t[0][j]*r1 + met->Tair_t[1][j]*r2;
      met->Pair[j] = met->Pair_t[0][j]*r1 + met->Pair_t[1][j]*r2;
      met->rain[j] = met->rain_t[0][j]*r1 + met->rain_t[1][j]*r2;
      met->RH[j] = met->RH_t[0][j]*r1 + met->RH_t[1][j]*r2;
      met->cloud[j] = met->cloud_t[0][j]*r1 + met->cloud_t[1][j]*r2;
      
      /* Place bounds on rain, humidity and cloud variables */
       if (met->cloud[j]<0.0) 
	 met->cloud[j]=0.0;
       if (met->cloud[j]>1.0)
	 met->cloud[j]=1.0;
       if (met->rain[j]<0.0) 
	 met->rain[j]=0.0;
       if (met->RH[j]<0.0)
	 met->RH[j]=0.0;
       if (met->RH[j]>100.0)
	 met->RH[j]=100.0;
    }
} // End of updateMetData



void AllocateMet(propT *prop, gridT *grid, metT **met , int myproc){
  /* Allocates memory to the meteorological structure array on the SUNTANS grid points*/
  int j, k;
  int Nc = grid->Nc;
  int nt = 2;
  
  if(VERBOSE>3 && myproc==0) printf("Allocating met structure...\n");
  *met = (metT *)SunMalloc(sizeof(metT),"AllocateMet");
  
  (*met)->z_Uwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->z_Vwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->z_Tair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->z_RH = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  
  (*met)->Uwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Vwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Tair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Pair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->rain = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->RH = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->cloud = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  
  (*met)->Hs = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Hl = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Hlw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Hsw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->tau_x = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->tau_y = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->ustar = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Tstar = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->qstar = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->EP = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Htmp = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->xtmp = (REAL *)SunMalloc(19*sizeof(REAL),"updateAirSeaFluxes"); 

  (*met)->Uwind_t = (REAL **)SunMalloc(nt*sizeof(REAL *),"AllocateMet");
  (*met)->Vwind_t = (REAL **)SunMalloc(nt*sizeof(REAL *),"AllocateMet");
  (*met)->Tair_t = (REAL **)SunMalloc(nt*sizeof(REAL *),"AllocateMet");
  (*met)->Pair_t = (REAL **)SunMalloc(nt*sizeof(REAL *),"AllocateMet");
  (*met)->rain_t = (REAL **)SunMalloc(nt*sizeof(REAL *),"AllocateMet");
  (*met)->RH_t = (REAL **)SunMalloc(nt*sizeof(REAL *),"AllocateMet");
  (*met)->cloud_t = (REAL **)SunMalloc(nt*sizeof(REAL *),"AllocateMet");
  for(j=0;j<nt;j++){
      (*met)->Uwind_t[j] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");   
      (*met)->Vwind_t[j] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->Tair_t[j] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->Pair_t[j] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->rain_t[j] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->RH_t[j] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->cloud_t[j] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  }
  
  for(k=0;k<Nc;k++){
      (*met)->z_Uwind[k] = 0.0;
      (*met)->z_Vwind[k] = 0.0;
      (*met)->z_Tair[k] = 0.0;
      (*met)->z_RH[k] = 0.0;
      
      (*met)->Uwind[k] = 0.0;
      (*met)->Vwind[k] = 0.0;
      (*met)->Tair[k] = 0.0;
      (*met)->Pair[k] = 0.0;
      (*met)->rain[k] = 0.0;
      (*met)->RH[k] = 0.0;
      (*met)->cloud[k] = 0.0;
      
      (*met)->Hs[k] = 0.0;
      (*met)->Hl[k] = 0.0;
      (*met)->Hlw[k] = 0.0;
      (*met)->Hsw[k] = 0.0;
      (*met)->tau_x[k] = 0.0;
      (*met)->tau_y[k] = 0.0;
      (*met)->ustar[k] = 0.0;
      (*met)->Tstar[k] = 0.0;
      (*met)->qstar[k] = 0.0;
      (*met)->EP[k] = 0.0;
      (*met)->Htmp[k] = 0.0;
      
      for(j=0;j<nt;j++){
	  (*met)->Uwind_t[j][k] = 0.0;
	  (*met)->Vwind_t[j][k] = 0.0;
	  (*met)->Tair_t[j][k] = 0.0;
	  (*met)->Pair_t[j][k] = 0.0;
	  (*met)->rain_t[j][k] = 0.0;
	  (*met)->RH_t[j][k] = 0.0;
	  (*met)->cloud_t[j][k] = 0.0;
      }  
  }
  for(j=0;j<19;j++){
    (*met)->xtmp[j]=0.0;
  }
} // End of function
  
void AllocateMetIn(propT *prop, gridT *grid, metinT **metin, int myproc){
  /* Allocates memory to the meteorological input structure array*/
  int j, k, retval;
  size_t NUwind;
  size_t NVwind;
  size_t NTair;
  size_t NPair;
  size_t Nrain;
  size_t NRH;
  size_t Ncloud;
  size_t nt;
  int Nc = grid->Nc;
  
  if(VERBOSE>3 && myproc==0) printf("Allocating metin structure...\n");
  *metin = (metinT *)SunMalloc(sizeof(metinT),"AllocateMetIn");
  
  
  /* Scalars */
#ifdef USENETCDF  
  (*metin)->NUwind = returndimlen(prop->metncid,"NUwind");
  (*metin)->NVwind = returndimlen(prop->metncid,"NVwind");
  (*metin)->NTair = returndimlen(prop->metncid,"NTair");
  (*metin)->NPair = returndimlen(prop->metncid,"NPair");
  (*metin)->Nrain = returndimlen(prop->metncid,"Nrain");
  (*metin)->NRH = returndimlen(prop->metncid,"NRH");
  (*metin)->Ncloud = returndimlen(prop->metncid,"Ncloud");
  (*metin)->nt = returndimlen(prop->metncid,"nt");
#else
  (*metin)->NUwind = 1;
  (*metin)->NVwind = 1;
  (*metin)->NTair = 1;
  (*metin)->NPair = 1;
  (*metin)->Nrain = 1;
  (*metin)->NRH = 1;
  (*metin)->Ncloud = 1;
  (*metin)->nt = 1;
#endif
  (*metin)->t0 = -999;
  (*metin)->t1 = -999;
  
 
  
  NUwind = (*metin)->NUwind;
  NVwind = (*metin)->NVwind;
  NTair = (*metin)->NTair;
  NPair = (*metin)->NPair;
  Nrain = (*metin)->Nrain;
  NRH = (*metin)->NRH;
  Ncloud = (*metin)->Ncloud;
  nt = (*metin)->nt;
  
  if (VERBOSE>3){
      printf("NUwind = %d\n",NUwind);
      printf("NVwind = %d\n",NVwind); 
      printf("NTair = %d\n",NTair);
      printf("NPair = %d\n",NPair);
      printf("Nrain = %d\n",Nrain);
      printf("NRH = %d\n",NRH);
      printf("Ncloud = %d\n",Ncloud);
      printf("nt = %d\n",nt);
      printf("Nc = %d\n",Nc);
  }
  /* Allocate the coordinate vectors*/
  printf("Allocating coordinates...\n");
  (*metin)->x_Uwind = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_Vwind = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_Tair = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_Pair = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_rain = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_RH = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_cloud = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");
  
  (*metin)->y_Uwind = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_Vwind = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_Tair = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_Pair = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_rain = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_RH = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_cloud = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");

  (*metin)->z_Uwind = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->z_Vwind = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->z_Tair = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
  (*metin)->z_RH = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn"); 
  
  (*metin)->time = (REAL *)SunMalloc(nt*sizeof(REAL),"AllocateMetIn");
  
  /* Allocate the 2-D grid based weights */
  (*metin)->WUwind = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WVwind = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WTair = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WPair = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Wrain = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WRH = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Wcloud = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  for(j=0;j<Nc;j++){
      (*metin)->WUwind[j] = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
      (*metin)->WVwind[j] = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
      (*metin)->WTair[j] = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
      (*metin)->WPair[j] = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
      (*metin)->Wrain[j] = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
      (*metin)->WRH[j] = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
      (*metin)->Wcloud[j] = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");
  }
  
  /* Allocate the 2-D variable data (2 time steps)*/
  (*metin)->Uwind = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Vwind = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Tair = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Pair = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocateMetIn");
  (*metin)->rain = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocateMetIn");
  (*metin)->RH = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocateMetIn");
  (*metin)->cloud = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocateMetIn");
  for(j=0;j<2;j++){
      (*metin)->Uwind[j] = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");   
      (*metin)->Vwind[j] = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
      (*metin)->Tair[j] = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
      (*metin)->Pair[j] = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
      (*metin)->rain[j] = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
      (*metin)->RH[j] = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
      (*metin)->cloud[j] = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");
  }
//   (*metin)->Uwind = (REAL **)SunMalloc(NUwind*sizeof(REAL *),"AllocateMetIn");
//   (*metin)->Vwind = (REAL **)SunMalloc(NVwind*sizeof(REAL *),"AllocateMetIn");
//   (*metin)->Tair = (REAL **)SunMalloc(NTair*sizeof(REAL *),"AllocateMetIn");
//   (*metin)->Pair = (REAL **)SunMalloc(NPair*sizeof(REAL *),"AllocateMetIn");
//   (*metin)->rain = (REAL **)SunMalloc(Nrain*sizeof(REAL *),"AllocateMetIn");
//   (*metin)->RH = (REAL **)SunMalloc(NRH*sizeof(REAL *),"AllocateMetIn");
//   (*metin)->cloud = (REAL **)SunMalloc(Ncloud*sizeof(REAL *),"AllocateMetIn");
//   for(j=0;j<NUwind;j++){
//       (*metin)->Uwind[j] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocateMetIn");}
//   for(j=0;j<NVwind;j++){    
//       (*metin)->Vwind[j] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocateMetIn");}
//   for(j=0;j<NTair;j++){
//       (*metin)->Tair[j] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocateMetIn");}
//   for(j=0;j<NPair;j++){
//       (*metin)->Pair[j] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocateMetIn");}
//   for(j=0;j<Nrain;j++){
//       (*metin)->rain[j] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocateMetIn");}
//   for(j=0;j<NRH;j++){
//       (*metin)->RH[j] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocateMetIn");}
//   for(j=0;j<Ncloud;j++){
//       (*metin)->cloud[j] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocateMetIn");}
  
  /* Initialises all of the input meteorological arrays with zeros*/ 
  // Need to allocate variable by variable as the lengths are different
  if(VERBOSE>2 && myproc==0) printf("Uwind, nj = %d, Nc = %d...\n",(*metin)->NUwind,grid->Nc);
  for(j=0;j<(*metin)->NUwind;j++){
      (*metin)->x_Uwind[j]=0.0;
      (*metin)->y_Uwind[j]=0.0;
      (*metin)->z_Uwind[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WUwind[k][j]=0.0;
      }
       for(k=0;k<2;k++){
	 (*metin)->Uwind[k][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("Vwind, nj = %d, Nc = %d...\n",(*metin)->NVwind,grid->Nc);
  for(j=0;j<(*metin)->NVwind;j++){
      (*metin)->x_Vwind[j]=0.0;
      (*metin)->y_Vwind[j]=0.0;
      (*metin)->z_Vwind[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WVwind[k][j]=0.0;
      }
       for(k=0;k<2;k++){
	 (*metin)->Vwind[k][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("Tair, nj = %d, Nc = %d...\n",(*metin)->NTair,grid->Nc);
  for(j=0;j<(*metin)->NTair;j++){
      (*metin)->x_Tair[j]=0.0;
      (*metin)->y_Tair[j]=0.0;
      (*metin)->z_Tair[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WTair[k][j]=0.0;
      }
       for(k=0;k<2;k++){
	 (*metin)->Tair[k][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("Pair, nj = %d, Nc = %d...\n",(*metin)->NPair,grid->Nc);
  for(j=0;j<(*metin)->NPair;j++){
      (*metin)->x_Pair[j]=0.0;
      (*metin)->y_Pair[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WPair[k][j]=0.0;
      }
       for(k=0;k<2;k++){
	 (*metin)->Pair[k][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("rain, nj = %d, Nc = %d...\n",(*metin)->Nrain,grid->Nc);
  for(j=0;j<(*metin)->Nrain;j++){
      (*metin)->x_rain[j]=0.0;
      (*metin)->y_rain[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->Wrain[k][j]=0.0;
      }
       for(k=0;k<2;k++){
	 (*metin)->rain[k][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("RH, nj = %d, Nc = %d...\n",(*metin)->NRH,grid->Nc);
  for(j=0;j<(*metin)->NRH;j++){
      (*metin)->x_RH[j]=0.0;
      (*metin)->y_RH[j]=0.0;
      (*metin)->z_RH[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WRH[k][j]=0.0;
      }
       for(k=0;k<2;k++){
	 (*metin)->RH[k][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("cloud, nj = %d, Nc = %d...\n",(*metin)->Ncloud,grid->Nc);
  for(j=0;j<(*metin)->Ncloud;j++){
      (*metin)->x_cloud[j]=0.0;
      (*metin)->y_cloud[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->Wcloud[k][j]=0.0;
      }
       for(k=0;k<2;k++){
	 (*metin)->cloud[k][j]=0.0;
      }
  }
  if(VERBOSE>2 && myproc==0) printf("time, nt = %d ...\n",(*metin)->nt);
  for(j=0;j<nt;j++){
      (*metin)->time[j]=0.0;
  }
}


int getTimeRec(REAL nctime, REAL *time, int nt){
   /* Retuns the index of the first preceding time step in the vector time*/
   int j;
   
   for(j=0;j<nt;j++){
      if (time[j]>=nctime)
	return j-1;
   }
   return nt;
}

REAL getToffSet(char *basetime, char *starttime){
    /* Returns the time offset in days between two time strings - starttime and basetime
     * 
     * The time string format is: yyyymmdd.HHMMSS (15 characters)
     * Uses the time.h libraries
     */
	
    //char *strptime(const char *buf, const char *format, struct tm *tm) 
    //time_t mktime ( struct tm * timeptr ); 
    struct tm tm0; 
    struct tm tm1;
    time_t t0, t1;

    strptime(basetime,"%Y%m%d.%H%M%S",&tm0);
    strptime(starttime,"%Y%m%d.%H%M%S",&tm1);
    
    t0 = mktime(&tm0);
    t1 = mktime(&tm1);
    return difftime(t1,t0)/86400.0;

}
/* These functions can be moved to util.c later*/
void calcInterpWeights(gridT *grid, propT *prop, REAL *xo, REAL *yo, int Ns, REAL **klambda, int myproc){
    /* Calculates the interpolation weights for all grid points based on "Ns" interpolants
     * at cooridinates (xo, yo)
     */
    
    int j,i, jj,ii;
    int Nc = grid->Nc;
    REAL sumgamma, dist, tmp;
    REAL *gamma;
    REAL **C, **Ctmp;
    // Allocate the arrays
    if(prop->varmodel==0){
      gamma = (REAL *)SunMalloc(Ns*sizeof(REAL),"CalcInterpWeights");
    }else{
	gamma = (REAL *)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	C = (REAL **)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	Ctmp = (REAL **)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	for (j=0;j<Ns+1;j++){
	  C[j] = (REAL *)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	  Ctmp[j] = (REAL *)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	}
    }
    
    if(prop->varmodel==0){ // Inverse distance weighting
      if(VERBOSE>1 && myproc==0) printf("Calculating interpolation weights using inverse distance weighting...\n");     
      for(i=0;i<Nc;i++){
	sumgamma=0.0;
	for(j=0;j<Ns;j++){
	    dist = pow(grid->xv[i]-xo[j],2) + pow(grid->yv[i]-yo[j],2);
	    gamma[j] = 1.0/dist;
	    sumgamma += gamma[j];
	    //printf("dist = %f, sum = %f\n",dist,sumgamma);
	}
	for(j=0;j<Ns;j++){
	    klambda[i][j] = gamma[j]/sumgamma;
	    //printf("weight = %f\n",klambda[i][j]);
	}
      }
    
    SunFree(gamma,Ns,"CalcInterpWeights");
    
    }else{ // kriging
	if(VERBOSE>1 && myproc==0)  printf("Calculating interpolation weights using kriging...\n");
	
	// Construct the LHS Matrix C
	for(i=0;i<Ns+1;i++){
	  for(j=0;j<Ns+1;j++){
	    C[i][j]=1.0;
	  }
	}
	for(i=0;i<Ns;i++){
	  //C[i][i]=0.0;
	  C[i][i]=semivariogram(prop->varmodel, prop->nugget, prop->sill, prop->range, 0.0);
	  for(j=i+1;j<Ns;j++){
	    dist = sqrt(  pow(xo[i]-xo[j],2) + pow(yo[i]-yo[j],2) );
	    C[i][j] = semivariogram(prop->varmodel, prop->nugget, prop->sill, prop->range, dist);
	    C[j][i]=C[i][j];
	  }
	}
	C[Ns][Ns]=0.0;
	
// 	printf("C[i][j]:\n");
// 	for(i=0;i<Ns+1;i++){
// 	  for(j=0;j<Ns+1;j++){
// 	    printf("%1.6f ",C[i][j]);
// 	  }
// 	  printf("\n");
// 	}
	
	// Loop through each model grid point and calculate the  weights
	for(i=0;i<Nc;i++){
	  for(j=0;j<Ns;j++){
	    dist = sqrt( pow(grid->xv[i]-xo[j],2) + pow(grid->yv[i]-yo[j],2) );
	    gamma[j] = semivariogram(prop->varmodel, prop->nugget, prop->sill, prop->range, dist);
	  }
	  gamma[Ns]=1.0;
	  
	  // Solve the linear system
	  for(ii=0;ii<Ns+1;ii++){
	    for(jj=0;jj<Ns+1;jj++){
	      Ctmp[ii][jj]=C[ii][jj];
	    }
	  }
	  linsolve(Ctmp,gamma,Ns+1);
	  
	  // Write to the weights array
	  for(jj=0;jj<Ns;jj++){
	    klambda[i][jj] = gamma[jj];
	  }
	  
	  // Check the weights
	  if(VERBOSE>3 && myproc==0) {
	    sumgamma=0.0;
	    printf("W[j]:\n");
	    for(jj=0;jj<Ns;jj++){ // don't include the last point
	      sumgamma+=gamma[jj];
	      printf("%3.6f ",gamma[jj]);
	    }
	    printf("\nSumW = %f (should equal 1.000)\n",sumgamma);
	  }
	} 
	
      // Free up the arrays
      SunFree(gamma,Ns+1,"CalcInterpWeights");
      for (j=0;j<Ns+1;j++){
	SunFree(C[j],Ns+1,"CalcInterpWeights");
	SunFree(Ctmp[j],Ns+1,"CalcInterpWeights");
      }
    }// end of kriging   
} // End of calcInterpWeights

static REAL semivariogram(int varmodel, REAL nugget, REAL sill, REAL range, REAL D){
  
  REAL tmp;
 // Calculates the semivariogram function 
 if (varmodel==1){
   // Spherical model
   if(D > range){
     return sill;
   }else{
      tmp = D/range;
      return (nugget + (sill - nugget) * (1.5*tmp - 0.5 * pow(tmp,3))); 
   }
   return 0.0;
 }  
}

void weightInterpArray(REAL **D, REAL **klambda, int Nc, int Ns, int nt, REAL **Dout){
  /* Perform weighted interpolation on a field D [2d-array]*/
  int i,j, k;
  for(k=0;k<nt;k++){
    for(i=0;i<Nc;i++){
	Dout[k][i] = 0.0;
	for(j=0;j<Ns;j++){
	    Dout[k][i] += klambda[i][j] * D[k][j];
	}
    }
  }
}

void weightInterpField(REAL *D, REAL **klambda, int Nc, int Ns, REAL *Dout){
  /* Perform weighted interpolation on a field D [vector]*/
  int i,j;
  
  for(i=0;i<Nc;i++){
      Dout[i] = 0.0;
      for(j=0;j<Ns;j++){
	  Dout[i] += klambda[i][j] * D[j];
      }
  }
}

void linsolve(REAL **A, REAL *b, int N){
  
  /* Solves a linear system of equations A.x=b
   * 
   * A is a square matrix and b is a vector.
   * The solution x is returned in vector b
   * 
   * Reference:
   * 	Golub and Van Loan, "Matrix Computations", 1999, Ch 3
   */

  int i,j,k;
  REAL sumi;
  
  // Algorithm to find LU decomp - See page 99
  for(k=0;k<N-1;k++){
    for(i=k+1;i<N;i++){
      A[i][k] = A[i][k]/A[k][k];
      for(j=k+1;j<N;j++){
	A[i][j] = A[i][j] - A[k][j]*A[i][k];
      }
    }
  }

  // Solve L.y=b via forward substitution (Alg 3.1.1);
  b[0] = b[0];
  for(i=1;i<N;i++){
    sumi=0.0;
    for(j=0;j<i;j++){
      sumi = sumi + A[i][j]*b[j];
    }
    b[i]=b[i]-sumi;
  }
  
  //Solve U.x=y via backward substitution (Alg 3.1.2)
  
  b[N-1] = b[N-1]/A[N-1][N-1];
  for(i=N-2;i>=0;i--){
    sumi=0.0;
    for(j=i+1;j<N;j++){
      sumi = sumi + A[i][j]*b[j];
    }
    b[i] = (b[i] - sumi)/A[i][i];
  }  
} // End of linsolve

void updateAirSeaFluxes(propT *prop, gridT *grid, physT *phys, metT *met,REAL **T){
 /*
  * Main routine for calculating the air-sea heat and salt fluxes
  *
  * Computed terms are stored in the met structure array
  */ 
  
  int j, ktop, iptr, n;
  int Nc = grid->Nc;
  REAL *x=met->xtmp; // pointer to vector passed to cor30
  REAL Umag; // Wind Speed magnitude
  // Constant flux coefficients
  REAL Cd = prop->Cda; // Drag coefficient
  REAL Ch = prop->Ch; // Stanton Number
  REAL Ce = prop->Ce; //Dalton Number
  REAL cp = 4186.0; // Specific heat of water
  REAL Lv = 2.50e6; // Latent heat of vaporization
  REAL rhoa = 1.20; // Density of air	
  

  /* Terms in the COARE3.0 input vector
   REAL u=x[0]; //wind speed (m/s]  at height zu [m]
   REAL us=x[1]; //surface current speed in the wind direction [m/s]
   REAL ts=x[2]; //bulk water temperature [C] if jcool=1, interface water T if jcool=0  
   REAL t=x[3]; //bulk air temperature [C], height zt
   REAL Qs=x[4]; //bulk water spec hum [kg/kg] if jcool=1, ...
   REAL Q=x[5]; //bulk air spec hum [kg/kg], height zq
   REAL Rs=x[6]; //downward solar flux [W/m^2]
   REAL Rl=x[7]; //downard IR flux [W/m^2]
   REAL rain=x[8]; //rain rate [mm/hr]
   REAL zi=x[9]; //PBL depth [m]
   REAL P=x[10]; //Atmos surface pressure [mb]
   REAL zu=x[11]; //wind speed measurement height [m]
   REAL zt=x[12]; //air T measurement height [m]
   REAL zq=x[13]; //air q measurement height [m]
   REAL lat=x[14]; //latitude [deg, N=+]
   REAL jcool=x[15]; //implement cool calculation skin switch, 0=no, 1=yes
   REAL jwave=x[16]; //implement wave dependent roughness model
   REAL twave=x[17]; //wave period [s]
   REAL hwave=x[18]; //wave height [m]
  */
  
  /* The goal is fill up vector x   */ 
  //x = (REAL *)SunMalloc(19*sizeof(REAL),"updateAirSeaFluxes");

  
  //if(myproc==0) printf(" j, Hs, Hl, tau, Hlw, Hsw\n"); 
  for(j=0;j<Nc;j++){
// for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
//    j = grid->cellp[iptr];  
    ktop = grid->ctop[j];
    // Wind speed
    Umag = sqrt( pow( (met->Uwind[j]-phys->uc[j][ktop]) ,2) + pow( (met->Vwind[j]-phys->vc[j][ktop]),2) );
    //Umag = sqrt( pow(met->Uwind[j],2) + pow(met->Vwind[j],2) );
    x[0] = Umag;

    // Surface current speed in wind direction
    // This is the projection of the water velocity vector onto the wind velocity vector
    x[1] = 0.0; 
    //x[1] = fabs(phys->uc[j][ktop]*met->Uwind[j]/Umag + phys->vc[j][ktop]*met->Vwind[j]/Umag); 

    // Water temperature
    x[2] = T[j][ktop];
    
    // Air temperature
    x[3] = met->Tair[j];
    
    // Water specific humidty
    x[4] = qsat(T[j][0], met->Pair[j]);
    
    // Air specific humidity
    x[5] = specifichumidity(met->RH[j],met->Tair[j],met->Pair[j]);
    
    // Longwave radiation
    met->Hlw[j] = longwave(met->Tair[j],T[j][ktop],met->cloud[j]);
    x[6] = met->Hlw[j];
    
    // Shortwave radiation
    met->Hsw[j] = shortwave(prop->nctime/86400.0,prop->latitude,met->cloud[j]);
    x[7] = met->Hsw[j];
    
    //rain [mm/hr] (rain heat flux is not included at the moment)
    x[8] = met->rain[j]*3600;
    
    //Air pressure [mb]
    x[10] = met->Pair[j];
    
    //wind speed height
    x[11] = met->z_Uwind[j];
    
    //air temp height
    x[12] = met->z_Tair[j];
    
    //humidity measurement height
    x[13] = met->z_RH[j];
    
     //Set some constant values
     x[9] = 600.0; 
     x[14] = prop->latitude;
     x[15] = 0.0;
     x[16] = 0.0;
     x[17] = 0.0;
     x[18] = 0.0;
    //printf(" %d, x[0]:%6.6f, x[1]:%6.6f, x[2]:%6.6f, x[3]:%6.6f, x[4]:%6.6f, x[5]:%6.6f\n",j,x[0],x[1],x[2],x[3],x[4],x[5]);
    if (prop->metmodel==2){
      /* Calculate the actual fluxes */
      cor30a(x);
      /* Output variable back into the input array 
      x[0] = hsb;
      x[1] = hlb;
      x[2] = tau;
      x[3] = zo;
      x[4] = zot;
      x[5] = zoq;
      x[6] = rhoa;
      x[7] = usr;
      x[8] = tsr;
      x[9] = qsr;
      x[10] = dter;
      x[11] = dqer;
      x[12] = tkt;
      x[13] = RF;
      x[14] = ut; // Speed including gustiness
      x[15] = Cd;
      x[16] = Ch;
      x[17] = Ce;
      x[18] = ug;
      */
      // Output the fluxes to the met structure array
      met->Hs[j] = -x[0];//Note the change of sign
      met->Hl[j] = -x[1];
      met->ustar[j] = x[7];//These are not used at present
      met->Tstar[j] = x[8];
      met->qstar[j] = x[9];
    
      /* Calculate the wind stress components
      * tau_x = rhoa * Cd * S * (Ucurrent - Uwind) : Fairall et al, 1996
      */
      // *** I think this is double counting if the surface currents have already been accounted for above***
      //met->tau_x[j] = x[6] * x[15] * x[14] * (phys->uc[j][ktop] - met->Uwind[j]);
      //met->tau_y[j] = x[6] * x[15] * x[14] * (phys->vc[j][ktop] - met->Vwind[j]);
      //met->tau_x[j] = x[6] * x[15] * x[14] * (met->Uwind[j] - phys->uc[j][ktop]);
      //met->tau_y[j] = x[6] * x[15] * x[14] * (met->Vwind[j] - phys->vc[j][ktop]);
      // No gust speed in stress term
      met->tau_x[j] = x[6] * x[15] * Umag * (met->Uwind[j] - phys->uc[j][ktop]);
      met->tau_y[j] = x[6] * x[15] * Umag * (met->Vwind[j] - phys->vc[j][ktop]);

      //No surface current dependence 
      //met->tau_x[j] = x[6] * x[15] * x[14] * met->Uwind[j];
      //met->tau_y[j] = x[6] * x[15] * x[14] * met->Vwind[j];
      //met->tau_x[j] = 1.2 * 0.0011 * x[14] * met->Uwind[j];
      //met->tau_y[j] = 1.2 * 0.0011 * x[14] * met->Vwind[j];
      //printf("%10.6f, %10.6f, %10.6f, %10.6f\n",x[6],x[15],x[14],met->Vwind[j]);

    }else if(prop->metmodel==3){// Compute fluxes with constant parameters
      met->Hs[j] = - rhoa * cp * Ch * Umag * (x[2] - x[3]);
      met->Hl[j] = - rhoa * Lv * Ce * Umag * (x[4] - x[5]);
      met->tau_x[j] = rhoa * Cd * Umag * (met->Uwind[j] - phys->uc[j][ktop]); 
      met->tau_y[j] = rhoa * Cd * Umag * (met->Vwind[j] - phys->vc[j][ktop]);
    }
    // Check for nans and dump the inputs
    for(n=0;n<19;n++){
    	if(x[n]!=x[n]){
	   printf("Error in COARE3.0 Algorithm at j = %d, x[%d] = nan.\n",j,n);
	   printf("Uwind[%d] = %6.10f, z_Uwind = %6.10f m\n",j,met->Uwind[j],met->z_Uwind[j]);
	   printf("Vwind[%d] = %6.10f, z_Vwind = %6.10f m\n",j,met->Vwind[j],met->z_Vwind[j]);
	   printf("Tair[%d] = %6.10f, z_Tair = %6.10f m\n",j,met->Tair[j],met->z_Tair[j]);
	   printf("Pair[%d] = %6.10f\n",j,met->Pair[j]);
	   printf("rain[%d] = %6.10f\n",j,met->rain[j]);
	   printf("RH[%d] = %6.10f, z_RH = %6.10f m\n",j,met->RH[j],met->z_RH[j]);
	   printf("cloud[%d] = %6.10f\n",j,met->cloud[j]);
	   MPI_Finalize();
	   exit(EXIT_FAILURE);
	}
    }
  }
  //SunFree(x,19,"updateAirSeaFluxes");
} // End updateAirFluxes

static REAL specifichumidity(REAL RH, REAL Ta, REAL Pair){
 /* 
  * Convert relative humidity (%) to specific humidity (kg/kg 
  */
 
 REAL cff;
 /*
  * Compute air saturation vapor pressure (mb), using Teten formula.
  */
  cff=(1.0007+3.46E-6*Pair)*6.1121*exp(17.502*Ta/(240.97+Ta));
  
  /*
  *  Compute specific humidity at Saturation, Qair (kg/kg).
  */
  cff=cff*RH/100;                    // Vapor pres (mb)
  return (0.62197*(cff/(Pair-0.378*cff))); // Spec hum (kg/kg)
} // End specifichumidity

static REAL qsat(REAL Tw, REAL Pair){
 /*
  *Compute water saturation vapor pressure (mb), using Teten formula.
  */
 REAL cff;
 
 cff=(1.0007+3.46E-6*Pair)*6.1121* exp(17.502*Tw/(240.97+Tw));

  //  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
  cff=cff*0.98;
 
  return (0.62197*(cff/(Pair-0.378*cff)));
} // End qsat


static REAL longwave(REAL Ta, REAL Tw, REAL C_cloud){
  /* Calculate net longwave radiation into water
   *
   * Ref: Martin and McCutcheon, "Hydrodynamics and Transport for Water
   * Quality Modeling", 1999
   */
  
  // Constants
  const REAL T_ref = 273.16;             // conversion from C to K
  const REAL sigma = 5.67e-8;         // Boltzmann constant (W m^{-2} K^{-4})
  const REAL alpha_0 = 0.937e-5;         // Proportionality constants ([alpha_0]= K^{-2}, all others dimensionless)
  const REAL alpha_LW = 0.17;            // LW cloud cover fraction coefficient
  const REAL r_LW = 0.03;                // Fraction of longwave radiation reflected at surface
  const REAL epsilon_w = 0.97;           // emissivity of water
  
  REAL H_LE, epsilon_a, H_LW;
  
  // Emitted Long Wave Radiation
  H_LE = epsilon_w*sigma*pow(Tw + T_ref,4);
  
  //  Incoming Long Wave Radiation
  
  //emissivity of air
  epsilon_a = alpha_0*(1+alpha_LW*C_cloud)*pow(Ta + T_ref,2);
  
  H_LW = epsilon_a*sigma*(1-r_LW)*pow(Ta+T_ref,4);
  
  return  (H_LW - H_LE);

} // End longwave

REAL shortwave(REAL time, REAL Lat,REAL C_cloud){
  /*
   *Compute solar radiation flux using the Gill, 1982 formulae 
   *
   */
  
  const REAL S=1368.0;  //[W m-2]
  const REAL albedo = 0.06;
  const REAL R = 0.76;
  REAL omega1, omega0;
  REAL delta, singamma, Qsc, Hsw;

  omega1 = 2*PI/1.0; // Diurnal cycle
  omega0 = 2*PI/365.25; // Annual cycle


  delta = 23.5*PI/180 * cos(omega0*time - 2.95);
  singamma = sin(delta)*sin(PI*Lat/180) - cos(delta)*cos(PI*Lat/180)*cos(omega1*time);

  // Clear sky radiation
  if(singamma >= 0.0){
      Qsc = R*S*singamma;
  }else{
      Qsc = 0;
  }
  // Corrected radiation
  Hsw = (1.0-albedo) * (1.0 - 0.65 * pow(C_cloud,2)) * Qsc;

  return Hsw;
  
} // End shortwave

static void cor30a(REAL *y){
  /* 
 * Air-Sea fluxes based on bulk flux formulation
 *
 * References:
 *	Fairall et al, 1996, JGR
 *  	Fairall et al, 2003, Journal of Climate
 *
 * See:
 *	http://coaps.fsu.edu/COARE/flux_algor/
 *
 * Adapted from cor30.m matlab function
 */
  
  /* Get variables from the input vector */
  REAL u=y[0]; //wind speed (m/s]  at height zu [m]
  REAL us=y[1]; //surface current speed in the wind direction [m/s]
  REAL ts=y[2]; //bulk water temperature [C] if jcool=1, interface water T if jcool=0  
  REAL t=y[3]; //bulk air temperature [C], height zt
  REAL Qs=y[4]; //bulk water spec hum [kg/kg] if jcool=1, ...
  REAL Q=y[5]; //bulk air spec hum [kg/kg], height zq
  REAL Rs=y[6]; //downward solar flux [W/m^2]
  REAL Rl=y[7]; //downard IR flux [W/m^2]
  REAL rain=y[8]; //rain rate [mm/hr]
  REAL zi=y[9]; //PBL depth [m]
  REAL P=y[10]; //Atmos surface pressure [mb]
  REAL zu=y[11]; //wind speed measurement height [m]
  REAL zt=y[12]; //air T measurement height [m]
  REAL zq=y[13]; //air q measurement height [m]
  REAL lat=y[14]; //latitude [deg, N=+]
  REAL jcool=y[15]; //implement cool calculation skin switch, 0=no, 1=yes
  REAL jwave=y[16]; //implement wave dependent roughness model
  REAL twave=y[17]; //wave period [s]
  REAL hwave=y[18]; //wave height [m]
  
  int nits =3;
  int i;
  
  REAL Beta, von, fdg, tdk, grav;
  REAL Rgas, Le, cpa, cpv, rhoa, visa;
  REAL Al, be, cpw, rhow, visw, tcw, bigc, wetc;
  REAL lwave, cwave, Rns, Rnl;
  REAL dt, du, dq, ta, ug, dter, dqer;
  REAL u10, usr, zo10, Cd10, Ch10, Ct10, zot10;
  REAL Cd, Ct, CC, Ribcu, Ribu, zetu, L10, tsr, qsr, tkt, charn, rr, Bf;
  REAL hsb, hlb, qout, dels, qcol, alq;
  REAL xlamx, dwat, dtmp, alfac, RF, wbar, hl_webb, zet;
  REAL zo, zot, zoq, L, ut, tau, Ch, Ce, Chn_10, Cdn_10, Cen_10;
  
  //printf("y[13]:%6.6f, y[14]:%6.6f, y[15]:%6.6f, y[16]:%6.6f, y[17]:%6.6f, y[18]:%6.6f\n",y[13],y[14],y[15],y[16],y[17],y[18]);
  
/**********   set constants *************/
  Beta=1.2;
  von=0.4;
  fdg=1.00;
  tdk=273.16;
  grav=9.81;
  /*************  air constants ************/
  Rgas=287.1;
  Le=(2.501-.00237*ts)*1e6;
  cpa=1004.67;
  cpv=cpa*(1+0.84*Q);
  rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*Q));
  visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t);
   /************  cool skin constants  *******/
   Al=2.1e-5*pow(ts+3.2,0.79);
   be=0.026;
   cpw=4000;
   rhow=1022;
   visw=1e-6;
   tcw=0.6;
   bigc=16*grav*cpw*pow(rhow*visw,3.0)/(tcw*tcw*rhoa*rhoa);
   wetc=0.622*Le*Qs/(Rgas*pow(ts+tdk,2.0));
   
   zo = 0.0;
   zot=0.0;
   zoq=0.0;
   L=0.0;
   /***************   wave parameters  *********/
   lwave=grav/2/PI*pow(twave,2);
   cwave=grav/2/PI*twave;
   
   /**************  compute aux stuff *******/
   Rns=Rs*.945;
   Rnl=0.97*(5.67e-8*pow(ts-0.3*jcool+tdk,4)-Rl);
   
   /***************   Begin bulk loop *******/
   
   /***************  first guess ************/
   du=u-us;
   dt=ts-t-.0098*zt;
   dq=Qs-Q;
   ta=t+tdk;
   ug=.5;
   dter=0.3; 
   dqer=wetc*dter;
   ut=sqrt(du*du+ug*ug);
   u10=ut*log(10/1e-4)/log(zu/1e-4);
   usr=.035*u10;
   zo10=0.011*usr*usr/grav+0.11*visa/usr;
   Cd10=pow(von/log(10/zo10),2);
   Ch10=0.00115;
   Ct10=Ch10/sqrt(Cd10);
   zot10=10/exp(von/Ct10);
   Cd=pow(von/log(zu/zo10),2);
   Ct=von/log(zt/zot10);
   CC=von*Ct/Cd;
   Ribcu=-zu/zi/.004/pow(Beta,3);
   Ribu=-grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/pow(ut,2);

   if(Ribu<0){
     zetu=CC*Ribu/(1+Ribu/Ribcu);
   }else{
     zetu=CC*Ribu*(1+27/9*Ribu/CC);
   }
   
   L10=zu/zetu;
   if (zetu>50){
    nits=1;
   }
   usr=ut*von/(log(zu/zo10) - psiu_30(zu/L10));
   tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_30(zt/L10));
   qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-psit_30(zq/L10));
   
   tkt=.001;
   
   charn=0.011;
   if (ut>10.0){
     charn=0.011+(ut-10)/(18-10)*(0.018-0.011);
   }
   if(ut>18.0){
     charn=0.018;
   }
   
   /***************  bulk loop ************/
  for(i=0;i<nits;i++){
     zet=von*grav*zu/ta*(tsr*(1+0.61*Q)+.61*ta*qsr)/(usr*usr)/(1+0.61*Q);
     //Hard-wire it without waves
     zo=charn*usr*usr/grav+0.11*visa/usr;

//       if((int)jwave==0){
//           zo=charn*usr*usr/grav+0.11*visa/usr;
//       }
//       if((int)jwave==1){
//           zo=50.0/2.0/PI*lwave*pow(usr/cwave,4.5)+0.11*visa/usr;
//       } // Oost et al
//       if((int)jwave==2){
//           zo=1200.0*hwave*pow(hwave/lwave,4.5)+0.11*visa/usr;
//       } // Taylor and Yelland
      
      rr=zo*usr/visa;
      L=zu/zet;
      zoq=Min(1.15e-4,5.5e-5/pow(rr,0.6));
      zot=zoq;
      usr=ut*von/(log(zu/zo)-psiu_30(zu/L));
      tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psit_30(zt/L));
      qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psit_30(zq/L));
      Bf=-grav/ta*usr*(tsr+.61*ta*qsr);
      
      if (Bf>0){
	ug=Beta*pow(Bf*zi,0.333);
      }else{
	ug=.2;
      }
      
      ut=sqrt(du*du+ug*ug);
      Rnl=0.97*(5.67e-8*pow(ts-dter*jcool+tdk,4)-Rl);
      hsb=-rhoa*cpa*usr*tsr;
      hlb=-rhoa*Le*usr*qsr;
      qout=Rnl+hsb+hlb;
      dels=Rns*(.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4))); 	// Eq.16 Shortwave
      qcol=qout-dels;
      alq=Al*qcol+be*hlb*cpw/Le;					// Eq. 7 Buoy flux water

     if (alq>0){
     	xlamx=6/pow(1+pow(bigc*alq/pow(usr,4),0.75),0.333);				// Eq 13 Saunders
        tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr);			//Eq.11 Sub. thk

     }else{
       xlamx=6.0;
       tkt=Min(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr));			// Eq.11 Sub. thk
     }
     
      dter=qcol*tkt/tcw; //  Eq.12 Cool skin
      dqer=wetc*dter;
     
  }//bulk iter loop
  
  tau=rhoa*usr*usr*du/ut;  //stress
  hsb=-rhoa*cpa*usr*tsr;
  hlb=-rhoa*Le*usr*qsr;
  
  
  /****************   rain heat flux ********/
  dwat=2.11e-5*pow((t+tdk)/tdk,1.94); // water vapour diffusivity
  dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); 	//heat diffusivity
  alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp));      	// wet bulb factor
  RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600;
  
  /****************   Webb et al. correection  ************/
  wbar=1.61*hlb/Le/(1+1.61*Q)/rhoa+hsb/rhoa/cpa/ta;//formulation in hlb already includes webb
  hl_webb=rhoa*wbar*Q*Le;
  
  /**************   compute transfer coeffs relative to ut @meas. ht **********/
  Cd=tau/rhoa/ut/max(.1,du);
  Ch=-usr*tsr/ut/(dt-dter*jcool);
  Ce=-usr*qsr/(dq-dqer*jcool)/ut;
  
  /************  10-m neutral coeff realtive to ut ********/
  Cdn_10=von*von/log(10/zo)/log(10/zo);
  Chn_10=von*von*fdg/log(10/zo)/log(10/zot);
  Cen_10=von*von*fdg/log(10/zo)/log(10/zoq);
  
  /* Output variable back into the input array */
  y[0] = hsb;
  y[1] = hlb;
  y[2] = tau;
  y[3] = zo;
  y[4] = zot;
  y[5] = zoq;
  y[6] = rhoa;
  y[7] = usr;
  y[8] = tsr;
  y[9] = qsr;
  y[10] = dter;
  y[11] = dqer;
  y[12] = tkt;
  y[13] = RF;
  y[14] = ut;
  y[15] = Cd;
  y[16] = Ch;
  y[17] = Ce;
  y[18] = ug;
  
  //printf("y[13]:%6.6f, y[14]:%6.6f, y[15]:%6.6f, y[16]:%6.6f, y[17]:%6.6f, y[18]:%6.6f\n",y[13],y[14],y[15],y[16],y[17],y[18]);
  //y=[hsb hlb tau zo zot zoq L usr tsr qsr dter dqer tkt RF wbar Cd Ch Ce Cdn_10 Chn_10 Cen_10 ug ];
  //   1   2   3   4  5   6  7  8   9  10   11   12  13  14  15  16 17 18    19      20    21  22
} // End cor30a

static REAL psiu_30(REAL zet){
  
  REAL x, psik, psic, f, c,  psi;

  x=pow(1.0-15.0*zet,0.25);
  psik=2.0*log((1.0+x)/2.0)+log((1+x*x)/2.0)-2*atan(x)+2.0*atan(1.0);
  x=pow(1-10.15*zet,0.3333);
  psic=1.5*log((1.0+x+x*x)/3.0)-sqrt(3.0)*atan((1.0+2.0*x)/sqrt(3.0))+4.0*atan(1.0)/sqrt(3.0);
  f=zet*zet/(1.0+zet*zet);
  psi=(1.0-f)*psik+f*psic;                                               
  
  if(zet>0){
    c=Min(50.0,0.35*zet);
    psi=-(pow(1.0+1.0*zet,1.0)+.667*(zet-14.28)/exp(c)+8.525);
  }
  return psi;
  
} // End psiu_30

static REAL psit_30(REAL zet){
 
  REAL x, psik, psic, f, c,  psi;
  
  x=pow(1.0-15.0*zet,0.5);
  psik=2.0*log((1+x)/2.0);
  x=pow(1-34.15*zet,0.3333);
  psic=1.5*log((1.0+x+x*x)/3.0)-sqrt(3.0)*atan((1.0+2.0*x)/sqrt(3.0))+4.0*atan(1.0)/sqrt(3.0);
  f=zet*zet/(1+zet*zet);
  psi=(1.0-f)*psik+f*psic;  
  
  if(zet>0){
    c=Min(50,0.35*zet);
    psi=-(pow(1+2/3*zet,1.5)+0.6667*(zet-14.28)/exp(c)+8.525);
  }
  return psi;
} // End psit_30
