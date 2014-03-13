/*
 * All functions related to the transport of age tracers
 * ----------------------------------------------------
 *  
 */


#include "age.h"
#include "memory.h"

/*
 * Private functions
 */
static void AllocateAgeVariables(gridT *grid, ageT **age, propT *prop);
static void InitializeAgeVariables(gridT *grid, propT *prop, int myproc);

/*
 * Function: AllocateAgeVariables()
 * ---------------------------------------
 * Allocate memory to the age structure
 *
 */

static void AllocateAgeVariables(gridT *grid, ageT **age, propT *prop){

  int i,j,jptr,Nc=grid->Nc;

  // allocate  structure
  *age = (ageT *)SunMalloc(sizeof(ageT),"AllocateAgeVariables");

  (*age)->agec = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAgeVariables");
  (*age)->agealpha = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAgeVariables");
  (*age)->Cn_Ac = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAgeVariables");
  (*age)->Cn_Aa = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAgeVariables");

  // for each cell allocate memory for the number of layers at that location
  for(i=0;i<Nc;i++) {
      (*age)->agec[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAgeVariables");
      (*age)->agealpha[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAgeVariables");
      (*age)->Cn_Ac[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAgeVariables");
      (*age)->Cn_Aa[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAgeVariables");
  }

  // allocate boundary variables 
  (*age)->boundary_age = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocateAgeVariables");
  (*age)->boundary_agealpha = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocateAgeVariables");

  // allocate over vertical layers
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
    j=grid->edgep[jptr];
    (*age)->boundary_age[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAgeVariables");
    (*age)->boundary_agealpha[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAgeVariables");
  }

}//End Function


/*
 * Function: InitializeAgeVariables()
 * ---------------------------------------
 * Initialize values into the (global) age structure
 *
 */

static void InitializeAgeVariables(gridT *grid, propT *prop, int myproc){

    int i,k,Nc=grid->Nc;
    REAL *ncscratch;
    int Nci, Nki, T0;

    // Read the netcdf properties and allocate a scratch array for reading in the data
    if (prop->readinitialnc>0){
	ReadInitialNCcoord(prop,grid,&Nci,&Nki,&T0,myproc);
	ncscratch = (REAL *)SunMalloc(Nki*Nci*sizeof(REAL),"InitializeAgeVariables");
    }

    for(i=0;i<Nc;i++) {
	for(k=0;k<grid->Nk[i];k++) {
	    age->agec[i][k]=0;
	    age->agealpha[i][k]=0;
	}
    }

    // Initialise the age arrays (netcdf only)
    if (prop->calcage && prop->readinitialnc)
       ReturnAgeNC(prop,grid,ncscratch,Nci,Nki,T0,myproc);

}//End function

/*
 * Function: UpdateAge()
 * ---------------------------------------
 * Update the age concentration quantities
 *
 */
void UpdateAge(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
    int i, ib, iptr, j, jptr, k;

    //Allocate arrays
    if(prop->n==prop->nstart+1){
    	AllocateAgeVariables(grid,&age,prop);
	InitializeAgeVariables(grid, prop, myproc);
    }

    // Specify age at boundaries for use in updatescalars. 
    //Type-2 -set value to 1 
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
        j = grid->edgep[jptr];
        ib = grid->grad[2*j];
        for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
          //age->boundary_tmp[jptr-grid->edgedist[2]][k]=age->agec[ib][k];
	  age->boundary_age[jptr-grid->edgedist[2]][k]=1;
    }
    //Type-3 -set value to 0 
    for(jptr=grid->edgedist[3];jptr<grid->edgedist[5];jptr++) {
        j = grid->edgep[jptr];
        ib = grid->grad[2*j];
        for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
	  age->boundary_age[jptr-grid->edgedist[3]][k]=0;
    }



    //printf("Updating agec...\n");
    //printf("prop->rtime = %f\n",prop->rtime);
    UpdateScalars(grid,phys,prop,phys->wnew,age->agec,age->boundary_age,age->Cn_Ac,prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,0,prop->TVDsalt);
//    UpdateScalars(grid,phys,prop,phys->wnew,phys->agec,phys->boundary_age,phys->Cn_Ac,prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,phys->uold,phys->wtmp,NULL,NULL,0,0,comm,myproc,0,prop->TVDsalt);

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      for(k=grid->ctop[i];k<grid->Nk[i];k++){
//	 age->agec[i][k] = age->agec[i][k]*prop->dt; 
	 age->agec[i][k] = age->agec[i][k]; 
      }
    }

    ISendRecvCellData3D(age->agec,grid,myproc,comm);

    for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
        j = grid->edgep[jptr];
        ib = grid->grad[2*j];

        for(k=grid->ctop[ib];k<grid->Nk[ib];k++) 
          //age->boundary_tmp[jptr-grid->edgedist[2]][k]=age->agealpha[ib][k];
	  age->boundary_agealpha[jptr-grid->edgedist[2]][k]=0.0;
    }

    //printf("Updating agealpha...\n");
    //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    //  i = grid->cellp[iptr];
    //  for(k=grid->ctop[i];k<grid->Nk[i];k++){
    //     phys->uold[i][k] = phys->agec[i][k]; 
    //     phys->wtmp[i][k] = 0;
    //  }
    //}


    //UpdateScalars(grid,phys,prop,phys->wnew,phys->agealpha,phys->boundary_agealpha,phys->Cn_Aa, prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,phys->uold,phys->wtmp,NULL,NULL,0,0,comm,myproc,0,prop->TVDsalt);
    UpdateScalars(grid,phys,prop,phys->wnew,age->agealpha,age->boundary_agealpha,age->Cn_Aa, prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,0,prop->TVDsalt);


    // Alpha parameter source term
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      for(k=grid->ctop[i];k<grid->Nk[i];k++){
         //These are the source terms for the alpha parameter
         //phys->wtmp[i][k] = 0;
         //phys->uold[i][k] = phys->agec[i][k]*prop->dt;
         age->agealpha[i][k] += age->agec[i][k]*prop->dt; 
         //phys->agealpha[i][k] = phys->agealpha[i][k] + phys->agec[i][k]*prop->rtime; 
      }
    }

    ISendRecvCellData3D(age->agealpha,grid,myproc,comm);
    
    //printf("Done\n");
    
}



