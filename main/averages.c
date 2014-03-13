/*
* Averaging functions
* -------------------
* All functions to calculate and output average quantities 
* 
*/ 

#include "averages.h"

/***********************************************
* Private functions
***********************************************/
/*
 * Function: AllocateAverageVariables()
 * ------------------------------------
 * Allocate memory to the average variable arrays
 *
 */
void AllocateAverageVariables(gridT *grid, averageT **average, propT *prop)
{
  int flag=0, i, j, jptr, ib, Nc=grid->Nc, Ne=grid->Ne, Np=grid->Np, nf, k;


  prop->avgctr=0;
  prop->avgtimectr=0;
  prop->avgfilectr=0;

  // allocate averageical structure
  *average = (averageT *)SunMalloc(sizeof(averageT),"AllocateAverageVariables");

  // Allocate 3D arrays
  (*average)->uc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->vc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->w = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->nu_v = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->s = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->rho = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");

  // Edge variables
  (*average)->U_F = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->s_F = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->T_F = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");

  if(prop->calcage)
      (*average)->agemean = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");

  // cell-centered averageical variables in plan (no vertical direction)
  (*average)->h = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
  (*average)->s_dz = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
  (*average)->T_dz = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");

  if(prop->metmodel>0){
    (*average)->Uwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Vwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Tair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Pair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->rain = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->RH = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->cloud = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");

    (*average)->Hs = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Hl = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Hlw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Hsw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->tau_x = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->tau_y = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->EP = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
  }
  
  // for each cell allocate memory for the number of layers at that location
  for(i=0;i<Nc;i++) {
      (*average)->uc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->vc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->w[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateAverageVariables");
      (*average)->s[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->rho[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->nu_v[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      if(prop->calcage)
         (*average)->agemean[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");     	

  }
  
  // Allocate edge variables
  for(j=0;j<Ne;j++){
      (*average)->U_F[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->s_F[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->T_F[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");

  }
  // Netcdf write variables
  (*average)->tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"AllocateAverageVariables");
  (*average)->tmpvarE = (REAL *)SunMalloc(grid->Ne*grid->Nkmax*sizeof(REAL),"AllocateAverageVariables");
  (*average)->tmpvarW = (REAL *)SunMalloc(grid->Nc*(grid->Nkmax+1)*sizeof(REAL),"AllocateAverageVariables");
 
}//End of function

/*
 * Function: ZeroAverageVariables() 
 * ---------------------------------
 * Zero the average arrays
 */
void ZeroAverageVariables(gridT *grid, averageT *average, propT *prop){

    int i,j,k,Nc=grid->Nc,Ne=grid->Ne;

  for(i=0;i<Nc;i++) {
    average->w[i][grid->Nk[i]]=0;
    for(k=0;k<grid->Nk[i];k++) {
      average->uc[i][k]=0;
      average->vc[i][k]=0;
      average->w[i][k]=0;
      average->s[i][k]=0;
      average->T[i][k]=0;
      average->rho[i][k]=0;
      average->nu_v[i][k]=0;

      if(prop->calcage){
      	average->agemean[i][k]=0;
      }
    }
    // 2D cell-centred variables
    average->h[i]=0;
    average->s_dz[i]=0;
    average->T_dz[i]=0;
    if(prop->metmodel>0){
	average->Uwind[i]=0;
	average->Vwind[i]=0;
	average->Tair[i]=0;
	average->Pair[i]=0;
	average->rain[i]=0;
	average->RH[i]=0;
	average->cloud[i]=0;
	average->Hs[i]=0;
	average->Hl[i]=0;
	average->Hlw[i]=0;
	average->Hsw[i]=0;
	average->tau_x[i]=0;
	average->tau_y[i]=0;
	average->EP[i]=0;
    }
  }

  for(j=0;j<Ne;j++){
      for(k=0;k<grid->Nke[j];k++){
	  average->U_F[j][k]=0;
	  average->s_F[j][k]=0;
	  average->T_F[j][k]=0;
      }
  }

}//End of function

/*
 * Function: UpdateAverageVariables() 
 * ---------------------------------
 * Updates the average arrays by summing the last time step
 *
 */
void UpdateAverageVariables(gridT *grid, averageT *average, physT *phys, metT *met, propT *prop,MPI_Comm comm, int myproc){

    int i,iptr,j,jptr,k,Nc=grid->Nc,Ne=grid->Ne;
    int nc1, nc2;
    REAL flx, sdz,Tdz,theta=prop->theta;

//  for(i=0;i<Nc;i++) {
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    average->w[i][grid->Nk[i]]+=phys->w[i][grid->Nk[i]];
    sdz=0;
    Tdz=0;
    //for(k=0;k<grid->Nk[i];k++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      average->uc[i][k]+=phys->uc[i][k];
      average->vc[i][k]+=phys->vc[i][k];
      average->w[i][k]+=phys->w[i][k];
      average->s[i][k]+=phys->s[i][k];
      average->T[i][k]+=phys->T[i][k];
      average->rho[i][k]+=phys->rho[i][k];
      average->nu_v[i][k]+=phys->nu_tv[i][k];

      if(prop->calcage){
        if(age->agec[i][k]>1e-10){
	    // Calculate the mean online here
	    average->agemean[i][k]+=age->agealpha[i][k]/age->agec[i][k];
	}else{
	    average->agemean[i][k]+=0;
        }
      }
      sdz+=phys->s[i][k]*grid->dzz[i][k];
      Tdz+=phys->T[i][k]*grid->dzz[i][k];
    }
    // 2D cell-centred variables
    average->h[i]+=phys->h[i];
    //average->s_dz[i]+=sdz;
    //average->T_dz[i]+=Tdz;
    //average->h[i]=phys->h[i];
    average->s_dz[i]=sdz; //Instantaneous values (used for scalar budget)
    average->T_dz[i]=Tdz; // Instantaneous values 
    if(prop->metmodel>0 && prop->metmodel<4){
	average->Uwind[i]+=met->Uwind[i];
	average->Vwind[i]+=met->Vwind[i];
	average->Tair[i]+=met->Tair[i];
	average->Pair[i]+=met->Pair[i];
	average->rain[i]+=met->rain[i];
	average->RH[i]+=met->RH[i];
	average->cloud[i]+=met->cloud[i];
	average->Hs[i]+=met->Hs[i];
	average->Hl[i]+=met->Hl[i];
	average->Hlw[i]+=met->Hlw[i];
	average->Hsw[i]+=met->Hsw[i];
	average->tau_x[i]+=met->tau_x[i];
	average->tau_y[i]+=met->tau_y[i];
	average->EP[i]+=met->EP[i];
    }
  }

//  for(j=0;j<Ne;j++){
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
//     for(k=0;k<grid->Nke[j];k++){
      for(k=grid->etop[j];k<grid->Nke[j];k++){
          //flx = phys->u[j][k]*grid->dzf[j][k]*grid->df[j]; 
	  // Flux needs to be consistent with continuity equation
	  flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k])*grid->dzf[j][k]*grid->df[j]; 
	  
	  average->U_F[j][k] += flx; 
      }
  }
  
  /*
   * Compute salinity and temperature fluxes using TVD scheme
   */
  if(prop->TVD ){

    //Salt
    // Compute the scalar on the vertical faces (for horiz. advection)
    HorizontalFaceScalars(grid,phys,prop,phys->s,phys->boundary_s,prop->TVDsalt,comm,myproc); 
//      for(j=0;j<Ne;j++){
//   	for(k=0;k<grid->Nke[j];k++) {
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      for(k=grid->etop[j];k<grid->Nke[j];k++){
	  //flx = phys->u[j][k]*grid->dzf[j][k]*grid->df[j]; 
	  //See equation 85 in SUNTANS paper
	  flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k])*grid->dzf[j][k]*grid->df[j]; 
	  //if(phys->u[j][k]>0)
	  if(phys->utmp2[j][k]>0)
	    average->s_F[j][k]+=phys->SfHp[j][k] * flx;
	  else
	    average->s_F[j][k]+=phys->SfHm[j][k] * flx;
	}
      }

     //Temperature
    HorizontalFaceScalars(grid,phys,prop,phys->T,phys->boundary_T,prop->TVDtemp,comm,myproc); 
//      for(j=0;j<Ne;j++){
//	for(k=0;k<grid->Nke[j];k++) {
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      for(k=grid->etop[j];k<grid->Nke[j];k++){
	  //flx = phys->u[j][k]*grid->dzf[j][k]*grid->df[j]; 
	  flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k])*grid->dzf[j][k]*grid->df[j]; 
	  //if(phys->u[j][k]>0)
	  if(phys->utmp2[j][k]>0)
	    average->T_F[j][k]+=phys->SfHp[j][k] * flx;
	  else
	    average->T_F[j][k]+=phys->SfHm[j][k] * flx;
	}
      }

  }else{ // No TVD
  
  }//End flux calculation
 /*
   * Compute salinity and temperature fluxes using central-difference
  */

    //Salt and temp
/*
//      for(j=0;j<Ne;j++){
//	for(k=0;k<grid->Nke[j];k++) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
       j = grid->edgep[jptr]; 
       for(k=grid->etop[j];k<grid->Nke[j];k++){

	  flx = phys->u[j][k]*grid->dzf[j][k]*grid->df[j]; 

	  nc1 = grid->grad[2*j];
	  nc2 = grid->grad[2*j+1];
	  if(nc1==-1) nc1=nc2;
	  if(nc2==-1) nc2=nc1;

	  average->s_F[j][k]+= 0.5*(phys->s[nc1][k]+phys->s[nc2][k])*flx;
	  average->T_F[j][k]+= 0.5*(phys->T[nc1][k]+phys->T[nc2][k])*flx;
	}
      }
*/


}//End of function

/*
 * Function: ComputeAverageVariables() 
 * ---------------------------------
 * Computes the average by dividing by the number averaging steps "ntaverage"
 * 
 */
void ComputeAverageVariables(gridT *grid, averageT *average, physT *phys, metT *met, int ntaverage, propT *prop){

    int i,j,k,Nc=grid->Nc,Ne=grid->Ne;
    REAL nt = 1.0/((REAL)ntaverage);

  for(i=0;i<Nc;i++) {
    average->w[i][grid->Nk[i]]+=phys->w[i][grid->Nk[i]];
    for(k=0;k<grid->Nk[i];k++) {
      average->uc[i][k] *= nt;
      average->vc[i][k] *= nt;
      average->w[i][k] *= nt;
      average->s[i][k] *= nt;
      average->T[i][k] *= nt;
      average->rho[i][k] *= nt;
      average->nu_v[i][k] *= nt;

      if(prop->calcage){
	    average->agemean[i][k] *= nt;
      }
    }
    // 2D cell-centred variables
    average->h[i] *= nt;
    //average->s_dz[i] *= nt;
    //average->T_dz[i] *= nt;
    if(prop->metmodel>0){
	average->Uwind[i] *= nt;
	average->Vwind[i] *= nt;
	average->Tair[i] *= nt;
	average->Pair[i] *= nt;
	average->rain[i] *= nt;
	average->RH[i] *= nt;
	average->cloud[i] *= nt;
	average->Hs[i] *= nt;
	average->Hl[i] *= nt;
	average->Hlw[i] *= nt;
	average->Hsw[i] *= nt;
	average->tau_x[i] *= nt;
	average->tau_y[i] *= nt;
	average->EP[i] *= nt;
    }
  }

  for(j=0;j<Ne;j++){
      for(k=0;k<grid->Nke[j];k++){
	  average->U_F[j][k] *= nt;
	  average->s_F[j][k] *= nt;
	  average->T_F[j][k] *= nt;
      }
  }

}//End of function

/*
 * Function: SendRecvAverages()
 * ----------------------------
 *  Communicate the average values amongst processors
 *  Only really necessary for writing
 *
 */
void SendRecvAverages(propT *prop, gridT *grid, averageT *average, MPI_Comm comm, int myproc){
    // Communicate 2D variables
    ISendRecvCellData2D(average->h,grid,myproc,comm);
    ISendRecvCellData2D(average->s_dz,grid,myproc,comm);
    ISendRecvCellData2D(average->T_dz,grid,myproc,comm);
    if(prop->metmodel>0){
	ISendRecvCellData2D(average->Uwind,grid,myproc,comm);
	ISendRecvCellData2D(average->Vwind,grid,myproc,comm);
	ISendRecvCellData2D(average->Tair,grid,myproc,comm);
	ISendRecvCellData2D(average->Pair,grid,myproc,comm);
	ISendRecvCellData2D(average->rain,grid,myproc,comm);
	ISendRecvCellData2D(average->RH,grid,myproc,comm);
	ISendRecvCellData2D(average->cloud,grid,myproc,comm);
	ISendRecvCellData2D(average->Hs,grid,myproc,comm);
	ISendRecvCellData2D(average->Hl,grid,myproc,comm);
	ISendRecvCellData2D(average->Hlw,grid,myproc,comm);
	ISendRecvCellData2D(average->Hsw,grid,myproc,comm);
	ISendRecvCellData2D(average->tau_x,grid,myproc,comm);
	ISendRecvCellData2D(average->tau_y,grid,myproc,comm);
	ISendRecvCellData2D(average->EP,grid,myproc,comm);
    }

    // Communicate the 3D cell data
    ISendRecvCellData3D(average->uc,grid,myproc,comm);
    ISendRecvCellData3D(average->vc,grid,myproc,comm);
    ISendRecvCellData3D(average->s,grid,myproc,comm);
    ISendRecvCellData3D(average->T,grid,myproc,comm);
    ISendRecvCellData3D(average->rho,grid,myproc,comm);
    ISendRecvCellData3D(average->nu_v,grid,myproc,comm);
    if(prop->calcage)
	ISendRecvCellData3D(average->agemean,grid,myproc,comm);

    ISendRecvWData(average->w,grid,myproc,comm);

    // Communicate the 3D edge data
    ISendRecvEdgeData3D(average->U_F,grid,myproc,comm);
    ISendRecvEdgeData3D(average->s_F,grid,myproc,comm);
    ISendRecvEdgeData3D(average->T_F,grid,myproc,comm);
     
}//End of function
