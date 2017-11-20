/*
 * File: sediment.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file contains physically-based functions for cohesive sediment calculation.
 * Cohesive suspended sediments calculation is done in updatescalars when prop->sedi==1
 * Here sediment settling velocity, erosion rate and consolidation in different layers are calculated
 * 
 *
 */
#include "suntans.h"
#include "grid.h"
#include "phys.h"
#include "sendrecv.h"
#include "timer.h"
#include "initialization.h"
#include "boundaries.h"
#include "sediments.h"
#include "util.h"
#include "tvd.h"
#include "mympi.h"
#include "scalars.h"
#include "physio.h"

/*
 * Function: ReadSediProperties
 * Usage: ReadSediProperties(grid,phys,prop,myproc);
 * ----------------------------------------------------
 * Based on sedi.dat, load in the important parameters for 
 * cohesive sediment transport model
 *
 */
void ReadSediProperties(int myproc) {    
  int m;
  char str[BUFFERLENGTH];
  sediments->Nlayer = MPI_GetValue(DATAFILE,"Nlayer","ReadSediProperties",myproc);
  sediments->Nsize = MPI_GetValue(DATAFILE,"Nsize","ReadSediProperties",myproc);   
  sediments->WSconstant = MPI_GetValue(DATAFILE,"WSconstant","ReadSediProperties",myproc);
  //sediments->SETsediment= MPI_GetValue(DATAFILE,"SETsediment","ReadSediProperties",myproc);
  sediments->readSediment= MPI_GetValue(DATAFILE,"readSediment","ReadSediProperties",myproc);
  sediments->bedInterval= MPI_GetValue(DATAFILE,"bedInterval","ReadSediProperties",myproc);
  sediments->bedComplex= MPI_GetValue(DATAFILE,"bedComplex","ReadSediProperties",myproc);
  sediments->ParabolKappa= MPI_GetValue(DATAFILE,"ParabolKappa","ReadSediProperties",myproc);
  sediments->TBMAX= MPI_GetValue(DATAFILE,"TBMAX","ReadSediProperties",myproc);
  if(sediments->bedComplex==1){
    printf("because bedComplex==1, so set bedInterval=1 automatically");
    sediments->bedInterval=1;
  }
  
  // condition check
  /*
    if(sediments->SETsediment==1 && sediments->Nsize<=3 && sediments->Nlayer<=5){
    printf("Nsize = %d, Nlayer= %d, but SETsediment==1 which means Nsize>3 or Nlayer>5. You should set SETsediment as 0.\n",sediments->Nsize, sediments->Nlayer );
    MPI_Finalize();
    exit(EXIT_FAILURE);
    } else if(sediments->SETsediment==0 && (sediments->Nsize>3 || sediments->Nlayer>5)){
    printf("Nsize = %d, Nlayer= %d, but SETsediment==0 which means Nsize<=3 or Nlayer<=5. You should set SETsediment as 1.\n",sediments->Nsize, sediments->Nlayer );
    MPI_Finalize();
    exit(EXIT_FAILURE);  
    }
  */
  
  // condition check
  if(sediments->readSediment==1 && sediments->Nsize>3){
    printf("Nsize = %d>1, but readSediment==1 which means Nsize==1. You should set readSediment as 0 or Nsize==1.\n",sediments->Nsize);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  //allocate spaces for all sediments properties
  sediments->Ds=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Ws0=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Gsedi=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Prt=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Consolid=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  sediments->E0=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  sediments->Taue=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  sediments->Taud=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  sediments->Drydensity=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  sediments->Thickness=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  sediments->Softhard=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  /*sediments->Bedmudratio=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");*/
  //if(sediments->SETsediment==0){
  for(m=1;m<=sediments->Nsize;m++) {
    sprintf(str,"Ds%d",m);
    sediments->Ds[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"Ws0%d",m);
    if(sediments->WSconstant==1)
      sediments->Ws0[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    else
      sediments->Ws0[m-1]=0.0;
    sprintf(str,"Gsedi%d",m);
    sediments->Gsedi[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc); 
    sprintf(str,"Prt%d",m);
    sediments->Prt[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
  }
  for(m=1;m<=sediments->Nlayer;m++) {
    sprintf(str,"Consolid%d",m);
    sediments->Consolid[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"E0%d",m);
    sediments->E0[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"Taue%d",m);
    sediments->Taue[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"Taud%d",m);
    sediments->Taud[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"Drydensity%d",m);
    sediments->Drydensity[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"Thickness%d",m);
    sediments->Thickness[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"softhard%d",m);
    sediments->Softhard[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    //sprintf(str,"Bedmudratio%d",m);
    //sediments->Bedmudratio[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
  }
  //}
  sediments->Chind=MPI_GetValue(DATAFILE,"Chind","ReadSediProperties",myproc);
  sediments->Cfloc=MPI_GetValue(DATAFILE,"Cfloc","ReadSediProperties",myproc);
}

/*
 * Function: AllocateSediment
 * Usage: allocate space for sediment variables
 * ----------------------------------------------------
 * Based on the value from ReadSedimentProperties
 *
 */
void AllocateSediment(gridT *grid, int myproc) { 
  int i,j,jptr,k;

  sediments->SediC = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->SediCbed = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  //sediments->Erosion = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");
  //sediments->Erosion_old = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");   
  sediments->boundary_sediC = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->Ws = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->Deposition = (REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");
  //sediments->Deposition_old = (REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");
  
  if(sediments->TBMAX==1){
    sediments->Seditb = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
    sediments->Seditbmax = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  }
  for(i=0;i<sediments->Nsize;i++){
    sediments->SediC[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    sediments->SediCbed[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables"); 
    sediments->Deposition[i]=(REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
    //Deposition_old[i]=(REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
    // boundary part  
    sediments->boundary_sediC[i] = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocateSediment");
    sediments->Ws[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    //sediments->Erosion[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    //sediments->Erosion_old[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    for(j=0;j<grid->Nc;j++){
      sediments->SediC[i][j] = (REAL *)SunMalloc(grid->Nk[j]*sizeof(REAL), "AllocateSediVariables");
      sediments->SediCbed[i][j] = (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables"); 
      //sediments->Erosion[i][j] = (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables"); 
      //sediments->Erosion_old[i][j] = (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables"); 
      // allocate boundary value
      for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
        k=grid->edgep[jptr];  
        sediments->boundary_sediC[i][jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[k]*sizeof(REAL), "AllocateSediVariables");
      }
      sediments->Ws[i][j] = (REAL *)SunMalloc((grid->Nk[j]+1)*sizeof(REAL), "AllocateSediVariables"); 
    }
  }

  //Layermass = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  sediments->Layerthickness = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  sediments->SediKappa_tv = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  sediments->Wnewsedi= (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  //sediments->Woldsedi= (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  sediments->Erosiontotal= (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  sediments->Erosiontotal_old= (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");

  for(i=0;i<grid->Nc;i++){
    //Layermass[i]=(REAL *)SunMalloc(Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Layerthickness[i]=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Erosiontotal[i] = (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Wnewsedi[i]= (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL), "AllocateSediVariables");
    sediments->SediKappa_tv[i]= (REAL *)SunMalloc((grid->Nk[i])*sizeof(REAL), "AllocateSediVariables");
    //Woldsedi[i]= (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL), "AllocateSediVariables");
    sediments->Erosiontotal_old[i] = (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
  }
  /*if(sediments->SETsediment==1){
    sediments->Ds=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "AllocateSediVariables");
    sediments->Ws0=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "AllocateSediVariables");
    sediments->Gsedi=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "AllocateSediVariables");
    sediments->Prt=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "AllocateSediVariables");
    sediments->Consolid=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->E0=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Taue=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Taud=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Drydensity=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Thickness=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Softhard=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    //sediments->Bedmudratio=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
  }*/
}

/*
 * Function: InitializeSediment
 * Usage: give initial value for sediment variables
 * ----------------------------------------------------
 * Based on the value from ReadSedimentProperties
 * here assume Erosion and Deposition initial value is 0, may have error
 */
void InitializeSediment(gridT *grid, physT *phys, propT *prop,  int myproc) { 
  int i,j,k,ne;
  REAL *stmp,z;
  FILE *InitSedimentFID;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  // give 0 for SediC boundary_sediC Ws
  for(i=0;i<sediments->Nsize;i++) {
    for(j=0;j<grid->Nc;j++){
      for(k=0;k<grid->Nk[j];k++) {
        sediments->SediC[i][j][k]=0;        
        sediments->Ws[i][j][k]=0;
      }
      sediments->Ws[i][j][grid->Nk[j]]=0;
    }

    for(j=0;j<(grid->edgedist[5]-grid->edgedist[2]);j++) {
      ne=grid->edgep[j+grid->edgedist[2]]; 
      for(k=0;k<grid->Nke[ne];k++)
        sediments->boundary_sediC[i][j][k]=0; // boundary_sediC will be set using BoundaryScalars function  
    }
  }
  
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i]+1;k++){
      sediments->Wnewsedi[i][k]=0;
      //Woldsedi[i][k]=0;
      if(k!=grid->Nk[i])
        sediments->SediKappa_tv[i][k]=0;
    }
  }

  // give 0 for SediCbed
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++){
      sediments->Deposition[i][j]=0;      
      //sediments->Deposition_old[i][j]=0;
      for(k=0;k<sediments->Nlayer;k++){
        sediments->SediCbed[i][j][k]=0;
        //sediments->Erosion[i][j][k]=0;
        //sediments->Erosion_old[i][j][k]=0;
      }
    }

  //give 0 for Erosion, Deposition, Layermass
  for(i=0;i<grid->Nc;i++) {
    for(j=0;j<sediments->Nlayer;j++){
      //Layermass[i][j]=0;
      sediments->Layerthickness[i][j]=0;
      sediments->Erosiontotal[i][j]=0;
      sediments->Erosiontotal_old[i][j]=0;
    }
  }
  
  // get Nsize and Nlayer sediment properties when Nsize>3 or Nlayer>5, otherwise they are set in ReadSediProperties
  /*if(SETsediment==1){
    SetSediment(grid,myproc);
  }*/
  
  // set first tb
  if(sediments->TBMAX==1)
    for(i=0;i<grid->Nc;i++) {
      sediments->Seditb[i]=0;
      sediments->Seditbmax[i]=0;
    }
  
  // give layer mass for each cell
  for(i=0;i<grid->Nc;i++){
    for(j=0;j<sediments->Nlayer;j++){
      //sediments->Layermass[i][j]=sediments->Drydensity[j]*sediments->Thickness[j]*grid->Ac[i];
      sediments->Layerthickness[i][j]=sediments->Thickness[j];//*sediments->Bedmudratio[j]; // consider there are not all mud in each layer
    }
  }

  // read or use function to set initial condition for SediC
  if(sediments->readSediment) { //it means Nsize==1
    MPI_GetFile(filename,DATAFILE,"InitSedimentFile","InitializeSediment",myproc);
    InitSedimentFID = MPI_FOpen(filename,"r","InitializeSediment",myproc);
    stmp = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializeSediment");
    if(fread(stmp,sizeof(REAL),grid->Nkmax,InitSedimentFID) != grid->Nkmax)
      printf("Error reading stmp second\n");
    fclose(InitSedimentFID);    
    
    for(i=0;i<grid->Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        sediments->SediC[0][i][k]=stmp[k];
    
    SunFree(stmp,grid->Nkmax,"InitializeSediment");
  } else {
    for(j=0;j<sediments->Nsize;j++){
      for(i=0;i<grid->Nc;i++) {
        z = 0;
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          z-=grid->dz[k]/2;
          sediments->SediC[j][i][k]=ReturnSediment(grid->xv[i],grid->yv[i],z,j); // add new functions in Initialization.c
          z-=grid->dz[k]/2;
        }
      }
    }
  }

  // give the initial value for SediCbed
  for(j=0;j<sediments->Nsize;j++)
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<sediments->Nlayer;k++){
        sediments->SediCbed[j][i][k]=sediments->Drydensity[k]*ReturnBedSedimentRatio(grid->xv[i],grid->yv[i],k,j,sediments->Nsize);
      }

  // calculate initial settling velocity 
  if(sediments->WSconstant==1){
    for(i=0;i<sediments->Nsize;i++)
      for(j=0;j<grid->Nc;j++)
        for(k=0;k<grid->Nk[j]+1;k++)
          sediments->Ws[i][j][k]=sediments->Ws0[i];
  } else {
    SettlingVelocity(grid,phys,prop,myproc);
  }
  // calculate initial Deposition & deposition
  CalculateErosion(grid,phys,prop,myproc);
  CalculateDeposition(grid,phys,myproc);
}

/*
 * Function: SetSediment
 * Usage: give value for sediment properties when 
 * Nsize>3 or Nlayer>5
 * ----------------------------------------------------
 * You can set any value you like
 *
 */
/*void SetSediment(gridT *grid, int myproc) { 
  int i;
  // Properties for each fraction of sediment, you can set any value
  for(i=0;i<Nsize;i++){
    Ds[i]=0.0002*i;
    Ws0[i]=0.002*i;
    Gsedi[i]=2.65;
    Prt[i]=1;
  }
  
  // Properties for each bed layer 
  for(i=0;i<Nlayer;i++){
    Consolid[i]=0.0002;
    E0[i]=0.01;
    Taue[i]=0.1;
    Taud[i]=0.1;
    Drydensity[i]=530000;
    Thickness[i]=0.5*i;
    Softhard[i]=0;
    //Bedmudratio[i]=0.8;
  }
}*/

/*
 * Function: FreeSediment
 * Usage: free space for all the variables
 * ----------------------------------------------------
 * Basic sunfree function
 *
 */
void FreeSediment(gridT *grid, int myproc) {
  int i,j;

  for(i=0;i<sediments->Nsize;i++){
    for(j=0;j<grid->Nc;j++){
      free(sediments->SediC[i][j]);
      free(sediments->Ws[i][j]);
      free(sediments->SediCbed[i][j]); 
      //free(sediments->Erosion[i][j]);
      //free(sediments->Erosion_old[i][j]);
    }
    free(sediments->SediC[i]);
    free(sediments->SediCbed[i]);
    free(sediments->Ws[i]);
    //free(sediments->Erosion[i]);
    //free(sediments->Erosion_old);
    free(sediments->Deposition[i]);
    //free(sediments->Deposition_old[i]);
  }
  free(sediments->SediC);
  free(sediments->SediCbed);
  //free(sediments->Erosion);
  //free(sediments->Erosion_old);
  free(sediments->Ws);
  free(sediments->Deposition);
  //free(sediments->Deposition_old);
  for(i=0;i<grid->Nc;i++){
    //free(sediments->Layermass[i]);
    free(sediments->Wnewsedi[i]);
    //free(sediments->Woldsedi[i]);
    free(sediments->Layerthickness[i]);
    free(sediments->SediKappa_tv[i]);
    free(sediments->Erosiontotal[i]);
    free(sediments->Erosiontotal_old[i]);
  }
  free(sediments->SediKappa_tv);
  free(sediments->Erosiontotal);
  free(sediments->Erosiontotal_old);
  free(sediments->Wnewsedi);
  //free(sediments->Woldsedi);
  //free(sediments->Layermass);
  free(sediments->Layerthickness);
  free(sediments->Ds);
  free(sediments->Ws0);
  free(sediments->Gsedi);
  free(sediments->Prt);
  free(sediments->Consolid);
  free(sediments->E0);
  free(sediments->Taue);
  free(sediments->Taud);
  free(sediments->Drydensity);
  free(sediments->Thickness);
  free(sediments->Softhard);
  //free(sediments->Bedmudratio);
  if(sediments->TBMAX==1){
    free(sediments->Seditb);
    free(sediments->Seditbmax);
  }
}

/*
 * Function: SettlingVelocity
 * this function still need to be tested
 * Usage: calculate flocculation and hindered settling velocity
 * ----------------------------------------------------
 * Richarson and Zaki(1954) equation
 * only work when WSconstant=0
 *
 */
void SettlingVelocity(gridT *grid, physT *phys, propT *prop,  int myproc) { 
  int i,j,k;
  REAL tmp,sum,Cgel=1800000,wsn=2.0;

  // roughly estimate by Stokes Law
  for(i=0;i<sediments->Nsize;i++)
    sediments->Ws0[i] = (sediments->Gsedi[i]-1)*prop->grav*pow(sediments->Ds[i], 2)/(18*prop->nu);
  // coefficient for Richarson and Zaki(1954) equation
  //Cgel = 1800000;
  //wsn = 2.0;

  for(j=0;j<grid->Nc;j++){  
    for(k=0;k<grid->Nk[j]+1;k++){   
      sum=0;
      if(k==0){
        for(i=0;i<sediments->Nsize;i++)
          sum+=sediments->SediC[i][j][k];
      } else if(k==grid->Nk[j]){
        for(i=0;i<sediments->Nsize;i++)
          sum+=sediments->SediC[i][j][k-1];
      } else {
        for(i=0;i<sediments->Nsize;i++)
          sum+=0.5*(sediments->SediC[i][j][k-1]+sediments->SediC[i][j][k+1]);
      }
      //for(i=0;i<sediments->Nsize;i++)
      //sum+=0.5*(sediments->SediC[i][j][k-1]+sediments->SediC[i][j][k+1]);
      for(i=0;i<sediments->Nsize;i++){
        if(sum < sediments->Cfloc)
          sediments->Ws[i][j][k] = sediments->Ws0[i];
        else if(sum>= sediments->Cfloc && sum < sediments->Chind)
          sediments->Ws[i][j][k] =sediments-> Ws0[i]*sum/sediments->Gsedi[i]/1000/1000;  ///??????????
        else if(sum>= sediments->Chind && sum < Cgel){
          tmp = Min(1,sum/Cgel);
          tmp = 1-tmp;
          sediments->Ws[i][j][k] = sediments->Ws0[i]*sediments->Chind/sediments->Gsedi[i]/1000/1000*pow(tmp,wsn);
        } else {
          sediments->Ws[i][j][k] = 0;
          //sediments->SediC[i][j][k] = Cgel;
        }
      }
    }
    /*for(i=0;i<Nsize;i++){
      sediments->Ws[i][j][grid->Nk[j]] = sediments->Ws[i][j][grid->Nk[i]-1];
      sediments->Ws[i][j][0]=Ws[i][j][1];
      }*/
  }
}

/*
 * Function: CalculateErosion
 * Usage: calculate Erosion Rate for suspended sediment transport
 *        calculate Deposition for n time step
 * ----------------------------------------------------
 * for two types, hard and soft
 * coefficients based on Mike21
 */
void CalculateErosion(gridT *grid, physT *phys, propT *prop, int myproc) { 
  int j,k;
  REAL taub,utmp,taubtmp1,taubtmp2,em,alpha,ratio, depotmp, nettmp, erosionmax,ds90;
  // assume ks=3*ds90
  ds90=MPI_GetValue(DATAFILE,"Ds90","CalculationErosion",myproc);
  // coefficient for hard erosion
  em=1.0;
  // coefficient for soft erosion
  alpha=4.2;  // 4.2 to 25.6 according to MIKE21

  for(j=0;j<grid->Nc;j++){
    // calculate taub
    //taub = prop->CdB*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*(pow(phys->uc[j][grid->Nk[j]-1], 2)+pow(phys->vc[j][grid->Nk[j]-1], 2));
    utmp=sqrt(pow(phys->uc[j][grid->Nk[j]-1], 2)+pow(phys->vc[j][grid->Nk[j]-1], 2));
    taub=pow(log(10*0.5*grid->dzz[j][grid->Nk[j]-1]/ds90)/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);

    if((phys->h[j]+grid->dv[j])/grid->dv[j]<0.001)
      taub=0;
    if(sediments->TBMAX==1){
      sediments->Seditb[j]=taub;
    if(sediments->Seditbmax[j]<=taub)
      sediments->Seditbmax[j]=taub;
    }
    erosionmax=0;
    for(k=0;k<sediments->Nlayer;k++){
      erosionmax=sediments->Layerthickness[j][k]*sediments->Drydensity[k]/prop->dt;
      sediments->Erosiontotal_old[j][k]=sediments->Erosiontotal[j][k];
      sediments->Erosiontotal[j][k]=0;
      if(taub>sediments->Taue[k]){
        if(sediments->Softhard[k]==0){
	  // soft erosion
          sediments->Erosiontotal[j][k]=sediments->E0[k]*exp(alpha*(taub-sediments->Taue[k]));
        } else {
          // hard erosion
          ratio=Max(taub/sediments->Taue[k]-1,0);
          sediments->Erosiontotal[j][k]=sediments->E0[k]*pow(ratio,em);
        }
        if(sediments->Erosiontotal[j][k]>erosionmax)
          sediments->Erosiontotal[j][k]=erosionmax;
      }     
    }
    if(sediments->bedComplex==1){
      nettmp=0;
      //assume last layer cannot be flushed away
      for(k=sediments->Nlayer-2;k>0;k--){
        nettmp=sediments->Layerthickness[j][k]*sediments->Drydensity[k]/prop->dt-sediments->Erosiontotal[j][k]+sediments->Erosiontotal[j][k+1]+sediments->Consolid[k-1]-sediments->Consolid[k]; //here can use theta method.
        if(nettmp<0)
          sediments->Erosiontotal[j][k]+=nettmp;
      }
    }
  }
  
  // calculate deposition for each fraction
  if(sediments->bedComplex==1){
    nettmp=0;
    for(j=0;j<grid->Nc;j++){
      depotmp=0;
      for(k=0;k<sediments->Nsize;k++){
        /*if(grid->Nk[j]-1>0){
	// here deposition is explicit n time steps while Erosion is n+1 time steps
	sediments->Deposition[k][j]=(1.5*sediments->SediC[k][j][grid->Nk[j]-1]-0.5*sediments->SediC[k][j][grid->Nk[j]-2])*(1.5*sediments->Ws[k][j][grid->Nk[j]-1]-0.5*sediments->Ws[k][j][grid->Nk[j]-2]);
        } else {
	sediments->Deposition[k][j]=sediments->SediC[k][j][grid->Nk[j]-1]*sediments->Ws[k][j][grid->Nk[j]-1];
        }*/
        depotmp+=sediments->Deposition[k][j];
      }
      if(sediments->Nlayer>1)
        nettmp=sediments->Layerthickness[j][0]*sediments->Drydensity[0]/prop->dt-sediments->Erosiontotal[j][0]+sediments->Erosiontotal[j][1]+depotmp-sediments->Consolid[1];
      if(nettmp<0)   
        sediments->Erosiontotal[j][0]+=nettmp;
    }
  }
}

/*
 * Function: CalculateDeposition
 * Usage: calculate deposition
 * ----------------------------------------------------
 * based on settling velocity and SediConcentration
 */
void CalculateDeposition(gridT *grid, physT *phys, int myproc) {
  int j,k;
  for(j=0;j<grid->Nc;j++){
    for(k=0;k<sediments->Nsize;k++){
      if(grid->Nk[j]+1-1>0){
        // here deposition is explicit n time steps while Erosion is n+1 time steps
        sediments->Deposition[k][j]=(1.5*sediments->SediC[k][j][grid->Nk[j]-1]-0.5*sediments->SediC[k][j][grid->Nk[j]-2])*(sediments->Ws[k][j][grid->Nk[j]]);
      } else {
        sediments->Deposition[k][j]=sediments->SediC[k][j][grid->Nk[j]-1]*sediments->Ws[k][j][grid->Nk[j]];
      }
    }
  }
}

/*
 * Function: BedChange
 * Usage: calculate Cbed and Thickness change
 * ----------------------------------------------------
 * update by intervals, calculate Cbed fraction to 
 * get erosion for each fraction
 * assume dry density is constant for each layer
 * assume the last layer cannot be flushed away 
 * when bedComplex=1, Bedchange will be calculated every time step
 * there may some error in extreme condition
 */
void BedChange(gridT *grid, physT *phys, propT *prop, int myproc) {
  int i,j,k;
  REAL depotmp, thicktmp;

  // update Deposition first
  CalculateDeposition(grid,phys,myproc);

  // update thickness
  // assume sum(Cbed) is constant=Drydensity
  for(j=0;j<grid->Nc;j++){
    depotmp=0;
    for(k=0;k<sediments->Nsize;k++)
      depotmp+=sediments->Deposition[k][j];
    // update thickness and Cbed
    // for the first layer assume Deposition can only go into 1st layer
    thicktmp=sediments->Layerthickness[j][0];
    if(sediments->Nlayer>1){
      // here can use theta method
      sediments->Layerthickness[j][0]+=sediments->bedInterval*prop->dt/sediments->Drydensity[0]*(-sediments->Erosiontotal[j][0]+sediments->Erosiontotal[j][1]+depotmp-sediments->Consolid[0]);
    } else {
      // if there is one layer, no consolidation and erosion[j][1];
      // because here we use deposition n+1 step, so there may be some error to make Layerthickness<0
      sediments->Layerthickness[j][0]+=sediments->bedInterval*prop->dt/sediments->Drydensity[0]*(-sediments->Erosiontotal[j][0]+depotmp); 
    }
    
    if(sediments->Layerthickness[j][0]<=0){
      sediments->Layerthickness[j][0]=0;
      for(k=0;k<sediments->Nsize;k++)
	sediments->SediCbed[k][j][0]=0;
    } else {
      for(k=0;k<sediments->Nsize;k++){
	if(sediments->Nlayer>1){
	  sediments->SediCbed[k][j][0]=(sediments->SediCbed[k][j][0]*thicktmp+sediments->bedInterval*prop->dt*(-sediments->SediCbed[k][j][0]/sediments->Drydensity[0]*(sediments->Consolid[0]+sediments->Erosiontotal[j][0])+sediments->SediCbed[k][j][1]/sediments->Drydensity[1]*sediments->Erosiontotal[j][1]+sediments->Deposition[k][j]))/sediments->Layerthickness[j][0];
	} else {
	  sediments->SediCbed[k][j][0]=(sediments->SediCbed[k][j][0]*thicktmp+sediments->bedInterval*prop->dt*(-sediments->SediCbed[k][j][0]/sediments->Drydensity[0]*(sediments->Erosiontotal[j][0])+sediments->Deposition[k][j]))/sediments->Layerthickness[j][0];
	}
      }
    }
    
    for(k=1;k<sediments->Nlayer-1;k++){
      thicktmp=sediments->Layerthickness[j][k];
      sediments->Layerthickness[j][k]+=sediments->bedInterval*prop->dt/sediments->Drydensity[k]*(-sediments->Erosiontotal[j][k]+sediments->Erosiontotal[j][k+1]+sediments->Consolid[k-1]-sediments->Consolid[k]);
      if(sediments->Layerthickness[j][k]<=0){
        sediments->Layerthickness[j][k]=0;
        for(i=0;i<sediments->Nsize;i++)
          sediments->SediCbed[i][j][k]=0;
      } else {
        for(i=0;i<sediments->Nsize;i++){
          sediments->SediCbed[i][j][k]=(sediments->SediCbed[i][j][k]*thicktmp+sediments->bedInterval*prop->dt*(-sediments->SediCbed[i][j][k]/sediments->Drydensity[k]*(sediments->Consolid[k]+sediments->Erosiontotal[j][k])+sediments->SediCbed[i][j][k+1]/sediments->Drydensity[k+1]*sediments->Erosiontotal[j][k+1]+sediments->Consolid[k-1]))/sediments->Layerthickness[j][k];
        }
      }            
    }
  }
}

/*
 * Function: SedimentSource
 * Usage: dT/dt + u dot grad T = d/dz ( kappaT dT/dz) + A + B*T
 *--------------------------------------------------------------
 * same use as Heatsource, add new parameter Nosize to represent the No. of Nsize for fraction
 *
 */
void SedimentSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop,int Nosize,REAL theta) {
  int i, k, layertop;
  REAL erosion, erosion_old;

  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i]-1;k++) {
      B[i][k]=0;
      A[i][k]=0;
    }
    //advection
    B[i][grid->Nk[i]-1]=sediments->Ws[Nosize][i][grid->Nk[i]]/grid->dzz[i][k]; 
    //erosion
    layertop=0;
    erosion=0;
    erosion_old=0;
    if(sediments->Nlayer>1){
      while(sediments->Layerthickness[i][layertop]==0 && layertop<sediments->Nlayer){
	erosion+=sediments->Erosiontotal[i][layertop];
	erosion_old+=sediments->Erosiontotal_old[i][layertop];
	layertop++;
      }
    }
    erosion+=sediments->Erosiontotal[i][layertop];
    erosion_old+=sediments->Erosiontotal_old[i][layertop];
    A[i][grid->Nk[i]-1]=((1-theta)*erosion_old+theta*erosion)*sediments->SediCbed[Nosize][i][layertop]/sediments->Drydensity[layertop]/grid->dzzold[i][grid->Nk[i]-1];
  }
}

/*
 * Function: SedimentVerticalVelocity
 * Usage: wnewsedi=phys->wnew-ws, woldsedi=phys->wtmp2-ws, ws is explicit
 *--------------------------------------------------------------
 * provide the vertical velocity field for updatescalar function
 *
 */
void SedimentVerticalVelocity(gridT *grid, physT *phys,int Nosize,int symbol, int myproc) {
  int i,k;

  if(symbol==1)
    for(i=0;i<grid->Nc;i++) 
      for(k=0;k<grid->Nk[i]+1;k++) {
        sediments->Wnewsedi[i][k]=phys->wnew[i][k]-sediments->Ws[Nosize][i][k];
        phys->wtmp2[i][k]=phys->wtmp2[i][k]-sediments->Ws[Nosize][i][k];
      }
  else
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i]+1;k++)
        phys->wtmp2[i][k]=phys->wtmp2[i][k]+sediments->Ws[Nosize][i][k];
}

/*
 * Function: OpenSediFiles
 * read sediment output file names and open sediment files
 *--------------------------------------------------------
 *
 */
void OpenSediFiles(propT *prop, int myproc) {
  int i;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  sediments->SedimentFID=(FILE **)SunMalloc(sediments->Nsize*sizeof(FILE *), "OpenSediFiles"); 
  for(i=0;i<sediments->Nsize;i++){
    sprintf(str,"Sediment%dFile",i+1);
    MPI_GetFile(filename,DATAFILE,str,"OpenSediFiles",myproc);

    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    sediments->SedimentFID[i] = MPI_FOpen(str,"w","OpenFiles",myproc); 
  }
  MPI_GetFile(filename,DATAFILE,"LayerFile","OpenSediFiles",myproc);
  if(prop->mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  sediments->LayerthickFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  
  if(sediments->TBMAX==1) {
    MPI_GetFile(filename,DATAFILE,"tbFile","OpenSediFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    sediments->SeditbFID = MPI_FOpen(str,"w","OpenFiles",myproc);

    MPI_GetFile(filename,DATAFILE,"tbmaxFile","OpenSediFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    sediments->SeditbmaxFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  }
}

/*
 * Function: OutputSediment
 * output sediment data for each time step
 *--------------------------------------------------------------
 *
 */
void OutputSediment(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, int blowup, MPI_Comm comm) {
  int i, j, jptr, k, nwritten, nosize, nolayer;
  char str[BUFFERLENGTH];
  FILE *ofile;
  REAL thicktmp;

  if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) {
    
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      sprintf(str,"Error outputting SSC data for size class %d of %d.\n",nosize+1,sediments->Nsize);
      Write3DData(sediments->SediC[nosize],phys->htmp,prop->mergeArrays,sediments->SedimentFID[nosize],
		  str,grid,numprocs,myproc,comm);
    }
    
    for(i=0;i<grid->Nc;i++){
      phys->htmp[i]=0;
      for(nolayer=0;nolayer<sediments->Nlayer;nolayer++)
	phys->htmp[i]+=sediments->Layerthickness[i][nolayer];
    }
    Write2DData(phys->htmp,prop->mergeArrays,sediments->LayerthickFID,"Error outputting bed thickness data!\n",
    		grid,numprocs,myproc,comm);
    
    if(sediments->TBMAX==1) {
      Write2DData(sediments->Seditb,prop->mergeArrays,sediments->SeditbFID,"Error outputting bed shear stress data!\n",
		  grid,numprocs,myproc,comm);
    }
  }
  
  if(prop->n==prop->nsteps+prop->nstart) {
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      fclose(sediments->SedimentFID[nosize]);
    }
    fclose(sediments->LayerthickFID);
    
    if(sediments->TBMAX==1) {
      fclose(sediments->SeditbFID);

      Write2DData(sediments->Seditbmax,prop->mergeArrays,sediments->SeditbmaxFID,"Error outputting max bed shear stress data!\n",
		  grid,numprocs,myproc,comm);    
      fclose(sediments->SeditbmaxFID);
    }
  }
}

/*
 * Function: CalculateSediDiffusivity
 * calculate sediment diffusivity 
 *--------------------------------------------------------------
 * may use phys->kappa_tv directly from my25 function or 
 * use parabolic diffusivity
 */
void CalculateSediDiffusivity(gridT *grid, physT *phys,int Nosize,int myproc) {  
  int ii, kk;
  REAL z;
  
  if(sediments->ParabolKappa==1){
    for(ii=0;ii<grid->Nc;ii++){
      z=grid->dv[ii]+phys->h[ii]-0.5*grid->dzz[ii][0];
      for(kk=0;kk<grid->Nk[ii];kk++){
	if(phys->CdB[ii]!=-1)
	  sediments->SediKappa_tv[ii][kk]=z*(1-z/(grid->dv[ii]+phys->h[ii]))*phys->uc[ii][grid->Nk[ii]-1]*sqrt(phys->CdB[ii])*0.41/sediments->Prt[Nosize];
	if(kk!=grid->Nk[ii]-1)
	  z-=0.5*(grid->dzz[ii][kk]+grid->dzz[ii][kk+1]);
      }
    }
  } else {
    for(ii=0;ii<grid->Nc;ii++)
      for(kk=0;kk<grid->Nk[ii];kk++)
	sediments->SediKappa_tv[ii][kk]=phys->kappa_tv[ii][kk]/sediments->Prt[Nosize];
  }
}

/*
 * Function: ISendRecvSediBedData3D
 * Usage: ISendRecvSediBedData3D(SediCbed or Thickness,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the 3D cell data for sediment bed back and forth between
 * processors using nonblocking sends/recvs.
 *
 */
void ISendRecvSediBedData3D(REAL **celldata, gridT *grid, int nlayer,int myproc,MPI_Comm comm) {
  int k, n, nstart, neigh, neighproc;
  REAL t0=Timer();
  
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    
    nstart=0;
    for(n=0;n<grid->num_cells_send[neigh];n++) {
      for(k=0;k<nlayer;k++) 
        grid->send[neigh][nstart+k]=celldata[grid->cell_send[neigh][n]][k];
      nstart+=nlayer;
    }
    
    MPI_Isend((void *)(grid->send[neigh]),grid->total_cells_send[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[neigh])); 
  }
  
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Irecv((void *)(grid->recv[neigh]),grid->total_cells_recv[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[grid->Nneighs+neigh]));
  }
  MPI_Waitall(2*grid->Nneighs,grid->request,grid->status);
  
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<grid->num_cells_recv[neigh];n++) {
      for(k=0;k<nlayer;k++) 
        celldata[grid->cell_recv[neigh][n]][k]=grid->recv[neigh][nstart+k];
      nstart+=nlayer;
    }
  }
  // t_comm+=Timer()-t0;
}

/*
 * Function: ComputeSediments
 * Usage: ComputeSediments(grid,phys,prop, myproc, numproc, blowup, comm);
 * ----------------------------------------------------
 * This function is the main function for sediment model part
 * include all the calculation for sediment transport
 * and called by phys.c
 *
*/
void ComputeSediments(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int k;
  if(prop->n==1+prop->nstart){
    sediments=(sedimentsT *)SunMalloc(sizeof(sedimentsT),"ComputeSediments");
    // allocate and initialize all the sediment variables
    ReadSediProperties(myproc); 
    OpenSediFiles(prop,myproc); 
    AllocateSediment(grid,myproc);  
    InitializeSediment(grid,phys,prop,myproc);
    BoundarySediment(grid,phys,prop);
  }
  
  // calculate n+theta Erosion for boundary 
  CalculateErosion(grid,phys,prop,myproc);
  // calculate n+1 Sediment concentration field
  for(k=0;k<sediments->Nsize;k++){
    
    SedimentSource(phys->wtmp,phys->uold,grid,phys,prop,k,prop->theta);
    
    SedimentVerticalVelocity(grid,phys,k,1,myproc);
    
    CalculateSediDiffusivity(grid,phys,k,myproc);
    
    UpdateScalars(grid,phys,prop,sediments->Wnewsedi,sediments->SediC[k],sediments->boundary_sediC[k],phys->Cn_T,0,0,sediments->SediKappa_tv,prop->theta,phys->uold,phys->wtmp,NULL,NULL,0,0,comm,myproc,0,prop->TVDtemp);
    SedimentVerticalVelocity(grid,phys,k,-1,myproc);
    ISendRecvCellData3D(sediments->SediC[k],grid,myproc,comm);
  }          
  
  if(sediments->WSconstant==0)
    SettlingVelocity(grid,phys,prop,myproc);
  if(prop->n%sediments->bedInterval==0 && sediments->bedInterval>0)
    BedChange(grid,phys,prop,myproc);   
  // get the boundary value for the next time step
  BoundarySediment(grid,phys,prop);
  // output sediment results
  OutputSediment(grid,phys,prop,myproc,numprocs,blowup,comm);
  // free space
  //if(prop->n==prop->nstart+prop->nsteps)
  //FreeSediment(grid,sediments,myproc);
}

