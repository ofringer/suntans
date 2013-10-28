/*
* NetCDF IO dummy functions
* -------------------
* All functions either generic or specific involved in netcdf io should go here
* 
*/ 

//#include "mynetcdf-nonetcdf.h"
#include "mynetcdf.h"

#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "met.h"
#include "boundaries.h"
#include "averages.h"


int MPI_NCOpen(char *file, int perms, char *caller, int myproc){
    return -1;
}

int MPI_NCClose(int ncid){
    return -1;
}

int getTimeRec(REAL nctime, REAL *time, int nt){
    return -1;
}

size_t returndimlen(int ncid, char *dimname){
    return 0;
}

size_t returndimlenBC(int ncid, char *dimname){
    return 0;
}

int getTimeRecBnd(REAL nctime, REAL *time, int nt){
    return -1;
}

int getICtime(propT *prop, int Nt, int myproc){
    return -1;
}


void InitialiseOutputNC(propT *prop, gridT *grid, physT *phys, metT *met, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set outputNetcdf = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void InitialiseOutputNCugrid(propT *prop, gridT *grid, physT *phys, metT *met, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set outputNetcdf = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void WriteOuputNC(propT *prop, gridT *grid, physT *phys, metT *met,int blowup, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set outputNetcdf = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void InitialiseAverageNCugrid(propT *prop, gridT *grid, averageT *average, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set calcaverge = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void WriteAverageNC(propT *prop, gridT *grid, averageT *average, physT *phys, metT *met, int blowup, MPI_Comm comm, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set calcaverage = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReadMetNCcoord(propT *prop, gridT *grid, metinT *metin,int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set metmodel = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReadMetNC(propT *prop, gridT *grid, metinT *metin,int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set metmodel = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReadBndNCcoord(int ncid, propT *prop, gridT *grid, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set netcdfBdy = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReadBdyNC(propT *prop, gridT *grid, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set netcdfBdy = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}
/*
void UpdateBdyNC(propT *prop, gridT *grid, int myproc,MPI_Comm comm){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set netcdfBdy = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}
*/
void ReadInitialNCcoord(propT *prop, gridT *grid, int *Nci, int *Nki, int *T0, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set readinitialnc = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReturnFreeSurfaceNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int T0, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set readinitialnc = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReturnTemperatureNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set readinitialnc= 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReturnSalinityNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set readinitialnc = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void ReturnAgeNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){

  if(myproc==0) printf("Error: NetCDF Libraries required. Set readinitialnc = 0\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

