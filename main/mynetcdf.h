/*
 * NetCDF IO header file
 *
 *
 */

#ifndef _mynetcdf_h
#define _mynetcdf_h

#ifndef NONETCDF
#include "netcdf.h"
#else
//Netcdf globals
#define NC_NOWRITE 0
#define NC_NETCDF4 0
#endif

#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "met.h"
#include "boundaries.h"
#include "averages.h"
#include "age.h"

/* Netcdf error */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); MPI_Finalize(); exit(ERRCODE);}
  
/* Public functions */
int getTimeRec(REAL nctime, REAL *time, int nt);
size_t returndimlen(int ncid, char *dimname);

void InitialiseOutputNC(propT *prop, gridT *grid, physT *phys, metT *met, int myproc);
void InitialiseOutputNCugrid(propT *prop, gridT *grid, physT *phys, metT *met, int myproc);
void WriteOutputNC(propT *prop, gridT *grid, physT *phys, metT *met,int blowup, int myproc);
void WriteOutputNCmerge(propT *prop, gridT *grid, physT *phys, metT *met, int blowup, int numprocs, int myproc, MPI_Comm comm);
void InitialiseAverageNCugrid(propT *prop, gridT *grid, averageT *average, int myproc);
void WriteAverageNC(propT *prop, gridT *grid, averageT *average, physT *phys, metT *met, int blowup, MPI_Comm comm, int myproc);
void WriteAverageNCmerge(propT *prop, gridT *grid, averageT *average, physT *phys, metT *met, int blowup, int numprocs, MPI_Comm comm, int myproc);
void ReadMetNCcoord(propT *prop, gridT *grid, metinT *metin,int myproc);
void ReadMetNC(propT *prop, gridT *grid, metinT *metin,int myproc);

void ReadBndNCcoord(int ncid, propT *prop, gridT *grid, int myproc, MPI_Comm comm);
void ReadBdyNC(propT *prop, gridT *grid, int myproc, MPI_Comm comm);
//void UpdateBdyNC(propT *prop, gridT *grid, int myproc,MPI_Comm comm);
size_t returndimlenBC(int ncid, char *dimname);
int getTimeRecBnd(REAL nctime, REAL *time, int nt);
void ReadInitialNCcoord(propT *prop, gridT *grid, int *Nci, int *Nki, int *T0, int myproc);
int getICtime(propT *prop, int Nt, int myproc);
void ReturnLatitudeNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int myproc);
void ReturnFreeSurfaceNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int T0, int myproc);
void ReturnTemperatureNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc);
void ReturnSalinityNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc);
void ReturnVelocityNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc);
void ReturnAgeNC(propT *prop, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc);
int MPI_NCOpen(char *file, int perms, char *caller, int myproc);
int MPI_NCClose(int ncid);
#endif
