/*
 * SUNTANS meteorological data
 * 
*/ 
#ifndef _met_h
#define _met_h

#include "suntans.h"
#include "fileio.h"
#include "memory.h"
#include "grid.h"
#include "phys.h"
#include "util.h"

/* Structure array for meteorological input data*/
typedef struct _metinT {
  
  // Number of points for each variable
  size_t NUwind;
  size_t NVwind;
  size_t NTair;
  size_t NPair;
  size_t Nrain;
  size_t NRH;
  size_t Ncloud;
  
  
  // Coordinate vectors for each variable
  REAL *x_Uwind;
  REAL *y_Uwind;
  REAL *x_Vwind;
  REAL *y_Vwind;
  REAL *x_Tair;
  REAL *y_Tair;
  REAL *x_Pair;
  REAL *y_Pair;
  REAL *x_rain;
  REAL *y_rain;
  REAL *x_RH;
  REAL *y_RH;
  REAL *x_cloud;
  REAL *y_cloud;
  
  // Height vectors
  REAL *z_Uwind;
  REAL *z_Vwind;
  REAL *z_Tair;
  REAL *z_RH;
  
  // Time records
  size_t nt; //Time
  REAL *time;
  
  // Time record locators
  int t0;
  int t1;
  
  // Interpolation Weights
  REAL **WUwind;
  REAL **WVwind;
  REAL **WTair;
  REAL **WPair;
  REAL **Wrain;
  REAL **WRH;
  REAL **Wcloud;
  
  // The actual input data
  REAL **Uwind;
  REAL **Vwind;
  REAL **Tair;
  REAL **Pair;
  REAL **rain;
  REAL **RH;
  REAL **cloud;
  
} metinT;

/* Meteorological data on SUNTANS grid*/

typedef struct _metT {
  
  // Height vectors
  REAL *z_Uwind;
  REAL *z_Vwind;
  REAL *z_Tair;
  REAL *z_RH;
  
  // Data on suntans grid centres (model time step)
  REAL *Uwind;
  REAL *Vwind;
  REAL *Tair;
  REAL *Pair;
  REAL *rain;
  REAL *RH;
  REAL *cloud;
  
  // Computed variables at grid centres
  REAL *Hs;
  REAL *Hl;
  REAL *Hlw;
  REAL *Hsw;
  REAL *tau_x;
  REAL *tau_y;
  REAL *ustar;
  REAL *Tstar;
  REAL *qstar;
  REAL *EP;
  REAL *Htmp;
  REAL *xtmp;  

  // Data on suntans grid centres (two time steps)
  REAL **Uwind_t;
  REAL **Vwind_t;
  REAL **Tair_t;
  REAL **Pair_t;
  REAL **rain_t;
  REAL **RH_t;
  REAL **cloud_t;
} metT;

/* Public functions*/
#ifdef USENETCDF
  void InitialiseOutputNC(propT *prop, gridT *grid, physT *phys, metT *met, int myproc);
  void WriteOuputNC(propT *prop, gridT *grid, physT *phys, metT *met,int blowup, int myproc);
  void InitialiseMetFields(propT *prop, gridT *grid, metinT *metin, metT *met, int myproc);
  void updateMetData(propT *prop, gridT *grid, metinT *metin, metT *met, int myproc);
  void AllocateMet(propT *prop, gridT *grid, metT **met, int myproc);
  void AllocateMetIn(propT *prop, gridT *grid, metinT **metin, int myproc);
#endif
void updateAirSeaFluxes(propT *prop, gridT *grid, physT *phys, metT *met,REAL **T);
REAL shortwave(REAL time, REAL lat,REAL C_cloud);
#endif
