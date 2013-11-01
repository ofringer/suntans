/*
 * Averages header file
 *
 */

#ifndef _averages_h
#define _averages_h

#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "met.h"
#include "tvd.h"
#include "age.h"

/*
 * Main averages variable struct.
 */
typedef struct _averageT {
  
  REAL **uc;
  REAL **vc;
  REAL **w;
  REAL **s;
  REAL **T;
  REAL **rho;
  REAL **nu_v;

  REAL *h;

  // Flux variables
  REAL **U_F; // volume flux
  REAL **s_F; //Salt flux @ face
  REAL ** T_F; // Temperature flux @ faceL

  // Depth-integrated T/S (for budgets)
  REAL *s_dz;
  REAL *T_dz;

  // Age variables
  REAL **agemean;

  // Atmospheric flux variables
  REAL *Uwind;
  REAL *Vwind;
  REAL *Tair;
  REAL *Pair;
  REAL *rain;
  REAL *RH;
  REAL *cloud;
  
  REAL *Hs;
  REAL *Hl;
  REAL *Hlw;
  REAL *Hsw;
  REAL *tau_x;
  REAL *tau_y;
  REAL *EP;
 
  // Variables for netcdf write
  REAL *tmpvar;
  REAL *tmpvarW;
  REAL *tmpvarE;

} averageT;

/* *** Public Functions *** */
void AllocateAverageVariables(gridT *grid, averageT **average, propT *prop);
void ZeroAverageVariables(gridT *grid, averageT *average, propT *prop);
void UpdateAverageVariables(gridT *grid, averageT *average, physT *phys, metT *met, propT *prop, MPI_Comm comm, int myproc);
void ComputeAverageVariables(gridT *grid, averageT *average, physT *phys, metT *met, int netaverage, propT *prop);
void SendRecvAverages(propT *prop, gridT *grid, averageT *average, MPI_Comm comm, int myproc);

#endif
