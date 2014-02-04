/*
 * File: phys.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for phys.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _phys_h
#define _phys_h

#include "suntans.h"
#include "grid.h"
#include "fileio.h"

/*
 * Enumerated type definitions
 *
 */

// effectively flags for use in interpolation schemes
typedef enum _interpolation {
  nRT1, nRT2,
  tRT1, tRT2,
  QUAD, PEROT, LSQ
} interpolation;

/*
 * Main physical variable struct.
 *
 */
typedef struct _physT {
  REAL **u;
  REAL **uc;
  REAL **vc;
  REAL **wc, **wf;

  /*  new variables for nodal and tangential velocities */
  // definitions follow from Wang et al 2011
  // nRT1[Np][Nk][numpcneighs] so that for each node there is a 
  // value for each cell neighbor (non-unique to a node) and 
  // note that the number of cell neighbors varies based on node
  REAL ***nRT1u;
  REAL ***nRT1v;
  // nRT2[Np][Nk] has a unique value for each node since it is 
  // the area-weigted average of all the nRT1 values 
  // around the node
  REAL **nRT2u;
  REAL **nRT2v;
  // tRT1[Ne][Nk] is the area-weighted average of the nRT1 values
  // for the neighboring cells of an edge
  REAL **tRT1;
  // tRT2[Ne][Nk] is the area-weighted average of the nRT2 values
  // for each node at the end of the edge
  REAL **tRT2;

  REAL **uold;
  REAL **vold;
  REAL *D;
  REAL **w;
  REAL **q;
  REAL **qc;
  REAL **s;
  REAL **T;
  REAL **s0;
  REAL **rho;
  REAL *h;
  REAL *hcorr;
  unsigned char *active;

  REAL **boundary_u;
  REAL **boundary_v;
  REAL **boundary_w;
  REAL **boundary_s;
  REAL **boundary_T;
  REAL **boundary_tmp;
  REAL **boundary_rho;
  REAL *boundary_h;
  REAL *boundary_flag;

  REAL **nu_tv;
  REAL **kappa_tv;
  REAL **nu_lax;
  REAL *tau_T;
  REAL *tau_B;
  REAL *CdT;
  REAL *CdB;
  REAL **qT;
  REAL **lT;

  REAL mass;
  REAL mass0;
  REAL volume;
  REAL volume0;
  REAL Ep;
  REAL Ep0;
  REAL Ek;
  REAL Eflux1;
  REAL Eflux2;
  REAL Eflux3;
  REAL Eflux4;
  REAL smin;
  REAL smax;

  REAL *htmp;
  REAL *hold;
  REAL *htmp2;
  REAL *htmp3;
  REAL *hcoef;
  REAL *hfcoef;
  REAL **stmp;
  REAL **stmp2;
  REAL **stmp3;
  REAL **utmp;
  REAL **utmp2;
  REAL **ut;
  REAL **Cn_R;
  REAL **Cn_T;
  REAL **Cn_U;
  REAL **Cn_U2; //AB3
  REAL **Cn_W;
  REAL **Cn_W2; //AB3
  REAL **Cn_q;
  REAL **Cn_l;
  REAL **wnew;
  REAL **wtmp;
  REAL **wtmp2;
  REAL **qtmp;
  //GLS Turbulence variables
  REAL **TP;
  REAL **TB;
  REAL **TD;

  REAL *ap;
  REAL *am;
  REAL *bp;
  REAL *bm;

  REAL *a;
  REAL *b;
  REAL *c;
  REAL *d;

  // Horizontal facial scalar
  REAL **SfHp;
  REAL **SfHm;

  //Define variables for TVD schemes
  REAL *Cp;
  REAL *Cm;
  REAL *rp;
  REAL *rm;
  REAL *wp;
  REAL *wm;
  REAL **gradSx;
  REAL **gradSy; 
  
  //Variables for heat flux model
  REAL *Tsurf;
  REAL *dT;
  REAL **Ttmp;

  // Variables for netcdf write
  REAL *tmpvar;
  REAL *tmpvarW;
  REAL *tmpvarE;

  // Least squares fit arrays
  REAL **A;
  REAL **Apr;
  REAL **AT;
  REAL *bpr;

} physT;

/*
 * Main property struct.
 *
 */
typedef struct _propT {
  REAL dt, Cmax, rtime, amp, omega, flux, timescale, theta0, theta, thetaM, 
       thetaS, thetaB, nu, nu_H, tau_T, z0T, CdT, z0B, CdB, CdW, relax, epsilon, qepsilon, resnorm, 
       dzsmall, beta, kappa_s, kappa_sH, gamma, kappa_T, kappa_TH, grav, Coriolis_f, CmaxU, CmaxW, 
       laxWendroff_Vertical, latitude;
  int ntout, ntoutStore, ntprog, nsteps, nstart, n, ntconserve, nonhydrostatic, cgsolver, maxiters, 
      qmaxiters, hprecond, qprecond, volcheck, masscheck, nonlinear, linearFS, newcells, wetdry, sponge_distance, 
    sponge_decay, thetaramptime, readSalinity, readTemperature, turbmodel, 
    TVD, horiTVD, vertTVD, TVDsalt, TVDtemp, TVDturb, laxWendroff, stairstep, AB, TVDmomentum, conserveMomentum,
    mergeArrays, computeSediments;
  FILE *FreeSurfaceFID, *HorizontalVelocityFID, *VerticalVelocityFID, *SalinityFID, *BGSalinityFID, 
       *InitSalinityFID, *InitTemperatureFID, *TemperatureFID, *PressureFID, *VerticalGridFID, *ConserveFID,    
       *StoreFID, *StartFID, *EddyViscosityFID, *ScalarDiffusivityFID; 
  interpolation interp; int prettyplot;
  int metmodel,  varmodel, outputNetcdf,  metncid, netcdfBdy, netcdfBdyFileID, readinitialnc, initialNCfileID, calcage, calcaverage;
  int outputNetcdfFileID, averageNetcdfFileID;
  REAL nctime, toffSet;
  int nctimectr, avgtimectr, avgctr, ntaverage, nstepsperncfile, ncfilectr;
  REAL nugget, sill, range, Lsw, Cda, Ce, Ch;
  char  starttime[15], basetime[15]; 
} propT;


/*
 * Public function declarations.
 *
 */
void Solve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void AllocatePhysicalVariables(gridT *grid, physT **phys, propT *prop);
void FreePhysicalVariables(gridT *grid, physT *phys, propT *prop);
void InitializePhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void InitializeVerticalGrid(gridT **grid,int myproc);
void ReadProperties(propT **prop, gridT *grid, int myproc);
void SetDragCoefficients(gridT *grid, physT *phys, propT *prop);
REAL DepthFromDZ(gridT *grid, physT *phys, int i, int kind);
REAL InterpToFace(int j, int k, REAL **phi, REAL **u, gridT *grid);
inline void ComputeUC(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc, interpolation interp) ;
void UpdateDZ(gridT *grid, physT *phys, propT *prop, int option);
void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void SetDensity(gridT *grid, physT *phys, propT *prop);

#endif
