/*
 * Header file for phys.c
 *
 * $Id: phys.h,v 1.2 2002-11-03 20:36:33 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:19:17  fringer
 * Initial revision
 *
 *
 */
#ifndef _phys_h
#define _phys_h

#include "suntans.h"
#include "grid.h"
#include "fileio.h"

/*
 * Main physical variable struct.
 *
 */
typedef struct _physT {
  REAL **u;
  REAL *D;
  REAL **w;
  REAL **wf;
  REAL **q;
  REAL **s;
  REAL *h;

  REAL **nu_tv;
  REAL *tau_T;
  REAL *tau_B;
  REAL *CdT;
  REAL *CdB;

  REAL mass;
  REAL mass0;
  REAL volume;
  REAL volume0;
  REAL Ep;
  REAL Ep0;
  REAL smin;
  REAL smax;

  REAL *htmp;
  REAL *hold;
  REAL **stmp;
  REAL **utmp;
  REAL **utmp2;
  REAL **wtmp;

  REAL *ap;
  REAL *am;
  REAL *bp;
  REAL *bm;


  REAL *a;
  REAL *b;
  REAL *c;
  REAL *d;

} physT;

/*
 * Main property struct.
 *
 */
typedef struct _propT {
  REAL dt, rtime, amp, omega, flux, theta, thetaAB, 
    thetaFS, nu, tau_T, CdT, CdB, relax, epsilon, dzsmall;
  int ntout, ntprog, nsteps, n, ntconserve, maxiters, volcheck, masscheck;
  FILE *FreeSurfaceFID, *HorizontalVelocityFID, *VerticalVelocityFID,
    *SalinityFID, *VerticalGridFID, *ConserveFID;

} propT;

/*
 * Public function declarations.
 *
 */
void Solve(gridT *grid, physT *phys, int myproc, int numprocs, MPI_Comm comm);
void AllocatePhysicalVariables(gridT *grid, physT **phys);
void FreePhysicalVariables(gridT *grid, physT *phys);
void InitializePhysicalVariables(gridT *grid, physT *phys);
void InitializeVerticalGrid(gridT **grid);

#endif
