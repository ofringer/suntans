/*
 * Header file for phys.c
 *
 * $Id: phys.h,v 1.6 2003-05-12 00:19:03 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2003/04/29 00:19:56  fringer
 * Added BGSalinityFID
 *
 * Revision 1.4  2002/12/01 10:42:46  fringer
 * Added initial density vector s0 as well as Flux terms to compute internal
 * wave energy fluxes at specified faces.
 *
 * Revision 1.3  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
 * Revision 1.2  2002/11/03 20:36:33  fringer
 * Added parameters to check for mass and volume conservation
 *
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
  REAL **s0;
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
  REAL Eflux1;
  REAL Eflux2;
  REAL Eflux3;
  REAL Eflux4;
  REAL smin;
  REAL smax;

  REAL *htmp;
  REAL *hold;
  REAL **stmp;
  REAL **utmp;
  REAL **ut;
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
  REAL dt, Cmax, rtime, amp, omega, flux, theta, thetaAB, 
    thetaFS, nu, tau_T, CdT, CdB, relax, epsilon, dzsmall, beta;
  int ntout, ntprog, nsteps, n, ntconserve, maxiters, volcheck, masscheck,
    kriging_cov, nonlinear;
  FILE *FreeSurfaceFID, *HorizontalVelocityFID, *VerticalVelocityFID,
    *SalinityFID, *BGSalinityFID, *VerticalGridFID, *ConserveFID;

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
