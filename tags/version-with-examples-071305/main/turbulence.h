/*
 * File: turbulence.h
 * Description: Header file for turbulence.c
 *
 * $Id: turbulence.h,v 1.2 2004-09-16 21:33:01 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2004/09/13 04:14:55  fringer
 * Header file for turbulence.c.
 *
 *
 */
#ifndef _turbulence_h
#define _turbulence_h

// Von Karman's Constant
#define KAPPA_VK 0.42

void my25(gridT *grid, physT *phys, propT *prop, REAL **q, REAL **l, REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc);

#endif
