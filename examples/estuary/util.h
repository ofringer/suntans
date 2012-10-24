/*
 * Header file for util.c
 *
 * $Id: util.h,v 1.1 2005-10-31 05:59:10 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2004/09/22 06:31:59  fringer
 * Added the Min function.
 *
 * Revision 1.4  2004/06/15 18:36:33  fringer
 * Changed definition of Max to static.
 *
 * Revision 1.3  2004/05/29 20:25:02  fringer
 *  Revision before converting to CVS.
 *
 * Revision 1.2  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
 * Revision 1.1  2002/11/03 00:24:01  fringer
 * Initial revision
 *
 *
 */
#ifndef _util_h
#define _util_h

#include "grid.h"
#include "suntans.h"

void Sort(int *a, int *v, int N);
void ReOrderIntArray(int *a, int *order, int *tmp, int N, int Num);
void ReOrderRealArray(REAL *a, int *order, REAL *tmp, int N, int Num);
int *ReSize(int *a, int N);
REAL *ReSizeReal(REAL *a, int N);
int IsMember(int i, int *points, int numpoints);
void FindNearest(int *points, REAL *x, REAL *y, int N, int np, REAL xi, REAL yi);
void Interp(REAL *x, REAL *y, REAL *z, int N, REAL *xi, REAL *yi, REAL *zi, int Ni);
void TriSolve(REAL *a, REAL *b, REAL *c, REAL *d, REAL *u, int N);
int IsNan(REAL x);
REAL UpWind(REAL u, REAL dz1, REAL dz2);
void Copy(REAL **from, REAL **to, gridT *grid);
static REAL Max(REAL x1, REAL x2);
REAL Min(REAL x, REAL y);

#endif
