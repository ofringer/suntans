/*
 * Header file for util.c
 *
 * $Id: util.h,v 1.2 2002-11-05 01:31:17 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:24:01  fringer
 * Initial revision
 *
 *
 */
#ifndef _util_h
#define _util_h

#include "suntans.h"

void Sort(int *a, int *v, int N);
void ReOrderIntArray(int *a, int *order, int *tmp, int N, int Num);
void ReOrderRealArray(REAL *a, int *order, REAL *tmp, int N, int Num);
int *ReSize(int *a, int N);
int IsMember(int i, int *points, int numpoints);
void FindNearest(int *points, REAL *x, REAL *y, int N, int np, REAL xi, REAL yi);
void Interp(REAL *x, REAL *y, REAL *z, int N, REAL *xi, REAL *yi, REAL *zi, int Ni);
void TriSolve(REAL *a, REAL *b, REAL *c, REAL *d, REAL *u, int N);
int IsNan(REAL x);
REAL UpWind(REAL u, REAL dz1, REAL dz2);

#endif
