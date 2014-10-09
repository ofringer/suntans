/*
 * File: util.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for util.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _util_h
#define _util_h

#include "grid.h"
#include "suntans.h"

enum Type 
{
  DOUBLE,
  INT
};

void Sort(int *a, int *v, int N);
void ReOrderIntArray(int *a, int *order, int *tmp, int N, int Num, int *nfaces, int *grad, int maxfaces);
void ReOrderRealArray(REAL *a, int *order, REAL *tmp, int N, int Num, int *nfaces, int *grad, int maxfaces);
int *ReSize(int *a, int N);
int IsMember(int i, int *points, int numpoints);
int FindNearest(int *points, REAL *x, REAL *y, int N, int np, REAL xi, REAL yi);
void Interp(REAL *x, REAL *y, REAL *z, int N, REAL *xi, REAL *yi, REAL *zi, int Ni, int maxFaces);
void TriSolve(REAL *a, REAL *b, REAL *c, REAL *d, REAL *u, int N);
int IsNan(REAL x);
REAL UpWind(REAL u, REAL dz1, REAL dz2);
void Copy(REAL **from, REAL **to, gridT *grid);
REAL Max(REAL x1, REAL x2);
REAL Min(REAL x, REAL y);
void ComputeGradient(REAL **gradient, REAL **phi, gridT *grid, int direction);
void PrintVectorToFile(enum Type etype, void *Vector, int M, char *filename, int myproc);
int max(int a, int b);
int SharedListValue(int *list1, int *list2, int listsize);
REAL QuadInterp(REAL x, REAL x0, REAL x1, REAL x2, REAL y0, REAL y1, REAL y2);
REAL getToffSet(char starttime[15], char basetime[15]);
#endif
