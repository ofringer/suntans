/*
 * File: kriging.c
 * ---------------
 *
 * $Id: kriging.c,v 1.3 2004-05-29 20:25:02 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2003/06/10 02:23:27  fringer
 * Added functions to precompute inverse of kriging matrices that are
 * stored during the first step of computation in phys.c
 *
 * Revision 1.1  2003/05/12 00:14:18  fringer
 * Initial revision
 *
 *
 */
#include<math.h>
#include "f2c.h"
#include "clapack.h"
#include "suntans.h"

void kriging(REAL xd, REAL yd, REAL *fd, REAL *x, REAL *y, REAL **f, int N0, 
	     int p, int N, REAL *Ai);
REAL *inversekriging(REAL *x, REAL *y, int p, int N);
void product(REAL *A, REAL *x, int N);
void inverse(float *A, int N);

static REAL K(float x,int p) {
  return pow(fabs(x),p);
}

void kriging(REAL xd, REAL yd, REAL *fd, REAL *x, REAL *y, REAL **f, int N0, 
	     int p, int N, REAL *Ai) {
  int i, j, *ipiv, N1=N+1, M=1;
  REAL *lam = (REAL *)SunMalloc((N+1)*sizeof(REAL),"kriging");

  for(i=0;i<N;i++)
    lam[i]=K(sqrt(pow(x[i]-xd,2)+pow(y[i]-yd,2)),p);
  lam[N]=1;

  product(Ai,lam,N+1);

  for(j=0;j<N0;j++) {
    fd[j]=0;
    for(i=0;i<N;i++)
      fd[j]+=lam[i]*f[j][i];
  }

  SunFree(lam,(N+1)*sizeof(float),"kriging");
}

REAL *inversekriging(REAL *x, REAL *y, int p, int N) {
  int i, j;
  float *A = (float *)SunMalloc((N+1)*(N+2)/2*sizeof(float),"inversekriging");
  REAL *Ai = (REAL *)SunMalloc((N+1)*(N+2)/2*sizeof(REAL),"inversekriging");

  for(i=0;i<(N+1)*(N+2)/2;i++)
    A[i]=0;

  for(i=0;i<N;i++) {
    for(j=i+1;j<N;j++)
      A[i + j*(j+1)/2]=K(sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)),p);
    A[i + j*(j+1)/2]=1;
  }
  A[(N+1)*(N+2)/2-1]=0;

  inverse(A,N+1);

  for(i=0;i<(N+1)*(N+2)/2;i++)
    Ai[i]=A[i];

  SunFree(A,(N+1)*(N+2)/2*sizeof(float),"inversekriging");
  return Ai;
}

void inverse(float *A, int N) {
  int i, j, info=0;
  int *ipiv = (int *)SunMalloc(N*sizeof(int),"inverse");
  float *X = (float *)SunMalloc(N*N*sizeof(float),"inverse");

  // X is the identity, since to find the inverse
  // of A that is the solution of AX=I.
  for(i=0;i<N*N;i++)
    X[i]=0;
  for(i=0;i<N;i++)
    X[i*(N+1)]=1;

  // Use sspsv to solve and place inv(A) into X.
  sspsv_("U",(integer *)(&N),(integer *)(&N),(real *)A,(integer *)ipiv,
	 (real *)X,(integer *)(&N),(integer *)(&info));

  // Now replace A with inv(A)
  for(j=0;j<N;j++)
    for(i=0;i<=j;i++) 
      A[i+j*(j+1)/2]=X[i+j*N];

  SunFree(X,N*N*sizeof(float),"inverse");
  SunFree(ipiv,N*sizeof(int),"inverse");
}
  
void product(REAL *A, REAL *x, int N) {
  int i, j, ind;
  REAL *y = (REAL *)SunMalloc(N*sizeof(REAL),"product");

  for(i=0;i<N;i++) {
    y[i]=0;
    for(j=i;j<N;j++) 
      y[i]+=A[i + j*(j+1)/2]*x[j];
    for(j=0;j<i;j++) 
      y[i]+=A[j + i*(i+1)/2]*x[j];
  }

  for(i=0;i<N;i++) 
    x[i]=y[i];

  SunFree(y,N*sizeof(REAL),"product");
}

