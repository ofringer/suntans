/*
 * File: kriging.c
 * ---------------
 *
 * $Id: kriging.c,v 1.1 2003-05-12 00:14:18 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#include<math.h>
#include "f2c.h"
#include "clapack.h"
#include "suntans.h"

void kriging(REAL xd, REAL yd, REAL *fd, REAL *x, REAL *y, REAL *uf, int p, int N, int *status);

static REAL K(float x,int p) {
  return pow(fabs(x),p);
}

void kriging(REAL xd, REAL yd, REAL *fd, REAL *x, REAL *y, REAL *f, int p, int N, int *status) {
  int i, j, *ipiv, N1=N+1, M=1;
  float *lam, *A;

  ipiv = (int *)malloc((N+1)*sizeof(int));
  A = (float *)malloc((N+1)*(N+2)/2*sizeof(float));
  lam = (float *)malloc((N+1)*sizeof(float));

  for(i=0;i<(N+1)*(N+2)/2;i++)
    A[i]=0;

  for(i=0;i<N;i++) {
    for(j=i+1;j<N;j++)
      A[i + j*(j+1)/2]=K(sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)),p);
    A[i + j*(j+1)/2]=1;
  }
  A[(N+1)*(N+2)/2-1]=0;

  for(i=0;i<N;i++)
    lam[i]=K(sqrt(pow(x[i]-xd,2)+pow(y[i]-yd,2)),p);
  lam[N]=1;

  *status=0;
  sspsv_("U",(integer *)(&N1),(integer *)(&M),(real *)A,(integer *)ipiv,
	 (real *)lam,(integer *)(&N1),(integer *)(status));

  *fd=0;
  for(i=0;i<N;i++)
    (*fd)+=lam[i]*f[i];

  free(ipiv);
  free(A);
  free(lam);

}
