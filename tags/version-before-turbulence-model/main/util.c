/*
 * File: util.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 10/16/02
 * --------------------------------
 * This file contains utility functions for array operations.
 *
 * $Id: util.c,v 1.3 2004-05-29 20:25:02 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2003/04/26 14:23:39  fringer
 * Changed the large number to INFTY, which is defined in suntans.h
 *
 * Revision 1.1  2002/11/03 00:23:27  fringer
 * Initial revision
 *
 *
 */
#include "grid.h"
#include "util.h"

void Sort(int *a, int *v, int N)
{
  int i, j, temp, *tmp;

  tmp = (int *)malloc(N*sizeof(int));
  for(i=0;i<N;i++)
    tmp[i] = v[i];

  for(i=0;i<N;i++) {
    for(j=i+1;j<N;j++) 
      if(tmp[j]<tmp[i]) {
	temp = a[j];
	a[j] = a[i];
	a[i] = temp;
	temp = tmp[j];
	tmp[j] = tmp[i];
	tmp[i] = temp;
      }
  }

  free(tmp);
}
      
int *ReSize(int *a, int N)
{
  int i;
  int *new = (int *)malloc(N*sizeof(int));
  for(i=0;i<N;i++)
    new[i]=a[i];
  free(a);

  return new;
}

void ReOrderRealArray(REAL *a, int *order, REAL *tmp, int N, int Num)
{
  int n, j;

  for(n=0;n<N;n++)
    for(j=0;j<Num;j++) 
      tmp[n*Num+j]=a[n*Num+j];

  for(n=0;n<N;n++)
    for(j=0;j<Num;j++)
      a[n*Num+j]=tmp[order[n]*Num+j];
}

void ReOrderIntArray(int *a, int *order, int *tmp, int N, int Num)
{
  int n, j;

  for(n=0;n<N;n++)
    for(j=0;j<Num;j++)
      tmp[n*Num+j]=a[n*Num+j];
  for(n=0;n<N;n++)
    for(j=0;j<Num;j++)
      a[n*Num+j]=tmp[order[n]*Num+j];
}

int IsMember(int i, int *points, int numpoints)
{
  int j;
  for(j=0;j<numpoints;j++)
    if(i==points[j]) return j;
  return -1;
}    

void Interp(REAL *x, REAL *y, REAL *z, int N, REAL *xi, REAL *yi, REAL *zi, int Ni)
{
  int j, n, numpoints=3, *points=(int *)malloc(numpoints*sizeof(int));
  REAL r;

  for(n=0;n<Ni;n++) {
    FindNearest(points,x,y,N,numpoints,xi[n],yi[n]);
    zi[n]=0;
    r=0;
    for(j=0;j<numpoints;j++) 
      r += sqrt(x[points[j]]*x[points[j]]+0*y[points[j]]*y[points[j]]);
    for(j=0;j<numpoints;j++)
      zi[n] += z[points[j]]*sqrt(x[points[j]]*x[points[j]]+0*y[points[j]]*y[points[j]])/r;
  }

  free(points);
}

void FindNearest(int *points, REAL *x, REAL *y, int N, int np, REAL xi, REAL yi)
{
  int i, n;
  REAL dist, d;
  for(n=0;n<np;n++) points[n]=-1;

  for(n=0;n<np;n++) {
    dist=INFTY;
    for(i=0;i<N;i++) {
      d = (x[i]-xi)*(x[i]-xi)+(y[i]-yi)*(y[i]-yi);
      if(d<dist & IsMember(i,points,np)==-1) {
	dist=d;
	points[n]=i;
      }
    }
  }
}


void TriSolve(REAL *a, REAL *b, REAL *c, REAL *d, REAL *u, int N)
{
  int k;

  //  printf("In trisolve with N = %d\n",N);

  for(k=1;k<N;k++) {
    //    if(b[k-1]==0) printf("b[%d]=0\n",k-1);
    b[k]-=a[k]*c[k-1]/b[k-1];
    d[k]-=a[k]*d[k-1]/b[k-1];
  }

  u[N-1]=d[N-1]/b[N-1];

  for(k=N-2;k>=0;k--)
    u[k] = d[k]/b[k]-c[k]*u[k+1]/b[k];
}
  
int IsNan(REAL x) 
{
  if(x!=x)
    return 1;
  return 0;
}

static REAL Max(REAL x1, REAL x2)
{
  if(x2>x1)
    return x2;
  return x1;
}

REAL UpWind(REAL u, REAL dz1, REAL dz2)
{
  REAL fluxheight;

  if(fabs(u)>0)
    if(u>0)
      fluxheight=dz2;
    else
      fluxheight=dz1;
  else
    fluxheight=Max(dz1,dz2);
  
  return fluxheight;
}

void Copy(REAL **from, REAL **to, gridT *grid) {
  int i, k;

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      to[i][k]=from[i][k];
}

