/*
 * File: turbulence.c
 * Description:  Contains the Mellor-Yamad level 2.5 turbulence model.
 *
 * $Id: turbulence.c,v 1.1 2004-09-13 04:14:36 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#include "math.h"
#include "phys.h"
#include "grid.h"
#include "turbulence.h"

void my25(gridT *grid, physT *phys, propT *prop, REAL **q, REAL **l, REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT) {
  int i, k;
  REAL fab;

  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++) {
      nuT[i][k]=0;
      kappaT[i][k]=0;
    }
  }

  if(prop->n==1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++) {
	Cn_q[i][k]=0;
	Cn_l[i][k]=0;
      }
  } else
    fab=1.5;
}
