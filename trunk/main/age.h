/*
 * Age tracer related header file 
 */
#ifndef _age_h
#define _age_h

#include "phys.h"
#include "grid.h"

// Age structure
typedef struct _ageT {
  REAL **agec;
  REAL **agealpha;
  REAL **Cn_Ac;
  REAL **Cn_Aa;

  REAL **boundary_age;
  REAL **boundary_agealpha;

} ageT;

//Make the age struct global
ageT *age;

/*
 * Public function declarations.
 *
 */

void UpdateAge(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);
#endif
