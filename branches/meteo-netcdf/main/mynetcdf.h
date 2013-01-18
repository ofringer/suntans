/*
 * NetCDF IO header file
 *
 *
 */

#include <netcdf.h>
#include "suntans.h"

/* Netcdf error */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
  
/* Public functions */
void nc_read_3D(int ncid, char *vname, size_t *start, size_t *count, REAL ***tmparray);
void nc_read_2D(int ncid, char *vname, size_t *start, size_t *count, REAL **tmparray, int myproc);
