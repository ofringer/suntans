/*
 * File: diffusion.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for diffusion.c
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _diffusion_h
#define _diffusion_h

#include "grid.h"
#include "phys.h"

void LaxWendroff(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);

#endif
