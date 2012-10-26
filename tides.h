/*
 * File: tides.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for tides.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _tides_h
#define _tides_h

void SetTideComponents(gridT *grid, int myproc);

int numtides;
REAL **u_amp, **v_amp, **h_amp, **u_phase, **v_phase, **h_phase, *omegas;

#endif

