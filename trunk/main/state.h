/*
 * File: state.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for state.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _state_h
#define _state_h

#include "grid.h"
#include "phys.h"

/*
 * Function: StateEquation
 * Usage: rho = StateEquation(prop,s,T,p);
 * ---------------------------------------
 * Returns the density as a function of temperature, salinity, and
 * pressure.
 *
 */
REAL StateEquation(const propT *prop, const REAL s, const REAL T, const REAL p);

#endif
