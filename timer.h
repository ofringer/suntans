/*
 * File: timer.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for timer.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _timer_h
#define _timer_h

#include "suntans.h"

// Global variables for timing
REAL t_start, t_source, t_predictor, t_nonhydro, t_turb, t_transport, t_io, t_comm,
  t_check, t_tictoc;

/*
 * Function: Timer
 * Usage: printf("Time = %f\n",Timer()-t0);
 * ----------------------------------------
 * Returns the time in seconds and uses the timer
 * defined in timer.h.
 *
 */
extern REAL Timer(void);
inline void Tic(void);
inline REAL Toc(void);

#endif
