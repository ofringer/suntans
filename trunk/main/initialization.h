/*
 * File: initialization.h
 * Description: Header file for initialization.c
 *
 * $Id: initialization.h,v 1.1 2003-04-29 00:23:33 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#ifndef _initialization_h
#define _initialization_h

REAL ReturnDepth(REAL x, REAL y);
REAL ReturnFreeSurface(REAL x, REAL y);
REAL ReturnSalinity(REAL x, REAL y, REAL z);
REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z);

#endif
