/*
 * File: initialization.h
 * Description: Header file for initialization.c
 *
 * $Id: initialization.h,v 1.2 2003-05-12 00:19:51 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2003/04/29 00:23:33  fringer
 * Initial revision
 *
 *
 */
#ifndef _initialization_h
#define _initialization_h

REAL ReturnDepth(REAL x, REAL y);
REAL ReturnFreeSurface(REAL x, REAL y, REAL d);
REAL ReturnSalinity(REAL x, REAL y, REAL z);
REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z);

#endif
