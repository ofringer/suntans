/*
 * File: memory.h
 * Description: Header file for memory.c
 *
 * $Id: memory.h,v 1.2 2004-05-29 20:25:02 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2003/04/29 00:10:59  fringer
 * Initial revision
 *
 *
 */
#ifndef _memory_h
#define _memory_h

int TotSpace, VerboseMemory;

/*
 * Function: SunMalloc
 * Usage: ptr=(int *)SunMalloc(N*sizeof(int),"Function");
 * ------------------------------------------------------
 * Same as the malloc function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void *SunMalloc(const int bytes, const char *function);

/*
 * Function: SunFree
 * Usage: SunFree(ptr,bytes,"Function");
 * ------------------------------------------
 * Same as the free function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void SunFree(void *ptr, const int bytes, const char *function);

#endif
