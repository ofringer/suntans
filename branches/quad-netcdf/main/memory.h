/*
 * File: memory.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for memory.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _memory_h
#define _memory_h

#include "suntans.h"

unsigned TotSpace;
int VerboseMemory;
char oldAllocFunction[BUFFERLENGTH], oldFreeFunction[BUFFERLENGTH];

/*
 * Function: SunMalloc
 * Usage: ptr=(int *)SunMalloc(N*sizeof(int),"Function");
 * ------------------------------------------------------
 * Same as the malloc function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void *SunMalloc(const unsigned bytes, const char *function);

/*
 * Function: SunFree
 * Usage: SunFree(ptr,bytes,"Function");
 * ------------------------------------------
 * Same as the free function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void SunFree(void *ptr, const unsigned bytes, const char *function);

#endif
