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
