/*
 * File: memory.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains memory allocation and deallocation functions
 * that keep track of the total memory used and freed throughout the
 * program with the global variable TotSpace, which contains the
 * total space used in bytes.  If the global varialbe VerboseMemory
 * is set to 1, then memory statistics will be printed as memory
 * is allocated and freed.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include<stdlib.h>
#include<stdio.h>
#include "memory.h"

/*
 * Function: SunMalloc
 * Usage: ptr=(int *)SunMalloc(N*sizeof(int),"Function");
 * ------------------------------------------------------
 * Same as the malloc function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void *SunMalloc(const int bytes, const char *function) {
  void *ptr = malloc(bytes);

  //VerboseMemory=1;

  if(TotSpace==-1)
    TotSpace=0;

  if(ptr==NULL) {
    printf("Error.  Out of memory!\n");
    printf("Total memory: %d, attempted to allocate: %d in function %s\n",
	   TotSpace,bytes,function);
    exit(1);
  } else {
    TotSpace+=bytes;
    if(VerboseMemory)
      printf("Allocated %d, Total: %d (%s)\n",
	     bytes,TotSpace,function);
    return ptr;
  }
}

/*
 * Function: SunFree
 * Usage: SunFree(ptr,bytes,"Function");
 * ------------------------------------------
 * Same as the free function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void SunFree(void *ptr, const int bytes, const char *function) {
  if(ptr==NULL) {
    printf("Error!  Attempting to free a NULL pointer in funciton %s\n",function);
    exit(1);
  } else {
    free(ptr);
    TotSpace-=bytes;
    if(VerboseMemory)
      printf("Freed %d, Total: %d (%s)\n",
	     bytes,TotSpace,function);
  }
}
