/*
 * File: memory.c
 * Description: Contains memory allocation and deallocation functions
 * that keep track of the total memory used and freed throughout the
 * program with the global variable TotSpace, which contains the
 * total space used in bytes.  If the global varialbe VerboseMemory
 * is set to 1, then memory statistics will be printed as memory
 * is allocated and freed.
 *
 * $Id: memory.c,v 1.3 2004-05-29 20:25:02 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2003/06/10 02:22:36  fringer
 * Updated to print out memory stats if VerboseMemory==1.
 *
 * Revision 1.1  2003/04/29 00:10:32  fringer
 * Initial revision
 *
 *
 */
#include<stdlib.h>
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
