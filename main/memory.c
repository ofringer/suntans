/*
 * File: memory.c
 * Description: Contains memory allocation and deallocation functions
 * that keep track of the total memory used and freed throughout the
 * program with the global variable TotSpace, which contains the
 * total space used in bytes.  If the global varialbe VerboseMemory
 * is set to 1, then memory statistics will be printed as memory
 * is allocated and freed.
 *
 * $Id: memory.c,v 1.1 2003-04-29 00:10:32 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#include<stdlib.h>
#include "memory.h"

/*
 * Function: MyMalloc
 * Usage: ptr=(int *)MyMalloc(N*sizeof(int));
 * ------------------------------------------
 * Same as the malloc function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void *MyMalloc(const int bytes) {
  void *ptr = malloc(bytes);

  if(TotSpace==-1)
    TotSpace=0;

  if(ptr==NULL) {
    printf("Error.  Out of memory!\n");
    printf("Total memory: %d, attempted to allocate: %d\n",
	   TotSpace,bytes);
    exit(1);
  } else {
    TotSpace+=bytes;
    if(VerboseMemory)
      printf("Allocated %d, Total: %d\n",
	     bytes,TotSpace);
    return ptr;
  }
}

/*
 * Function: MyFree
 * Usage: MyFree(ptr);
 * ------------------------------------------
 * Same as the free function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void MyFree(void *ptr, int bytes) {
  if(ptr==NULL) {
    printf("Error!  Attempting to free a NULL pointer!\n");
    exit(1);
  } else {
    free(ptr);
    TotSpace-=bytes;
    if(VerboseMemory)
      printf("Freed %d, Total: %d\n",
	     bytes,TotSpace);
  }
}
