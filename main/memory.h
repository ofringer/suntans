/*
 * File: memory.h
 * Description: Header file for memory.c
 *
 * $Id: memory.h,v 1.1 2003-04-29 00:10:59 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#ifndef _memory_h
#define _memory_h

int TotSpace, VerboseMemory;

/*
 * Function: MyMalloc
 * Usage: ptr=(int *)MyMalloc(N*sizeof(int));
 * ------------------------------------------
 * Same as the malloc function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void *MyMalloc(const int bytes);

/*
 * Function: MyFree
 * Usage: MyFree(ptr);
 * ------------------------------------------
 * Same as the free function in stdlib.h, but this
 * one keeps track of the total memory with the global
 * variable TotSpace.
 *
 */
void MyFree(void *ptr, int bytes);

#endif
