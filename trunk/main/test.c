/*
 * file: test.c
 * ------------
 *
 * $Id: test.c,v 1.1 2003-07-31 09:01:09 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */

#include "mpi.h"

int main(void) {
  
  int *a = (int *)MyMalloc(20);
  MyFree(a,20);
}
