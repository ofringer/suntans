/*
 * file: test.c
 * ------------
 *
 * $Id: test.c,v 1.3 2004-06-23 05:29:15 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2004/05/29 20:25:02  fringer
 *  Revision before converting to CVS.
 *
 * Revision 1.1  2003/07/31 09:01:09  fringer
 * Initial revision
 *
 *
 */

#include "mpi.h"

int main(void) {
  
  int *a = (int *)MyMalloc(20);
  MyFree(a,20);

}
