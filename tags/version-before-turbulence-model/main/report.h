/*
 * Header file for report.c
 *
 * $Id: report.h,v 1.2 2004-05-29 20:25:02 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:23:06  fringer
 * Initial revision
 *
 *
 */
#ifndef _report_h
#define _report_h

#include<mpi.h>
#include "grid.h"

void ReportPartition(gridT *maingrid, 
		     gridT *localgrid, 
		     int myproc, 
		     MPI_Comm comm);
void ReportConnectivity(gridT *grid, 
			gridT *maingrid, 
			int myproc);
void ReportBoundaryDistributions(gridT *maingrid, 
				 gridT *localgrid, 
				 int myproc);
void ParseFlags(int argc, char *argv[], int myproc);
void Usage(char *str);

#endif
