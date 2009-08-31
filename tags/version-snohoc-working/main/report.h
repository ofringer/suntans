/*
 * File: report.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for report.c.
 * 
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _report_h
#define _report_h

#include "mympi.h"
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
