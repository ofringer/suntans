/*
 * File: sendrecv.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for sendrecv.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _sendrecv_h
#define _sendrecv_h

#include "grid.h"
#include "mympi.h"

void AllocateTransferArrays(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void FreeTransferArrays(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
void ISendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm);
void ISendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm);
void ISendRecvWData(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm);
void ISendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, MPI_Comm comm);
void CheckCommunicateCells(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void CheckCommunicateEdges(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);

#endif

