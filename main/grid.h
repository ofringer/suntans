/*
 * File: grid.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for grid.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _grid_h
#define _grid_h

#include "suntans.h"
#include "fileio.h"
#include "mympi.h"

#define MAXBCTYPES 4
#define MAXMARKS 7
#define NUMDGHIST 10

/*
 * Main grid struct.
 *
 */
typedef struct _gridT {
  REAL *xp;
  REAL *yp;
  REAL *xv;
  REAL *yv;
  REAL *xe;
  REAL *ye;
  REAL *dv;
  REAL *dz;
  REAL **dzz;
  REAL **dzf;
  REAL *dzfB;
  REAL **dzzold;
  REAL *dztop;
  int *face;
  int *edges;
  int *cells;
  int *neigh;
  int *eneigh;
  int *grad;
  int *gradf;
  int *mark;
  int *normal;

  int *xadj;
  int *adjncy;
  int *vtxdist;
  int *cellbnddist;
  int *edgebnddist;

  int *celldist;
  int *edgedist;
  int *cellp;
  int *edgep;
  int **cell_send;
  int **cell_recv;
  int **edge_send;
  int **edge_recv;
  int *num_cells_send;
  int *num_cells_recv;
  int *num_cells_send_first;
  int *num_cells_recv_first;
  int *num_edges_send;
  int *num_edges_recv;

  int maxtosend, maxtorecv;
  MPI_Status *status;
  MPI_Request *request;
  REAL **recv;
  REAL **send;
  int *total_cells_send;
  int *total_cells_recv;
  int *total_cells_sendW;
  int *total_cells_recvW;
  int *total_edges_send;
  int *total_edges_recv;

  int *mnptr;
  int *eptr;
  int *part;
  int *order;
  int *epart;
  int *vwgt;

  REAL *df;
  REAL *dg;
  REAL *def;
  REAL *Ac;
  REAL *n1;
  REAL *n2;
  REAL *xi;

  int *numneighs;
  int **neighs;

  int Nneighs;
  int *myneighs;

  int *Nk;
  int *ctop;
  int *ctopold;
  int *etop;
  int *etopold;
  int *Nke;
  int *Nkc;
  int Nkmax;
  int Nc;
  int Nc_comp;
  int Ne;
  int Ne_comp;
  int Np;
  int Nge;

} gridT;

/*
 * Public function declarations.
 *
 */
void GetGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void Partition(gridT *maingrid, gridT **localgrid, MPI_Comm comm);
void SendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, 
			MPI_Comm comm);
void ISendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, 
			 MPI_Comm comm);
void SendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, 
			MPI_Comm comm);
void ISendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, 
			MPI_Comm comm);
void SendRecvWData(REAL **celldata, gridT *grid, int myproc, 
		   MPI_Comm comm);
void ISendRecvWData(REAL **celldata, gridT *grid, int myproc, 
		    MPI_Comm comm);
void SendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, 
			MPI_Comm comm);
void ISendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, 
			 MPI_Comm comm);
void CheckCommunicateCells(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void CheckCommunicateEdges(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void InitMainGrid(gridT **grid, int Np, int Ne, int Nc);
void ReadMainGrid(gridT *grid, int myproc);
void GetDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
void CreateCellGraph(gridT *grid);
void CreateEdgeGraph(gridT *grid);
void Connectivity(gridT *grid, int myproc);
void ReadFileNames(int myproc);
int IsBoundaryCell(int mgptr, gridT *maingrid, int myproc);
void ReadGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void AllocateTransferArrays(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void FreeTransferArrays(gridT *grid, int myproc, int numprocs, MPI_Comm comm);

#endif
