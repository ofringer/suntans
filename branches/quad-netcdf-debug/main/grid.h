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

//Parmetis 2.0
//#include "parmetis.h"
//Parmetis 3.1
//#include "parmetislib.h"

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
  REAL *dzbot; 
  REAL dzsmall;
  int *face;
  int *edges;
  int *cells;
  int *neigh;
  int *eneigh;
  int *grad;
  int *gradf;
  int *mark;
  int *edge_id;
  int *normal;

  // new deinitions for nodal neighbor data structure
  int *localtoglobalpoints; //global index from local point index
  int *numppneighs;   // number of point neighbors to point
  int **ppneighs;     // index of points neighbors to a point
  int *numpeneighs;   // number of edge neighbors to point
  int **peneighs;     // index of edge neighbors to a point
  int *numpcneighs;   // number of cell neighbors to point
  int **pcneighs;     // index of cell neighbors to a point
  REAL **Actotal;     // total area of cell neighbor areas to point
                      // note that this is over all the cells and 
                      // the layer too
// eventually want to get these here and not use global variables
//  // variables to define how the halo works
//  boundaryselection g_halolist[2] ={EDGE, NODE}; 
//  int g_halolistsize = 2;

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
  int *Nkp;
  int *Nke;
  int *Nkc;
  int Nkmax;
  int Nc;
  int Nc_comp;
  int Ne;
  int Ne_comp;
  int Np;
  int Nge;
  int *nfaces;
  int maxfaces;
  int stairstep;
  int fixdzz;
  REAL smoothbot;
} gridT;

// enums used to set type of boundary selection we will use
typedef enum _boundaryselection{
  EDGE,
  NODE
} boundaryselection;


/*
 * Public function declarations.
 *
 */
void GetGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void Partition(gridT *maingrid, gridT **localgrid, MPI_Comm comm);
void InitMainGrid(gridT **grid, int Np, int Ne, int Nc, int myproc);
void ReadMainGrid(gridT *grid, int myproc);
void GetDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
void CreateCellGraph(gridT *grid);
void CreateEdgeGraph(gridT *grid);
void Connectivity(gridT *grid, int myproc);
inline int IsBoundaryCell(int mgptr, gridT *maingrid, int myproc);
void AllocateTransferArrays(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void FreeTransferArrays(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
REAL GetArea(REAL *xt, REAL *yt, int Nf);
void InitLocalGrid(gridT **grid);

#endif
