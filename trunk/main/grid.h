/*
 * Header file for grid.c
 *
 * $Id: grid.h,v 1.4 2003-04-29 00:21:25 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2003/04/21 20:26:35  fringer
 * Working version before addition of ghost cells for kriging.
 *
 * Revision 1.2  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
 * Revision 1.1  2002/11/03 00:18:19  fringer
 * Initial revision
 *
 *
 */
#ifndef _grid_h
#define _grid_h

#include "parmetis.h"
#include "suntans.h"
#include "fileio.h"
#include "mympi.h"

#define MAXBCTYPES 4
#define MAXMARKS 7

/*
 * Main grid struct.
 *
 */
typedef struct _gridT {
  REAL *xp;
  REAL *yp;
  REAL *xv;
  REAL *yv;
  REAL *dv;
  REAL *dz;
  REAL **dzz;
  REAL **dzzold;
  REAL *dztop;
  int *face;
  int *edges;
  int *cells;
  int *neigh;
  int *eneigh;
  int *grad;
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
  int *num_edges_send;
  int *num_edges_recv;

  int *mnptr;
  int *eptr;
  int *part;
  int *order;
  int *epart;
  int *vwgt;

  REAL *df;
  REAL *dg;
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

  int Nnearestcells;
  int Nnearestedges;
  int **nearestcells;
  int **nearestedges;

} gridT;

/*
 * Public function declarations.
 *
 */
void GetGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void Partition(gridT *maingrid, gridT **localgrid, MPI_Comm comm);
void SendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm);
void SendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm);
void CheckCommunicateCells(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void CheckCommunicateEdges(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void InitMainGrid(gridT **grid, int Np, int Ne, int Nc);
//void InitMainGrid(gridT **grid);
void ReadMainGrid(gridT *grid);
void GetDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
void CreateCellGraph(gridT *grid);
void CreateEdgeGraph(gridT *grid);
void Connectivity(gridT *grid, int myproc);
void ReadFileNames(int myproc);
int IsBoundaryCell(int mgptr, gridT *maingrid, int myproc);

#endif
