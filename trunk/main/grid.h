/*
 * Header file for grid.c
 *
 * $Id: grid.h,v 1.7 2003-12-02 02:05:46 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.6  2003/05/12 00:20:25  fringer
 * Added ReadGrid and ISendRecvCellData2D, as well as xe and ye to
 * the gridT struct.
 *
 * Revision 1.5  2003/05/01 00:32:39  fringer
 * Added prototypes for first/all transfer of 2D/3D data.
 *
 * Revision 1.4  2003/04/29 00:21:25  fringer
 * The include line for parmetis.h has to come after the mympi.h
 * line since parmetis.h does not have the mpi.h header in it.
 * Changed prototype for Connectivity to use myproc to restrict printing on
 * proc 0.
 *
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
  REAL *xe;
  REAL *ye;
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

typedef enum {
  first, all
} transferType;

/*
 * Public function declarations.
 *
 */
void GetGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void Partition(gridT *maingrid, gridT **localgrid, MPI_Comm comm);
void SendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, 
			MPI_Comm comm, transferType type);
void ISendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, 
			MPI_Comm comm, transferType type);
void SendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, 
			MPI_Comm comm, transferType type);
void SendRecvWData(REAL **celldata, gridT *grid, int myproc, 
		   MPI_Comm comm, transferType type);
void SendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, 
			MPI_Comm comm);
void CheckCommunicateCells(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void CheckCommunicateEdges(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void InitMainGrid(gridT **grid, int Np, int Ne, int Nc);
//void InitMainGrid(gridT **grid);
void ReadMainGrid(gridT *grid, int myproc);
void GetDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
void CreateCellGraph(gridT *grid);
void CreateEdgeGraph(gridT *grid);
void Connectivity(gridT *grid, int myproc);
void ReadFileNames(int myproc);
int IsBoundaryCell(int mgptr, gridT *maingrid, int myproc);
void ReadGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm);

#endif
