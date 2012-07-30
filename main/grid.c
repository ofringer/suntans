/*
 * File: grid.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains grid-based functions.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "grid.h"
#include "partition.h"
#include "util.h"
#include "initialization.h"
#include "memory.h"
#include "triangulate.h"
#include "report.h"
#include "timer.h"

#define VTXDISTMAX 100

/*
 * Private function declarations.
 */
static void InitLocalGrid(gridT **grid);
static void VertGrid(gridT *maingrid, gridT **localgrid, MPI_Comm comm);
void Topology(gridT **maingrid, gridT **localgrid, int myproc, int numprocs);
static int GetNumCells(gridT *grid, int proc);
static void TransferData(gridT *maingrid, gridT **localgrid, int myproc);
static int GetNumEdges(gridT *grid);
static void Geometry(gridT *maingrid, gridT **grid, int myproc);
static REAL GetCircumcircleRadius(REAL *xt, REAL *yt, int Nf);
static REAL GetArea(REAL *xt, REAL *yt, int Nf);
static void EdgeMarkers(gridT *maingrid, gridT **localgrid, int myproc);
static void ReOrder(gridT *grid);
static int IsCellNeighborProc(int nc, gridT *maingrid, gridT *localgrid, 
			      int myproc, int neighproc);
static int IsEdgeNeighborProc(int ne, gridT *maingrid, gridT *localgrid, 
			      int myproc, int neighproc);
static void MakePointers(gridT *maingrid, gridT **localgrid, int myproc, MPI_Comm comm);
static void ResortBoundaries(gridT *localgrid, int myproc);
static void InterpDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static void FreeGrid(gridT *grid, int numprocs);
static void OutputData(gridT *maingrid, gridT *grid, int myproc, int numprocs);
static void CreateFaceArray(int *grad, int *gradf, int *neigh, int *face, int Nc, int Ne);
static void CreateNormalArray(int *grad, int *face, int *normal, int Nc);
static void ReadDepth(gridT *grid, int myproc);
static int CorrectVoronoi(gridT *grid, int myproc);
static int CorrectAngles(gridT *grid, int myproc);
static void VoronoiStats(gridT *grid);
static void FixDZZ(gridT *grid, REAL maxdepth, int Nkmax, int fixdzz, int myproc);
static int GetNk(REAL *dz, REAL localdepth, int Nkmax);

/************************************************************************/
/*                                                                      */
/*                       Public Functions                               */
/*                                                                      */
/************************************************************************/

/*
 * Function: GetGrid
 * Usage: GetGrid(grid,myproc,numprocs,comm);
 * ------------------------------------------
 * Initialize the unstructured grids on each processor.
 *
 */
void GetGrid(gridT **localgrid, int myproc, int numprocs, MPI_Comm comm)
{
  int Np, Ne, Nc, numcorr;
  gridT *maingrid;

  ReadFileNames(myproc);

  if(!TRIANGULATE) {
    Np = MPI_GetSize(POINTSFILE,"GetGrid",myproc);
    Ne = MPI_GetSize(EDGEFILE,"GetGrid",myproc);
    Nc = MPI_GetSize(CELLSFILE,"GetGrid",myproc);

    // Every processor will know about data read in from
    // triangle as well as the depth.
    if(myproc==0 && VERBOSE>0) printf("Initializing Main Grid...\n");
    InitMainGrid(&maingrid,Np,Ne,Nc);

    if(myproc==0 && VERBOSE>0) printf("Reading Grid...\n");
    ReadMainGrid(maingrid,myproc);
  } else {
    if(myproc==0 && VERBOSE>0) printf("Triangulating the point set...\n");
    if(!GetTriangulation(&maingrid,myproc)) {
      printf("Error computing triangulation!\n");
      MPI_Finalize();
    exit(EXIT_FAILURE);
    }
  }

  // Adjust Voronoi points for obtuse triangles (CorrectVoronoi==-1) 
  // or for small Voronoi edge lengths (CorrectVoronoi==+1).
  if(myproc==0 && VERBOSE>1) {
    printf("Voronoi statistics:\n");
    VoronoiStats(maingrid);
  }
  if((int)MPI_GetValue(DATAFILE,"CorrectVoronoi","ReadMainGrid",0)>0) {
    numcorr=CorrectVoronoi(maingrid,myproc);
    if(myproc==0 && VERBOSE>1) {
      if(numcorr) {
	printf("Voronoi statistics after correction:\n");
	VoronoiStats(maingrid);
      } else 
	printf("No Voronoi points were corrected.\n");
    }
  }
  if((int)MPI_GetValue(DATAFILE,"CorrectVoronoi","ReadMainGrid",0)<0) {
    numcorr=CorrectAngles(maingrid,myproc);
    if(myproc==0 && VERBOSE>1) {
      if(numcorr) {
	printf("Voronoi statistics after correction:\n");
	VoronoiStats(maingrid);
      } else 
	printf("No Voronoi points were corrected.\n");
    }
  }

  if(myproc==0 && VERBOSE>0) printf("Getting the depth for graph weights...\n");
  GetDepth(maingrid,myproc,numprocs,comm);

  if(myproc==0 && VERBOSE>0) printf("Computing Graph...\n");
  CreateCellGraph(maingrid);

  if(myproc==0 && VERBOSE>0) printf("Computing Connectivity...\n");
  Connectivity(maingrid,myproc);

  if(myproc==0 && VERBOSE>0) printf("Partitioning...\n");
  Partition(maingrid,localgrid,comm);

  //CheckCommunicateCells(maingrid,*localgrid,myproc,comm);
  //  CheckCommunicateEdges(maingrid,*localgrid,myproc,comm);
  //  SendRecvCellData2D((*localgrid)->dv,*localgrid,myproc,comm);

  OutputData(maingrid,*localgrid,myproc,numprocs);
  //  FreeGrid(maingrid,numprocs);
}

void Partition(gridT *maingrid, gridT **localgrid, MPI_Comm comm)
{
  int j;
  int myproc, numprocs;

  InitLocalGrid(localgrid);

  MPI_Comm_size(comm, &numprocs);
  MPI_Comm_rank(comm, &myproc);

  GetPartitioning(maingrid,localgrid,myproc,numprocs,comm);
  Topology(&maingrid,localgrid,myproc,numprocs);

  if(myproc==0 && VERBOSE>2) printf("Transferring data...\n");
  TransferData(maingrid,localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("Creating edge markers...\n");
  EdgeMarkers(maingrid,localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("Computing edge and voronoi distances and areas...\n");
  Geometry(maingrid,localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("Vert grid...\n");
  VertGrid(maingrid,localgrid,comm);

  if(myproc==0 && VERBOSE>2) printf("Reordering...\n");
  //  ReOrder(*localgrid);

  if(myproc==0 && VERBOSE>2) printf("Making pointers...\n");
  MakePointers(maingrid,localgrid,myproc,comm);

  // NEED THIS?
  //  ResortBoundaries(*localgrid,myproc);

  if(VERBOSE>3) ReportConnectivity(*localgrid,maingrid,myproc);
  if(VERBOSE>2) ReportPartition(maingrid,*localgrid,myproc,comm);
}

/*
 * Function: SendRecvCellData2D
 * Usage: SendRecvCellData2D(grid->h,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the cell data back and forth between
 * processors.  
 *
 */
void SendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, 
			MPI_Comm comm)
{
  int n, neigh, neighproc, *num_send, *num_recv;
  REAL **recv, **send;
  MPI_Status status;

  num_send=grid->num_cells_send;
  num_recv=grid->num_cells_recv;
    
  recv = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvCellData2D");
  send = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvCellData2D");

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    send[neigh] = (REAL *)SunMalloc(num_send[neigh]*sizeof(REAL),"SendRecvCellData2D");
    recv[neigh] = (REAL *)SunMalloc(num_recv[neigh]*sizeof(REAL),"SendRecvCellData2D");

    for(n=0;n<num_send[neigh];n++)
      send[neigh][n]=celldata[grid->cell_send[neigh][n]];

    MPI_Send((void *)(send[neigh]),num_send[neigh],
	     MPI_DOUBLE,neighproc,1,comm); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Recv((void *)(recv[neigh]),num_recv[neigh],
	     MPI_DOUBLE,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    for(n=0;n<num_recv[neigh];n++)
      celldata[grid->cell_recv[neigh][n]]=recv[neigh][n];
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    SunFree(send[neigh],num_send[neigh]*sizeof(REAL),"SendRecvCellData2D");
    SunFree(recv[neigh],num_recv[neigh]*sizeof(REAL),"SendRecvCellData2D");
  }
  SunFree(send,grid->Nneighs*sizeof(REAL *),"SendRecvCellData2D");
  SunFree(recv,grid->Nneighs*sizeof(REAL *),"SendRecvCellData2D");
}

/*
 * Function: SendRecvCellData3D
 * Usage: SendRecvCellData3D(grid->s,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the 3D cell data back and forth between
 * processors.  
 *
 */
void SendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, 
			MPI_Comm comm)
{
  int k, n, nstart, neigh, neighproc, *Nsend, *Nrecv, *num_send, *num_recv;
  REAL **recv, **send;
  MPI_Status status;

  num_send=grid->num_cells_send;
  num_recv=grid->num_cells_recv;

  recv = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvCellData3D");
  send = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvCellData3D");
  Nsend = (int *)SunMalloc(grid->Nneighs*sizeof(int),"SendRecvCellData3D");
  Nrecv = (int *)SunMalloc(grid->Nneighs*sizeof(int),"SendRecvCellData3D");

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    Nsend[neigh] = 0;
    Nrecv[neigh] = 0;
    for(n=0;n<num_send[neigh];n++) 
      Nsend[neigh]+=grid->Nk[grid->cell_send[neigh][n]];
    for(n=0;n<num_recv[neigh];n++) 
      Nrecv[neigh]+=grid->Nk[grid->cell_recv[neigh][n]];

    send[neigh] = (REAL *)SunMalloc(Nsend[neigh]*sizeof(REAL),"SendRecvCellData3D");
    recv[neigh] = (REAL *)SunMalloc(Nrecv[neigh]*sizeof(REAL),"SendRecvCellData3D");

    nstart=0;
    for(n=0;n<num_send[neigh];n++) {
      for(k=0;k<grid->Nk[grid->cell_send[neigh][n]];k++) 
	send[neigh][nstart+k]=celldata[grid->cell_send[neigh][n]][k];
      nstart+=grid->Nk[grid->cell_send[neigh][n]];
    }

    MPI_Send((void *)(send[neigh]),Nsend[neigh],MPI_DOUBLE,neighproc,1,comm); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Recv((void *)(recv[neigh]),Nrecv[neigh],MPI_DOUBLE,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<num_recv[neigh];n++) {
      for(k=0;k<grid->Nk[grid->cell_recv[neigh][n]];k++) 
	celldata[grid->cell_recv[neigh][n]][k]=recv[neigh][nstart+k];
      nstart+=grid->Nk[grid->cell_recv[neigh][n]];
    }
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    SunFree(send[neigh],Nsend[neigh]*sizeof(REAL),"SendRecvCellData3D");
    SunFree(recv[neigh],Nrecv[neigh]*sizeof(REAL),"SendRecvCellData3D");
  }
  SunFree(send,grid->Nneighs*sizeof(REAL *),"SendRecvCellData3D");
  SunFree(recv,grid->Nneighs*sizeof(REAL *),"SendRecvCellData3D");
  SunFree(Nsend,grid->Nneighs*sizeof(int),"SendRecvCellData3D");
  SunFree(Nrecv,grid->Nneighs*sizeof(int),"SendRecvCellData3D");
}

/*
 * Function: SendRecvWData
 * Usage: SendRecvWData(grid->w,grid,myproc,comm);
 * -----------------------------------------------
 * This function will transfer the 3D w data back and forth between
 * processors.  
 *
 */
void SendRecvWData(REAL **celldata, gridT *grid, int myproc, 
		   MPI_Comm comm)
{
  int k, n, nstart, neigh, neighproc, *Nsend, *Nrecv, *num_send, *num_recv;
  REAL **recv, **send;
  MPI_Status status;

  num_send=grid->num_cells_send;
  num_recv=grid->num_cells_recv;

  recv = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvWData");
  send = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvWData");
  Nsend = (int *)SunMalloc(grid->Nneighs*sizeof(int),"SendRecvWData");
  Nrecv = (int *)SunMalloc(grid->Nneighs*sizeof(int),"SendRecvWData");

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    Nsend[neigh] = 0;
    Nrecv[neigh] = 0;
    for(n=0;n<num_send[neigh];n++) 
      Nsend[neigh]+=(1+grid->Nk[grid->cell_send[neigh][n]]);
    for(n=0;n<num_recv[neigh];n++) 
      Nrecv[neigh]+=(1+grid->Nk[grid->cell_recv[neigh][n]]);

    send[neigh] = (REAL *)SunMalloc(Nsend[neigh]*sizeof(REAL),"SendRecvWData");
    recv[neigh] = (REAL *)SunMalloc(Nrecv[neigh]*sizeof(REAL),"SendRecvWData");

    nstart=0;
    for(n=0;n<num_send[neigh];n++) {
      for(k=0;k<(1+grid->Nk[grid->cell_send[neigh][n]]);k++) 
	send[neigh][nstart+k]=celldata[grid->cell_send[neigh][n]][k];
      nstart+=(1+grid->Nk[grid->cell_send[neigh][n]]);
    }

    MPI_Send((void *)(send[neigh]),Nsend[neigh],MPI_DOUBLE,neighproc,1,comm); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Recv((void *)(recv[neigh]),Nrecv[neigh],MPI_DOUBLE,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<num_recv[neigh];n++) {
      for(k=0;k<(1+grid->Nk[grid->cell_recv[neigh][n]]);k++) 
	celldata[grid->cell_recv[neigh][n]][k]=recv[neigh][nstart+k];
      nstart+=(1+grid->Nk[grid->cell_recv[neigh][n]]);
    }
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    SunFree(send[neigh],Nsend[neigh]*sizeof(REAL),"SendRecvWData");
    SunFree(recv[neigh],Nrecv[neigh]*sizeof(REAL),"SendRecvWData");
  }
  SunFree(send,grid->Nneighs*sizeof(REAL *),"SendRecvWData");
  SunFree(recv,grid->Nneighs*sizeof(REAL *),"SendRecvWData");
  SunFree(Nsend,grid->Nneighs*sizeof(int),"SendRecvWData");
  SunFree(Nrecv,grid->Nneighs*sizeof(int),"SendRecvWData");
}

/*
 * Function: SendRecvEdgeData3D
 * Usage: SendRecvEdgeData3D(grid->u,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the 3D edge data back and forth between
 * processors.  
 *
 */
void SendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, 
			MPI_Comm comm)
{
  int k, n, nstart, neigh, neighproc, *Nsend, *Nrecv, *num_send, *num_recv;
  REAL **recv, **send;
  MPI_Status status;

  num_send=grid->num_edges_send;
  num_recv=grid->num_edges_recv;

  recv = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvEdgeData3D");
  send = (REAL **)SunMalloc(grid->Nneighs*sizeof(REAL *),"SendRecvEdgeData3D");
  Nsend = (int *)SunMalloc(grid->Nneighs*sizeof(int),"SendRecvEdgeData3D");
  Nrecv = (int *)SunMalloc(grid->Nneighs*sizeof(int),"SendRecvEdgeData3D");

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    Nsend[neigh] = 0;
    Nrecv[neigh] = 0;
    for(n=0;n<num_send[neigh];n++) 
      Nsend[neigh]+=grid->Nke[grid->edge_send[neigh][n]];
    for(n=0;n<num_recv[neigh];n++) 
      Nrecv[neigh]+=grid->Nke[grid->edge_recv[neigh][n]];

    send[neigh] = (REAL *)SunMalloc(Nsend[neigh]*sizeof(REAL),"SendRecvEdgeData3D");
    recv[neigh] = (REAL *)SunMalloc(Nrecv[neigh]*sizeof(REAL),"SendRecvEdgeData3D");

    nstart=0;
    for(n=0;n<num_send[neigh];n++) {
      for(k=0;k<grid->Nke[grid->edge_send[neigh][n]];k++) 
	send[neigh][nstart+k]=edgedata[grid->edge_send[neigh][n]][k];
      nstart+=grid->Nke[grid->edge_send[neigh][n]];
    }

    MPI_Send((void *)(send[neigh]),Nsend[neigh],MPI_DOUBLE,neighproc,1,comm); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Recv((void *)(recv[neigh]),Nrecv[neigh],MPI_DOUBLE,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<num_recv[neigh];n++) {
      for(k=0;k<grid->Nke[grid->edge_recv[neigh][n]];k++) {
	edgedata[grid->edge_recv[neigh][n]][k]=recv[neigh][nstart+k];
      }
      nstart+=grid->Nke[grid->edge_recv[neigh][n]];
    }
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    SunFree(send[neigh],Nsend[neigh]*sizeof(REAL),"SendRecvEdgeData3D");
    SunFree(recv[neigh],Nrecv[neigh]*sizeof(REAL),"SendRecvEdgeData3D");
  }
  SunFree(send,grid->Nneighs*sizeof(REAL *),"SendRecvEdgeData3D");
  SunFree(recv,grid->Nneighs*sizeof(REAL *),"SendRecvEdgeData3D");
  SunFree(Nsend,grid->Nneighs*sizeof(int),"SendRecvEdgeData3D");
  SunFree(Nrecv,grid->Nneighs*sizeof(int),"SendRecvEdgeData3D");
}

void CheckCommunicateCells(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm)
{
  int n, neigh, neighproc, noncontig, flag;
  int **recv, **send;
  MPI_Status status;

  recv = (int **)SunMalloc(maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateCells");
  send = (int **)SunMalloc(maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateCells");

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];
    send[neigh] = (int *)SunMalloc(localgrid->num_cells_send[neigh]*sizeof(int),"CheckCommunicateCells");
    recv[neigh] = (int *)SunMalloc(localgrid->num_cells_recv[neigh]*sizeof(int),"CheckCommunicateCells");
    for(n=0;n<localgrid->num_cells_send[neigh];n++)
      send[neigh][n]=localgrid->mnptr[localgrid->cell_send[neigh][n]];

    MPI_Send((void *)(send[neigh]),localgrid->num_cells_send[neigh],
	     MPI_INT,neighproc,1,comm); 
  }

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];
    MPI_Recv((void *)(recv[neigh]),localgrid->num_cells_recv[neigh],
	     MPI_INT,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  noncontig=0;
  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    for(n=0;n<localgrid->num_cells_recv[neigh];n++)
      if(localgrid->mnptr[localgrid->cell_recv[neigh][n]] !=
	 recv[neigh][n]) {
	printf("Warning!  Non-contiguous cell data! %d <-> %d\n",myproc,maingrid->neighs[myproc][neigh]);
	noncontig=1;
	break;
      }
  }

  if(noncontig) {
    for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
      flag=1;
      for(n=0;n<localgrid->num_cells_recv[neigh];n++)
	if(IsMember(localgrid->mnptr[localgrid->cell_recv[neigh][n]],recv[neigh],
		    localgrid->num_cells_recv[neigh])==-1) {
	  flag=0;
	  break;
	}
      if(!flag)
	printf("Warning!  Incorrect cell data! %d <-> %d\n",myproc,maingrid->neighs[myproc][neigh]);
      else
	printf("Cell data is non-contiguous but correctly transferred. %d <-> %d\n",
	       myproc,maingrid->neighs[myproc][neigh]);
    }

    for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
      printf("%d <-> %d, Sent: %d, Received: %d\n",myproc,maingrid->neighs[myproc][neigh],
	     localgrid->num_cells_send[neigh],localgrid->num_cells_recv[neigh]);
      printf("Sent: ");
      for(n=0;n<localgrid->num_cells_send[neigh];n++) printf("%d ",send[neigh][n]);
      printf("\n");
      printf("Rec: ");
      for(n=0;n<localgrid->num_cells_recv[neigh];n++) printf("%d ",recv[neigh][n]);
      printf("\n");
      printf("Loc: ");
      for(n=0;n<localgrid->num_cells_recv[neigh];n++) 
	printf("%d ",localgrid->mnptr[localgrid->cell_recv[neigh][n]]);
      printf("\n");
    }
  } else 
    printf("Processor %d, Cell data transferred correctly.\n", myproc);

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    SunFree(send[neigh],localgrid->num_cells_send[neigh]*sizeof(int),"CheckCommunicateCells");
    SunFree(recv[neigh],localgrid->num_cells_recv[neigh]*sizeof(int),"CheckCommunicateCells");
  }
  SunFree(send,maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateCells");
  SunFree(recv,maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateCells");
}

void CheckCommunicateEdges(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm)
{
  int n, neigh, neighproc, noncontig, flag;
  int **recv, **send;
  MPI_Status status;

  recv = (int **)SunMalloc(maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateEdges");
  send = (int **)SunMalloc(maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateEdges");

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];

    send[neigh] = (int *)SunMalloc(localgrid->num_edges_send[neigh]*sizeof(int),"CheckCommunicateEdges");
    recv[neigh] = (int *)SunMalloc(localgrid->num_edges_recv[neigh]*sizeof(int),"CheckCommunicateEdges");

    for(n=0;n<localgrid->num_edges_send[neigh];n++)
      send[neigh][n]=localgrid->eptr[localgrid->edge_send[neigh][n]];

    MPI_Send((void *)(send[neigh]),localgrid->num_edges_send[neigh],
	     MPI_INT,neighproc,1,comm); 
  }
  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];
    MPI_Recv((void *)(recv[neigh]),localgrid->num_edges_recv[neigh],
	     MPI_INT,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);


  noncontig=0;
  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    for(n=0;n<localgrid->num_edges_recv[neigh];n++)
      if(localgrid->eptr[localgrid->edge_recv[neigh][n]] !=
	 recv[neigh][n]) {
	printf("Warning!  Non-contiguous edge data! %d <-> %d\n",myproc,maingrid->neighs[myproc][neigh]);
	noncontig=1;
	break;
      }
  }

  if(noncontig) {
    for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
      flag=1;
      for(n=0;n<localgrid->num_edges_recv[neigh];n++)
	if(IsMember(localgrid->eptr[localgrid->edge_recv[neigh][n]],recv[neigh],
		    localgrid->num_edges_recv[neigh])==-1) {
	  flag=0;
	  break;
	}
      if(!flag)
	printf("Warning!  Incorrect edge data! %d <-> %d\n",myproc,maingrid->neighs[myproc][neigh]);
      else
	printf("Edge data is non-contiguous but correctly transferred. %d <-> %d\n",
	       myproc,maingrid->neighs[myproc][neigh]);
    }

    for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
      printf("%d <-> %d, Sent: %d, Received: %d\n",myproc,maingrid->neighs[myproc][neigh],
	     localgrid->num_edges_send[neigh],localgrid->num_edges_recv[neigh]);
      printf("Sent: ");
      for(n=0;n<localgrid->num_edges_send[neigh];n++) printf("%d ",send[neigh][n]);
      printf("\n");
      printf("Rec: ");
      for(n=0;n<localgrid->num_edges_recv[neigh];n++) printf("%d ",recv[neigh][n]);
      printf("\n");
      printf("Loc: ");
      for(n=0;n<localgrid->num_edges_recv[neigh];n++) 
	printf("%d ",localgrid->eptr[localgrid->edge_recv[neigh][n]]);
      printf("\n");
    }
  } else
    printf("Processor %d, Edge data transferred correctly.\n",myproc);

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    SunFree(send[neigh],localgrid->num_edges_send[neigh]*sizeof(int),"CheckCommunicateEdges");
    SunFree(recv[neigh],localgrid->num_edges_recv[neigh]*sizeof(int),"CheckCommunicateEdges");
  }
  SunFree(send,maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateEdges");
  SunFree(recv,maingrid->numneighs[myproc]*sizeof(int *),"CheckCommunicateEdges");
}
  
/*
 * Function: InitMainGrid
 * Usage: InitMainGrid(&grid);
 * ---------------------------
 * Initialize the main grid struct.
 *
 */
void InitMainGrid(gridT **grid, int Np, int Ne, int Nc)
{
  *grid = (gridT *)SunMalloc(sizeof(gridT),"InitMainGrid");
  
  // Number of cells
  //  (*grid)->Nc = getsize(CELLSFILE);
  (*grid)->Nc = Nc;
  // Number of edges
  //(*grid)->Ne = getsize(EDGEFILE);
  (*grid)->Ne = Ne;
  // Number of points defining the vertices of the polygons
  //(*grid)->Np = getsize(POINTSFILE);
  (*grid)->Np = Np;

  // (x,y) coordinates of vertices
  (*grid)->xp = (REAL *)SunMalloc((*grid)->Np*sizeof(REAL),"InitMainGrid");
  (*grid)->yp = (REAL *)SunMalloc((*grid)->Np*sizeof(REAL),"InitMainGrid");
  // (x,y) coordinates of voronoi points
  (*grid)->xv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"InitMainGrid");
  (*grid)->yv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"InitMainGrid");

  // Pointers to xp,yp coordinates that define two endpoints of faces (0<edges<Np)
  (*grid)->edges = (int *)SunMalloc(NUMEDGECOLUMNS*(*grid)->Ne*sizeof(int),"InitMainGrid");
  // Pointers to xv,yv coordinates that define two endpoints of voronoi edges (0<grad<Np)
  (*grid)->grad = (int *)SunMalloc(2*(*grid)->Ne*sizeof(int),"InitMainGrid");
  // gradf contains the face number (0<gradf<NFACES) of the given cell that corresponds to
  // the given edge.
  (*grid)->gradf = (int *)SunMalloc(2*(*grid)->Ne*sizeof(int),"InitMainGrid");
  // Pointers to xp,yp coordinates of vertices that make up polygons (0<cells<Np)
  (*grid)->cells = (int *)SunMalloc(NFACES*(*grid)->Nc*sizeof(int),"InitMainGrid");
  // Pointers to neighboring cells (0<neigh<Nc)
  (*grid)->neigh = (int *)SunMalloc(NFACES*(*grid)->Nc*sizeof(int),"InitMainGrid");
  // Dot product of unique normal with outward normal
  (*grid)->normal = (int *)SunMalloc(NFACES*(*grid)->Ne*sizeof(int),"InitMainGrid");
  // Indices of voronoi edges to cells
  (*grid)->grad = (int *)SunMalloc(2*(*grid)->Ne*sizeof(int),"InitMainGrid");
  // Indices of pointers to faces of each cell
  (*grid)->face = (int *)SunMalloc(NFACES*(*grid)->Nc*sizeof(int),"InitMainGrid");
  // Indices to edges for momentum control volume
  (*grid)->eneigh = (int *)SunMalloc(2*(NFACES-1)*(*grid)->Ne*sizeof(int),"InitMainGrid");

  // Depth at Voronoi points
  (*grid)->dv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"InitMainGrid");
  // Weights for partitioning cell graph
  (*grid)->vwgt = (int *)SunMalloc((*grid)->Nc*sizeof(int),"InitMainGrid");

  // Assigned cell partition
  (*grid)->part = (int *)SunMalloc((*grid)->Nc*sizeof(int),"InitMainGrid");
  // For ordering and making sure boundary edge data is contiguous.
  (*grid)->order = (int *)SunMalloc((*grid)->Ne*sizeof(int),"InitMainGrid");
  // Edge markers
  (*grid)->mark = (int *)SunMalloc((*grid)->Ne*sizeof(int),"InitMainGrid");
  // Stores the indices to the start and end points of the adjncy array
  (*grid)->xadj = (int *)SunMalloc(((*grid)->Nc+1)*sizeof(int),"InitMainGrid");
  (*grid)->vtxdist = (int *)SunMalloc(VTXDISTMAX*sizeof(int),"InitMainGrid");
}

/*
 * Function: ReadMainGrid
 * Usage: ReadmainGrid(grid,myproc);
 * ---------------------------------
 * Read in the cell data. If the CORRECTVORONOI
 * flag is set to 1 then the voronoi points are set to be the centroids of
 * the cells.
 *
 * The files contain the data as follows:
 * POINTSFILE (x,y) list of points corresponding to the vertices of the cells, which also contains
 *            a third column with markers to those points.
 * CELLSFILE voronoi_x voronoi_y cell_pt1 cell_pt2 cell_pt3 neigh_pt1 neigh_pt2 neigh_pt3
 * EDGEFILE list of indices to points in POINTSFILE (always 2 columns + edge markers = 3 columns)
 *
 */
void ReadMainGrid(gridT *grid, int myproc)
{
  int j, n, nei, nf;
  char str[BUFFERLENGTH];
  FILE *ifile;

  ifile = MPI_FOpen(POINTSFILE,"r","ReadMainGrid",myproc);
  for(n=0;n<grid->Np;n++) {
    grid->xp[n]=getfield(ifile,str);
    grid->yp[n]=getfield(ifile,str);
    getfield(ifile,str);
  }
  fclose(ifile);
  
  ifile = MPI_FOpen(EDGEFILE,"r","ReadMainGrid",myproc);
  for(n=0;n<grid->Ne;n++) {
    for(j=0;j<NUMEDGECOLUMNS-1;j++) 
      grid->edges[NUMEDGECOLUMNS*n+j]=(int)getfield(ifile,str);
    grid->mark[n]=(int)getfield(ifile,str);
    for(j=0;j<2;j++) 
      grid->grad[2*n+j]=(int)getfield(ifile,str);
  }
  fclose(ifile);

  ifile = MPI_FOpen(CELLSFILE,"r","ReadMainGrid",myproc);
  for(n=0;n<grid->Nc;n++) {
    grid->xv[n] = getfield(ifile,str);
    grid->yv[n] = getfield(ifile,str);
    for(nf=0;nf<NFACES;nf++)
      grid->cells[n*NFACES+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<NFACES;nf++) {
      if((nei=(int)getfield(ifile,str))!=-1)
	grid->neigh[n*NFACES+nf]=nei;
      else
	grid->neigh[n*NFACES+nf]=-1;
    }
  }
  
}

/*
 * Function: ReadFileNames
 * Usage: ReadFileNames(myproc);
 * -----------------------------
 * Reads the names of the files containing the grid data
 * from the file defined by GRIDDATAFILELIST in suntans.h
 * 
 */
void ReadFileNames(int myproc)
{
  MPI_GetFile(POINTSFILE,DATAFILE,"points","OpenFiles",myproc);
  MPI_GetFile(EDGEFILE,DATAFILE,"edges","OpenFiles",myproc);
  MPI_GetFile(CELLSFILE,DATAFILE,"cells","OpenFiles",myproc);
  MPI_GetFile(INPUTDEPTHFILE,DATAFILE,"depth","OpenFiles",myproc);
  MPI_GetFile(CELLCENTEREDFILE,DATAFILE,"celldata","OpenFiles",myproc);
  MPI_GetFile(EDGECENTEREDFILE,DATAFILE,"edgedata","OpenFiles",myproc);
  MPI_GetFile(VERTSPACEFILE,DATAFILE,"vertspace","OpenFiles",myproc);
  MPI_GetFile(TOPOLOGYFILE,DATAFILE,"topology","OpenFiles",myproc);
}

static void ReadDepth(gridT *grid, int myproc) {
  int n;
  char str[BUFFERLENGTH];
  FILE *fid;
  sprintf(str,"%s-voro",INPUTDEPTHFILE);

  fid = MPI_FOpen(str,"r","ReadDepth",myproc);
  for(n=0;n<grid->Nc;n++) {
    getfield(fid,str);
    getfield(fid,str);
    grid->dv[n]=getfield(fid,str);
  }
  fclose(fid);
}
  
void GetDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm)
{
  int n, maxgridweight=1000, IntDepth, Nkmax, stairstep, fixdzz, kount=0;
  REAL mindepth, maxdepth, maxdepth0, minimum_depth, *dz;

  Nkmax = MPI_GetValue(DATAFILE,"Nkmax","GetDepth",myproc);
  stairstep = MPI_GetValue(DATAFILE,"stairstep","GetDepth",myproc);

  dz = (REAL *)SunMalloc(Nkmax*sizeof(REAL),"GetDepth");

  maxdepth=0.0;
  mindepth=INFTY;
  IntDepth=(int)MPI_GetValue(DATAFILE,"IntDepth","GetDepth",myproc);
  minimum_depth=(REAL)MPI_GetValue(DATAFILE,"minimum_depth","GetDepth",myproc);
  //fixdzz=(REAL)MPI_GetValue(DATAFILE,"fixdzz","GetDepth",myproc);
  grid->dzsmall = (REAL)MPI_GetValue(DATAFILE,"dzsmall","FixDZZ",myproc);
  
  if(IntDepth==1) 
    InterpDepth(grid,myproc,numprocs,comm);
  else if(IntDepth==2) 
    ReadDepth(grid,myproc);
  else {
    for(n=0;n<grid->Nc;n++) {
      grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
    }
  }

  for(n=0;n<grid->Nc;n++) 
    if(grid->dv[n]>maxdepth)
      maxdepth = grid->dv[n];
  maxdepth0 = maxdepth;

  GetDZ(dz,maxdepth,maxdepth,Nkmax,myproc);  

  //if(!stairstep && fixdzz) 
  //  FixDZZ(grid,maxdepth,Nkmax,fixdzz,myproc);

  if(minimum_depth!=0) {
    if(minimum_depth>0) {
      printf("Setting minimum depth to %f\n",minimum_depth);
      for(n=0;n<grid->Nc;n++) {
	if(grid->dv[n]<minimum_depth) {
	  grid->dv[n]=minimum_depth;
	  kount++;
	}
      }
      if(VERBOSE>0 && myproc==0 && kount>0) 
	printf("Increased the depth to the minimum set value of %.2f %d times.\n",minimum_depth,kount);
    } else if(minimum_depth<0) {
      for(n=0;n<grid->Nc;n++) {
	if(grid->dv[n]<dz[0]) {
	  grid->dv[n]=dz[0];
	  kount++;
	}
      }
      if(VERBOSE>0 && myproc==0 && kount>0) 
	printf("Increased the depth to the minimum set value of dz[0]=%.2f %d times.\n",dz[0],kount);
    }
  }
  
  for(n=0;n<grid->Nc;n++) {
    if(grid->dv[n]>maxdepth)
      maxdepth = grid->dv[n];
    if(grid->dv[n]<mindepth)
      mindepth = grid->dv[n];
  }

  if(maxdepth!=maxdepth0 && myproc==0) printf("WARNING!!!!!!!!!  Maximum depth changed after fixing vertical spacing...\n");

  if(Nkmax>1) {
    if(mindepth!=maxdepth) 
      for(n=0;n<grid->Nc;n++) {
	grid->vwgt[n]=(int)(maxgridweight*(float)GetNk(dz,grid->dv[n],Nkmax)/(float)Nkmax);
	//	grid->vwgt[n] = (int)(maxgridweight*(grid->dv[n]-mindepth)/(maxdepth-mindepth));
    } else
      for(n=0;n<grid->Nc;n++) 
	grid->vwgt[n] = maxgridweight;
  } else {
      for(n=0;n<grid->Nc;n++) 
	grid->vwgt[n] = Nkmax;
  }

  if(VERBOSE>3) 
    for(n=0;n<grid->Nc;n++) 
      printf("grid->vwgt[%d]=%d\n",n,grid->vwgt[n]);

  SunFree(dz,Nkmax*sizeof(REAL),"GetDepth");
}

void CreateCellGraph(gridT *grid)
{
  int n, j, Nge;

  Nge = 0;
  for(n=0;n<grid->Nc;n++) 
    for(j=0;j<NFACES;j++)
      if(grid->neigh[NFACES*n+j]!=-1) Nge++;
  grid->Nge = Nge/2;
  grid->adjncy = (int *)SunMalloc(Nge*sizeof(int),"CreateCellGraph");

  Nge=0;
  grid->xadj[0]=Nge;
  for(n=0;n<grid->Nc;n++) { 
    for(j=0;j<NFACES;j++) 
      if(grid->neigh[NFACES*n+j]!=-1) {
	grid->adjncy[Nge++]=grid->neigh[NFACES*n+j];
      }
    grid->xadj[n+1]=Nge;
    }
  if(VERBOSE>3) {
    printf("xadj: ");
    for(n=0;n<grid->Nc+1;n++) 
      printf("%d ",grid->xadj[n]);
    printf("\n");
    printf("adjncy: ");
    for(j=0;j<Nge;j++) 
      printf("%d ",grid->adjncy[j]);
    printf("\n");
  }
}

void CreateEdgeGraph(gridT *grid)
{
  int n, j, Nge;

  Nge = 0;
  for(n=0;n<grid->Ne;n++) 
    for(j=0;j<2*(NFACES-1);j++)
      if(grid->eneigh[2*(NFACES-1)*n+j]!=-1) Nge++;
  grid->Nge = Nge/2;
  grid->adjncy = (int *)SunMalloc(Nge*sizeof(int),"CreateEdgeGraph");

  Nge=0;
  grid->xadj[0]=Nge;
  for(n=0;n<grid->Ne;n++) { 
    for(j=0;j<2*(NFACES-1);j++) 
      if(grid->eneigh[2*(NFACES-1)*n+j]!=-1) {
	grid->adjncy[Nge++]=grid->eneigh[2*(NFACES-1)*n+j];
      }
    grid->xadj[n+1]=Nge;
    }
  if(VERBOSE>3) {
    printf("xadj: ");
    for(n=0;n<grid->Ne+1;n++) 
      printf("%d ",grid->xadj[n]);
    printf("\n");
    printf("adjncy: ");
    for(j=0;j<Nge;j++) 
      printf("%d ",grid->adjncy[j]);
    printf("\n");
  }
}

static void CreateFaceArray(int *grad, int *gradf, int *neigh, int *face, int Nc, int Ne)
{
  int n, j, nf, nc, nc1, nc2, nb;

  for(n=0;n<Nc;n++)
    for(nf=0;nf<NFACES;nf++)
      face[n*NFACES+nf]=-1;

  for(n=0;n<Ne;n++) {
    nc1 = grad[2*n];
    nc2 = grad[2*n+1];
    if(nc1!=-1 && nc2!=-1)
      for(nf=0;nf<NFACES;nf++) {
	nb = neigh[nc1*NFACES+nf];
	if(nb==nc2) {
	  face[nc1*NFACES+nf]=n;
	  gradf[2*n] = nf;
	}
	nb = neigh[nc2*NFACES+nf];
	if(nb==nc1) {
	  face[nc2*NFACES+nf]=n;
	  gradf[2*n+1] = nf;
	}
      }
  }

  for(n=0;n<Ne;n++)
    for(j=0;j<2;j++) {
      nc = grad[2*n+j];
      if(nc != -1)
	if(IsMember(n,&(face[nc*NFACES]),NFACES)==-1)
	  for(nf=0;nf<NFACES;nf++) 
	    if(face[nc*NFACES+nf] == -1) {
	      face[nc*NFACES+nf]=n;
	      gradf[2*n+j]=nf;
	      break;
	    }
    }
}

static void CreateNormalArray(int *grad, int *face, int *normal, int Nc)
{
  int n, nf;
  for(n=0;n<Nc;n++){
    for(nf=0;nf<NFACES;nf++)
      if(n==grad[2*face[NFACES*n+nf]])
	normal[n*NFACES+nf]=-1;
      else
	normal[n*NFACES+nf]=1;
  }
}

void Connectivity(gridT *grid, int myproc)
{
  int n, nf, ng, ne;

  /* Create the face array, which contains indexes to the edges of
     each cell */
  if(myproc==0 && VERBOSE>2) printf("Creating face array and normals\n");
  CreateFaceArray(grid->grad,grid->gradf,grid->neigh,grid->face,grid->Nc,grid->Ne);
  CreateNormalArray(grid->grad,grid->face,grid->normal,grid->Nc);

  /* Create the edge connectivity array eneigh, which points to the 2(NFACES-1)
     neighbors that comprise the edges for the momentum control volume */
  if(myproc==0 && VERBOSE>2) printf("Creating edge connectivity array...\n");
  for(n=0;n<grid->Ne;n++) {
    for(ne=0;ne<2*(NFACES-1);ne++)
      grid->eneigh[2*(NFACES-1)*n+ne]=-1;
    ne = 0;
    for(ng=0;ng<2;ng++) 
      if(grid->grad[2*n+ng]!=-1)
	for(nf=0;nf<NFACES;nf++) {
	  if(grid->face[NFACES*grid->grad[2*n+ng]+nf] != n) 
	    grid->eneigh[2*(NFACES-1)*n+ne++]=
	      grid->face[NFACES*grid->grad[2*n+ng]+nf];
	}
  }
}

/*
 * Function: IsBoundaryCell
 * Usage: IsBoundaryCell(i,maingrid,myproc);
 * -----------------------------------------
 * Determines whether or not the global cell index
 * i is a boundary cell. 
 * Returns:
 *   3 if this is a ghost cell.
 *   2 if this is an interproc boundary cell.
 *   1 if this is a non-interproc boundary cell 
 *     (not a ghost cell), i.e. free-surface specified 
 *      containing an edge with mark=3
 *   0 otherwise.
 *  -1 for error.
 *
 * If it is both a ghost cell and a boundary cell or
 * both a 
 *
 */
int IsBoundaryCell(int mgptr, gridT *maingrid, int myproc)
{
  int nf, nei;

  if(mgptr >= maingrid->Nc || mgptr < 0) {
    printf("Error in IsBoundaryCell: index out of bounds!\n");
    return -1;
  }
  if(maingrid->part[mgptr]==myproc) {
    for(nf=0;nf<NFACES;nf++) {
      if(maingrid->mark[maingrid->face[mgptr*NFACES+nf]]==3)
	return 1;
    }
    for(nf=0;nf<NFACES;nf++) {
      nei = maingrid->neigh[mgptr*NFACES+nf];
      if(nei != -1)
	if(maingrid->part[nei] != myproc)
	  return 2;
    }
    return 0;
  } else
    return 3;
}

/*
 * Function: OutputData
 * Usage: OutputData(maingrid,localgrid,myproc);
 * ---------------------------------------------
 * Outputs the required grid data.
 *
 */
static void OutputData(gridT *maingrid, gridT *grid, int myproc, int numprocs)
{
  int j, n, nf, neigh, Np=maingrid->Np, Nc=grid->Nc, Ne=grid->Ne;
  char str[BUFFERLENGTH];
  FILE *ofile;

  if(TRIANGULATE && myproc==0) {
    if(VERBOSE>2) printf("Outputting %s...\n",POINTSFILE);
    ofile = MPI_FOpen(POINTSFILE,"w","OutputData",myproc);
    for(j=0;j<Np;j++)
      fprintf(ofile,"%f %f 0\n",maingrid->xp[j],maingrid->yp[j]);
    fclose(ofile);
  
    ofile = MPI_FOpen(EDGEFILE,"w","OutputData",myproc);
    for(n=0;n<maingrid->Ne;n++) {
      for(j=0;j<NUMEDGECOLUMNS-1;j++) 
	fprintf(ofile,"%d ",maingrid->edges[NUMEDGECOLUMNS*n+j]);
      fprintf(ofile,"%d ",maingrid->mark[n]);
      for(j=0;j<2;j++) 
	fprintf(ofile,"%d ",maingrid->grad[2*n+j]);
      fprintf(ofile,"\n");
    }
    fclose(ofile);
    
    ofile = MPI_FOpen(CELLSFILE,"w","ReadMainGrid",myproc);
    for(n=0;n<maingrid->Nc;n++) {
      fprintf(ofile,"%f %f ",maingrid->xv[n],maingrid->yv[n]);
      for(nf=0;nf<NFACES;nf++)
	fprintf(ofile,"%d ",maingrid->cells[n*NFACES+nf]);
      for(nf=0;nf<NFACES;nf++) 
	fprintf(ofile,"%d ",maingrid->neigh[n*NFACES+nf]);
      fprintf(ofile,"\n");
    }
    fclose(ofile);
  }

  sprintf(str,"%s.%d",CELLSFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);
  ofile = MPI_FOpen(str,"w","OutputData",myproc);
  for(j=0;j<grid->Nc;j++) {
    fprintf(ofile,"%e %e ",grid->xv[j],grid->yv[j]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->cells[j*NFACES+nf]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->neigh[j*NFACES+nf]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  sprintf(str,"%s.%d",EDGEFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);

  ofile = MPI_FOpen(str,"w","OutputData",myproc);
  for(j=0;j<grid->Ne;j++) {
    for(nf=0;nf<2;nf++)
      fprintf(ofile,"%d ",grid->edges[j*NUMEDGECOLUMNS+nf]);
    fprintf(ofile,"%d ",grid->mark[j]);
    for(nf=0;nf<2;nf++)
      fprintf(ofile,"%d ",grid->grad[2*j+nf]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);

  ofile = MPI_FOpen(str,"w","OutputData",myproc);
  for(n=0;n<Nc;n++) {
    fprintf(ofile,"%e %e %e %e %d ",grid->xv[n],grid->yv[n],
	    grid->Ac[n],grid->dv[n],grid->Nk[n]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->face[NFACES*n+nf]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->neigh[NFACES*n+nf]);
    for(nf=0;nf<NFACES;nf++) 
      fprintf(ofile,"%d ",grid->normal[NFACES*n+nf]);
    for(nf=0;nf<NFACES;nf++) 
      fprintf(ofile,"%e ",grid->def[NFACES*n+nf]);
    // Removed this line since it is never reloaded...
    //fprintf(ofile,"%d ",grid->mnptr[n]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);

  ofile = MPI_FOpen(str,"w","OutputData",myproc);
  for(n=0;n<Ne;n++) {
    // Removed the two zeros because they are redundant
    // fprintf(ofile,"%e %e 0 0 %e %e %d %d %d %d ",
    fprintf(ofile,"%e %e %e %e %e %e %d %d %d %d %d %d %d ",
	    grid->df[n],grid->dg[n],grid->n1[n],grid->n2[n],grid->xe[n],grid->ye[n],
	    grid->Nke[n],grid->Nkc[n],grid->grad[2*n],grid->grad[2*n+1],
	    grid->gradf[2*n],grid->gradf[2*n+1],grid->mark[n]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  sprintf(str,"%s.%d",TOPOLOGYFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);

  ofile = MPI_FOpen(str,"w","OutputData",myproc);
  fprintf(ofile,"%d %d\n",numprocs,grid->Nneighs);
  for(neigh=0;neigh<grid->Nneighs;neigh++) 
    fprintf(ofile,"%d ",grid->myneighs[neigh]);
  fprintf(ofile,"\n");
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    fprintf(ofile,"%d %d %d %d\n",
	    grid->num_cells_send[neigh],
	    grid->num_cells_recv[neigh],
	    grid->num_edges_send[neigh],
	    grid->num_edges_recv[neigh]);
    for(n=0;n<grid->num_cells_send[neigh];n++)
      fprintf(ofile,"%d ",grid->cell_send[neigh][n]);
    fprintf(ofile,"\n");
    for(n=0;n<grid->num_cells_recv[neigh];n++)
      fprintf(ofile,"%d ",grid->cell_recv[neigh][n]);
    fprintf(ofile,"\n");
    for(n=0;n<grid->num_edges_send[neigh];n++)
      fprintf(ofile,"%d ",grid->edge_send[neigh][n]);
    fprintf(ofile,"\n");
    for(n=0;n<grid->num_edges_recv[neigh];n++)
      fprintf(ofile,"%d ",grid->edge_recv[neigh][n]);
    fprintf(ofile,"\n");
  }
  for(n=0;n<MAXBCTYPES-1;n++) 
    fprintf(ofile,"%d ",grid->celldist[n]);
  fprintf(ofile,"\n");
  for(n=0;n<MAXMARKS-1;n++) 
    fprintf(ofile,"%d ",grid->edgedist[n]);
  fprintf(ofile,"\n");
  for(n=0;n<grid->Nc;n++) 
    fprintf(ofile,"%d ",grid->cellp[n]);
  fprintf(ofile,"\n");
  for(n=0;n<grid->Ne;n++) 
    fprintf(ofile,"%d ",grid->edgep[n]);
  fprintf(ofile,"\n");
  fclose(ofile);

  if(myproc==0 && VERBOSE>2) printf("Outputting %s...\n",VERTSPACEFILE);
  if(myproc==0) {
    sprintf(str,"%s",VERTSPACEFILE);
    ofile = MPI_FOpen(str,"w","OutputData",myproc);
    for(n=0;n<grid->Nkmax;n++)
      fprintf(ofile,"%e\n",grid->dz[n]);
    fclose(ofile);
  }
}

/*
 * Function: ReadGrid
 * Usage: ReadGrid(&grid,myproc,numprocs,comm);
 * -------------------------------------------
 * Reads the partitioned grid data and allocates space
 * for the required arrays -- must have been
 * called with the right number of processors!
 *
 */
void ReadGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm) 
{
  int neigh, n, nf, Nkmax;
  char str[BUFFERLENGTH];
  FILE *ifile;

  ReadFileNames(myproc);

  InitLocalGrid(grid);

  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  (*grid)->Nc = MPI_GetSize(str,"ReadGrid",myproc);
  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  (*grid)->Ne = MPI_GetSize(str,"ReadGrid",myproc);

  /*
   * First read in the topology file
   *
   */
  // Here check to make sure you're reading in a topology file that
  // corresponds to the right number of processors. All processors
  // need to read in the 0 topo file to check this (rather than doing an mpi_send/recv
  sprintf(str,"%s.0",TOPOLOGYFILE);
  ifile = MPI_FOpen(str,"r","ReadGrid",myproc);
  if(numprocs!=((int)getfield(ifile,str))) {
    if(myproc==0) {
      printf("Error! Topology file(s) %s.*\n",TOPOLOGYFILE);
      printf("is/are not written for %d processors!\n",numprocs);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  sprintf(str,"%s.%d",TOPOLOGYFILE,myproc);
  if(VERBOSE>2) printf("Reading %s...\n",str);
  ifile = MPI_FOpen(str,"r","ReadGrid",myproc);

  numprocs=(int)getfield(ifile,str);
  (*grid)->Nneighs=(int)getfield(ifile,str);
  (*grid)->myneighs=(int *)SunMalloc((*grid)->Nneighs*sizeof(int),"ReadGrid");
  (*grid)->num_cells_send=(int *)SunMalloc((*grid)->Nneighs*sizeof(int),"ReadGrid");
  (*grid)->num_cells_recv=(int *)SunMalloc((*grid)->Nneighs*sizeof(int),"ReadGrid");
  (*grid)->num_edges_send=(int *)SunMalloc((*grid)->Nneighs*sizeof(int),"ReadGrid");
  (*grid)->num_edges_recv=(int *)SunMalloc((*grid)->Nneighs*sizeof(int),"ReadGrid");
  (*grid)->cell_send=(int **)SunMalloc((*grid)->Nneighs*sizeof(int *),"ReadGrid");
  (*grid)->cell_recv=(int **)SunMalloc((*grid)->Nneighs*sizeof(int *),"ReadGrid");
  (*grid)->edge_send=(int **)SunMalloc((*grid)->Nneighs*sizeof(int *),"ReadGrid");
  (*grid)->edge_recv=(int **)SunMalloc((*grid)->Nneighs*sizeof(int *),"ReadGrid");

  for(neigh=0;neigh<(*grid)->Nneighs;neigh++) 
    (*grid)->myneighs[neigh]=(int)getfield(ifile,str);

  for(neigh=0;neigh<(*grid)->Nneighs;neigh++) {

    (*grid)->num_cells_send[neigh]=(int)getfield(ifile,str);
    (*grid)->num_cells_recv[neigh]=(int)getfield(ifile,str);
    (*grid)->num_edges_send[neigh]=(int)getfield(ifile,str);
    (*grid)->num_edges_recv[neigh]=(int)getfield(ifile,str);

    (*grid)->cell_send[neigh]=(int *)SunMalloc((*grid)->num_cells_send[neigh]*sizeof(int),"ReadGrid");
    (*grid)->cell_recv[neigh]=(int *)SunMalloc((*grid)->num_cells_recv[neigh]*sizeof(int),"ReadGrid");
    (*grid)->edge_send[neigh]=(int *)SunMalloc((*grid)->num_edges_send[neigh]*sizeof(int),"ReadGrid");
    (*grid)->edge_recv[neigh]=(int *)SunMalloc((*grid)->num_edges_recv[neigh]*sizeof(int),"ReadGrid");

    for(n=0;n<(*grid)->num_cells_send[neigh];n++)
      (*grid)->cell_send[neigh][n]=(int)getfield(ifile,str);
    for(n=0;n<(*grid)->num_cells_recv[neigh];n++)
      (*grid)->cell_recv[neigh][n]=(int)getfield(ifile,str);
    for(n=0;n<(*grid)->num_edges_send[neigh];n++)
      (*grid)->edge_send[neigh][n]=(int)getfield(ifile,str);
    for(n=0;n<(*grid)->num_edges_recv[neigh];n++)
      (*grid)->edge_recv[neigh][n]=(int)getfield(ifile,str);
  }

  (*grid)->celldist = (int *)SunMalloc((MAXBCTYPES-1)*sizeof(int),"ReadGrid");
  (*grid)->edgedist = (int *)SunMalloc((MAXMARKS-1)*sizeof(int),"ReadGrid");
  (*grid)->cellp = (int *)SunMalloc((*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->edgep = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid");

  for(n=0;n<MAXBCTYPES-1;n++) 
    (*grid)->celldist[n]=(int)getfield(ifile,str);
  for(n=0;n<MAXMARKS-1;n++) 
    (*grid)->edgedist[n]=(int)getfield(ifile,str);
  for(n=0;n<(*grid)->Nc;n++) 
    (*grid)->cellp[n]=(int)getfield(ifile,str);
  for(n=0;n<(*grid)->Ne;n++) 
    (*grid)->edgep[n]=(int)getfield(ifile,str);

  /*
   * Now read in cell-centered data.dat
   *
   */
  (*grid)->xv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->yv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->dv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->Ac = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");

  (*grid)->Nk = (int *)SunMalloc((*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->neigh = (int *)SunMalloc(NFACES*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->face = (int *)SunMalloc(NFACES*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->normal = (int *)SunMalloc(NFACES*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->def = (REAL *)SunMalloc(NFACES*(*grid)->Nc*sizeof(REAL),"ReadGrid");

  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Reading %s...\n",str);

  ifile = MPI_FOpen(str,"r","ReadGrid",myproc);
  for(n=0;n<(*grid)->Nc;n++) {
    (*grid)->xv[n]=getfield(ifile,str);
    (*grid)->yv[n]=getfield(ifile,str);
    (*grid)->Ac[n]=getfield(ifile,str);
    (*grid)->dv[n]=getfield(ifile,str);
    (*grid)->Nk[n]=(int)getfield(ifile,str);
    for(nf=0;nf<NFACES;nf++)
      (*grid)->face[NFACES*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<NFACES;nf++)
      (*grid)->neigh[NFACES*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<NFACES;nf++)
      (*grid)->normal[NFACES*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<NFACES;nf++)
      (*grid)->def[NFACES*n+nf]=getfield(ifile,str);
  }
  fclose(ifile);

  /*
   * Now read in edge-centered data.dat
   *
   */
  (*grid)->df = (REAL *)SunMalloc((*grid)->Ne*sizeof(REAL),"ReadGrid");
  (*grid)->dg = (REAL *)SunMalloc((*grid)->Ne*sizeof(REAL),"ReadGrid");
  (*grid)->n1 = (REAL *)SunMalloc((*grid)->Ne*sizeof(REAL),"ReadGrid");
  (*grid)->n2 = (REAL *)SunMalloc((*grid)->Ne*sizeof(REAL),"ReadGrid");
  (*grid)->xe = (REAL *)SunMalloc((*grid)->Ne*sizeof(REAL),"ReadGrid");
  (*grid)->ye = (REAL *)SunMalloc((*grid)->Ne*sizeof(REAL),"ReadGrid");

  (*grid)->Nke = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid");
  (*grid)->Nkc = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid");
  (*grid)->grad = (int *)SunMalloc(2*(*grid)->Ne*sizeof(int),"ReadGrid");
  (*grid)->gradf = (int *)SunMalloc(2*(*grid)->Ne*sizeof(int),"ReadGrid");
  (*grid)->mark = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid");

  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Reading %s...\n",str);

  ifile = MPI_FOpen(str,"r","ReadGrid",myproc);
  for(n=0;n<(*grid)->Ne;n++) {
      (*grid)->df[n] = getfield(ifile,str);
      (*grid)->dg[n] = getfield(ifile,str);
      (*grid)->n1[n] = getfield(ifile,str);
      (*grid)->n2[n] = getfield(ifile,str);
      (*grid)->xe[n] = getfield(ifile,str);
      (*grid)->ye[n] = getfield(ifile,str);
      (*grid)->Nke[n] = (int)getfield(ifile,str);
      (*grid)->Nkc[n] = (int)getfield(ifile,str);
      (*grid)->grad[2*n] = (int)getfield(ifile,str);
      (*grid)->grad[2*n+1] = (int)getfield(ifile,str);
      (*grid)->gradf[2*n] = (int)getfield(ifile,str);
      (*grid)->gradf[2*n+1] = (int)getfield(ifile,str);
      (*grid)->mark[n] = (int)getfield(ifile,str);
  }
  fclose(ifile);

  /*
   * Now read in vertical grid spacing...
   *
   */
  (*grid)->Nkmax=MPI_GetValue(DATAFILE,"Nkmax","VertGrid",myproc);
  Nkmax=MPI_GetSize(VERTSPACEFILE,"ReadGrid",myproc);
  if((*grid)->Nkmax!=Nkmax) {
    printf("Error in reading in grid data!\n");
    printf("Length of %s: %d is not consistent with what is in %s: %d\n",
	   VERTSPACEFILE,Nkmax,DATAFILE,(*grid)->Nkmax);
    MPI_Finalize();
    exit(1);
  }
  (*grid)->dz = (REAL *)SunMalloc((*grid)->Nkmax*sizeof(REAL),"ReadGrid");

  if(myproc==0 && VERBOSE>2) printf("Reading %s...\n",VERTSPACEFILE);
  ifile = MPI_FOpen(VERTSPACEFILE,"r","ReadGrid",myproc);
  for(n=0;n<(*grid)->Nkmax;n++)
    (*grid)->dz[n]=getfield(ifile,str);
  fclose(ifile);

  // These are not read in but just initialized
  (*grid)->ctop = (int *)SunMalloc((*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->ctopold = (int *)SunMalloc((*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->etop = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid");
  (*grid)->etopold = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid");

  for(n=0;n<(*grid)->Nc;n++) 
    (*grid)->ctop[n]=0;
  for(n=0;n<(*grid)->Ne;n++) 
    (*grid)->etop[n]=0;
}

/************************************************************************/
/*                                                                      */
/*                       Private Functions                              */
/*                                                                      */
/************************************************************************/

/*
 * Function: InitLocalGrid
 * Usage: InitLocalGrid(&grid);
 * ----------------------------
 * Initialize the local grid struct.
 *
 */
static void InitLocalGrid(gridT **grid)
{
  *grid = (gridT *)SunMalloc(sizeof(gridT),"InitLocalGrid");
}

/*
 * Function: FreeGrid
 * Usage: FreeGrid(grid,numprocs);
 * -------------------------------
 * Free up memory associated with the grid.
 *
 */
static void FreeGrid(gridT *grid, int numprocs)
{
  int proc;

  SunFree(grid->xp,grid->Np*sizeof(REAL),"FreeGrid");
  SunFree(grid->yp,grid->Np*sizeof(REAL),"FreeGrid");
  SunFree(grid->xv,grid->Nc*sizeof(REAL),"FreeGrid");
  SunFree(grid->yv,grid->Nc*sizeof(REAL),"FreeGrid");
  SunFree(grid->dv,grid->Nc*sizeof(REAL),"FreeGrid");
  SunFree(grid->vwgt,grid->Nc*sizeof(REAL),"FreeGrid");

  SunFree(grid->edges,NUMEDGECOLUMNS*grid->Ne*sizeof(int),"FreeGrid");
  SunFree(grid->cells,NFACES*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->neigh,NFACES*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->eneigh,(2*NFACES-1)*grid->Ne*sizeof(int),"FreeGrid");
  SunFree(grid->part,grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->order,grid->Ne*sizeof(int),"FreeGrid");
  SunFree(grid->Nk,grid->Nc*sizeof(int),"FreeGrid");

  SunFree(grid->dz,grid->Nkmax*sizeof(REAL),"FreeGrid");
  SunFree(grid->face,NFACES*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->grad,2*grid->Ne*sizeof(int),"FreeGrid");
  SunFree(grid->mark,grid->Ne*sizeof(int),"FreeGrid");

  SunFree(grid->normal,NFACES*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->xadj,(grid->Nc+1)*sizeof(int),"FreeGrid");
  SunFree(grid->vtxdist,VTXDISTMAX,"FreeGrid");

  for(proc=0;proc<numprocs;proc++) 
    SunFree(grid->neighs[proc],numprocs*sizeof(int),"FreeGrid");
  SunFree(grid->neighs,numprocs*sizeof(int *),"FreeGrid");
  SunFree(grid->numneighs,numprocs*sizeof(int),"FreeGrid");

  SunFree(grid,sizeof(gridT),"FreeGrid");
}

static void VertGrid(gridT *maingrid, gridT **localgrid, MPI_Comm comm)
{
  int i, j, k, ne, myproc, numprocs, vertgridcorrect, stairstep;
  REAL dz0, dmin, dmax, dmaxtest;

  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myproc);

  maingrid->Nkmax = MPI_GetValue(DATAFILE,"Nkmax","VertGrid",myproc);
  vertgridcorrect = MPI_GetValue(DATAFILE,"vertgridcorrect","VertGrid",myproc);
  stairstep = MPI_GetValue(DATAFILE,"stairstep","VertGrid",myproc);

  if(myproc==0) {
    dmax = maingrid->dv[0];
    dmin = maingrid->dv[0];
    for(i=0;i<maingrid->Nc;i++) {
      if(maingrid->dv[i]>dmax)
	dmax=maingrid->dv[i];
      if(maingrid->dv[i]<dmin)
	dmin=maingrid->dv[i];
    }
    dz0 = dmax/(REAL)(maingrid->Nkmax);

    if(dmin < dz0 && vertgridcorrect) {
      dz0 = dmin;
      if(WARNING) {
	printf("Warning!\n");
	printf("Not enough vertical grid points to resolve the most shallow cell!\n");
	printf("Changing maximum number of vertical cells from %d to %d\n",
	       maingrid->Nkmax,(int)(dmax/dz0));
      }
      maingrid->Nkmax = 1+(int)(dmax/dz0);
    }
    maingrid->dz = (REAL *)SunMalloc(maingrid->Nkmax*sizeof(REAL),"VertGrid");

    GetDZ(maingrid->dz,dmax,dmax,maingrid->Nkmax,myproc);
    dmaxtest=0;
    for(k=0;k<maingrid->Nkmax;k++)
      dmaxtest+=maingrid->dz[k];
    for(i=0;i<maingrid->Nc;i++)
      if(fabs(maingrid->dv[i]-dmaxtest>SMALL) && WARNING) {
	printf("Warning...sum of grid spacings dz is less than depth! %e\n",maingrid->dv[i]-dmaxtest);
	break;
      }
  }
  MPI_Bcast(&(maingrid->Nkmax),1,MPI_INT,0,comm);
  if(myproc!=0)
    maingrid->dz = (REAL *)SunMalloc(maingrid->Nkmax*sizeof(REAL),"VertGrid");
  MPI_Bcast((void *)maingrid->dz,maingrid->Nkmax,MPI_DOUBLE,0,comm);
  MPI_Bcast(&dz0,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(&dmax,1,MPI_DOUBLE,0,comm);

  (*localgrid)->Nkmax = maingrid->Nkmax;
  (*localgrid)->dz = (REAL *)SunMalloc((*localgrid)->Nkmax*sizeof(REAL),"VertGrid");
  for(k=0;k<(*localgrid)->Nkmax;k++) 
    (*localgrid)->dz[k]=maingrid->dz[k];

  maingrid->Nk = (int *)SunMalloc(maingrid->Nc*sizeof(int),"VertGrid");
  (*localgrid)->dztop = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"VertGrid");
  (*localgrid)->dzbot = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"VertGrid");
  (*localgrid)->ctop = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),"VertGrid");
  (*localgrid)->ctopold = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),"VertGrid");
  (*localgrid)->etop = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),"VertGrid");
  (*localgrid)->etopold = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),"VertGrid");
  (*localgrid)->Nk = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),"VertGrid");
  (*localgrid)->Nke = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),"VertGrid");
  (*localgrid)->Nkc = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),"VertGrid");

  if((*localgrid)->Nkmax>1) {
    for(i=0;i<(*localgrid)->Nc;i++) 
      if((*localgrid)->dv[i]==dmax) 
	(*localgrid)->Nk[i] = (*localgrid)->Nkmax;
      else {
	dmaxtest=0;
	for(k=0;k<(*localgrid)->Nkmax;k++) {
	  dmaxtest+=(*localgrid)->dz[k];
	  (*localgrid)->Nk[i] = k+1;
	  if(dmaxtest>=(*localgrid)->dv[i]) 
	    break;
	}
      }
    for(i=0;i<maingrid->Nc;i++) 
      if(maingrid->dv[i]==dmax) 
	maingrid->Nk[i] = maingrid->Nkmax;
      else {
	dmaxtest=0;
	for(k=0;k<maingrid->Nkmax;k++) {
	  dmaxtest+=maingrid->dz[k];
	  maingrid->Nk[i] = k+1;
	  if(dmaxtest>=maingrid->dv[i]) 
	    break;
	}
      }
  } else {
    for(i=0;i<(*localgrid)->Nc;i++)
      (*localgrid)->Nk[i] = (*localgrid)->Nkmax;
    for(i=0;i<maingrid->Nc;i++)
      maingrid->Nk[i] = maingrid->Nkmax;
  }

  for(i=0;i<(*localgrid)->Nc;i++) {
    (*localgrid)->ctop[i] = 0;
    (*localgrid)->dztop[i] = (*localgrid)->dz[(*localgrid)->ctop[i]];
  }

  for(j=0;j<(*localgrid)->Ne;j++) {
    (*localgrid)->etop[j] = 0;
    (*localgrid)->etopold[j] = 0;
  }
  
  // To use stair-stepping do this (not partial-stepping).
  if(stairstep) 
    for(i=0;i<(*localgrid)->Nc;i++) {
      dmaxtest=0;
      for(k=0;k<(*localgrid)->Nk[i];k++)
	dmaxtest+=(*localgrid)->dz[k];
      (*localgrid)->dv[i]=dmaxtest;
    }

  for(j=0;j<(*localgrid)->Ne;j++) {
    ne = (*localgrid)->eptr[j];
    if(maingrid->grad[2*ne]==-1) {
      (*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
      (*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
    } else if(maingrid->grad[2*ne+1]==-1) {
      (*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne]];
      (*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne]];
    }
    else {
      if(maingrid->Nk[maingrid->grad[2*ne]]<
	 maingrid->Nk[maingrid->grad[2*ne+1]]) {
	(*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne]];
	(*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
      } else if(maingrid->Nk[maingrid->grad[2*ne]]>
		maingrid->Nk[maingrid->grad[2*ne+1]]) {
	(*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
	(*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne]];
      } else {
	(*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne]];
	(*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne]];
      }
    }
  }
}

/*
 * Function: ResortBoundaries
 * Usage: ResortBoundaries(grid,myproc);
 * -------------------------------------
 *
 */
static void ResortBoundaries(gridT *localgrid, int myproc)
{
  int neigh, n, *tmp;

  for(neigh=0;neigh<localgrid->Nneighs;neigh++) {
    tmp = (int *)SunMalloc(localgrid->num_cells_send[neigh]*sizeof(int),"ResortBoundaries");
    for(n=0;n<localgrid->num_cells_send[neigh];n++)
      tmp[n]=localgrid->mnptr[localgrid->cell_send[neigh][n]];
    Sort(localgrid->cell_send[neigh],tmp,localgrid->num_cells_send[neigh]);
    free(tmp);

    tmp = (int *)SunMalloc(localgrid->num_cells_recv[neigh]*sizeof(int),"ResortBoundaries");
    for(n=0;n<localgrid->num_cells_recv[neigh];n++)
      tmp[n]=localgrid->mnptr[localgrid->cell_recv[neigh][n]];
    Sort(localgrid->cell_recv[neigh],tmp,localgrid->num_cells_recv[neigh]);
    free(tmp);

    tmp = (int *)SunMalloc(localgrid->num_edges_send[neigh]*sizeof(int),"ResortBoundaries");
    for(n=0;n<localgrid->num_edges_send[neigh];n++)
      tmp[n]=localgrid->eptr[localgrid->edge_send[neigh][n]];
    Sort(localgrid->edge_send[neigh],tmp,localgrid->num_edges_send[neigh]);
    free(tmp);

    tmp = (int *)SunMalloc(localgrid->num_edges_recv[neigh]*sizeof(int),"ResortBoundaries");
    for(n=0;n<localgrid->num_edges_recv[neigh];n++)
      tmp[n]=localgrid->eptr[localgrid->edge_recv[neigh][n]];
    Sort(localgrid->edge_recv[neigh],tmp,localgrid->num_edges_recv[neigh]);
    free(tmp);
  }
}

static void MakePointers(gridT *maingrid, gridT **localgrid, int myproc, MPI_Comm comm)
{
  int i, n, nf, ne, neigh, neighproc, j, k, mark;
  int **cell_send, **cell_recv, **edge_send, **edge_recv;
  int *num_cells_send, *num_cells_recv,
    *num_edges_send, *num_edges_recv;
  int *cellp, *edgep, *celldist, *edgedist, *lcptr, *leptr;
  unsigned short *flagged;
  MPI_Status status;

  cellp = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),"MakePointers");
  edgep = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),"MakePointers");
  celldist = (int *)SunMalloc((MAXBCTYPES-1)*sizeof(int),"MakePointers");
  edgedist = (int *)SunMalloc((MAXMARKS-1)*sizeof(int),"MakePointers");
  lcptr = (int *)SunMalloc(maingrid->Nc*sizeof(int),"MakePointers");
  leptr = (int *)SunMalloc(maingrid->Ne*sizeof(int),"MakePointers");
  flagged = (unsigned short *)SunMalloc((*localgrid)->Ne*sizeof(int),"MakePointers");

  // Set up the global pointer arrays.
  for(i=0;i<(*localgrid)->Nc;i++)
    lcptr[(*localgrid)->mnptr[i]]=i;
  for(i=0;i<(*localgrid)->Ne;i++)
    leptr[(*localgrid)->eptr[i]]=i;

  k=0;
  celldist[0]=0;
  // Put computational cells first
  for(n=0;n<(*localgrid)->Nc;n++)
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==0 ||
       IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2)
	cellp[k++]=n;
  celldist[1]=k;
  // Followed by boundary cells
  for(n=0;n<(*localgrid)->Nc;n++)
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==1)
      cellp[k++]=n;
  celldist[2]=k;

  k=0;
  edgedist[0]=0;
  // Put computational edges first
  for(n=0;n<(*localgrid)->Ne;n++)
    if((*localgrid)->mark[n]==0 || (*localgrid)->mark[n]==5)
      edgep[k++]=n;
  edgedist[1]=k;
  // Followed by boundary cells
  for(mark=1;mark<MAXMARKS-2;mark++) {
    for(n=0;n<(*localgrid)->Ne;n++)
      if((*localgrid)->mark[n]==mark)
	edgep[k++]=n;
    edgedist[mark+1]=k;
  }

  (*localgrid)->cellp=cellp;
  (*localgrid)->edgep=edgep;
  (*localgrid)->celldist=celldist;
  (*localgrid)->edgedist=edgedist;

  /*
  for(bctype=0;bctype<MAXBCTYPES-2;bctype++) {
    printf("Cells of type %d: ",bctype);
    for(j=celldist[bctype];j<celldist[bctype+1];j++)
      printf("%d ",cellp[j]);
    printf("\n");
  }
  for(mark=0;mark<MAXMARKS-2;mark++) {
    printf("Edges of type %d: ",mark);
    for(j=edgedist[mark];j<edgedist[mark+1];j++)
      printf("%d ",edgep[j]);
    printf("\n");
  }
  */

  cell_send = (int **)SunMalloc((*localgrid)->Nneighs*sizeof(int *),"MakePointers");
  cell_recv = (int **)SunMalloc((*localgrid)->Nneighs*sizeof(int *),"MakePointers");
  edge_send = (int **)SunMalloc((*localgrid)->Nneighs*sizeof(int *),"MakePointers");
  edge_recv = (int **)SunMalloc((*localgrid)->Nneighs*sizeof(int *),"MakePointers");
  num_cells_send = (int *)SunMalloc((*localgrid)->Nneighs*sizeof(int),"MakePointers");
  num_cells_recv = (int *)SunMalloc((*localgrid)->Nneighs*sizeof(int),"MakePointers");
  num_edges_send = (int *)SunMalloc((*localgrid)->Nneighs*sizeof(int),"MakePointers");
  num_edges_recv = (int *)SunMalloc((*localgrid)->Nneighs*sizeof(int),"MakePointers");

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    num_cells_send[neigh]=0;
    num_cells_recv[neigh]=0;
    num_edges_send[neigh]=0;
    num_edges_recv[neigh]=0;
  }

  // Set up the pointers for receiving first and place the global
  // indices in the recv pointer array.  Place the cells that border
  // the local processor first followed by the others.  This will be
  // useful when transferring data for free-surface and nonhydrostatic
  // pressure iterations, since these are the only ones that are needed.
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    for(n=0;n<(*localgrid)->Nc;n++) 
      if(maingrid->part[(*localgrid)->mnptr[n]]==neighproc)
	num_cells_recv[neigh]++;
    cell_recv[neigh]=(int *)SunMalloc(num_cells_recv[neigh]*sizeof(int),"MakePointers");

    k=0;
    for(n=0;n<(*localgrid)->Nc;n++) 
      if(maingrid->part[(*localgrid)->mnptr[n]]==neighproc) {
	for(nf=0;nf<NFACES;nf++) {
	  ne=(*localgrid)->neigh[n*NFACES+nf];
	  if(ne!=-1)
	    if(maingrid->part[(*localgrid)->mnptr[ne]]==myproc) {
	      cell_recv[neigh][k++]=(*localgrid)->mnptr[n];
	      break;
	    }
	}
      }
  }
  if(VERBOSE>2) 
    for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++)
      printf("Proc: %d, neighbor %d, receiving %d\n",
	     myproc,(*localgrid)->myneighs[neigh],
	     num_cells_recv[neigh]);

  // Send out the number of cells that are being received and then
  // the actual global indices being sent.
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Send(&(num_cells_recv[neigh]),1,MPI_INT,neighproc,1,comm); 
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Recv(&(num_cells_send[neigh]),1,MPI_INT,neighproc,1,comm,&status);
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Send(cell_recv[neigh],num_cells_recv[neigh],MPI_INT,neighproc,1,comm); 
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    cell_send[neigh]=(int *)SunMalloc(num_cells_send[neigh]*sizeof(int),"MakePointers");
    MPI_Recv(cell_send[neigh],num_cells_send[neigh],MPI_INT,neighproc,1,comm,&status);
  }

  // Now set the indices to point to the local grid.
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    for(j=0;j<num_cells_send[neigh];j++) 
      cell_send[neigh][j]=lcptr[cell_send[neigh][j]];
    for(j=0;j<num_cells_recv[neigh];j++)
      cell_recv[neigh][j]=lcptr[cell_recv[neigh][j]];
  }

  // Now do the edges.  The edges that correspond to the cells that
  // are sent are added to the edge-send/recv arrays.
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    for(j=0;j<(*localgrid)->Ne;j++)
      flagged[j]=0;
    for(i=0;i<num_cells_send[neigh];i++)
      for(nf=0;nf<NFACES;nf++)
	if(!flagged[(*localgrid)->face[NFACES*cell_send[neigh][i]+nf]]) {
	  flagged[(*localgrid)->face[NFACES*cell_send[neigh][i]+nf]]=1;
	  num_edges_send[neigh]++;
	}
    edge_send[neigh]=(int *)SunMalloc(num_edges_send[neigh]*sizeof(int),"MakePointers");
    for(j=0;j<(*localgrid)->Ne;j++)
      flagged[j]=0;
    k=0;
    for(i=0;i<num_cells_send[neigh];i++)
      for(nf=0;nf<NFACES;nf++)
	if(!flagged[(*localgrid)->face[NFACES*cell_send[neigh][i]+nf]]) {
	  flagged[(*localgrid)->face[NFACES*cell_send[neigh][i]+nf]]=1;
	  edge_send[neigh][k++]=
	    (*localgrid)->eptr[(*localgrid)->face[NFACES*cell_send[neigh][i]+nf]];
	}
    MPI_Send(&(num_edges_send[neigh]),1,MPI_INT,neighproc,1,comm); 
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Recv(&(num_edges_recv[neigh]),1,MPI_INT,neighproc,1,comm,&status);
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Send(edge_send[neigh],num_edges_send[neigh],MPI_INT,neighproc,1,comm); 
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    edge_recv[neigh] = (int *)SunMalloc(num_edges_recv[neigh]*sizeof(int),"MakePointers");
    MPI_Recv(edge_recv[neigh],num_edges_recv[neigh],MPI_INT,neighproc,1,comm,&status);
  }

  // Now convert back to local indices
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    for(j=0;j<num_edges_send[neigh];j++)
      edge_send[neigh][j]=leptr[edge_send[neigh][j]];
    for(j=0;j<num_edges_recv[neigh];j++) 
      edge_recv[neigh][j]=leptr[edge_recv[neigh][j]];
  }

  if(VERBOSE>3) {
    printf("CELLS: \n");
    for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
      neighproc=(*localgrid)->myneighs[neigh];
      printf("Proc (%d<-->%d): S/R= %d/%d\n",
	     myproc,neighproc,num_cells_send[neigh],num_cells_recv[neigh]);
    }

    for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
      printf("Send %d-->%d: ",myproc,(*localgrid)->myneighs[neigh]);
      for(j=0;j<num_cells_send[neigh];j++)
	printf("%d ",cell_send[neigh][j]);
      printf("\n");
      printf("Recv %d<--%d: ",myproc,(*localgrid)->myneighs[neigh]);
      for(j=0;j<num_cells_recv[neigh];j++)
	printf("%d ",cell_recv[neigh][j]);
      printf("\n");
    }    
    
    printf("EDGES: \n");
    for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
      neighproc=(*localgrid)->myneighs[neigh];
      printf("Proc (%d<-->%d): S/R= %d/%d\n",
	     myproc,neighproc,num_edges_send[neigh],num_edges_recv[neigh]);
    }

    for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
      printf("Send %d-->%d: ",myproc,(*localgrid)->myneighs[neigh]);
      for(j=0;j<num_edges_send[neigh];j++)
	printf("%d ",edge_send[neigh][j]);
      printf("\n");
      printf("Recv %d<--%d: ",myproc,(*localgrid)->myneighs[neigh]);
      for(j=0;j<num_edges_recv[neigh];j++)
	printf("%d ",edge_recv[neigh][j]);
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  (*localgrid)->cell_send=cell_send;
  (*localgrid)->cell_recv=cell_recv;
  (*localgrid)->edge_send=edge_send;
  (*localgrid)->edge_recv=edge_recv;
  (*localgrid)->num_cells_send=num_cells_send;
  (*localgrid)->num_cells_recv=num_cells_recv;
  (*localgrid)->num_edges_send=num_edges_send;
  (*localgrid)->num_edges_recv=num_edges_recv;

  free(lcptr);
  free(leptr);
}

static void ReOrder(gridT *grid) 
{
  int n, nf, numflag, options[8], Nc = grid->Nc, Ne=grid->Ne;
  int *corder, *corderp, *eorder, *eorderp;
  REAL *tmp;

  corder = (int *)SunMalloc(Nc*sizeof(int),"ReOrder");
  corderp = (int *)SunMalloc(Nc*sizeof(int),"ReOrder");
  eorder = (int *)SunMalloc(Ne*sizeof(int),"ReOrder");
  eorderp = (int *)SunMalloc(Ne*sizeof(int),"ReOrder");
  tmp = (REAL *)SunMalloc(2*(NFACES-1)*Ne*sizeof(REAL),"ReOrder");

  for(n=0;n<Nc;n++) {
    corder[n]=n;
    corderp[n]=n;
  }
  for(n=0;n<Ne;n++) {
    eorder[n]=n;
    eorderp[n]=n;
  }

  numflag=0;
  options[0]=0;

  /*
  grid->xadj = (int *)SunMalloc((Nc+1)*sizeof(int),"ReOrder");
  CreateCellGraph(grid);
  METIS_NodeND(&Nc,grid->xadj,grid->adjncy,&numflag,options,corder,corderp);
  free(grid->xadj);
  free(grid->adjncy);
  */

  /*
  grid->xadj = (int *)SunMalloc((Ne+1)*sizeof(int),"ReOrder");
  CreateEdgeGraph(grid);
  METIS_NodeND(&Ne,grid->xadj,grid->adjncy,&numflag,options,eorder,eorderp);
  free(grid->xadj);
  free(grid->adjncy);
  */

  // Reorder the data corresponding to cells with the corder array
  ReOrderRealArray(grid->Ac,corder,tmp,grid->Nc,1);
  ReOrderRealArray(grid->xv,corder,tmp,grid->Nc,1);
  ReOrderRealArray(grid->yv,corder,tmp,grid->Nc,1);
  ReOrderRealArray(grid->dv,corder,tmp,grid->Nc,1);
  ReOrderIntArray(grid->cells,corder,(int *)tmp,grid->Nc,NFACES);
  ReOrderIntArray(grid->face,corder,(int *)tmp,grid->Nc,NFACES);
  ReOrderIntArray(grid->normal,corder,(int *)tmp,grid->Nc,NFACES);
  ReOrderIntArray(grid->neigh,corder,(int *)tmp,grid->Nc,NFACES);
  ReOrderIntArray(grid->vwgt,corder,(int *)tmp,grid->Nc,1);
  ReOrderIntArray(grid->mnptr,corder,(int *)tmp,grid->Nc,1);
  ReOrderIntArray(grid->Nk,corder,(int *)tmp,grid->Nc,1);

  //Reorder the data corresponding to edges with the eorder array
  ReOrderRealArray(grid->df,eorder,tmp,grid->Ne,1);
  ReOrderRealArray(grid->dg,eorder,tmp,grid->Ne,1);
  ReOrderRealArray(grid->n1,eorder,tmp,grid->Ne,1);
  ReOrderRealArray(grid->n2,eorder,tmp,grid->Ne,1);
  ReOrderRealArray(grid->xi,eorder,tmp,grid->Ne,2*(NFACES-1));
  ReOrderIntArray(grid->grad,eorder,(int *)tmp,grid->Ne,2);
  ReOrderIntArray(grid->eneigh,eorder,(int *)tmp,grid->Ne,2*(NFACES-1));
  ReOrderIntArray(grid->edges,eorder,(int *)tmp,grid->Ne,NUMEDGECOLUMNS);
  ReOrderIntArray(grid->mark,eorder,(int *)tmp,grid->Ne,1);
  ReOrderIntArray(grid->eptr,eorder,(int *)tmp,grid->Ne,1);
  ReOrderIntArray(grid->Nke,eorder,(int *)tmp,grid->Ne,1);
  ReOrderIntArray(grid->Nkc,eorder,(int *)tmp,grid->Ne,1);

  // Now adjust the pointers to point to the new locations
  // Face and eneigh arrays point to edges.  We need to use the createfacearray
  // again because it guarantees consistency for the CG solver.  This is
  // because reordering does not necessarily guarantee that the ordering of
  // the faces in the face array are in the same as the corresponding neighbors
  // in the neigh array.
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<NFACES;nf++) 
      tmp[n*NFACES+nf]=grid->face[n*NFACES+nf];
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<NFACES;nf++) 
      grid->face[n*NFACES+nf]=eorderp[(int)(tmp[n*NFACES+nf])];

  for(n=0;n<Ne;n++)
    for(nf=0;nf<2*(NFACES-1);nf++) 
      tmp[2*(NFACES-1)*n+nf]=grid->eneigh[2*(NFACES-1)*n+nf];
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2*(NFACES-1);nf++)
      if((int)(tmp[n*2*(NFACES-1)+nf])!=-1)
	grid->eneigh[n*2*(NFACES-1)+nf]=eorderp[(int)(tmp[n*2*(NFACES-1)+nf])];
      else
	grid->eneigh[n*2*(NFACES-1)+nf]=-1;

  // Grad and neigh arrays point to cells
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2;nf++)
      tmp[2*n+nf]=grid->grad[2*n+nf];
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2;nf++)
      if((int)(tmp[2*n+nf])!=-1)
	grid->grad[2*n+nf]=corderp[(int)(tmp[2*n+nf])];
      else
	grid->grad[2*n+nf]=-1;

  for(n=0;n<Nc;n++)
    for(nf=0;nf<NFACES;nf++)
      tmp[n*NFACES+nf]=grid->neigh[n*NFACES+nf];
  for(n=0;n<Nc;n++)
    for(nf=0;nf<NFACES;nf++)
      if((int)(tmp[n*NFACES+nf])!=-1)
	grid->neigh[n*NFACES+nf]=corderp[(int)(tmp[n*NFACES+nf])];
      else
  	grid->neigh[n*NFACES+nf]=-1;

  free(corder);
  free(corderp);
  free(eorder);
  free(eorderp);
  free(tmp);
}

static void EdgeMarkers(gridT *maingrid, gridT **localgrid, int myproc)
{
  int n, ne, nf;

  for(n=0;n<(*localgrid)->Nc;n++) 
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==3) 
      for(nf=0;nf<NFACES;nf++) {
	ne = (*localgrid)->face[n*NFACES+nf];
	if(!(*localgrid)->mark[ne]) 
	  (*localgrid)->mark[ne]=6;
      }
  for(n=0;n<(*localgrid)->Nc;n++) 
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2) 
      for(nf=0;nf<NFACES;nf++) {
	ne = (*localgrid)->face[n*NFACES+nf];
	if((*localgrid)->mark[ne]==6) 
	  (*localgrid)->mark[ne]=5;
      }
}

/*
 * Function: GetNumCells
 * Usage: GetNumCells(grid,myproc);
 * --------------------------------
 * Returns the number of cells (including ghost points)
 * belonging to a particular processor that result from
 * the partitioning.
 *
 */
int GetNumCells(gridT *grid, int proc)
{
  int j, nf, flag, ncells, ne;

  ncells=0;
  for(j=0;j<grid->Nc;j++) {
    flag = 0;
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->neigh[j*NFACES+nf];
      if(ne!=-1)
	if(grid->part[ne]==proc) {
	flag=1;
      }
    }
    if(grid->part[j]==proc || flag==1) {
      ncells++;
    }
  }
  return ncells;
}

/*
 * Function: GetNumEdges
 * Usage: GetNumEdges(grid);
 * -------------------------
 * Returns the number of edges (including ghost edges)
 * belonging to a particular processor that result from
 * the partitioning.
 *
 */
int GetNumEdges(gridT *grid)
{
  int j, nf, nedges, ne;

  unsigned short *flagged = (unsigned short *)SunMalloc(grid->Nc*sizeof(unsigned short),"GetNumEdges");

  for(j=0;j<grid->Nc;j++) 
    flagged[j]=0;

  nedges=0;
  for(j=0;j<grid->Nc;j++) {
    flagged[j]=1;
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->neigh[j*NFACES+nf];
      if(ne == -1)
	nedges++;
      else if(!flagged[ne])
	nedges++;
    }
  }
  free(flagged);

  return nedges;
}

int IsCellNeighborProc(int nc, gridT *maingrid, gridT *localgrid, 
		       int myproc, int neighproc)
{
  int nf, ne;

  if(maingrid->part[localgrid->mnptr[nc]]!=myproc &&
     maingrid->part[localgrid->mnptr[nc]]!=neighproc)
    return 0;

  if(maingrid->part[localgrid->mnptr[nc]]==neighproc)
    return 1;

  for(nf=0;nf<NFACES;nf++) {
    ne = localgrid->neigh[nc*NFACES+nf];
    if(ne != -1)
      if(maingrid->part[localgrid->mnptr[ne]] == neighproc)
	return 1;
  }
  return 0;
}

int IsEdgeNeighborProc(int ne, gridT *maingrid, gridT *localgrid, 
		       int myproc, int neighproc)
{
  int j, cnp, nc;

  for(j=0;j<2;j++) {
    nc = localgrid->grad[2*ne+j];
    if(nc != -1) {
      cnp = IsCellNeighborProc(nc,maingrid,localgrid,myproc,neighproc);
      if(cnp)
	return 1;
    }
  }
  return 0;
}

/*
 * Function: Topology
 * Usage: Topology(maingrid,localgrid,myproc,numprocs);
 * ----------------------------------------------------
 * Compute the number of neighbors for each processor and
 * determine the processor ids of each neighbor.  
 *
 */
void Topology(gridT **maingrid, gridT **localgrid, int myproc, int numprocs)
{
  int nf, j, proc, neigh1, neigh2, loc;

  (*maingrid)->numneighs=(int *)SunMalloc(numprocs*sizeof(int),"Topology");
  (*maingrid)->neighs=(int **)SunMalloc(numprocs*sizeof(int *),"Topology");

  // Initialize the arrays
  for(proc=0;proc<numprocs;proc++) {
    (*maingrid)->numneighs[proc]=0;
    (*maingrid)->neighs[proc]=(int *)SunMalloc(numprocs*sizeof(int),"Topology");
    for(j=0;j<numprocs;j++) 
      (*maingrid)->neighs[proc][j]=-1;
  }

  // Count the number of neighbors each processor grid has and
  // also count the number of interproc edges between each
  // neighbor pair.  
  for(j=0;j<(*maingrid)->Nc;j++)
    for(nf=0;nf<NFACES;nf++)
      if((*maingrid)->neigh[NFACES*j+nf] != -1)
	if((*maingrid)->part[j] != (*maingrid)->part[(*maingrid)->neigh[NFACES*j+nf]]) {
	  neigh1=(*maingrid)->part[j];
	  neigh2=(*maingrid)->part[(*maingrid)->neigh[NFACES*j+nf]];
	  loc = IsMember(neigh2,(*maingrid)->neighs[neigh1],numprocs);
	  if(loc<0) {
	    (*maingrid)->neighs[neigh1][(*maingrid)->numneighs[neigh1]]=neigh2;
	    (*maingrid)->numneighs[neigh1]++;
	  }
	}

  (*localgrid)->Nneighs=(*maingrid)->numneighs[myproc];
  (*localgrid)->myneighs=(int *)SunMalloc((*localgrid)->Nneighs*sizeof(int),"Topology");

  //  printf("Proc %d, neighs: ",myproc);
  for(j=0;j<(*localgrid)->Nneighs;j++) {
    (*localgrid)->myneighs[j]=(*maingrid)->neighs[myproc][j];
    //    printf("%d ",(*localgrid)->myneighs[j]);
  }
  //  printf("\n");
}

/*
 * Function: TransferData
 * Usage: TransferData(maingrid,localgrid,myproc);
 * -----------------------------------------------
 *
 * This function transfer cell data from the main grid onto the local 
 * grid and creates the global pointer that allows the identification of
 * the global cell index given a local index.  Because the main grid 
 * is stored on every processor, no communcation is needed. The computational
 * cells are printed first, followed by the boundary cells.
 *
 * The local grid will have new ghost cells due to the cut edges.
 * There is one ghost cell for each cut edge.  The local cells
 * array points to indices in the main xp,yp arrays.  But the local
 * neigh pointer points to local cell indices.
 * lcptr points from a cell in the main grid to one in the local grid.
 * mnptr points from a cell in the local grid to one in the main grid
 * and is stored on every local grid.
 *
 */
static void TransferData(gridT *maingrid, gridT **localgrid, int myproc)
{
  int i, j, k, n, nc, nf, ne, ng, flag, mgptr, *lcptr, *leptr, bctype, iface;
  unsigned short *flagged = 
    (unsigned short *)SunMalloc(maingrid->Ne*sizeof(unsigned short),"TransferData");

  lcptr = (int *)SunMalloc(maingrid->Nc*sizeof(int),"TransferData");
  leptr = (int *)SunMalloc(maingrid->Ne*sizeof(int),"TransferData");

  (*localgrid)->Nc = GetNumCells(maingrid,myproc);
  (*localgrid)->mnptr = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),"TransferData");
  (*localgrid)->vwgt = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),"TransferData");
  (*localgrid)->cells = (int *)SunMalloc(NFACES*(*localgrid)->Nc*sizeof(int),"TransferData");
  (*localgrid)->xv = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"TransferData");
  (*localgrid)->yv = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"TransferData");
  (*localgrid)->dv = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"TransferData");
  (*localgrid)->neigh = (int *)SunMalloc(NFACES*(*localgrid)->Nc*sizeof(int),"TransferData");
  (*localgrid)->normal = (int *)SunMalloc(NFACES*(*localgrid)->Nc*sizeof(int),"TransferData");

  for(j=0;j<maingrid->Nc;j++) 
    lcptr[j]=-1;
  for(j=0;j<maingrid->Ne;j++) 
    leptr[j]=-1;

  k=0;
  for(bctype=0;bctype<MAXBCTYPES;bctype++)
    for(j=0;j<maingrid->Nc;j++) 
      if(IsBoundaryCell(j,maingrid,myproc)==bctype) {
	flag=0;
	for(nf=0;nf<NFACES;nf++)
	  if(maingrid->neigh[j*NFACES+nf]!=-1)
	    if(maingrid->part[maingrid->neigh[j*NFACES+nf]]==myproc) {
	      flag=1;
	      break;
	    }
	if(flag==1 || maingrid->part[j]==myproc) {
	  lcptr[j]=k;
	  (*localgrid)->mnptr[k]=j;
	  (*localgrid)->xv[k]=maingrid->xv[j];
	  (*localgrid)->yv[k]=maingrid->yv[j];
	  (*localgrid)->dv[k]=maingrid->dv[j];
	  (*localgrid)->vwgt[k]=maingrid->vwgt[j];
	  for(nf=0;nf<NFACES;nf++)
	    (*localgrid)->cells[k*NFACES+nf]=maingrid->cells[j*NFACES+nf];
	  k++;
	}
      }
  /*
  for(i=0;i<(*localgrid)->Nc;i++) {
    j = (*localgrid)->mnptr[i];
    lcptr[j]=i;
    (*localgrid)->xv[i]=maingrid->xv[j];
    (*localgrid)->yv[i]=maingrid->yv[j];
    (*localgrid)->dv[i]=maingrid->dv[j];
    (*localgrid)->vwgt[i]=maingrid->vwgt[j];
    for(nf=0;nf<NFACES;nf++)
      (*localgrid)->cells[i*NFACES+nf]=maingrid->cells[j*NFACES+nf];
  }
  */
  for(j=0;j<(*localgrid)->Nc;j++) 
    for(nf=0;nf<NFACES;nf++) {
      mgptr = maingrid->neigh[(*localgrid)->mnptr[j]*NFACES+nf];
      if(mgptr>=0 && IsMember(mgptr,(*localgrid)->mnptr,(*localgrid)->Nc)>=0)
	(*localgrid)->neigh[j*NFACES+nf]=lcptr[mgptr];
      else
	(*localgrid)->neigh[j*NFACES+nf]=-1;
    }

  (*localgrid)->Ne = GetNumEdges(*localgrid);
  (*localgrid)->edges = (int *)SunMalloc(NUMEDGECOLUMNS*(*localgrid)->Ne*sizeof(int),"TransferData");
  (*localgrid)->mark = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),"TransferData");
  (*localgrid)->eptr = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),"TransferData");
  (*localgrid)->eneigh = (int *)SunMalloc(2*(NFACES-1)*(*localgrid)->Ne*sizeof(int),"TransferData");
  (*localgrid)->face = (int *)SunMalloc(NFACES*(*localgrid)->Nc*sizeof(int),"TransferData");
  (*localgrid)->grad = (int *)SunMalloc(2*(*localgrid)->Ne*sizeof(int),"TransferData");
  (*localgrid)->gradf = (int *)SunMalloc(2*(*localgrid)->Ne*sizeof(int),"TransferData");

  for(j=0;j<maingrid->Ne;j++) 
    flagged[j]=0;
  
  k=0;
  for(i=0;i<(*localgrid)->Nc;i++) 
    for(nf=0;nf<NFACES;nf++) {
      iface = maingrid->face[(*localgrid)->mnptr[i]*NFACES+nf];
      if(!flagged[iface]) {
	flagged[iface]=1;
	for(j=0;j<NUMEDGECOLUMNS;j++)
	  (*localgrid)->edges[k*NUMEDGECOLUMNS+j]=maingrid->edges[iface*NUMEDGECOLUMNS+j];
	(*localgrid)->mark[k]=maingrid->mark[iface];
	leptr[iface]=k;
	(*localgrid)->eptr[k++]=iface;
      }
    }

/*
  for(i=0;i<(*localgrid)->Nc;i++) 
    for(nf=0;nf<NFACES;nf++) 
      (*localgrid)->face[i*NFACES+nf]=leptr[maingrid->face[(*localgrid)->mnptr[i]*NFACES+nf]];
*/
  for(i=0;i<(*localgrid)->Nc;i++) 
    for(n=0;n<(*localgrid)->Ne;n++)
      for(j=0;j<2;j++) {
	(*localgrid)->grad[2*n+j]=-1;
	nc = maingrid->grad[2*(*localgrid)->eptr[n]+j];
	if(nc != -1)
	  (*localgrid)->grad[2*n+j]=lcptr[nc];
      }

  CreateFaceArray((*localgrid)->grad,(*localgrid)->gradf,(*localgrid)->neigh,(*localgrid)->face,
		  (*localgrid)->Nc,(*localgrid)->Ne);
  CreateNormalArray((*localgrid)->grad,(*localgrid)->face,(*localgrid)->normal,(*localgrid)->Nc);

  for(n=0;n<(*localgrid)->Ne;n++) {
    for(ne=0;ne<2*(NFACES-1);ne++)
      (*localgrid)->eneigh[2*(NFACES-1)*n+ne]=-1;
    ne = 0;
    for(ng=0;ng<2;ng++) 
      if((*localgrid)->grad[2*n+ng]!=-1)
	for(nf=0;nf<NFACES;nf++) {
	  if((*localgrid)->face[NFACES*(*localgrid)->grad[2*n+ng]+nf] != n) 
	    (*localgrid)->eneigh[2*(NFACES-1)*n+ne++]=
	      (*localgrid)->face[NFACES*(*localgrid)->grad[2*n+ng]+nf];
	}
  }

  free(flagged);
  free(lcptr);
  free(leptr);
}

static void Geometry(gridT *maingrid, gridT **grid, int myproc)
{
  int n, nf, ne, k, j, Nc=(*grid)->Nc, Ne=(*grid)->Ne, p1, p2;
  REAL xt[NFACES], yt[NFACES], xc, yc, den, R0, tx, ty, tmag, xdott;
  
  (*grid)->Ac = (REAL *)SunMalloc(Nc*sizeof(REAL),"Geometry");
  (*grid)->df = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->dg = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->def = (REAL *)SunMalloc(NFACES*Nc*sizeof(REAL),"Geometry");
  (*grid)->n1 = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->n2 = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->xe = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->ye = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->xi = (REAL *)SunMalloc(2*(NFACES-1)*Ne*sizeof(REAL),"Geometry");

  /* Compute the area of each cell.*/
  if(myproc==0 && VERBOSE>2) printf("Computing areas...\n");
  for(n=0;n<Nc;n++) {
    for(nf=0;nf<NFACES;nf++) {
      xt[nf]=maingrid->xp[(*grid)->cells[n*NFACES+nf]];
      yt[nf]=maingrid->yp[(*grid)->cells[n*NFACES+nf]];
    }
    (*grid)->Ac[n] = GetArea(xt,yt,NFACES);
  }
  
  /* Compute the length of each edge */
  if(myproc==0 && VERBOSE>2) printf("Computing lengths...\n");
  for(n=0;n<Ne;n++) 
    (*grid)->df[n]=
      sqrt(pow(maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n]]-
	       maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n+1]],2)+
	   pow(maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n]]-
	       maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n+1]],2));

  // Compute the centers of each edge, which are defined by the intersection
  // of the Voronoi Edge with the Delaunay Edge.  Note that this point is
  // not the midpoint of the edge when one of the neighboring triangles is obtuse 
  // and it has been corrected.
  for(n=0;n<Ne;n++) {
    p1 = (*grid)->edges[NUMEDGECOLUMNS*n];
    p2 = (*grid)->edges[NUMEDGECOLUMNS*n+1];

    tx = maingrid->xp[p2]-maingrid->xp[p1];
    ty = maingrid->yp[p2]-maingrid->yp[p1];
    tmag = sqrt(tx*tx+ty*ty);
    tx = tx/tmag;
    ty = ty/tmag;
    xdott = ((*grid)->xv[(*grid)->grad[2*n]]-maingrid->xp[p1])*tx+
      ((*grid)->yv[(*grid)->grad[2*n]]-maingrid->yp[p1])*ty;
    (*grid)->xe[n] = maingrid->xp[p1]+xdott*tx;
    (*grid)->ye[n] = maingrid->yp[p1]+xdott*ty;

    // This is the midpoint of the edge and is not used.
    /*
    (*grid)->xe[n] = 0.5*(maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n]]+
			  maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n+1]]);
    (*grid)->ye[n] = 0.5*(maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n]]+
			  maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n+1]]);
    */
  }

  /* Compute the normal distances between Voronoi points to compute the 
     gradients and then output the lengths to a data file. Also, compute 
     n1 and n2 which make up the normal vector components. */
  if(myproc==0 && VERBOSE>2) printf("Computing n1, n2, and dg..\n");
  for(n=0;n<Ne;n++) {
    if((*grid)->grad[2*n]!=-1 && (*grid)->grad[2*n+1]!=-1) {
      (*grid)->n1[n] = (*grid)->xv[(*grid)->grad[2*n]]-(*grid)->xv[(*grid)->grad[2*n+1]];
      (*grid)->n2[n] = (*grid)->yv[(*grid)->grad[2*n]]-(*grid)->yv[(*grid)->grad[2*n+1]];
      (*grid)->dg[n] = sqrt(pow((*grid)->n1[n],2)+pow((*grid)->n2[n],2));
      (*grid)->n1[n] = (*grid)->n1[n]/(*grid)->dg[n];
      (*grid)->n2[n] = (*grid)->n2[n]/(*grid)->dg[n];
      if(((*grid)->xv[(*grid)->grad[2*n]]==
	  (*grid)->xv[(*grid)->grad[2*n+1]]) &&
         ((*grid)->yv[(*grid)->grad[2*n]]==
	  (*grid)->yv[(*grid)->grad[2*n+1]])) {
	printf("Coincident Voronoi points on edge %d (%d,%d)!\n",n,(*grid)->grad[2*n],(*grid)->grad[2*n+1]);
      }
    }
    else {
      xc = (*grid)->xe[n];
      yc = (*grid)->ye[n];      
      if((*grid)->grad[2*n]==-1) {
        (*grid)->n1[n] = xc-(*grid)->xv[(*grid)->grad[2*n+1]];
        (*grid)->n2[n] = yc-(*grid)->yv[(*grid)->grad[2*n+1]];
	(*grid)->dg[n] = sqrt(pow((*grid)->n1[n],2)+
			      pow((*grid)->n2[n],2));
	(*grid)->n1[n] = (*grid)->n1[n]/(*grid)->dg[n];
	(*grid)->n2[n] = (*grid)->n2[n]/(*grid)->dg[n];
	(*grid)->dg[n] = 2.0*(*grid)->dg[n];
      } else {
        (*grid)->n1[n] = (*grid)->xv[(*grid)->grad[2*n]]-xc;
        (*grid)->n2[n] = (*grid)->yv[(*grid)->grad[2*n]]-yc;
	(*grid)->dg[n] = sqrt(pow((*grid)->n1[n],2)+
			      pow((*grid)->n2[n],2));
	(*grid)->n1[n] = (*grid)->n1[n]/(*grid)->dg[n];
	(*grid)->n2[n] = (*grid)->n2[n]/(*grid)->dg[n];
	(*grid)->dg[n] = 2.0*(*grid)->dg[n];
      }
    }
  }

  // Compute the length of the perpendicular bisector.  Note that this is not necessarily
  // the distance to the edge midpoint unless both triangles have not been corrected.  
  // Note that def points outward, so if def<0 then uc and vc will be 
  // computed correctly for obtuse triangles and advection of momentum is conservative.  
  // However, the method is unstable for obtuse triangles and correction must be used.  
  for(n=0;n<Nc;n++) {
    for(nf=0;nf<NFACES;nf++) {
      ne = (*grid)->face[n*NFACES+nf];
      (*grid)->def[n*NFACES+nf] = 
      	-(((*grid)->xv[n]-maingrid->xp[(*grid)->edges[ne*NUMEDGECOLUMNS]])*(*grid)->n1[ne]+
	  ((*grid)->yv[n]-maingrid->yp[(*grid)->edges[ne*NUMEDGECOLUMNS]])*(*grid)->n2[ne])*
	(*grid)->normal[n*NFACES+nf];

      // Distance to the edge midpoint. Not used.
      /*
      (*grid)->def[n*NFACES+nf]=sqrt(pow((*grid)->xv[n]-(*grid)->xe[ne],2)+
				     pow((*grid)->yv[n]-(*grid)->ye[ne],2));
      */
    }
  }
    /*
    for(nf=0;nf<NFACES;nf++) {
      xt[nf]=maingrid->xp[(*grid)->cells[n*NFACES+nf]];
      yt[nf]=maingrid->yp[(*grid)->cells[n*NFACES+nf]];
    }
    // Radius of the circumcircle
    R0 = GetCircumcircleRadius(xt,yt,NFACES);
    for(nf=0;nf<NFACES;nf++) {
      (*grid)->def[n*NFACES+nf]=sqrt(R0*R0-pow((*grid)->df[(*grid)->face[n*NFACES+nf]]/2,2));
      if(IsNan((*grid)->def[n*NFACES+nf]) || (*grid)->def[n*NFACES+nf]==0) {
	printf("Corrected at x=%f\n",(*grid)->xv[n]);
	(*grid)->def[n*NFACES+nf]=(*grid)->dg[(*grid)->face[n*NFACES+nf]]/2;
      }
    }
  }
    */

  /* Now compute the coefficients that make up the tangents to compute advection */
  if(myproc==0 && VERBOSE>2) printf("Computing xi coefficients...\n");
  for(n=0;n<Ne;n++) 
    for(nf=0;nf<2*(NFACES-1);nf++) 
      (*grid)->xi[2*(NFACES-1)*n+nf]=0;
  for(n=0;n<Ne;n++) 
    if((*grid)->n1[n]!=0 && (*grid)->n2[n]!=0) {
      for(nf=0;nf<2;nf++) {
	j = (*grid)->eneigh[2*(NFACES-1)*n+2*nf];
	k = (*grid)->eneigh[2*(NFACES-1)*n+2*nf+1];
	if(k!=-1 && j!=-1) {
	  den = sqrt(pow((*grid)->n1[k]*(*grid)->n1[n]+
			 (*grid)->n2[k]*(*grid)->n2[n],2)-
		     2.0*((*grid)->n1[j]*(*grid)->n1[n]+
			  (*grid)->n2[j]*(*grid)->n2[n])*
		     ((*grid)->n1[j]*(*grid)->n1[k]+
		      (*grid)->n2[j]*(*grid)->n2[k])*
		     ((*grid)->n1[k]*(*grid)->n1[n]+
		      (*grid)->n2[k]*(*grid)->n2[n])+
		     pow((*grid)->n1[j]*(*grid)->n1[n]+
			 (*grid)->n2[j]*(*grid)->n2[n],2));
	  (*grid)->xi[2*(NFACES-1)*n+2*nf]=
	    ((*grid)->n1[k]*(*grid)->n1[n]+(*grid)->n2[k]*(*grid)->n2[n])/den;
	  (*grid)->xi[2*(NFACES-1)*n+2*nf+1]=
	    -((*grid)->n1[j]*(*grid)->n1[n]+(*grid)->n2[j]*(*grid)->n2[n])/den;
	}
      }
    }
}

REAL GetCircumcircleRadius(REAL *xt, REAL *yt, int Nf)
{
  int nf;
  REAL R0 = 1;
  for(nf=0;nf<Nf;nf++)
    if(nf==Nf-1)
      R0*=sqrt(pow(xt[Nf-1]-xt[0],2)+pow(yt[Nf-1]-yt[0],2));
    else
      R0*=sqrt(pow(xt[nf+1]-xt[nf],2)+pow(yt[nf+1]-yt[nf],2));
  R0/=(4*GetArea(xt,yt,Nf));
  
  return R0;
}
  
REAL GetArea(REAL *xt, REAL *yt, int Nf)
{
  REAL b,r1,r2,h,l,xt2[NFACES],yt2[NFACES],area;

  if(Nf==3) {
    b=sqrt(pow(xt[1]-xt[0],2)+
	   pow(yt[1]-yt[0],2));
    r1=sqrt(pow(xt[2]-xt[0],2)+
	   pow(yt[2]-yt[0],2));
    r2=sqrt(pow(xt[2]-xt[1],2)+
	   pow(yt[2]-yt[1],2));
    l = ((xt[0]-xt[1])*(xt[2]-xt[1])+(yt[0]-yt[1])*(yt[2]-yt[1]))/b;
    h = sqrt(r2*r2-l*l);
    return 0.5*b*h;
  }
  else 
    if(Nf==4) {
      area=GetArea(xt,yt,3);
      xt2[0]=xt[2];
      xt2[1]=xt[3];
      xt2[2]=xt[0];
      yt2[0]=yt[2];
      yt2[1]=yt[3];
      yt2[2]=yt[0];
      area+=GetArea(xt2,yt2,3);
      return area;
    }
  return 0;
}

static void InterpDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm)
{
  int n, Nd, proc, nstart, ncount, scaledepth;
  REAL *xd, *yd, *d, scaledepthfactor;
  char str[BUFFERLENGTH];
  FILE *ifile, *ofile;
  MPI_Status status;

  scaledepth=(int)MPI_GetValue(DATAFILE,"scaledepth","InterpDepth",myproc);
  scaledepthfactor=MPI_GetValue(DATAFILE,"scaledepthfactor","InterpDepth",myproc);

  Nd = MPI_GetSize(INPUTDEPTHFILE,"InterpDepth",myproc);
  xd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpDepth");
  yd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpDepth");
  d = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpDepth");

  ifile = MPI_FOpen(INPUTDEPTHFILE,"r","InterpDepth",myproc);
  for(n=0;n<Nd;n++) {
    xd[n]=getfield(ifile,str);
    yd[n]=getfield(ifile,str);
    d[n]=fabs(getfield(ifile,str));
  }

  nstart = myproc*floor(grid->Nc/numprocs);
  if(myproc==numprocs-1 && grid->Nc%numprocs)
    ncount = grid->Nc - nstart;
  else
    ncount = floor(grid->Nc/numprocs);

  Interp(xd,yd,d,Nd,&(grid->xv[nstart]),
	 &(grid->yv[nstart]),&(grid->dv[nstart]),ncount);

  if(scaledepth)
    for(n=nstart;n<nstart+ncount;n++)
      grid->dv[n]*=scaledepthfactor;

  if(myproc!=0) 
    MPI_Send((void *)(&(grid->dv[nstart])),ncount,MPI_DOUBLE,0,1,comm); 
  else {
    for(proc=1;proc<numprocs;proc++) {
      nstart = proc*floor(grid->Nc/numprocs);
      if(proc==numprocs-1 && grid->Nc%numprocs)
	ncount = grid->Nc - nstart;
      else
	ncount = floor(grid->Nc/numprocs);

      MPI_Recv((void *)(&(grid->dv[nstart])),ncount,MPI_DOUBLE,proc,1,comm,&status);
    }
  }
  MPI_Bcast((void *)grid->dv,grid->Nc,MPI_DOUBLE,0,comm);

  if(myproc==0) {
    sprintf(str,"%s-voro",INPUTDEPTHFILE);
    ofile = MPI_FOpen(str,"w","InterpDepth",myproc);
    for(n=0;n<grid->Nc;n++) 
      fprintf(ofile,"%f %f %f\n",grid->xv[n],grid->yv[n],grid->dv[n]);
    fclose(ofile);
  }

  free(xd);
  free(yd);
  free(d);
}

static int CorrectVoronoi(gridT *grid, int myproc)
{
  int n, nf, nc1, nc2, numcorr=0;
  REAL xc, yc, xv1, xv2, yv1, yv2, xc1, xc2, yc1, yc2, dg, dg0;
  REAL VoronoiRatio=MPI_GetValue(DATAFILE,"VoronoiRatio","CorrectVoronoi",0);

  for(n=0;n<grid->Ne;n++) {
    nc1 = grid->grad[2*n];
    nc2 = grid->grad[2*n+1];
    
    if(nc1 != -1 && nc2 != -1) {
      xv1 = grid->xv[nc1];
      xv2 = grid->xv[nc2];
      yv1 = grid->yv[nc1];
      yv2 = grid->yv[nc2];
      xc1 = 0;
      xc2 = 0;
      yc1 = 0;
      yc2 = 0;
      for(nf=0;nf<NFACES;nf++) {
	xc1 += grid->xp[grid->cells[nc1*NFACES+nf]]/NFACES;
	xc2 += grid->xp[grid->cells[nc2*NFACES+nf]]/NFACES;
	yc1 += grid->yp[grid->cells[nc1*NFACES+nf]]/NFACES;
	yc2 += grid->yp[grid->cells[nc2*NFACES+nf]]/NFACES;
      }
      xc = 0.5*(xc1+xc2);
      yc = 0.5*(yc1+yc2);
      dg0 = sqrt(pow(xc2-xc1,2)+pow(yc2-yc1,2));
      dg = sqrt(pow(xv2-xv1,2)+pow(yv2-yv1,2));
      if(dg < VoronoiRatio*dg0) {
	if(VERBOSE>3) printf("Correcting Voronoi points %d and %d.\n",nc1,nc2);
	grid->xv[nc1]=xc+VoronoiRatio*(xc1-xc);
	grid->xv[nc2]=xc+VoronoiRatio*(xc2-xc);
	grid->yv[nc1]=yc+VoronoiRatio*(yc1-yc);
	grid->yv[nc2]=yc+VoronoiRatio*(yc2-yc);
	numcorr++;
      }
    } else {
      xc = 0.5*(grid->xp[grid->edges[NUMEDGECOLUMNS*n]]+
		grid->xp[grid->edges[NUMEDGECOLUMNS*n+1]]);
      yc = 0.5*(grid->yp[grid->edges[NUMEDGECOLUMNS*n]]+
		grid->yp[grid->edges[NUMEDGECOLUMNS*n+1]]);
      if(nc1 == -1) {
	xv2 = grid->xv[nc2];
	yv2 = grid->yv[nc2];
	xc2 = 0;
	yc2 = 0;
	for(nf=0;nf<NFACES;nf++) {
	  xc2 += grid->xp[grid->cells[nc2*NFACES+nf]]/NFACES;
	  yc2 += grid->yp[grid->cells[nc2*NFACES+nf]]/NFACES;
	}
	dg0 = sqrt(pow(xc2-xc,2)+pow(yc2-yc,2));
	dg = sqrt(pow(xv2-xc,2)+pow(yv2-yc,2));
	if(dg < VoronoiRatio*dg0) {
	  if(VERBOSE>3) printf("Correcting Voronoi point %d.\n",nc2);
	  grid->xv[nc2]=xc+VoronoiRatio*(xc2-xc);
	  grid->yv[nc2]=yc+VoronoiRatio*(yc2-yc);
	  numcorr++;
	}
      } else {
	xv1 = grid->xv[nc1];
	yv1 = grid->yv[nc1];
	xc1 = 0;
	yc1 = 0;
	for(nf=0;nf<NFACES;nf++) {
	  xc1 += grid->xp[grid->cells[nc1*NFACES+nf]]/NFACES;
	  yc1 += grid->yp[grid->cells[nc1*NFACES+nf]]/NFACES;
	}
	dg0 = sqrt(pow(xc1-xc,2)+pow(yc1-yc,2));
	dg = sqrt(pow(xv1-xc,2)+pow(yv1-yc,2));
	if(dg < VoronoiRatio*dg0) {
	  if(VERBOSE>3) printf("Correcting Voronoi point %d.\n",nc1);
	  grid->xv[nc1]=xc+VoronoiRatio*(xc1-xc);
	  grid->yv[nc1]=yc+VoronoiRatio*(yc1-yc);
	  numcorr++;
	}
      }
    }
  }
  if(numcorr>0 && VERBOSE>1 && myproc==0)
    printf("Corrected %d of %d edges with dg < %.2f dg0 (%.2f%%).\n",
	   numcorr,grid->Ne,VoronoiRatio,(REAL)numcorr/(REAL)grid->Ne*100.0);
  
  return numcorr;
}

static int CorrectAngles(gridT *grid, int myproc)
{
  int n, nf, nc1, nc2, numcorr;
  REAL xc, yc, xv1, xv2, yv1, yv2, xc1, xc2, xc3, yc1, yc2, yc3, dg, dg0, ang1, ang2,
    dot1, dot2, dot3, mag1, mag2, mag3, cosVoronoiRatio;
  REAL VoronoiRatio=MPI_GetValue(DATAFILE,"VoronoiRatio","CorrectVoronoi",0);
  cosVoronoiRatio = cos(VoronoiRatio*PI/180.0);

  numcorr=0;
  for(n=0;n<grid->Nc;n++) {
    xc1 = grid->xp[grid->cells[n*NFACES]];
    xc2 = grid->xp[grid->cells[n*NFACES+1]];
    xc3 = grid->xp[grid->cells[n*NFACES+2]];
    yc1 = grid->yp[grid->cells[n*NFACES]];
    yc2 = grid->yp[grid->cells[n*NFACES+1]];
    yc3 = grid->yp[grid->cells[n*NFACES+2]];
    mag1 = sqrt(pow(xc2-xc1,2)+pow(yc2-yc1,2));
    mag2 = sqrt(pow(xc3-xc1,2)+pow(yc3-yc1,2));
    mag3 = sqrt(pow(xc3-xc2,2)+pow(yc3-yc2,2));
    dot1 = (mag1*mag1+mag2*mag2-mag3*mag3)/(2*mag1*mag2);
    dot2 = (mag2*mag2+mag3*mag3-mag1*mag1)/(2*mag2*mag3);
    dot3 = (mag3*mag3+mag1*mag1-mag2*mag2)/(2*mag3*mag1);

    if(dot1<=cosVoronoiRatio ||
       dot2<=cosVoronoiRatio ||
       dot3<=cosVoronoiRatio) {
      numcorr++;
      grid->xv[n] = (xc1+xc2+xc3)/NFACES;
      grid->yv[n] = (yc1+yc2+yc3)/NFACES;
    }
  }
  if(numcorr>0 && VERBOSE>1 && myproc==0)
    printf("Corrected %d of %d cells with angles > %.1f degrees (%.2f%%).\n",
	   numcorr,grid->Nc,VoronoiRatio,(REAL)numcorr/(REAL)grid->Nc*100.0);

  return numcorr;
}

static void VoronoiStats(gridT *grid) {
  int i, n, nc1, nc2, numdg, dghist;
  REAL *dg, dgmin, dgmax, dgmean, dgstd, deldg, histmin, histmax;
  dg = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"VoronoiStats");

  numdg=0;
  for(n=0;n<grid->Ne;n++) {
    nc1=grid->grad[2*n];
    nc2=grid->grad[2*n+1];
    if(nc2>=0) {
      dg[numdg] = sqrt(pow(grid->xv[nc1]-grid->xv[nc2],2)+
		       pow(grid->yv[nc2]-grid->yv[nc1],2));
      numdg++;
    }
  }

  dgmin=INFTY;
  dgmax=dgmean=0;
  for(n=0;n<numdg;n++) {
    if(dg[n]<dgmin)
      dgmin=dg[n];
    if(dg[n]>dgmax)
      dgmax=dg[n];
    dgmean+=dg[n];
  }
  dgmean/=(REAL)numdg;

  dgstd=0;
  for(n=0;n<numdg;n++) 
    dgstd+=pow(dg[n]-dgmean,2);
  dgstd=sqrt(dgstd/(REAL)numdg);

  printf("\tMinimum distance: %.2e\n",dgmin);
  printf("\tMaximum distance: %.2e\n",dgmax);
  printf("\tMean distance: %.2e\n",dgmean);
  printf("\tStandard deviation: %.2e\n",dgstd);

  if(VERBOSE>2) {
    printf("\tVoronoi histogram:\n");
    printf("\t\t%.2e\t0\n",dgmin);
    deldg = (dgmax-dgmin)/(REAL)(NUMDGHIST-1);
    for(i=0;i<NUMDGHIST-1;i++) {
      dghist=0;
      histmin = dgmin+deldg*(REAL)i;
      histmax = histmin+deldg;
      for(n=0;n<numdg;n++) {
	if((dg[n]>=histmin &&
	    dg[n]<histmax) || 
	   (i==NUMDGHIST-2 && dg[n]==histmax))
	  dghist++;
      }
      printf("\t\t%.2e\t%d\n",histmax,dghist);
    }
  }

  SunFree(dg,grid->Ne*sizeof(REAL),"VoronoiStats");
}

/*
 * Function: AllocateTransferArrays
 * Usage: AllocateTransferArrays(&grid,myproc,numprocs,comm);
 * ----------------------------------------------------------
 * Allocate memory for the arrays used to send back and forth interprocessor data.
 *
 */
void AllocateTransferArrays(gridT **grid, int myproc, int numprocs, MPI_Comm comm) {
  int n, neigh, neighproc, maxtosend, maxtorecv;
  (*grid)->status = (MPI_Status *)SunMalloc(2*(*grid)->Nneighs*sizeof(MPI_Status),"AllocateTransferrays");
  (*grid)->request = (MPI_Request *)SunMalloc(2*(*grid)->Nneighs*sizeof(MPI_Request),"AllocateTransferrays");

  (*grid)->recv = (REAL **)SunMalloc((*grid)->Nneighs*sizeof(REAL *),"AllocateTransferArrays");
  (*grid)->send = (REAL **)SunMalloc((*grid)->Nneighs*sizeof(REAL *),"AllocateTransferArrays");
  (*grid)->total_cells_send = (int *)SunMalloc((*grid)->Nneighs*sizeof(int),"AllocateTransferArrays");
  (*grid)->total_cells_recv = (int *)SunMalloc((*grid)->Nneighs*sizeof(int),"AllocateTransferArrays");
  (*grid)->total_cells_sendW = (int *)SunMalloc((*grid)->Nneighs*sizeof(int),"AllocateTransferArrays");
  (*grid)->total_cells_recvW = (int *)SunMalloc((*grid)->Nneighs*sizeof(int),"AllocateTransferArrays");
  (*grid)->total_edges_send = (int *)SunMalloc((*grid)->Nneighs*sizeof(int),"AllocateTransferArrays");
  (*grid)->total_edges_recv = (int *)SunMalloc((*grid)->Nneighs*sizeof(int),"AllocateTransferArrays");

  for(neigh=0;neigh<(*grid)->Nneighs;neigh++) {
    neighproc = (*grid)->myneighs[neigh];

    (*grid)->total_cells_send[neigh] = 0;
    (*grid)->total_cells_recv[neigh] = 0;
    (*grid)->total_cells_sendW[neigh] = 0;
    (*grid)->total_cells_recvW[neigh] = 0;
    (*grid)->total_edges_send[neigh] = 0;
    (*grid)->total_edges_recv[neigh] = 0;
    for(n=0;n<(*grid)->num_cells_send[neigh];n++) 
      (*grid)->total_cells_send[neigh]+=(*grid)->Nk[(*grid)->cell_send[neigh][n]];
    for(n=0;n<(*grid)->num_cells_recv[neigh];n++) 
      (*grid)->total_cells_recv[neigh]+=(*grid)->Nk[(*grid)->cell_recv[neigh][n]];
    for(n=0;n<(*grid)->num_cells_send[neigh];n++) 
      (*grid)->total_cells_sendW[neigh]+=(1+(*grid)->Nk[(*grid)->cell_send[neigh][n]]);
    for(n=0;n<(*grid)->num_cells_recv[neigh];n++) 
      (*grid)->total_cells_recvW[neigh]+=(1+(*grid)->Nk[(*grid)->cell_recv[neigh][n]]);
    for(n=0;n<(*grid)->num_edges_send[neigh];n++) 
      (*grid)->total_edges_send[neigh]+=(*grid)->Nke[(*grid)->edge_send[neigh][n]];
    for(n=0;n<(*grid)->num_edges_recv[neigh];n++) 
      (*grid)->total_edges_recv[neigh]+=(*grid)->Nke[(*grid)->edge_recv[neigh][n]];

    if((*grid)->total_cells_sendW[neigh]>(*grid)->total_edges_send[neigh])
      (*grid)->maxtosend=(*grid)->total_cells_sendW[neigh];
    else
      (*grid)->maxtosend=(*grid)->total_edges_send[neigh];
    if((*grid)->total_cells_recvW[neigh]>(*grid)->total_edges_recv[neigh])
      (*grid)->maxtorecv=(*grid)->total_cells_recvW[neigh];
    else
      (*grid)->maxtorecv=(*grid)->total_edges_recv[neigh];
    (*grid)->send[neigh] = (REAL *)SunMalloc((*grid)->maxtosend*sizeof(REAL),"AllocateTransferArrays");
    (*grid)->recv[neigh] = (REAL *)SunMalloc((*grid)->maxtorecv*sizeof(REAL),"AllocateTransferArrays");
  }
}

/*
 * Function: FreeTransferArrays
 * Usage: FreeTransferArrays(&grid,myproc,numprocs,comm);
 * ------------------------------------------------------
 * Free memory for the arrays used to send back and forth interprocessor data.
 *
 */
void FreeTransferArrays(gridT *grid, int myproc, int numprocs, MPI_Comm comm) {
  int neigh, neighproc, maxtosend, maxtorecv;

  SunFree(grid->status,2*grid->Nneighs*sizeof(MPI_Status),"FreeTransferArrays");
  SunFree(grid->request,2*grid->Nneighs*sizeof(MPI_Request),"FreeTransferArrays");

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    if(grid->total_cells_sendW[neigh]>grid->total_edges_send[neigh])
      maxtosend=grid->total_cells_sendW[neigh];
    else
      maxtosend=grid->total_edges_send[neigh];
    if(grid->total_cells_recvW[neigh]>grid->total_edges_recv[neigh])
      maxtorecv=grid->total_cells_recvW[neigh];
    else
      maxtorecv=grid->total_edges_recv[neigh];
    SunFree(grid->send[neigh],maxtosend*sizeof(REAL),"FreeTransferArrays");
    SunFree(grid->recv[neigh],maxtorecv*sizeof(REAL),"FreeTransferArrays");
  }
  SunFree(grid->total_cells_send,grid->Nneighs*sizeof(int),"FreeTransferArrays");
  SunFree(grid->total_cells_recv,grid->Nneighs*sizeof(int),"FreeTransferArrays");
  SunFree(grid->total_cells_sendW,grid->Nneighs*sizeof(int),"FreeTransferArrays");
  SunFree(grid->total_cells_recvW,grid->Nneighs*sizeof(int),"FreeTransferArrays");
  SunFree(grid->total_edges_send,grid->Nneighs*sizeof(int),"FreeTransferArrays");
  SunFree(grid->total_edges_recv,grid->Nneighs*sizeof(int),"FreeTransferArrays");
  SunFree(grid->send,grid->Nneighs*sizeof(REAL *),"FreeTransferArrays");
  SunFree(grid->recv,grid->Nneighs*sizeof(REAL *),"FreeTransferArrays");
}

/*
 * Function: ISendRecvCellData2D
 * Usage: ISendRecvCellData2D(grid->h,grid,myproc,comm);
 * -----------------------------------------------------
 * This function will transfer the cell data back and forth between
 * processors.  
 *
 * This is the non-blocking version of SendRecvCellData2D.
 *
 */
void ISendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, 
			 MPI_Comm comm)
{
  int n, neigh, neighproc;
  REAL t0=Timer();

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    for(n=0;n<grid->num_cells_send[neigh];n++)
      grid->send[neigh][n]=celldata[grid->cell_send[neigh][n]];

    MPI_Isend((void *)(grid->send[neigh]),grid->num_cells_send[neigh],
	     MPI_DOUBLE,neighproc,1,comm,&(grid->request[neigh])); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Irecv((void *)(grid->recv[neigh]),grid->num_cells_recv[neigh],
	     MPI_DOUBLE,neighproc,1,comm,&(grid->request[grid->Nneighs+neigh]));
  }
  MPI_Waitall(2*grid->Nneighs,grid->request,grid->status);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    for(n=0;n<grid->num_cells_recv[neigh];n++)
      celldata[grid->cell_recv[neigh][n]]=grid->recv[neigh][n];
  }
  t_comm+=Timer()-t0;
}

/*
 * Function: ISendRecvCellData3D
 * Usage: ISendRecvCellData3D(grid->s,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the 3D cell data back and forth between
 * processors using nonblocking sends/recvs.
 *
 */
void ISendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, 
			MPI_Comm comm)
{
  int k, n, nstart, neigh, neighproc;
  REAL t0=Timer();

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    nstart=0;
    for(n=0;n<grid->num_cells_send[neigh];n++) {
      for(k=0;k<grid->Nk[grid->cell_send[neigh][n]];k++) 
	grid->send[neigh][nstart+k]=celldata[grid->cell_send[neigh][n]][k];
      nstart+=grid->Nk[grid->cell_send[neigh][n]];
    }

    MPI_Isend((void *)(grid->send[neigh]),grid->total_cells_send[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[neigh])); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Irecv((void *)(grid->recv[neigh]),grid->total_cells_recv[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[grid->Nneighs+neigh]));
  }
  MPI_Waitall(2*grid->Nneighs,grid->request,grid->status);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<grid->num_cells_recv[neigh];n++) {
      for(k=0;k<grid->Nk[grid->cell_recv[neigh][n]];k++) 
	celldata[grid->cell_recv[neigh][n]][k]=grid->recv[neigh][nstart+k];
      nstart+=grid->Nk[grid->cell_recv[neigh][n]];
    }
  }
  t_comm+=Timer()-t0;
}

/*
 * Function: ISendRecvWData
 * Usage: ISendRecvWData(grid->w,grid,myproc,comm);
 * -----------------------------------------------
 * This function will transfer the 3D w data back and forth between
 * processors using nonblocking sends/recvs.
 *
 */
void ISendRecvWData(REAL **celldata, gridT *grid, int myproc, 
		   MPI_Comm comm)
{
  int k, n, nstart, neigh, neighproc;
  REAL t0=Timer();

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    nstart=0;
    for(n=0;n<grid->num_cells_send[neigh];n++) {
      for(k=0;k<(1+grid->Nk[grid->cell_send[neigh][n]]);k++) 
	grid->send[neigh][nstart+k]=celldata[grid->cell_send[neigh][n]][k];
      nstart+=(1+grid->Nk[grid->cell_send[neigh][n]]);
    }

    MPI_Isend((void *)(grid->send[neigh]),grid->total_cells_sendW[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[neigh])); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Irecv((void *)(grid->recv[neigh]),grid->total_cells_recvW[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[grid->Nneighs+neigh]));
  }
  MPI_Waitall(2*grid->Nneighs,grid->request,grid->status);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<grid->num_cells_recv[neigh];n++) {
      for(k=0;k<(1+grid->Nk[grid->cell_recv[neigh][n]]);k++) 
	celldata[grid->cell_recv[neigh][n]][k]=grid->recv[neigh][nstart+k];
      nstart+=(1+grid->Nk[grid->cell_recv[neigh][n]]);
    }
  }
  t_comm+=Timer()-t0;
}

/*
 * Function: ISendRecvEdgeData3D
 * Usage: ISendRecvEdgeData3D(grid->u,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the 3D edge data back and forth between
 * processors using onblocking sends/recvs.
 *
 */
void ISendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, 
			MPI_Comm comm)
{
  int k, n, nstart, neigh, neighproc;
  REAL t0=Timer();

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    nstart=0;
    for(n=0;n<grid->num_edges_send[neigh];n++) {
      for(k=0;k<grid->Nke[grid->edge_send[neigh][n]];k++) 
	grid->send[neigh][nstart+k]=edgedata[grid->edge_send[neigh][n]][k];
      nstart+=grid->Nke[grid->edge_send[neigh][n]];
    }

    MPI_Isend((void *)(grid->send[neigh]),grid->total_edges_send[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[neigh])); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Irecv((void *)(grid->recv[neigh]),grid->total_edges_recv[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[grid->Nneighs+neigh]));
  }
  MPI_Waitall(2*grid->Nneighs,grid->request,grid->status);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<grid->num_edges_recv[neigh];n++) {
      for(k=0;k<grid->Nke[grid->edge_recv[neigh][n]];k++) {
	edgedata[grid->edge_recv[neigh][n]][k]=grid->recv[neigh][nstart+k];
      }
      nstart+=grid->Nke[grid->edge_recv[neigh][n]];
    }
  }
  t_comm+=Timer()-t0;
}

/*
 * Function: FixDZZ
 * Usage: FixDZZ(grid,maxdepth,Nkmax,fixdzz,myproc);
 * -------------------------------------------------
 * Check to make sure that the deepest vertical grid spacing is not too small.
 * and if so then increase the depth.
 *
 */

static void FixDZZ(gridT *grid, REAL maxdepth, int Nkmax, int fixdzz, int myproc) {
  int i, k, kount=0, mindepth0;
  REAL z, dzz, dzsmall, *dz = (REAL *)SunMalloc(Nkmax*sizeof(REAL),"FixDZZ");

  dzsmall=(REAL)MPI_GetValue(DATAFILE,"dzsmall","FixDZZ",myproc);
  
  GetDZ(dz,maxdepth,maxdepth,Nkmax,myproc);  
  for(i=0;i<grid->Nc;i++) {
    z=0;
    for(k=0;k<Nkmax;k++) {
      if(z>-grid->dv[i] && z-dz[k]<-grid->dv[i]) {
	dzz=z+grid->dv[i];
	
	if(fixdzz>0) {
	  if(dzz<dz[k]*dzsmall && k>0) {
	    if(myproc==0 && VERBOSE>2) printf("Fixing small bottom dz of %.2e by increasing the depth from %.2f to %.2f\n",
					      dzz,grid->dv[i],-z+dz[k]*dzsmall);
	    kount++;
	    grid->dv[i]=-z+dz[k]*dzsmall;
	  }
	} else {
	  if(dzz<dzsmall && k>0) {
	    if(myproc==0 && VERBOSE>2) printf("Fixing small bottom dz of %.2e by increasing the depth from %.2f to %.2f\n",
					      dzz,grid->dv[i],-z+dzsmall);
	    kount++;
	    grid->dv[i]=-z+dzsmall;
	  }
	}	  
	break;
      }
      z-=dz[k];
    }
  }
  if(VERBOSE>0 && myproc==0 && kount!=0) {
    printf("Fixed %d bottom cells with heights < %.2e dz.\n",kount,dzsmall);
    printf("To eliminate this action set fixdzz to 0 in suntans.dat.\n");
  }

  SunFree(dz,Nkmax*sizeof(REAL),"FixDZZ");
}



static int GetNk(REAL *dz, REAL localdepth, int Nkmax) {
  int k;
  REAL z=0;

  for(k=0;k<Nkmax;k++) {
    z-=dz[k];
    if(z < -localdepth)
      break;
  }
  return k;
}
  
  

