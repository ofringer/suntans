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
 * Global variables for halo region for inner-processor boundaries
 */
//boundaryselection g_halolist[1] ={EDGE};// classic suntans
//const int g_halolistsize = 1;
//// for quadratic interpolation
boundaryselection g_halolist[2] ={NODE, NODE}; 
const int g_halolistsize = 2;

/*
 * Private function declarations.
 */
static void InitLocalGrid(gridT **grid);
static void VertGrid(gridT *maingrid, gridT **localgrid, MPI_Comm comm);
static void Topology(gridT **maingrid, gridT **localgrid, int myproc, int numprocs);
static int GetNumCells(gridT *grid, int proc);
static void TransferData(gridT *maingrid, gridT **localgrid, int myproc);
static int GetNumEdges(gridT *grid);
static int GetNumPoints(gridT *localgrid, gridT *maingrid);
static void GetProcPoints(gridT *localgrid, gridT *maingrid);
static void Geometry(gridT *maingrid, gridT **grid, int myproc);
static REAL GetCircumcircleRadius(REAL *xt, REAL *yt, int Nf);
static void EdgeMarkers(gridT *maingrid, gridT **localgrid, int myproc);
static void ReOrder(gridT *grid);
static int IsCellNeighborProc(int nc, gridT *maingrid, gridT *localgrid, 
    int myproc, int neighproc);
static int IsEdgeNeighborProc(int ne, gridT *maingrid, gridT *localgrid, 
    int myproc, int neighproc);
static int IsNodeNeighborProc(int np, gridT *maingrid, gridT *localgrid, 
    int myproc, int neighproc);
static int IsNeighborGlobal(int cell, gridT *maingrid, 
    gridT **localgrid, int myproc, boundaryselection *list, int listsize);
static int IsNeighborGlobalByEdge(int cell, gridT *maingrid, 
    gridT **localgrid, int myproc, boundaryselection *list, int listsize);
static int IsNeighborGlobalByNode(int cell, gridT *maingrid, 
    gridT **localgrid, int myproc, boundaryselection *list, int listsize);
static inline int IsNeighborLocal(int cell, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize);
static int IsNeighborLocalByEdge(int cell, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize);
static int IsNeighborLocalByNode(int cell, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize);
static int IsBoundaryCellByHalo(int mgptr, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize);
static int IsBoundaryCellByEdge(int mgptr, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize);
static int IsBoundaryCellByNode(int mgptr, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize);
static inline void  BoundaryTopologyByHalo(gridT *maingrid, int numprocs,
    boundaryselection *list, int listsize);
static void  BoundaryTopologyByEdge(int cellbase, int cell, gridT *maingrid, 
    int numprocs, boundaryselection *list, int listsize);
static void  BoundaryTopologyByNode(int cellbase, int cell, gridT *maingrid, 
    int numprocs, boundaryselection *list, int listsize);
static void MakePointers(gridT *maingrid, gridT **localgrid, int myproc, 
    MPI_Comm comm);
static void ResortBoundaries(gridT *localgrid, int myproc);
static void InterpDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static void FreeGrid(gridT *grid, int numprocs);
static void OutputData(gridT *maingrid, gridT *grid, int myproc, int numprocs);
static inline void CreateFaceArray(int *grad, int *gradf, int *neigh, int *face, int *nfaces, int maxfaces, int Nc, int Ne);
static inline void ReorderCellPoints(int *face, int *edges, int *cells, int *nfaces, int maxfaces, int Nc);
static inline void Reordergradf(int *face, int *grad, int *gradf, int *nfaces, int maxfaces, int Nc, int Ne);
static inline void CreateNodeArray(gridT *grid, int Np, int Ne, int Nc, int myproc);
static inline void CreateNormalArray(int *grad, int *face, int *normal, int *nfaces, int maxfaces, int Nc);
static inline void CreateArray(int *grad, int *gradf, int *neigh, int *face, 
    int Nc, int Ne);
static inline void CreateMomentumCV(gridT *grid);
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
  int Np, Ne, Nc, numcorr, maxFaces, i;
  gridT *maingrid;

  // read in all the filenames from the suntans.dat file
  ReadFileNames(myproc);

  //  compute the grid with triangle or read it in if in suntans format
  if(!TRIANGULATE) {
    // get the number of points, edges and cells from suntans.dat
    Np = MPI_GetSize(POINTSFILE,"GetGrid",myproc);
    Ne = MPI_GetSize(EDGEFILE,"GetGrid",myproc);
    Nc = MPI_GetSize(CELLSFILE,"GetGrid",myproc);
    maxFaces=(int)MPI_GetValue(DATAFILE,"maxFaces","GetGrid",0);
    
    // Every processor will know about data read in from
    // triangle as well as the depth.
    if(myproc==0 && VERBOSE>0) printf("Initializing Main Grid...\n");
    // allocate memory and define main grid structure members
    InitMainGrid(&maingrid,Np,Ne,Nc,myproc,maxFaces);

    if(myproc==0 && VERBOSE>0) printf("Reading Grid...\n");
    // read in cells, edges, points files (suntans grid format files)
    // note that this only affects points: xp,yp, edges: edges, mark, grad
    // cells: xv, yv, cells, neigh
    ReadMainGrid(maingrid,myproc);
  } 
  //  if we need to compute the grid with the triangle library
  else {
    // get the triangulation via the triangle libraries
    if(myproc==0 && VERBOSE>0) printf("Triangulating the point set...\n");
    if(!GetTriangulation(&maingrid,myproc)) {
      printf("Error computing triangulation!\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  // Adjust Voronoi points for obtuse triangles (CorrectVoronoi==-1) 
  // or for small Voronoi edge lengths (CorrectVoronoi==+1).
  // this seems like it should be wrapped into its own function
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

  // get graph weights
  if(myproc==0 && VERBOSE>0) printf("Getting the depth for graph weights...\n");
  Tic();
  GetDepth(maingrid,myproc,numprocs,comm);
  if(myproc==0 && VERBOSE>0) printf("... time used is %f\n", Toc());

  // functions to set up parallelization of the code and compute geometrical 
  // information needed for grid computations

  // makes the face and normal arrays as well as the control volume for momentum 
  // based on the dual
  if(myproc==0 && VERBOSE>0) printf("Computing Connectivity...\n");
  Tic();
  Connectivity(maingrid,myproc);

  if(myproc==0 && VERBOSE>0) printf("... time used is %f\n", Toc());

  if(myproc==0 && VERBOSE>0) printf("Partitioning...\n");
  Tic();
  // most important and complicated function to ensure parallelism for grid
  // takes the main grid and redistributes to local grids on each processor 
  Partition(maingrid,localgrid,comm);

  if(myproc==0 && VERBOSE>0) printf("... time used is %f\n", Toc());
  
  // functions to check to make sure all the communications are set up properly
  //  CheckCommunicateCells(maingrid,*localgrid,myproc,comm);
  //  CheckCommunicateEdges(maingrid,*localgrid,myproc,comm);
  //  SendRecvCellData2D((*localgrid)->dv,*localgrid,myproc,comm);
 
  if(myproc==0 && VERBOSE>0) printf("Outputing Data...\n");
  Tic();

  OutputData(maingrid,*localgrid,myproc,numprocs);
  if(myproc==0 && VERBOSE>0) printf("... time used is %f\n", Toc());
  //FreeGrid(maingrid,numprocs);
}

/*
 * Function: Partition
 * Usage: Partition(maingrid, localgrid, comm)
 * ------------------------------------------
 * Partitions the cells unto each processor where ParMETIS 
 * colors each cell according to its processor based on 
 * weighting each cell by the depth (so that shallow 
 * regions of the domain which are less computationally
 * expensive are weighted less) in GetPartitioning.  
 * Topology determines neighboring processors for each
 * processor and TransferData determines which cells must 
 * have their information passed between processors.
 *
 */
void Partition(gridT *maingrid, gridT **localgrid, MPI_Comm comm)
{
  int j;
  int myproc, numprocs;

  // allocate memory for the local grid structures
  InitLocalGrid(localgrid);
  // is this redundant from suntans.c StartMPI function?
  // - no, since numprocs and proc are not defined here in the function (local vars)
  MPI_Comm_size(comm, &numprocs);
  MPI_Comm_rank(comm, &myproc);

  if(myproc==0 && VERBOSE>0) printf("\tComputing Graph...\n");
  // cell graph generates Nge, adjncy, xadj members of maingrid
  // Useful to create cell graph compatible with ParMETIS to "color" each
  // cell based on its partitioning on a processor in Partition:GetPartitioning
  Tic();
  CreateCellGraph(maingrid);
  if(myproc==0 && VERBOSE>0) printf("\t... time used is %f\n", Toc());

  // use ParMETIS API get get parallelized partioning, this basically 
  // just assigns each cell to a particular processor number  
  // (color cells with mapping!)
  if(myproc==0 && VERBOSE>0) printf("\tComputing Partitioning...\n");
  Tic();
  GetPartitioning(maingrid,localgrid,myproc,numprocs,comm);
  if(myproc==0 && VERBOSE>0) printf("\t... time used is %f\n", Toc());

//  // output a list of each cells processor for debugging
//  PrintVectorToFile(INT, maingrid->part, maingrid->Nc, "CellProcessors.dat",0);

  // location of these functions seems a little funny, they could be "popped" out 
  // of this function up a level for a little more clarity

  // topology basically just computes the number of neighbors "Nneighs"
  // each processor has and gives the list of the neighboring processors
  // which is also stored on the localgrid, "myneighs"
  // now we have a complete description of all the neighboring cell processors
  // for cells sharing an edges on the processor boundary.
  // NOte that this does not allow for cells to complete the whole layer, which may 
  // also be on a separate processor (like a wedge).  Hence, an additional loop in 
  // Topology is needed to consider the cases where a node is connected to an 
  // interprocessor edge 
  // which has a different processor number than node on the edge considered.  
  // This way the full amount of topology information can be gained.  
  // It is likely that TransferData will also have to be changed to account for this
  // behavior as it is no longer as simple as saying that we have an edge and its
  // associated neighs (cells on either side).  We now have an edge connected to 
  // two nodes which may be shared between multiple cells and the computational 
  // complexity and time to do the gridding increases 
  // this seems to necessitate the idea of a nodepoints.dat.proc type file where we 
  // have the node number p in 0 <= p <= Np, number of cell neighbors to node, numpneighs, 
  // a list of the cell numbers associated with each of the neighbors to the node
  // as in Pniegh (point neighbor).  For simplicity it is aslo probably good to 
  // store the number of the neighboring nodes (sharing an edge with p), particularly
  // to speed up the interpolation when it is time to do it. Thus the number of entries
  // will be free form.  This structure will be needed to adequately compute the 
  // number of neighboring processors and processor info for Topology
  // doing the nodal neighbor approach as opposed to edge approach will increase the 
  // number of neighbors that must be transfered (as all cells touching a boundary
  // node are now considered important ghost cells)
  // Topology needs updated so that any cell touching a node is considered to allow 
  // flagging as a boundary node
  if(myproc==0 && VERBOSE>0) printf("\tComputing Topology...\n");
  Tic();
  Topology(&maingrid,localgrid,myproc,numprocs);
  if(myproc==0 && VERBOSE>0) printf("\t... time used is %f\n", Toc());

  if(myproc==0 && VERBOSE>2) printf("\tTransferring data...\n");
  // function to move data between main and local grid where boundary points were
  // defined with convention of shared edge with number and list of neighbors 
  // described in Topology
  Tic(); 
  TransferData(maingrid,localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("\t... time used is %f\n", Toc());

  if(myproc==0 && VERBOSE>2) printf("\tCreating edge markers...\n");
  // get edge marker information where we designate mark=5 and mark=6 for interproc
  // boundary ghost cells and mark=6 for pure ghost cells
  Tic();
  EdgeMarkers(maingrid,localgrid,myproc);
  if(myproc==0 && VERBOSE>2) printf("\t... time used is %f\n", Toc());

  // note that Vert grid was moved above geometry so that Nkp can
  // be used in geometry
  if(myproc==0 && VERBOSE>2) printf("\tVert grid...\n");
  // get vertical grid, notably Nke, Nkc to distinguish number of layers at edges, cells
  Tic();
  VertGrid(maingrid,localgrid,comm);
  if(myproc==0 && VERBOSE>2) printf("\t... time used is %f\n", Toc());

  if(myproc==0 && VERBOSE>2) 
    printf("\tComputing edge and voronoi distances and areas...\n");
  Tic();
  Geometry(maingrid,localgrid,myproc);
  if(myproc==0 && VERBOSE>2) printf("\t... time used is %f\n", Toc());


  //if(myproc==0 && VERBOSE>2) printf("\tReordering...\n");
  // function to reorder the grid for increased computational efficiency, this 
  // is not currently known to work
  //  ReOrder(*localgrid);

  if(myproc==0 && VERBOSE>2) printf("\tMaking pointers...\n");
  // compute cellp, edgep, lcptr, leptrs and also get the reference for local cell processor
  // interprocessor send-receives
  Tic();
  MakePointers(maingrid,localgrid,myproc,comm);
  if(myproc==0 && VERBOSE>2) printf("\t... time used is %f\n", Toc());

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
void InitMainGrid(gridT **grid, int Np, int Ne, int Nc, int myproc, int maxFaces)
{ int n;
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
 
  // Max number of edges of cells
  (*grid)->maxfaces=maxFaces;
 
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

  // parts added
  // allocate nfaces and maxfaces
  (*grid)->nfaces= (int *)SunMalloc((*grid)->Nc*sizeof(int),"InitMainGrid");
  // Load the faces' number of each cell and calculate the maximum nfaces
  // if -t not use READNFACES
  if(!TRIANGULATE && maxFaces!=3){ 
    ReadNfaces(*grid, myproc, maxFaces);
  }else{ 
    for(n=0;n<Nc;n++)
      (*grid)->nfaces[n]=3;
   }

  // Pointers to xp,yp coordinates of vertices that make up polygons (0<cells<Np)
  (*grid)->cells = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"InitMainGrid");

  // Pointers to neighboring cells (0<neigh<Nc)
  (*grid)->neigh = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"InitMainGrid");
  // Dot product of unique normal with outward normal
  (*grid)->normal = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"InitMainGrid");
  // Indices of voronoi edges to cells
  (*grid)->grad = (int *)SunMalloc(2*(*grid)->Ne*sizeof(int),"InitMainGrid");
  // Indices of pointers to faces of each cell
  (*grid)->face = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"InitMainGrid");
  // Indices to edges for momentum control volume
  (*grid)->eneigh = (int *)SunMalloc(2*((*grid)->maxfaces-1)*(*grid)->Ne*sizeof(int),"InitMainGrid");

  // allocae nodal neighbor array information (list of neighbors for each point)
  // note that since the number of entries for each node is variable, 
  // these are a list to Np lists
  (*grid)->numppneighs = (int *)SunMalloc((*grid)->Np*sizeof(int),"InitMainGrid");
  (*grid)->ppneighs = (int **)SunMalloc((*grid)->Np*sizeof(int*),"InitMainGrid");
  (*grid)->numpeneighs = (int *)SunMalloc((*grid)->Np*sizeof(int),"InitMainGrid");
  (*grid)->peneighs = (int **)SunMalloc((*grid)->Np*sizeof(int*),"InitMainGrid");
  (*grid)->numpcneighs = (int *)SunMalloc((*grid)->Np*sizeof(int),"InitMainGrid");
  (*grid)->pcneighs = (int **)SunMalloc((*grid)->Np*sizeof(int*),"InitMainGrid");

  // Depth at Voronoi points
  (*grid)->dv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"InitMainGrid");
  // Weights for partitioning cell graph
  (*grid)->vwgt = (int *)SunMalloc((*grid)->Nc*sizeof(int),"InitMainGrid");

  // Assigned cell partition
  // (this is the number of processor cell is assigned to)
  (*grid)->part = (int *)SunMalloc((*grid)->Nc*sizeof(int),"InitMainGrid");
  // For ordering and making sure boundary edge data is contiguous.
  (*grid)->order = (int *)SunMalloc((*grid)->Ne*sizeof(int),"InitMainGrid");
  // Edge markers
  (*grid)->mark = (int *)SunMalloc((*grid)->Ne*sizeof(int),"InitMainGrid");
  // Stores the indices to the start and end points of the adjncy array
  (*grid)->xadj = (int *)SunMalloc(((*grid)->Nc+1)*sizeof(int),"InitMainGrid");
  (*grid)->vtxdist = (int *)SunMalloc(VTXDISTMAX*sizeof(int),"InitMainGrid");
}

/* added functions
 * Function: ReadNfaces
 * Usage: ReadNfaces(grid,myproc,maxfaces);
 * ---------------------------------
 * Read the first column in the cell data.
 * Create the array of nfaces for each cell
 * Calculate the maximum nfaces
 */ 
void ReadNfaces(gridT *grid, int myproc, int maxFaces)
{
  int n, max, i;
  char str[BUFFERLENGTH];
  FILE *ifile;
  ifile=MPI_FOpen(CELLSFILE,"r","ReadNfaces",myproc);
  max=0;

  for(n=0;n<grid->Nc;n++){
    grid->nfaces[n]=(int)getfield(ifile,str);
    if(grid->nfaces[n]>maxFaces){
    printf("Max number of edges is %d of No. %d in Cell.dat is more than maxFaces %d in suntans.dat!",grid->nfaces[n],n,maxFaces);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    }   
    for(i=0;i<(2*grid->nfaces[n]+2);i++)
      getfield(ifile,str);
  }    
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
 * EDGEFILE list of indices to points in POINTSFILE (always 2 columns + edge marker + 2 pointers to
 * neighboring cells (grad) = 5 columns)
 *
 */
void ReadMainGrid(gridT *grid, int myproc)
{
  int j, n, nei, nf,k;
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
    //added part, the first column is the number of faces for each cell, not xv 
    if(getcolumn(CELLSFILE)>8)
      getfield(ifile,str);
    grid->xv[n] = getfield(ifile,str);
    grid->yv[n] = getfield(ifile,str);
    for(nf=0;nf<grid->nfaces[n];nf++)
      grid->cells[n*grid->maxfaces+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<grid->nfaces[n];nf++) {
      if((nei=(int)getfield(ifile,str))!=-1)
	grid->neigh[n*grid->maxfaces+nf]=nei;
      else
	grid->neigh[n*grid->maxfaces+nf]=-1;
    }
  }
}

/*
 * Function: ReadFileNames
 * Usage: ReadFileNames(myproc);
 * -----------------------------
 * Reads the names of the files containing the grid data
 * from the file defined by GRIDDATAFILELIST in suntans.h
 * This is the global list mapping names to suntans.dat
 * 
 */
void ReadFileNames(int myproc)
{
  MPI_GetFile(POINTSFILE,DATAFILE,"points","ReadFileNames",myproc);
  MPI_GetFile(EDGEFILE,DATAFILE,"edges","ReadFileNames",myproc);
  MPI_GetFile(CELLSFILE,DATAFILE,"cells","ReadFileNames",myproc);
  MPI_GetFile(NODEFILE,DATAFILE,"nodes","ReadFileNames",myproc);
  MPI_GetFile(INPUTDEPTHFILE,DATAFILE,"depth","ReadFileNames",myproc);
  MPI_GetFile(CELLCENTEREDFILE,DATAFILE,"celldata","ReadFileNames",myproc);
  MPI_GetFile(EDGECENTEREDFILE,DATAFILE,"edgedata","ReadFileNames",myproc);
  MPI_GetFile(VERTSPACEFILE,DATAFILE,"vertspace","ReadFileNames",myproc);
  MPI_GetFile(TOPOLOGYFILE,DATAFILE,"topology","ReadFileNames",myproc);
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
    // interpolate the depth if desired 
    InterpDepth(grid,myproc,numprocs,comm);
  else if(IntDepth==2) 
    // already have cell center values so just read these
    ReadDepth(grid,myproc);
  else {
    for(n=0;n<grid->Nc;n++) {
      // get from from the initialization.c file
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

/*
 * Function: CreateCellGraph
 * Usage: CreateCellGraph(grid)
 * ------------------------------------------
 * Initialize the serial CSR xadj and adjncy 
 * data structures to import graph structure into
 * ParMETIS for graph partitioning (See ParMETIS
 * API for definitions and more info).
 *
 */
void CreateCellGraph(gridT *grid)
{
  int n, j, Nge;

  // compute the number of ghost edges summing up over the entire graph
  Nge = 0;
  // looping over each cell
  for(n=0;n<grid->Nc;n++) {
    // and over each face
    for(j=0;j<grid->nfaces[n];j++) {
      // check to see if there is a real neighbor (not ghost)
      if(grid->neigh[grid->maxfaces*n+j]!=-1) Nge++;
    }
  }
  grid->Nge = Nge/2;
  grid->adjncy = (int *)SunMalloc(Nge*sizeof(int),"CreateCellGraph");

  Nge=0;
  grid->xadj[0]=Nge;
  for(n=0;n<grid->Nc;n++) { 
    for(j=0;j<grid->nfaces[n];j++) {
      // for each cell's face where the neighbor is real (non-ghost)
      if(grid->neigh[grid->maxfaces*n+j]!=-1) {
        // populate adjacency matrix with non-ghost cells where 
        // each neighbor for a cell is added together
        grid->adjncy[Nge++]=grid->neigh[grid->maxfaces*n+j];
      }
      // record the cummulative number of cells adjacent so that
      // adjacency graph is split up by cells
    grid->xadj[n+1]=Nge;
    }
  }
  // output adjacency information
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

/*
 * Function: CreateEdgeGraph
 * Usage: CreateEdgeGraph(grid)
 * ------------------------------------------
 * Initialize the serial CSR xadj and adjncy 
 * data structures to import graph structure into
 * ParMETIS for graph partitioning (See ParMETIS
 * API for definitions and more info).
 *
 */
void CreateEdgeGraph(gridT *grid)
{
  int grad1, grad2, enei, n, j, Nge;

  Nge = 0;
  for(n=0;n<grid->Ne;n++) 
    //corrected part for get the grad for each edge, then get the eneigh
    grad1=grid->grad[n*2];
    grad2=grid->grad[n*2+1];
    // enei:the number of eneigh for each edge 
    for(j=0;j<enei;j++)
      if(grid->eneigh[2*(grid->maxfaces-1)*n+j]!=-1) Nge++;
  grid->Nge = Nge/2;
  grid->adjncy = (int *)SunMalloc(Nge*sizeof(int),"CreateEdgeGraph");

  Nge=0;
  grid->xadj[0]=Nge;
  for(n=0;n<grid->Ne;n++) { 
    //corrected part for get the grad for each edge, then get the eneigh
    grad1=grid->grad[n*2];
    grad2=grid->grad[n*2+1];
    // enei:the number of eneigh for each edge 
    for(j=0;j<enei;j++) 
      if(grid->eneigh[2*(grid->maxfaces-1)*n+j]!=-1) {
	grid->adjncy[Nge++]=grid->eneigh[2*(grid->maxfaces-1)*n+j];
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

/*
 * Function: CreateFaceArray
 * Usage: CreateFaceArray(grad, gradf, neigh, face, Nc, Ne)
 * ------------------------------------------
 * CreateFaceArray information, called for both the local and
 * global grids.
 *
 */
static inline void CreateFaceArray(int *grad, int *gradf, int *neigh, int *face, int *nfaces, int maxfaces, int Nc, int Ne) // added two elements
{
  int n, j, nf, nc, nc1, nc2, nb;

  // initialize all the faces to ghost faces in each cell over all faces of a cell
  for(n=0;n<Nc;n++)
    for(nf=0;nf<nfaces[n];nf++)
      face[n*maxfaces+nf]=-1; 
  //for(n=0;n<Ne;n++)
  //  for(j=0;j<2;j++)
  //    gradf[2*n+j]=-1;

  // over each edge
  for(n=0;n<Ne;n++) {
    nc1 = grad[2*n];
    nc2 = grad[2*n+1];
    // for a purely computational edge (no ghost cells adjacent to face)
    if(nc1!=-1 && nc2!=-1)
      for(nf=0;nf<nfaces[nc1];nf++) {
        // get a face for neighboring cell adjacent to edge for first cell touching edge
        nb = neigh[nc1*maxfaces+nf];
	//printf("nb is %d \n", nb);
        if(nb==nc2) {
          // if the neighbor and the adjacent cell are the same then we 
          // can store the actual face number and also the pointer to the face normal to 
          // the edge
          face[nc1*maxfaces+nf]=n;
          gradf[2*n] = nf;
          //printf("face = %d",face[nc1*maxfaces+nf]);
        }
        // get a face for neighboring cell adjacent to edge for second cell touching edge
        nb = neigh[nc2*maxfaces+nf];
        if(nb==nc1) {
          face[nc2*maxfaces+nf]=n;
          gradf[2*n+1] = nf;
        }
      }
  }

  // consider possibility that one of the faces is a boundary face (ghost cell adjacent = -1)
  for(n=0;n<Ne;n++) {
    for(j=0;j<2;j++) {
      // consider a cell corresponding to voronoi edge (nc)
      nc = grad[2*n+j];
      if(nc != -1)  { // if it is NOT a ghost cell
        // if particular edge is not a face of the cell (assigning remaining -1 faces not covered in
        // first case above in first set of loops)
        if(IsMember(n,&(face[nc*maxfaces]),nfaces[nc])==-1) {
          // consider each face of the cell to see if it is a boundary cell (ghost = -1)
          for(nf=0;nf<nfaces[nc];nf++) {
            if(face[nc*maxfaces+nf] == -1) {
              // assign actual value for the faces
              face[nc*maxfaces+nf]=n;
	      //printf("face = %d",face[nc*maxfaces+nf]);
              gradf[2*n+j]=nf;
              break;
            }
          }
        }
      }
    }
  }

  // - seems like the possibility that gradf wouldn't get assigned for 
  // a cell with ghost edge- but 
  // probably it is not needed/never called-- this is somewhat unclear
}

/*
 * Function: Reordergradf
 * usage: reorder gradf according to face after reordercellpoints
 * --------------------------------------------------------------
 * make gradf cosistent with faces
 *
 */
static inline void Reordergradf(int *face, int *grad, int *gradf, int *nfaces, int maxfaces, int Nc, int Ne){
  int i, j, k, nc;
  for(i=0;i<Ne;i++){
     for(j=0;j<2;j++){
       nc=grad[i*2+j];
       if(nc!=-1){
         for(k=0;k<nfaces[nc];k++){
            if(i==face[nc*maxfaces+k])
              gradf[2*i+j]=k;
              break;
         }
       }
     }
  }
}

/*
 * Function: ReorderCellPoints
 * Usage: ReorderCellPoints(face, cells,  Nc)
 * ------------------------------------------
 * Reorder the list of points in cells so that they follow a defined
 * winding structure, called for both the local and
 * global grids.
 *
 *///added two elements
static inline void ReorderCellPoints(int *face, int *edges, int *cells, int *nfaces, int maxfaces, int Nc) {
  // arrange face values so that they corresond to structure for quadratic interpolation
  // based on organization of cells (i.e. we will change cells (pointers to the nodes)
  // to be consistent with the existing edge order for interpolation

  int n1[2], n2[2], e3[maxfaces]; int nc, nf, e1, e2, sharednode, i;
  // over each cell
  for(nc=0; nc<Nc; nc++) {

    // consider each face in a cycle
    for(nf=0;nf<nfaces[nc]-1;nf++) {
      //  consider pairs of touching faces in cycle
      e1 = face[nc*maxfaces+nf];
      e2 = face[nc*maxfaces+nf+1];
      // nf+1 nodal point will be shared node

      // get each node
      n1[0] = edges[e1*NUMEDGECOLUMNS  ];
      n1[1] = edges[e1*NUMEDGECOLUMNS+1];
      n2[0] = edges[e2*NUMEDGECOLUMNS  ];
      n2[1] = edges[e2*NUMEDGECOLUMNS+1];
      // compare nodes to find first value for cells
      sharednode = SharedListValue(n1, n2, 2);
      //if edge is not clockwise or counter-clockwise 
      //exchang with the next one
      if(sharednode==-1 && maxfaces>3){
         for(i=0;i<(nfaces[nc]-1-nf);i++)
           e3[i]=face[nc*maxfaces+nf+1+i];
         for(i=0;i<(nfaces[nc]-2-nf);i++)
           face[nc*maxfaces+nf+1+i]=e3[i+1];
         face[nc*maxfaces+nfaces[nc]-1]=e3[0];
         nf=nf-1;
      } else if(sharednode==-1 && maxfaces==3) {
        printf("Error as sharednode in %s not found %d\n", "ReorderCellPoints",maxfaces);     
      }
      // check for match
      if(sharednode != -1) 
        cells[maxfaces*nc + (nf+1)]  = sharednode;

    }

    // consider last pairs of touching faces in cycle
    e1 = face[nc*maxfaces+nfaces[nc]-1];
    e2 = face[nc*maxfaces  ];
    // get each node
    n1[0] = edges[e1*NUMEDGECOLUMNS  ];
    n1[1] = edges[e1*NUMEDGECOLUMNS+1];
    n2[0] = edges[e2*NUMEDGECOLUMNS  ];
    n2[1] = edges[e2*NUMEDGECOLUMNS+1];

    // 0th (1st) nodal point will be shared node
    sharednode = SharedListValue(n1, n2, 2);

    // check for match
    if(sharednode != -1)
      cells[maxfaces*nc  ]  = sharednode;
    else
      printf("Error as sharednode in %s not found\n", "ReorderCellPoints");  

  }

  //PJw now we have sorted our list of cell vertexes to correspond with the 
  // structure in the face (edge) array for the cell
}

/*
 * Function: CreateNodeArray
 * Usage: CreateNodeArray(grad, gradf, neigh, face, Nc, Ne)
 * ------------------------------------------
 * CreatNodeArray information so that all the nodal and cell 
 * neighboors for a nodal point can be easily accessed, 
 * called for both the global and local grids.
 *
 */
static inline void CreateNodeArray(gridT *grid, int Np, int Ne, int Nc, int myproc)
{
  //  printf("Nodal array stub with Np=%d Ne=%d\n", Np, Ne);

  if(myproc==0 && VERBOSE>2)  printf("\t\t\tCompute nodal node neighbors...\n");
  Tic();

  // these numbers reflect the assumption that a node will never be surrounded by
  // cells with angles less than or equal 30 deg.
  const int maxnodeneighs = 24; 
  const int maxedgeneighs = 24; //since the edges are intentionally counted twice (once per cell)
  const int maxcellneighs = 12;
  int in, ie, ic, nf, np1, np2, inn, inode;

  // assume that at most there will only be a maximum of edge neighbors to a node
  //  int tempppneighs[Np][maxnodeneighs]; int tempnumppneighs;
  //  int temppeneighs[Np][maxedgeneighs];
  //  int temppcneighs[Np][maxcellneighs];
  int **tempppneighs, tempnumppneighs, **temppeneighs, **temppcneighs;
  tempppneighs = (int **)malloc(Np*sizeof(int *));
  temppeneighs = (int **)malloc(Np*sizeof(int *));
  temppcneighs = (int **)malloc(Np*sizeof(int *));
  for(in=0;in<Np;in++) {
    tempppneighs[in] = (int *)malloc(maxnodeneighs*sizeof(int));
    temppeneighs[in] = (int *)malloc(maxedgeneighs*sizeof(int));
    temppcneighs[in] = (int *)malloc(maxcellneighs*sizeof(int));
  }
  // initialize the number of neighbors for each point
  for (in = 0; in < Np; in++) { 
    grid->numppneighs[in] = 0;
    grid->numpeneighs[in] = 0;
    grid->numpcneighs[in] = 0;
  }

  /* get cell neighbors */


  // loop over each cell and assign it as a neighbor to each of its three points
  for(ic = 0; ic < Nc; ic++) {
    // get each of its nodes and assign to the respective point
    for(in =0; in <grid->nfaces[ic]; in++) {
      inode = grid->cells[ic*grid->maxfaces+ in];
      temppcneighs[inode][grid->numpcneighs[inode]++] = ic;
    }
  }

  // transfer temporary values to datastructure
  for(in = 0; in < Np; in++) {
    // transfer values
    //    // take care of TRIANGLE's bogus points and just make them zero anyway
    if(grid->numpcneighs[in] == 0) {
      //        printf("Warning!!:  Fake point %d found\n", in);
      grid->pcneighs[in] = (int*)SunMalloc(sizeof(int),"CreateNodeArray");
      grid->pcneighs[in][0] = -1;
    }
    else {
      //        printf("Good!!:  Non-Fake point %d found\n", in);
      // allocate memory
      grid->pcneighs[in] = (int*)SunMalloc(grid->numpcneighs[in]*sizeof(int),"CreateNodeArray");
    }
    for(ic = 0; ic < grid->numpcneighs[in]; ic++) {
      grid->pcneighs[in][ic] = temppcneighs[in][ic];
    }
  }

  /* get edge neighbors (correlated to cell neighbors) */
  /* hence, size 2*numpcneighs[in]                     */
  /* as well as node neighbors                         */

  // from cell neighboors, get edge neighbors
  for(in = 0; in < Np; in++) {
    // over each cell neighbor
    for(ic = 0; ic < grid->numpcneighs[in]; ic++) {
      // get cell edges (faces)
      for(nf = 0; nf < grid->nfaces[grid->pcneighs[in][ic]]; nf++) {
        // get particlar face
        ie = grid->face[grid->maxfaces*grid->pcneighs[in][ic] + nf];
        // get face nodes
        np1 = grid->edges[NUMEDGECOLUMNS*ie  ];
        np2 = grid->edges[NUMEDGECOLUMNS*ie+1];
        // check to see if a node is shared with in
        if(in == np1) {
          // if so this is a shared edge- add it!
          temppeneighs[in][grid->numpeneighs[in]++] = ie;
          // if so this is a shared node - add it!
          tempppneighs[in][grid->numppneighs[in]++] = np2;
        }
        if(in == np2) {
          // if so this is a shared edge- add it!
          temppeneighs[in][grid->numpeneighs[in]++] = ie;
          // if so this is a shared node - add it!
          tempppneighs[in][grid->numppneighs[in]++] = np1;
        }
      }
    }
  }
  // transfer termporary values to datastructures
  for(in = 0; in < Np; in++) {
    // allocate memory 
    grid->peneighs[in] = (int*)SunMalloc(grid->numpeneighs[in]*sizeof(int),"CreateNodeArray");
    grid->ppneighs[in] = (int*)SunMalloc(grid->numppneighs[in]*sizeof(int),"CreateNodeArray");
    // transfer values
    for(ie = 0; ie < grid->numpeneighs[in]; ie++) {
      grid->peneighs[in][ie] = temppeneighs[in][ie];
    }
    for(ie = 0; ie < grid->numppneighs[in]; ie++) {
      grid->ppneighs[in][ie] = tempppneighs[in][ie];
    }
  }

  /* simplify nodal neighbors (for no redundancy) and transfer values */

  // for each node
  for(in = 0; in < Np; in++) {
    // allocate memory (note that we are going to 
    //1/2 the size because of redundancy reduction)
    //  printf("grid->numppneighs[in]/2=%d\n", grid->numppneighs[in]);
    //    grid->ppneighs[in] = (int*)SunMalloc(grid->numppneighs[in]/2*sizeof(int),"CreateNodeArray");
    // get list of temp list
    tempnumppneighs = grid->numppneighs[in];
    // reset counter for actual list
    grid->numppneighs[in] = 0;
    // for each node in temp list
    for(inn = 0; inn < tempnumppneighs; inn++) {
      // determine if it has already been counted
      inode = tempppneighs[in][inn];
      // if it isn't a member add it
      if(IsMember(inode, tempppneighs[in], grid->numppneighs[in]) == -1) {
        // reassign to the same temp array for sorting
        tempppneighs[in][grid->numppneighs[in]++] = inode;
      }
    }
    // allocate memory for limited nodal neighbors
    grid->ppneighs[in] = (int*)SunMalloc(grid->numppneighs[in]*sizeof(int),"CreateNodeArray");
    // assign values from temp
    for(inn = 0; inn < grid->numppneighs[in]; inn++) {
      grid->ppneighs[in][inn] = tempppneighs[in][inn];
    }
  }

  //  printf("Finished computing nodal array information\n");
  //printf("Freeing neighs...\n");
  for(in=0;in<Np;in++) 
    SunFree(tempppneighs[in],Np*sizeof(int *),"CreateNodeArray");
  //printf("Freeing edges...\n");
  for(in=0;in<Np;in++) 
    SunFree(temppeneighs[in],Np*sizeof(int *),"CreateNodeArray");
  printf("Freeing cells...\n");
  for(in=0;in<Np;in++) 
    SunFree(temppcneighs[in],Np*sizeof(int *),"CreateNodeArray");

  printf("Points\n");
  SunFree(tempppneighs,maxnodeneighs*sizeof(int),"CreateNodeArray");
  printf("Edges\n");
  SunFree(temppeneighs,maxedgeneighs*sizeof(int),"CreateNodeArray");
  printf("Cells\n");
  SunFree(temppcneighs,maxcellneighs*sizeof(int),"CreateNodeArray");
  printf("Done!\n");
}

/*
 * Function: CreateNormalArray
 * Usage: CreateNormalArray(grad, face, normal, Nc)
 * ------------------------------------------
 * CreatNormalArray information, called for both the local and
 * global grids
 *
 */
static inline void CreateNormalArray(int *grad, int *face, int *normal, int *nfaces, int maxfaces, int Nc)
{
  int n, nf;
  // for each cell
  for(n=0;n<Nc;n++){
    // over each face on the cell
    for(nf=0;nf<nfaces[n];nf++) {
      // for each face arbitrarily set normal to -1 if the first grad 
      // is the same cell as the current cell, otherwise 1
      if(n==grad[2*face[maxfaces*n+nf]])
        normal[n*maxfaces+nf]=-1;
      else
        normal[n*maxfaces+nf]=1;
    }
  }
}

/*
 * Function: Connectivity
 * Usage: Connectivity(grid, myproc)
 * ------------------------------------------
 * Determine how grid is connected.
 *
 */
void Connectivity(gridT *grid, int myproc)
{
  int n, nf, ng, ne;

  /* Create the face array, which contains indexes to the edges of
     each cell */
  if(myproc==0 && VERBOSE>2) printf("\tCreating face array and normals and nodal neighboors\n");
  // get pointers to the faces of each cell (grid->face[NFACES*Nc]) as well as 
  // pointers to the face number of given cell corresponding to given edge 
  // (grid->gradf[2*Ne])
  CreateFaceArray(grid->grad,grid->gradf,grid->neigh,grid->face,grid->nfaces,grid->maxfaces,grid->Nc,grid->Ne);

  // printf("faces is %d", grid->face[99]);
  // reorder the cell points so that they have structure for use in interpolation
  // now reordercellpoints may reorder face, so gradf should be reordered
  ReorderCellPoints(grid->face, grid->edges, grid->cells, grid->nfaces, grid->maxfaces, grid->Nc);

  // reorder gradf
  Reordergradf(grid->face, grid->grad, grid->gradf, grid->nfaces, grid->maxfaces, grid->Nc, grid->Ne);
  // create dot product of unique normal with outward normal
  // (grid->normal[NFACES*Ne]) (always +/- 1)
  CreateNormalArray(grid->grad,grid->face,grid->normal,grid->nfaces,grid->maxfaces,grid->Nc);

  // create nodal neighboor values for the global grid

  if(myproc==0 && VERBOSE>2) printf("\tCreating global node array...\n");
  Tic();
  CreateNodeArray(grid, grid->Np, grid->Ne, grid->Nc, myproc);
  if(myproc==0 && VERBOSE>2) printf("\t... time used is %f\n", Toc());
  
  /* Create the edge connectivity array eneigh, which points to the 2(NFACES-1)
     neighbors that comprise the edges for the momentum control volume */
  if(myproc==0 && VERBOSE>2) printf("\tCreating edge connectivity array...\n");
  // get the edges comprizing the momentum control volume grid->eneigh[4*Ne]
  CreateMomentumCV(grid);
 
}

/*
 * Function: CreateMomentumCV
 * Usage: CreateMomentumCV(grid)
 * ------------------------------------------
 * Determine the edges for the momentum control volume, 
 * grid->eneigh[2(maxfaces-1)
 *
 */
static inline void CreateMomentumCV(gridT *grid)
{
  int n, ne, ng, nf;
  // loop over all the edges
  for(n=0; n < grid->Ne; n++) {
    // for each edge of the momentum control volume
    for(ne = 0; ne < 2*(grid->maxfaces-1); ne++) {
      // initialize the edge to be ghost for the CV
      grid->eneigh[2*(grid->maxfaces-1)*n+ne] = -1;
    }
    ne = 0;
    // for each neighbor
    for(ng=0; ng < 2; ng++) {
      // if the neighbor is not a ghost
      if(grid->grad[2*n+ng]!=-1) {
        // for each face
        for(nf=0; nf < grid->nfaces[grid->grad[2*n+ng]]; nf++) {
          // determine if the face for each neighbor 
          // isn't the current edge
          if(grid->face[grid->maxfaces*grid->grad[2*n+ng]+nf]!=n) {
            // if it isn't, set the momentum control volume
            grid->eneigh[2*(grid->maxfaces-1)*n+ne++]=grid->face[grid->maxfaces*grid->grad[2*n+ng]+nf];
          }
        }
      }
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
 *   2 if this is an inner interproc boundary cell.
 *   1 if this is a non-interproc boundary cell 
 *     (not a ghost cell), i.e. free-surface specified 
 *      containing an edge with mark=3
 *   0 otherwise, (e.g., computational cells)
 *  -1 for error.
 *
 * If it is both a ghost cell and a boundary cell or
 * both a 
 *
 */
// may have to be changed to incorporate full row of 
// interprocessor cells for the interpolation, generalized
// halo region
inline int IsBoundaryCell(int mgptr, gridT *maingrid, int myproc)
{
  int nf, nei;

//  if(myproc==0)  printf("IsBoundaryCell j=%d\n",mgptr);
  // check for array out-of-bounds error
  if(mgptr >= maingrid->Nc || mgptr < 0) {
    printf("Error in IsBoundaryCell: index out of bounds!\n");
    return -1;
  }
  // if the partition of the pointer is on the same processor
  if(maingrid->part[mgptr]==myproc) {
    for(nf=0;nf<maingrid->nfaces[mgptr];nf++) {
      // we have a type 3 boundary condition on same proc
      // this indicates you have a non-interprocessor boundary
      // cell
      if(maingrid->mark[maingrid->face[mgptr*maingrid->maxfaces+nf]]==3)
        return 1;
    }
    // search through the halo and return the result
    // if we found an inner interproc boundary cell
    if(IsBoundaryCellByHalo(mgptr, maingrid, myproc, g_halolist, g_halolistsize)!=-1)
      return 2;
    // for the case where we have no type 3 BC on the cell or
    // cell neighbors that are on a different processor
    // this is probably for debugging
    return 0;
  } else
    //  the cell is a ghost cell to the local processor (not 
    //  on the set of cells on the local processor)
    return 3;
}

/*
 * Function: IsBoundaryCellByHalo
 * Usage: IsBoundaryCellByHalo(i,maingrid,myproc);
 * -----------------------------------------
 * Determines whether or not the global cell index
 * i is a boundary cell by searching through all
 * it's halo points. 
 * Returns:
 *   2 if this is an inner interproc boundary cell.
 *  -1 for no hit
 *
 */
static int IsBoundaryCellByHalo(int mgptr, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize)
{
  // if we don't have a good list, exit
  if(listsize < 1) {
    printf("Error!!: List for boundary neighbors needs to be larger!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }

  // select which function to call off the first item in the list
  switch(list[0]) {
    case EDGE:
      // original selection by edge
      return IsBoundaryCellByEdge(mgptr, maingrid, myproc, list, listsize-1);
      break;
    case NODE:
      // new selection based on nodal neighbors (more inclusive)
      return IsBoundaryCellByNode(mgptr, maingrid, myproc, list, listsize-1);
      break;
    default:
      printf("Error in IsBoundaryCellByHalo!!\n");
      break;
  }
}
/*
 * Function: IsBoundaryCellByEdge
 * Usage: IsBoundaryCellByEdge(i,maingrid,myproc);
 * -----------------------------------------
 * Determines whether or not the global cell index
 * i is a boundary cell. 
 * Returns:
 *   2 if this is an inner interproc boundary cell.
 *  -1 for no hit
 *
 */
static int IsBoundaryCellByEdge(int cell, gridT *maingrid, int myproc, 
    boundaryselection *list, int listsize)
{
  int nf, nei;

//  printf("IsBoundaryCellByEdge\n");
  // loop over all the faces of cell and get the neighbors
  for(nf=0;nf<maingrid->nfaces[cell];nf++) {
    nei = maingrid->neigh[cell*maingrid->maxfaces+nf];
    // if any neighbor is not a ghost cell check to see
    // if it is on the processor too, if not it is an 
    // interboundary processor cell
    if(nei != -1) {
      if(maingrid->part[nei] != myproc)
        return 2;
      // if we didn't get the hit we need to pass this downstream
      // to be checked by the next part of the halo, so pass the 
      // current cell to be evaluated (and how to evaluate it 
      // via the next entry on the halo list).
      /* now process other items on the list, if any */

      // if there are more items to be considered
      if(listsize > 0) {

        // select which function to call off the first item in the list
        switch(list[1]) {
          case EDGE:
            // original selection by edge
            // check to see if we have a hit, if so we can quit and return 2
            if(IsBoundaryCellByEdge(nei, maingrid, myproc, &list[1], listsize-1) !=-1 )
              return 2;
            break;
          case NODE:
            // new selection based on nodal neighbors (more inclusive)
            // check to see if we have a hit, if so we can quit and return 2
            if(IsBoundaryCellByNode(nei, maingrid, myproc, &list[1], listsize-1) != -1)
              return 2;
            break;
          default:
            printf("Error in BoundaryTopologyByEdge!!\n");
            break;
        }
      }
    }
  }
  // no hit
  return -1;
}

/*
 * Function: IsBoundaryCellByNode
 * Usage: IsBoundaryCellByNode(cell,maingrid,myproc);
 * -----------------------------------------
 * Determines whether or not the global cell index
 * i is a boundary cell. 
 * Returns:
 *   2 if this is an inner interproc boundary cell.
 *  -1 for no hit
 *
 */
static int IsBoundaryCellByNode(int cell, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize)
{
  int nf, nei, in, nn;

  //  printf("IsBoundaryCellByNode\n");
  // loop over all the nodes of cell and get the neighbors
  for(nf=0;nf<maingrid->nfaces[cell];nf++) {
    in = maingrid->cells[cell*maingrid->maxfaces+nf];
    // loop over all the nodal cell neighbors
    for(nn = 0; nn < maingrid->numpcneighs[in]; nn++) {
      // get the nodal cell neighbors
      nei = maingrid->pcneighs[in][nn];
      // if any neighbor is not a ghost cell check to see
      // if it is on the processor too, if not it is an 
      // interboundary processor cell
      if(nei != -1) {
        if(maingrid->part[nei] != myproc)
          return 2;
        // if we didn't get the hit we need to pass this downstream
        // to be checked by the next part of the halo, so pass the 
        // current cell to be evaluated (and how to evaluate it 
        // via the next entry on the halo list).
        /* now process other items on the list, if any */

        // if there are more items to be considered
        if(listsize > 0) {

          // select which function to call off the first item in the list
          switch(list[1]) {
            case EDGE:
              // original selection by edge
              // check to see if we have a hit, if so we can quit and return 2
              if(IsBoundaryCellByEdge(nei, maingrid, myproc, &list[1], listsize-1)!=-1)
                return 2;
              break;
            case NODE:
              // new selection based on nodal neighbors (more inclusive)
              // check to see if we have a hit, if so we can quit and return 2
              if(IsBoundaryCellByNode(nei, maingrid, myproc, &list[1], listsize-1)!=-1)
                return 2;
              break;
            default:
              printf("Error in BoundaryTopologyByNode!!\n");
              break;
          }
        }
      }
    }
  }
  // no hit
  return -1;
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
    
    ofile = MPI_FOpen(CELLSFILE,"w","OutputData",myproc);
    for(n=0;n<maingrid->Nc;n++) {
      // add a column in cells.dat
      fprintf(ofile,"%d %f %f ",maingrid->nfaces[n],maingrid->xv[n],maingrid->yv[n]);
      // corrected part
      fprintf(ofile,"%f %f ",maingrid->xv[n],maingrid->yv[n]);
      for(nf=0;nf<maingrid->nfaces[n];nf++)
	fprintf(ofile,"%d ",maingrid->cells[n*maingrid->maxfaces+nf]);
      for(nf=0;nf<maingrid->nfaces[n];nf++) 
	fprintf(ofile,"%d ",maingrid->neigh[n*maingrid->maxfaces+nf]);
      fprintf(ofile,"\n");
    }
    fclose(ofile);
  }
  
  sprintf(str,"%s.%d",CELLSFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str); 
  ofile = MPI_FOpen(str,"w","OutputData",myproc);

  for(j=0;j<grid->Nc;j++) {
    fprintf(ofile,"%d %e %e ",grid->nfaces[j],grid->xv[j],grid->yv[j]);
    for(nf=0;nf<grid->nfaces[j];nf++)
      fprintf(ofile,"%d ",grid->cells[j*grid->maxfaces+nf]);
    for(nf=0;nf<grid->nfaces[j];nf++)
      fprintf(ofile,"%d ",grid->neigh[j*grid->maxfaces+nf]);
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
    //corrected part //add another column for celldata.dat for nfaces 
    fprintf(ofile,"%d %e %e %e %e %d ",grid->nfaces[n],grid->xv[n],grid->yv[n],grid->Ac[n],grid->dv[n],grid->Nk[n]);
    for(nf=0;nf<grid->nfaces[n];nf++)
     //printf("grid->face is %d\n",grid->face[grid->maxfaces*n+nf]);
      fprintf(ofile,"%d ",grid->face[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++)
      fprintf(ofile,"%d ",grid->neigh[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++) 
      fprintf(ofile,"%d ",grid->normal[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++) 
      fprintf(ofile,"%e ",grid->def[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++) 
      fprintf(ofile,"%d ",grid->cells[grid->maxfaces*n+nf]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);

  ofile = MPI_FOpen(str,"w","OutputData",myproc);
  for(n=0;n<Ne;n++) {
    fprintf(ofile,"%e %e %e %e %e %e %d %d %d %d %d %d %d %d %d\n",
	    grid->df[n],grid->dg[n],grid->n1[n],grid->n2[n],grid->xe[n],grid->ye[n],
	    grid->Nke[n],grid->Nkc[n],grid->grad[2*n],grid->grad[2*n+1],
	    grid->gradf[2*n],grid->gradf[2*n+1],grid->mark[n], 
      grid->edges[n*NUMEDGECOLUMNS], grid->edges[n*NUMEDGECOLUMNS+1]);
  }
  fclose(ofile);

  // added nodal structure information
  sprintf(str,"%s.%d",NODEFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);

  ofile = MPI_FOpen(str,"w","OutputData",myproc);
//  // output the number of points to the domain
//  fprintf(ofile,"%d\n", grid->Np);
  // over each nodal point
  for(n=0;n<Np;n++) { // Np needs to be local, global right now
    // NOTE THAT n SHOULD BE SUBSTITUTED FOR THE GLOBAL ID FOR THE NODE 
    // FOR THE PARALLLELIZED IMPLEMENTATION probably via local to global
    // nodal pointer
    
    // print the nodal data structure information (will need conversion pointer)
    fprintf(ofile,"%d %e %e ", n, maingrid->xp[n], maingrid->yp[n]);

    // print the number of nodal neighbors
    fprintf(ofile,"%d ", grid->numppneighs[n]);
//    printf("Np = %d numppneighs= %d ppneighbors = { ", n, grid->numppneighs[n]);
    // now loop over all the nodal neighbors for each node
    for(neigh = 0; neigh < grid->numppneighs[n]; neigh++ ) {
      fprintf(ofile,"%d ", grid->ppneighs[n][neigh]);
//      printf("%d ", grid->ppneighs[n][neigh]);
    }
//    printf("}");

    // now list the edge neighbors
    fprintf(ofile,"%d ", grid->numpeneighs[n]);
    // now loop over all the nodal neighbors for each node
    for(neigh = 0; neigh < grid->numpeneighs[n]; neigh++ ) {
      fprintf(ofile,"%d ", grid->peneighs[n][neigh]);
//      printf("%d ", grid->ppneighs[n][neigh]);
    }
//    printf("}");

    // now list all the cell neighboors
//    printf(" numpcneighs= %d pcneighbors = { ", grid->numpcneighs[n]);
    fprintf(ofile,"%d ", grid->numpcneighs[n]);
    for(neigh = 0; neigh < grid->numpcneighs[n]; neigh++) {
      // now loop over all the cell neighbors
      fprintf(ofile,"%d ", grid->pcneighs[n][neigh]);
//      printf("%d ", grid->pcneighs[n][neigh]);
    }
    // now add the maximum number of layers for the node
    fprintf(ofile,"%d ", grid->Nkp[n]);
    // now add the total area in each layer to the end
    for(neigh = 0; neigh < grid->Nkp[n]; neigh++){
      fprintf(ofile,"%e ", grid->Actotal[n][neigh]);
    }
    // now end the line since we've outputed all the data for 
    // nodal neighbor information
    fprintf(ofile,"\n");
//    printf("}\n");
  }
  fclose(ofile);

//  // added nodal structure information -- DEBUG FOR LIST OF POINTS CHECK
//  sprintf(str,"%s.%d",NODEFILE,myproc);
//  if(VERBOSE>2) printf("Outputting %s...\n",str);
//
//  ofile = MPI_FOpen(str,"w","OutputData",myproc);
//  // over each nodal point
//  for(n=0;n<grid->Np;n++) {
//    // print the nodal point
//    fprintf(ofile,"%d ", grid->localpoints[n]);
//    // now end the line since we've outputed all the data for 
//    // nodal neighbor information
//    fprintf(ofile,"\n");
//  }
//  fclose(ofile);

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
  int neigh, n, nf, np, ne, nc, Nkmax;
  int Np;
  char str[BUFFERLENGTH];
  FILE *ifile;

  ReadFileNames(myproc);

  InitLocalGrid(grid);

  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  (*grid)->Nc = MPI_GetSize(str,"ReadGrid",myproc);
  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  (*grid)->Ne = MPI_GetSize(str,"ReadGrid",myproc);
  sprintf(str,"%s.%d",NODEFILE,myproc);
  (*grid)->Np = MPI_GetSize(str,"ReadGrid",myproc);
 
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
  (*grid)->nfaces = (int *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");	
  (*grid)->xv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->yv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->dv = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->Ac = (REAL *)SunMalloc((*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->Nk = (int *)SunMalloc((*grid)->Nc*sizeof(int),"ReadGrid");

  //read maxFaces and nfaces from celldata.dat first
  (*grid)->maxfaces=(int)MPI_GetValue(DATAFILE,"maxFaces","ReadGrid",0);
  (*grid)->neigh = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->face = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->normal = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->def = (REAL *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->cells = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(REAL),"ReadGrid");

  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Reading %s...\n",str);

  ifile = MPI_FOpen(str,"r","ReadGrid",myproc);
  for(n=0;n<(*grid)->Nc;n++) {
    (*grid)->nfaces[n]=(int)getfield(ifile,str);
    if((*grid)->nfaces[n]>(*grid)->maxfaces){
    printf("Max number of edges is %d of No. %d in Cell.dat is more than maxFaces %d in suntans.dat!",(*grid)->nfaces[n],n,(*grid)->maxfaces);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    }
    (*grid)->xv[n]=getfield(ifile,str);
    (*grid)->yv[n]=getfield(ifile,str);
    (*grid)->Ac[n]=getfield(ifile,str);
    (*grid)->dv[n]=getfield(ifile,str);
    (*grid)->Nk[n]=(int)getfield(ifile,str);
    for(nf=0;nf<(*grid)->nfaces[n];nf++)
      (*grid)->face[(*grid)->maxfaces*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<(*grid)->nfaces[n];nf++)
      (*grid)->neigh[(*grid)->maxfaces*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<(*grid)->nfaces[n];nf++)
      (*grid)->normal[(*grid)->maxfaces*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<(*grid)->nfaces[n];nf++)
      (*grid)->def[(*grid)->maxfaces*n+nf]=getfield(ifile,str);
    for(nf=0;nf<(*grid)->nfaces[n];nf++)
      (*grid)->cells[(*grid)->maxfaces*n+nf]=(int)getfield(ifile,str);
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
  (*grid)->edges = (int *)SunMalloc((*grid)->Ne*NUMEDGECOLUMNS*sizeof(int),"ReadGrid");

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
      (*grid)->edges[NUMEDGECOLUMNS*n] = (int)getfield(ifile,str);
      (*grid)->edges[NUMEDGECOLUMNS*n+1] = (int)getfield(ifile,str);
  }
  fclose(ifile);

  /* 
   * Now read in node data
   *
   */

  sprintf(str,"%s.%d",NODEFILE,myproc);
  if(VERBOSE>2) printf("Reading %s...\n",str);
  ifile = MPI_FOpen(str,"r","ReadGrid",myproc);

  // allocate memory
  Np = (*grid)->Np;
  (*grid)->xp = (REAL *)SunMalloc(Np*sizeof(REAL),"ReadGrid");
  (*grid)->yp = (REAL *)SunMalloc(Np*sizeof(REAL),"ReadGrid");
  (*grid)->localtoglobalpoints = (int*)SunMalloc(Np*sizeof(int),"ReadGrid");
  (*grid)->numppneighs = (int*)SunMalloc(Np*sizeof(int),"ReadGrid");
  (*grid)->ppneighs = (int**)SunMalloc(Np*sizeof(int*),"ReadGrid");
  (*grid)->numpeneighs = (int*)SunMalloc(Np*sizeof(int),"ReadGrid");
  (*grid)->peneighs = (int**)SunMalloc(Np*sizeof(int*),"ReadGrid");
  (*grid)->numpcneighs = (int*)SunMalloc(Np*sizeof(int),"ReadGrid");
  (*grid)->pcneighs = (int**)SunMalloc(Np*sizeof(int*),"ReadGrid");
  (*grid)->Nkp= (int*)SunMalloc(Np*sizeof(int),"ReadGrid");
  (*grid)->Actotal = (REAL **)SunMalloc(Np*sizeof(REAL*),"ReadGrid");

  // for each line in the nodal file
  for(n=0; n < Np; n++) {
    // get the pointer to the global value
    (*grid)->localtoglobalpoints[n] = (int)getfield(ifile,str);
    // get the coordinates of the point
    (*grid)->xp[n] = getfield(ifile,str);
    (*grid)->yp[n] = getfield(ifile,str);
    // now get the nodal neighbors
    // now read in the number of nodal neighbors
    (*grid)->numppneighs[n] = (int)getfield(ifile,str);
    // now allocate the memory
    (*grid)->ppneighs[n] = (int*)SunMalloc((*grid)->numppneighs[n]*sizeof(int),"ReadGrid");
    // now loop over all these and get the values for the nodal nodal neighbors
    for(np=0; np < (*grid)->numppneighs[n]; np++) {
      (*grid)->ppneighs[n][np] = (int)getfield(ifile,str);
    }
    // now get the number of edge neighbors
    // now read in the number of edge neighbors
    (*grid)->numpeneighs[n] = (int)getfield(ifile,str);
    // now allocate the memory
    (*grid)->peneighs[n] = (int*)SunMalloc((*grid)->numpeneighs[n]*sizeof(int),"ReadGrid");
    // now loop over all these and get the values for the nodal edge neighbors
    for(np=0; np < (*grid)->numpeneighs[n]; np++) {
      (*grid)->peneighs[n][np] = (int)getfield(ifile,str);
    }
    // now get the number of cell neighbors
    // now read in the number of cell neighbors
    (*grid)->numpcneighs[n] = (int)getfield(ifile,str);
    // now allocate the memory
    (*grid)->pcneighs[n] = (int*)SunMalloc((*grid)->numpcneighs[n]*sizeof(int),"ReadGrid");
    // now loop over all these and get the values for the nodal cell neighbors
    for(np=0; np < (*grid)->numpcneighs[n]; np++) {
      (*grid)->pcneighs[n][np] = (int)getfield(ifile,str);
    }
    // now get the number of vertical layers Nkp
    (*grid)->Nkp[n] = (int)getfield(ifile,str);
    // allocate memory for each layer
    (*grid)->Actotal[n] = (REAL *)SunMalloc((*grid)->Nkp[n]*sizeof(REAL),"ReadGrid");
    // now get the total area of shared cells at this point  (Actotal)
    for(np=0; np < (*grid)->Nkp[n]; np++) {
      (*grid)->Actotal[n][np] = (REAL)getfield(ifile,str);
    }
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
  SunFree(grid->cells,grid->maxfaces*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->neigh,grid->maxfaces*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->eneigh,2*(grid->maxfaces-1)*grid->Ne*sizeof(int),"FreeGrid");
  SunFree(grid->part,grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->order,grid->Ne*sizeof(int),"FreeGrid");
  SunFree(grid->Nk,grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->nfaces,grid->Nc*sizeof(int),"FreeGrid");

  // NEED TO ADD NEW VARIABLES TO FREE UP MEMORY
  SunFree(grid->dz,grid->Nkmax*sizeof(REAL),"FreeGrid");
  SunFree(grid->face,grid->maxfaces*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->grad,2*grid->Ne*sizeof(int),"FreeGrid");
  SunFree(grid->mark,grid->Ne*sizeof(int),"FreeGrid");

  SunFree(grid->normal,grid->maxfaces*grid->Nc*sizeof(int),"FreeGrid");
  SunFree(grid->xadj,(grid->Nc+1)*sizeof(int),"FreeGrid");
  SunFree(grid->vtxdist,VTXDISTMAX,"FreeGrid");

  for(proc=0;proc<numprocs;proc++) 
    SunFree(grid->neighs[proc],numprocs*sizeof(int),"FreeGrid");
  SunFree(grid->neighs,numprocs*sizeof(int *),"FreeGrid");
  SunFree(grid->numneighs,numprocs*sizeof(int),"FreeGrid");

  SunFree(grid,sizeof(gridT),"FreeGrid");
}

/*
 * Function: VertGrid
 * Usage:  VertGrid(maingrid, localgrid, comm)
 * -------------------------
 * Functon to setup vertical grid (z-levels).
 */
static void VertGrid(gridT *maingrid, gridT **localgrid, MPI_Comm comm)
{
  int i, j, k, ne, npc, myproc, numprocs, vertgridcorrect, stairstep, maxNk;
  REAL dz0, dmin, dmax, dmaxtest;

  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myproc);

  maingrid->Nkmax = MPI_GetValue(DATAFILE,"Nkmax","VertGrid",myproc);
  vertgridcorrect = MPI_GetValue(DATAFILE,"vertgridcorrect","VertGrid",myproc);
  stairstep = MPI_GetValue(DATAFILE,"stairstep","VertGrid",myproc);

  // compute dmin and dmax for all the points on the main grid
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

    // for the case where shallowest point cannot be resolved
    if(dmin < dz0 && vertgridcorrect) {
      dz0 = dmin;
      if(WARNING) {
        printf("Warning!\n");
        printf("Not enough vertical grid points to resolve the most shallow cell!\n");
        printf("Changing maximum number of vertical cells from %d to %d\n",
            maingrid->Nkmax,(int)(dmax/dz0));
      }
      // fix this so that we can resolve the shallowest part of the domain
      maingrid->Nkmax = 1+(int)(dmax/dz0);
    }
    // allocate memory for size of vertical cells
    maingrid->dz = (REAL *)SunMalloc(maingrid->Nkmax*sizeof(REAL),"VertGrid");

    // get vertical grid spacing from initialization file
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
  // communicate Nkmax to all the processors
  MPI_Bcast(&(maingrid->Nkmax),1,MPI_INT,0,comm);
  // allocate memory for size of vertical cells
  if(myproc!=0)
    maingrid->dz = (REAL *)SunMalloc(maingrid->Nkmax*sizeof(REAL),"VertGrid");
  // transfer dz, dz0, and dmax to other processors
  MPI_Bcast((void *)maingrid->dz,maingrid->Nkmax,MPI_DOUBLE,0,comm);
  MPI_Bcast(&dz0,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(&dmax,1,MPI_DOUBLE,0,comm);

  (*localgrid)->Nkmax = maingrid->Nkmax;
  (*localgrid)->dz = (REAL *)SunMalloc((*localgrid)->Nkmax*sizeof(REAL),"VertGrid");
  // transfer grid spacing from global grid on proc to local grid on proc
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
  (*localgrid)->Nkp = (int *)SunMalloc((*localgrid)->Np*sizeof(int),"VertGrid");

  // for each local cell set number of layers (3D problem)
  if((*localgrid)->Nkmax>1) {
    for(i=0;i<(*localgrid)->Nc;i++) 
      if((*localgrid)->dv[i]==dmax) 
        (*localgrid)->Nk[i] = (*localgrid)->Nkmax;
      else {
        dmaxtest=0;
        for(k=0;k<(*localgrid)->Nkmax;k++) {
          // add up the cumulative cell depth until it exceeds the actual depth
          dmaxtest+=(*localgrid)->dz[k];
          (*localgrid)->Nk[i] = k+1;
          if(dmaxtest>=(*localgrid)->dv[i]) 
            break;
        }
      }
    // for each global cell set number of layeres (similarly)
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
  } 
  // 2D problem
  else {
    for(i=0;i<(*localgrid)->Nc;i++)
      (*localgrid)->Nk[i] = (*localgrid)->Nkmax;
    for(i=0;i<maingrid->Nc;i++)
      maingrid->Nk[i] = maingrid->Nkmax;
  }

  // initialize ctop, dztop
  for(i=0;i<(*localgrid)->Nc;i++) {
    (*localgrid)->ctop[i] = 0;
    (*localgrid)->dztop[i] = (*localgrid)->dz[(*localgrid)->ctop[i]];
  }

  // initialize etop, etopold
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
      // depth adjusted for cell disctretization in vertical
      (*localgrid)->dv[i]=dmaxtest;
    }

  // over all the local grid edges 
  for(j=0;j<(*localgrid)->Ne;j++) {
    // get clogal grid pointers
    ne = (*localgrid)->eptr[j];
    // if the neighboring cell is a ghost cell
    if(maingrid->grad[2*ne]==-1) {
      // set Nke = Nke = Nk for cell not ghost
      (*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
      (*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
    }
    // if the other neighboring cell is a ghost cell
    else if(maingrid->grad[2*ne+1]==-1) {
      // set Nke = Nke = Nk for cell not ghost
      (*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne]];
      (*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne]];
    }
    // if we have a computational cell (not by ghost cell)
    else {
      // if we don't have the same number of layers for 
      // adjacent cells 
      if(maingrid->Nk[maingrid->grad[2*ne]]<
          maingrid->Nk[maingrid->grad[2*ne+1]]) {
        // (Nke is edge layers which is smaller
        // than Nkc which is cell layers)
        (*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne]];
        (*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
      } 
      // other case (switched orientation)
      else if(maingrid->Nk[maingrid->grad[2*ne]]>
          maingrid->Nk[maingrid->grad[2*ne+1]]) {
        (*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne+1]];
        (*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne]];
      } 
      else {
        // case where the cell layers are the same Nke=Nkc=Nk
        (*localgrid)->Nke[j]=maingrid->Nk[maingrid->grad[2*ne]];
        (*localgrid)->Nkc[j]=maingrid->Nk[maingrid->grad[2*ne]];
      }
    }
  }

  /* compute Nkp for the number of layers for each node*/
  // for each node
  for(i = 0; i < (*localgrid)->Np; i++) {
    // look over its cell neighbors and take the max of Nkc
    // for each of the cell neighbors
    
    maxNk = (*localgrid)->Nk[(*localgrid)->pcneighs[i][0]];
    for(npc = 1; npc < (*localgrid)->numpcneighs[i]; npc++) {
      // find the max
      maxNk = max(maxNk, (*localgrid)->Nk[(*localgrid)->pcneighs[i][npc]]);
    }
    // set value to max found (-1 if this is a bogus [TRIANGLE] point)
    (*localgrid)->Nkp[i] = maxNk;
    if((*localgrid)->numpcneighs[i]==0)
      (*localgrid)->Nkp[i]=0;
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

/*
 * Function: MakePointers
 * Usage:  MakePoniters(maingrid, localgrid, myproc, comm)
 * -------------------------
 * Functon to make pointers to cellp, edgep, celldist, edgedist, lcptr, leptr
 * for convenience in specification in boundaries.c and initialization.c etc
 *
 */
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

  // place these values on localgrid
  (*localgrid)->cellp=cellp;
  (*localgrid)->edgep=edgep;
  (*localgrid)->celldist=celldist;
  (*localgrid)->edgedist=edgedist;

  // debug info
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

  // initialize to 0
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

  // this is cool because all we have to do is make sure we have included
  // all the correct cells in the following loop in the cell_recv structure
  // and the rest of the code does the rest

  // for each of the neighboring processors
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    
    // get number of neighboring processor
    neighproc=(*localgrid)->myneighs[neigh];

    // for each cell on local grid
    for(n=0;n<(*localgrid)->Nc;n++) {
      // see if the processor number is the same as the neighboring processor
      if(maingrid->part[(*localgrid)->mnptr[n]]==neighproc) {
        // if so, add to count for cell that will need to be recieved for particular
        // processor
        num_cells_recv[neigh]++;
      }
    }
    
    // allocae memory now that we have count of cells to be recieved (just based on 
    // coloring from the partitioning and should be automatically updated if the 
    // part array is set to be based on nodes instead of edges)
    cell_recv[neigh]=(int *)SunMalloc(num_cells_recv[neigh]*sizeof(int),"MakePointers");

    k=0;

    // for each cell on the local grid
    for(n=0;n<(*localgrid)->Nc;n++) {
      // if the cell is on the same processor as the neighboring processor
      if(maingrid->part[(*localgrid)->mnptr[n]]==neighproc) {
        if(IsNeighborGlobal(n, maingrid, localgrid, myproc, g_halolist, g_halolistsize)) {
              cell_recv[neigh][k++]=(*localgrid)->mnptr[n];
        }
      }
    }
  }
  // debug information for processor passing
  if(VERBOSE>2) 
    for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++)
      printf("Proc: %d, neighbor %d, receiving %d\n",
	     myproc,(*localgrid)->myneighs[neigh],
	     num_cells_recv[neigh]);

  // Send out the number of cells that are being received and then
  // the actual global indices being sent.
 
  // send data
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    // over each of the neighboring processors
    neighproc=(*localgrid)->myneighs[neigh];
    // send the amount of data that must be transfered 
    MPI_Send(&(num_cells_recv[neigh]),1,MPI_INT,neighproc,1,comm); 
  }

  // recieve sent data
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Recv(&(num_cells_send[neigh]),1,MPI_INT,neighproc,1,comm,&status);
  }

  // send data
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Send(cell_recv[neigh],num_cells_recv[neigh],MPI_INT,neighproc,1,comm); 
  }

  // receive data
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    cell_send[neigh]=(int *)SunMalloc(num_cells_send[neigh]*sizeof(int),"MakePointers");
    MPI_Recv(cell_send[neigh],num_cells_send[neigh],MPI_INT,neighproc,1,comm,&status);
  }

  // Now set the indices to point to the local grid.
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    for(j=0;j<num_cells_send[neigh];j++) 
      // flipping from global to local coordinates
      cell_send[neigh][j]=lcptr[cell_send[neigh][j]];
    for(j=0;j<num_cells_recv[neigh];j++)
      // flipping from global to local coordinates
      cell_recv[neigh][j]=lcptr[cell_recv[neigh][j]];
  }

  // Now do the edges.  The edges that correspond to the cells that
  // are sent are added to the edge-send/recv arrays.

  // accumulate number of edges for each neighboring processor
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    // get the neighboring processor
    neighproc=(*localgrid)->myneighs[neigh];
    // for each edge initialize the flag
    for(j=0;j<(*localgrid)->Ne;j++)
      flagged[j]=0;
    // for each of the boundary cells 
    for(i=0;i<num_cells_send[neigh];i++)
      // for each of the faces 
      for(nf=0;nf<(*localgrid)->nfaces[cell_send[neigh][i]];nf++)
        // check to make sure the particular edge hasn't been flagged
        if(!flagged[(*localgrid)->face[(*localgrid)->maxfaces*cell_send[neigh][i]+nf]]) {
          flagged[(*localgrid)->face[(*localgrid)->maxfaces*cell_send[neigh][i]+nf]]=1;
          num_edges_send[neigh]++;
        }
    // allocate memory for flagged edges
    edge_send[neigh]=(int *)SunMalloc(num_edges_send[neigh]*sizeof(int),"MakePointers");
    // set flagged to 0
    for(j=0;j<(*localgrid)->Ne;j++)
      flagged[j]=0;
    k=0;
    // compute sent data (edge pointers) 
    for(i=0;i<num_cells_send[neigh];i++)
      for(nf=0;nf<(*localgrid)->nfaces[cell_send[neigh][i]];nf++)
        if(!flagged[(*localgrid)->face[(*localgrid)->maxfaces*cell_send[neigh][i]+nf]]) {
          flagged[(*localgrid)->face[(*localgrid)->maxfaces*cell_send[neigh][i]+nf]]=1;
          edge_send[neigh][k++]=
            (*localgrid)->eptr[(*localgrid)->face[(*localgrid)->maxfaces*cell_send[neigh][i]+nf]];
        }
    // send data (number of edges on neighbor)
    MPI_Send(&(num_edges_send[neigh]),1,MPI_INT,neighproc,1,comm); 
  }

  // recieve data
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Recv(&(num_edges_recv[neigh]),1,MPI_INT,neighproc,1,comm,&status);
  }

  // send data (edges data)
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    MPI_Send(edge_send[neigh],num_edges_send[neigh],MPI_INT,neighproc,1,comm); 
  }

  // receive data
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

  // debugging output
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
  // use lock to make sure all communication is completed before continuing
  MPI_Barrier(comm);

  // store all data just sent/received
  (*localgrid)->cell_send=cell_send;
  (*localgrid)->cell_recv=cell_recv;
  (*localgrid)->edge_send=edge_send;
  (*localgrid)->edge_recv=edge_recv;
  (*localgrid)->num_cells_send=num_cells_send;
  (*localgrid)->num_cells_recv=num_cells_recv;
  (*localgrid)->num_edges_send=num_edges_send;
  (*localgrid)->num_edges_recv=num_edges_recv;

  // free temporary arrays
  free(lcptr);
  free(leptr);
}

/*
 * Function: ReOrder
 * Usage:  ReOrder(grid)
 * -------------------------
 * Functon to reorder grid 
 */
static void ReOrder(gridT *grid) 
{
  int n, nf, numflag, options[8], Nc = grid->Nc, Ne=grid->Ne, grad1, grad2, enei;
  int *corder, *corderp, *eorder, *eorderp;
  REAL *tmp;

  corder = (int *)SunMalloc(Nc*sizeof(int),"ReOrder");
  corderp = (int *)SunMalloc(Nc*sizeof(int),"ReOrder");
  eorder = (int *)SunMalloc(Ne*sizeof(int),"ReOrder");
  eorderp = (int *)SunMalloc(Ne*sizeof(int),"ReOrder");
  tmp = (REAL *)SunMalloc(2*(grid->maxfaces-1)*Ne*sizeof(REAL),"ReOrder");

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
  ReOrderRealArray(grid->Ac,corder,tmp,grid->Nc,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderRealArray(grid->xv,corder,tmp,grid->Nc,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderRealArray(grid->yv,corder,tmp,grid->Nc,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderRealArray(grid->dv,corder,tmp,grid->Nc,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->cells,corder,(int *)tmp,grid->Nc,grid->maxfaces,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->face,corder,(int *)tmp,grid->Nc,grid->maxfaces,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->normal,corder,(int *)tmp,grid->Nc,grid->maxfaces,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->neigh,corder,(int *)tmp,grid->Nc,grid->maxfaces,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->vwgt,corder,(int *)tmp,grid->Nc,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->mnptr,corder,(int *)tmp,grid->Nc,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->Nk,corder,(int *)tmp,grid->Nc,1,grid->nfaces,grid->grad,grid->maxfaces);

  //Reorder the data corresponding to edges with the eorder array
  ReOrderRealArray(grid->df,eorder,tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderRealArray(grid->dg,eorder,tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderRealArray(grid->n1,eorder,tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderRealArray(grid->n2,eorder,tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderRealArray(grid->xi,eorder,tmp,grid->Ne,2*(grid->maxfaces-1),grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->grad,eorder,(int *)tmp,grid->Ne,2,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->eneigh,eorder,(int *)tmp,grid->Ne,2*(grid->maxfaces-1),grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->edges,eorder,(int *)tmp,grid->Ne,NUMEDGECOLUMNS,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->mark,eorder,(int *)tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->eptr,eorder,(int *)tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->Nke,eorder,(int *)tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);
  ReOrderIntArray(grid->Nkc,eorder,(int *)tmp,grid->Ne,1,grid->nfaces,grid->grad,grid->maxfaces);

  // Now adjust the pointers to point to the new locations
  // Face and eneigh arrays point to edges.  We need to use the createfacearray
  // again because it guarantees consistency for the CG solver.  This is
  // because reordering does not necessarily guarantee that the ordering of
  // the faces in the face array are in the same as the corresponding neighbors
  // in the neigh array.
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<grid->nfaces[n];nf++) 
      tmp[n*grid->maxfaces+nf]=grid->face[n*grid->maxfaces+nf];
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<grid->nfaces[n];nf++) 
      grid->face[n*grid->maxfaces+nf]=eorderp[(int)(tmp[n*grid->maxfaces+nf])];

  for(n=0;n<Ne;n++)
    grad1=grid->grad[2*n];
    grad2=grid->grad[2*n+1];
    enei=grad1+grad2-2;
    for(nf=0;nf=enei;nf++) 
      tmp[2*(grid->maxfaces-1)*n+nf]=grid->eneigh[2*(grid->maxfaces-1)*n+nf];
  for(n=0;n<Ne;n++)
    grad1=grid->grad[2*n];
    grad2=grid->grad[2*n+1];
    enei=grad1+grad2-2;
    for(nf=0;nf<enei;nf++)
      if((int)(tmp[n*2*(grid->maxfaces-1)+nf])!=-1)
	grid->eneigh[n*2*(grid->maxfaces-1)+nf]=eorderp[(int)(tmp[n*2*(grid->maxfaces-1)+nf])];
      else
	grid->eneigh[n*2*(grid->maxfaces-1)+nf]=-1;

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
    for(nf=0;nf<grid->nfaces[n];nf++)
      tmp[n*grid->maxfaces+nf]=grid->neigh[n*grid->maxfaces+nf];
  for(n=0;n<Nc;n++)
    for(nf=0;nf<grid->nfaces[n];nf++)
      if((int)(tmp[n*grid->maxfaces+nf])!=-1)
	grid->neigh[n*grid->maxfaces+nf]=corderp[(int)(tmp[n*grid->maxfaces+nf])];
      else
  	grid->neigh[n*grid->maxfaces+nf]=-1;

  free(corder);
  free(corderp);
  free(eorder);
  free(eorderp);
  free(tmp);
}

/*
 * Function: EdgeMarkers
 * Usage: EdgeMarkers(maingrid, localgrid, myproc)
 * --------------------------------
 * Compute interprocessor edge markers where 
 *    mark = 5    edge is a ghost cell on interprocessor
 *                boundary
 *    mark = 6    edge is ghost or computational but not
 *                an interprocessor boundary
 *
 */
// could be written to be more efficient, potentially missing some
// of the important ideas here
static void EdgeMarkers(gridT *maingrid, gridT **localgrid, int myproc)
{
  int n, ne, nf;

  // for the cells on the local grid
  for(n=0;n<(*localgrid)->Nc;n++) {
    // if the cell is a ghost cell
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==3) {
      // for each face on the ghost cell
      for(nf=0;nf<(*localgrid)->nfaces[n];nf++) {
        // get the edge information
        ne = (*localgrid)->face[n*(*localgrid)->maxfaces+nf];
        // if the edge is marked as ghost or computational
        if(!(*localgrid)->mark[ne]) 
          // we have type 6 boundary condition
          (*localgrid)->mark[ne]=6;
      }
    }
  }
  // for each cell on the local grid 
  for(n=0;n<(*localgrid)->Nc;n++) {
    // if the cell is an inter processor boundary cell
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2) {
      // for each face on the inter proc boundary cell
      for(nf=0;nf<(*localgrid)->nfaces[n];nf++) {
        // get local edge ponter
        ne = (*localgrid)->face[n*(*localgrid)->maxfaces+nf];
        // if the marked cell has an edge which is marked as
        // ghost or computational but we know it is an 
        // interproc boundary, mark as 5
        if((*localgrid)->mark[ne]==6) 
          (*localgrid)->mark[ne]=5;
      }
    }
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
  // for each cell on the grid
  for(j=0;j<grid->Nc;j++) {
    // determine if the given cell is local to the processor
    // by determining if any of it's neighbors are local to the
    // processor
    flag = IsNeighborLocal(j,grid,proc, g_halolist, g_halolistsize);
    // if the current cell is on the processor or if it's 
    // neighbor is on the processor (indicating cell j
    //PJw is a ghost point), then we can add this cell to the 
    // total number of cells on the processor
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

  //  initiallize all flags corresponding to local cell to 0
  for(j=0;j<grid->Nc;j++) 
    flagged[j]=0;

  nedges=0;
  // for each cell on the local grid
  for(j=0;j<grid->Nc;j++) {
    // mark that the cell has already been considered
    flagged[j]=1;
    // over each of the faces consider the neighbor to 
    // the edge
    for(nf=0;nf<grid->nfaces[j];nf++) {
      ne = grid->neigh[j*grid->maxfaces+nf];
      // if the neighbor to the edge is a ghost cell
      // add the edge to the count
      if(ne == -1)
        nedges++;
      // if the neighbor hasn't been flagged (add edge to count)
      else if(!flagged[ne])
        nedges++;
    }
  }
  free(flagged);

  return nedges;
}

/*
 * Function: GetNumPoints
 * Usage: GetNumPoints(grid);
 * -------------------------
 * Returns the number of points
 * belonging to a particular processor that result from
 * the partitioning.
 *
 */
int GetNumPoints(gridT *localgrid, gridT *maingrid)
{
  int j, nf, npoints, np1, np2;

  unsigned short *flagged = (unsigned short *)
    SunMalloc(maingrid->Np*sizeof(unsigned short),"GetNumPoints");

  //  initiallize all flags corresponding to local cell to 0
  for(j=0;j<maingrid->Np;j++) {
    flagged[j]=0;
  }
  npoints=0;

  // for each edge on the local localgrid
  for(j=0;j<localgrid->Ne;j++) {
    // get the number of points comprising the edge (since points are global)
    np1 = localgrid->edges[NUMEDGECOLUMNS*j];
    np2 = localgrid->edges[NUMEDGECOLUMNS*j+1];
    // mark that the points forming the edge have already been considered
    flagged[np1]=1; 
    flagged[np2]=1;
  }
  // sum up all the flagged values to get the total number of points
  for(j=0; j < maingrid->Np; j++) {
    if(flagged[j] == 1) {
      npoints += 1;
    }
  }

  free(flagged);

  return npoints;
}

// function
/*
 * Function: GetProcPoints
 * Usage: GetProcPoints(grid);
 * -------------------------
 * Returns a list of the points
 * belonging to a particular processor that result from
 * the partitioning.
 *
 */
void GetProcPoints(gridT *localgrid, gridT *maingrid)
{
  int j, nf, npoints, np1, np2;

  unsigned short *flagged = (unsigned short *)
    SunMalloc(maingrid->Np*sizeof(unsigned short),"GetNumPoints");

  // note that sum of all processors number of local points will be greater than the
  // actual number of points since points on the boundary will be counted more than once.
  // this also means that the nodal neighbors for these points will be bogus as only the 
  // non-halo points will make sense and have the same nodal data structure as found in 
  // global.  The most extreme boundary points will only have partial data compared to 
  // the global and should never be  called.  This is expected.  Thus, the most
  // rigorous way to set this up is so that the sum of all the processor points is the
  // global number of points. The definition gets looser for ELM anyway so this condition
  // isn't necessary are will only be minorly related to efficiency (at leaset it 
  // seems that way for now)

  // assume for now that we have only the points from the main grid
  localgrid->Np = maingrid->Np;
  //  just copy all the points
  for(j=0;j<maingrid->Np;j++) {
    flagged[j]=1;
  }

  // allocate the memory from the previous computation of the total number of 
  // points on the localgrid
  localgrid->localtoglobalpoints = (int *)SunMalloc(localgrid->Np*sizeof(int),"GetProcPoints");

  // make list with all the flagged values and check the total number of points
  // this may not be necessary just make sure to pass the total number of 
  // points
  npoints = 0;
  for(j=0; j < maingrid->Np; j++) {
    if(flagged[j] == 1) {
      localgrid->localtoglobalpoints[npoints++] = j;
    }
  }

  free(flagged);
}

/*
 * Function: IsCellNeighborProc
 * Usage: IsCellNeighborProc(nc, maingrid, localgrid, myproc, neighproc)
 * ----------------------------------------------------
 * Determine if the cell is neighboring processor neighproc
 *
 */
int IsCellNeighborProc(int nc, gridT *maingrid, gridT *localgrid, 
    int myproc, int neighproc)
{
  int nf, ne;

  // case where the cell is not on the local grid or neighboring grid
  // first part of condition filters out other cells to make it local
  if(maingrid->part[localgrid->mnptr[nc]]!=myproc &&
      maingrid->part[localgrid->mnptr[nc]]!=neighproc)
    return 0;

  // cell does neighbor neighboring proc
  if(maingrid->part[localgrid->mnptr[nc]]==neighproc)
    return 1;

  // consider each of the faces and get neighbors
  for(nf=0;nf<localgrid->nfaces[nc];nf++) {
    ne = localgrid->neigh[nc*localgrid->maxfaces+nf];
    if(ne != -1)
      // for non-ghost cell if neighboring processor is 
      // on neighproc
      if(maingrid->part[localgrid->mnptr[ne]] == neighproc)
        return 1;
  }
  return 0;
}

/*
 * Function: IsEdgeNeighborProc
 * Usage: IsEdgeNeighborProc(ne, maingrid, localgrid, myproc, neighproc)
 * ----------------------------------------------------
 * Determine if the edge is neighboring neighproc
 *
 */
int IsEdgeNeighborProc(int ne, gridT *maingrid, gridT *localgrid, 
    int myproc, int neighproc)
{
  int j, cnp, nc;

  for(j=0;j<2;j++) {
    nc = localgrid->grad[2*ne+j];
    // for both cells neighboring edge (not ghost)
    if(nc != -1) {
      cnp = IsCellNeighborProc(nc,maingrid,localgrid,myproc,neighproc);
      if(cnp)
        return 1;
    }
  }
  return 0;
}

/*
 * Function: IsNodeNeighborProc
 * Usage: IsNodeNeighborProc(np, maingrid, localgrid, myproc, neighproc)
 * ----------------------------------------------------
 * Determine if the node is neighboring neighproc.  
 *
 */
int IsNodeNeighborProc(int np, gridT *maingrid, gridT *localgrid, 
		       int myproc, int neighproc)
{
  printf("Header for IsNodeNeighborProc.  Note implemented yet!!!\n");
  
}

/*
 * Function: Topology
 * Usage: Topology(maingrid,localgrid,myproc,numprocs);
 * ----------------------------------------------------
 * Compute the number of neighbors for each processor and
 * determine the processor ids of each neighbor.  This function
 * effectively computes grid->numneighs, grid->neighs, 
 * localgrid->Nneighs, and localgrid->myneighs, which is used 
 * to set up passing of the interprocessor boundary values within 
 * MakePointers.
 *
 */
void Topology(gridT **maingrid, gridT **localgrid, int myproc, int numprocs)
{
  int proc, j;
  char filename[BUFFERLENGTH];
  //int nf, j, cfneigh, neigh1, neigh2, loc, pcell, pneigh, proc;

  // set up primary maingrid data structures for neighbors
  // this gives maingrid->numneighs=0 and maingrid->neighs=[-1,...,-1]
  (*maingrid)->numneighs=(int *)SunMalloc(numprocs*sizeof(int),"Topology");
  (*maingrid)->neighs=(int **)SunMalloc(numprocs*sizeof(int *),"Topology");

  // Initialize the arrays
  for(proc=0;proc<numprocs;proc++) {
    // initially there are no neighboring processors (as for the serial case)
    (*maingrid)->numneighs[proc]=0;
    // assume the worst that each processor neighbors this processor's data,
    // in the case of a perfect island
    (*maingrid)->neighs[proc]=(int *)SunMalloc(numprocs*sizeof(int),"Topology");
    // initialize all neighboring processors to be -1 (ghost) until otherwise noted
    // where proc is in the set 0...numprocs
    for(j=0;j<numprocs;j++) 
      (*maingrid)->neighs[proc][j]=-1;
  }

  BoundaryTopologyByHalo(*maingrid, numprocs, g_halolist, g_halolistsize);
  

  // transfer the number of neighbors information to the local grids
  // irrspective of method used to find them
  // maingrid->numneighs => localgrid->Nneighs and 
  // maingrid->neighs    => localgrid->myneighs
  (*localgrid)->Nneighs=(*maingrid)->numneighs[myproc];
  // allocate memory for the list of neighboring processors to the local grid
  // being conservative with memory this time (contrast with the global)
  (*localgrid)->myneighs=(int *)SunMalloc((*localgrid)->Nneighs*sizeof(int),"Topology");

    printf("Proc %d, neighs: ",myproc);
  // transfer the list of neighbors to the local grid
  for(j=0;j<(*localgrid)->Nneighs;j++) {
    (*localgrid)->myneighs[j]=(*maingrid)->neighs[myproc][j];
        printf("%d ",(*localgrid)->myneighs[j]);
  }
    printf("\n");
  
}

/*
 * Function: BoundaryTopology
 * Usage: BoundaryTopology(*maingrid, numprocs);
 * ----------------------------------------------------
 * Called by Topology and used to compute the number 
 * of neighbors for each processor and
 * determine the processor ids of each neighbor by using
 * shared edges to determine the boundary cells.  
 * This function effectively computes 
 * grid->numneighs, grid->neighs. 
 * These are both of size [numprocs][numprocs]
 * This function only shows how processor regions are
 * connected and says nothing about how cells are 
 * connected.  The criteria is selectable by edge or node.
 *
 */
static inline void  BoundaryTopologyByHalo(gridT *maingrid, int numprocs, 
    boundaryselection *list, int listsize)
{
  int j; 

  // if we don't have a good list, exit
  if(listsize < 1) {
    printf("Error!!: List for boundary neighbors needs to be larger!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }

  // over each cell for the whole domain
  for(j=0;j < maingrid->Nc; j++) {

    // select which function to call off the first item in the list
    switch(list[0]) {
      case EDGE:
        // original selection by edge
        BoundaryTopologyByEdge(j, j, maingrid, numprocs, list, listsize-1);
        break;
      case NODE:
        // new selection based on nodal neighbors (more inclusive)
        BoundaryTopologyByNode(j, j, maingrid, numprocs, list, listsize-1);
        break;
      default:
        printf("Error in BoundaryTopologyByHalo!!\n");
        break;
    }
  }
}

/*
 * Function: BoundaryTopologyByEdge
 * Usage: BoundaryTopologyByEdge(*maingrid, numprocs);
 * ----------------------------------------------------
 * Called by Topology and used to compute the number 
 * of neighbors for each processor and
 * determine the processor ids of each neighbor by using
 * shared edges to determine the boundary cells.  
 * This function effectively computes 
 * grid->numneighs, grid->neighs. 
 * These are both of size [numprocs][numprocs]
 * This function only shows how processor regions are
 * connected and says nothing about how cells are 
 * connected.  The criteria is for processors which 
 * share the same edge (i.e. cells on different processors
 * have the same edge).
 *
 */
static void  BoundaryTopologyByEdge(int cellbase, int cell, gridT *maingrid, 
    int numprocs, boundaryselection *list, int listsize)
{
  int nf, cfneigh, loc, pcell, pneigh;

  // if we don't have a good list, exit
  if(listsize < 0) {
    printf("Error!!: List for boundary edge neighbors needs to be larger!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }

  // consider each face (edge) on the cells 
  for(nf=0;nf < maingrid->nfaces[cell]; nf++) {
    // if the neighbor to the cell's face (cfneigh) is not a ghost 
    cfneigh = maingrid->neigh[maingrid->maxfaces*cell+nf];
    // check to see if the cell should be added to the boundary
    if(cfneigh != -1) {
      // if it is not a ghost cell, get the cell processors (compare with base)
      pcell = maingrid->part[cellbase];
      pneigh = maingrid->part[cfneigh];
      // if the cell proc (pcell) is not the same as it's neighbor's proc (pneigh)
      if(pcell != pneigh) {
        // loc will only show up if it has already been added to list, 
        // meaning pneigh has been recorded as a neighboor of pcell if it has, do 
        // nothing, otherwise add it to the list
        loc = IsMember(pneigh,maingrid->neighs[pcell],numprocs);
        if(loc<0) {
          // add to list where numneighs is effectively an incrementer for the end
          // of the list
          maingrid->neighs[pcell][maingrid->numneighs[pcell]]=pneigh;
          maingrid->numneighs[pcell]++;
        }
      }
      /* now process other items on the list, if any */

      // if there are more items to be considered
      if(listsize > 0) {

        // select which function to call off the first item in the list
        switch(list[1]) {
          case EDGE:
            // original selection by edge
            BoundaryTopologyByEdge(cellbase, cfneigh, maingrid, numprocs,
                &list[1], listsize-1);
            break;
          case NODE:
            // new selection based on nodal neighbors (more inclusive)
            BoundaryTopologyByNode(cellbase, cfneigh, maingrid, numprocs, 
                &list[1], listsize-1);
            break;
          default:
            printf("Error in BoundaryTopologyByEdge!!\n");
            break;
        }
      }
    }
  }
}

/*
 * Function: BoundaryTopologyByNode
 * Usage: BoundaryTopologyByNode(*maingrid, numprocs);
 * ----------------------------------------------------
 * Called by Topology and used to compute the number 
 * of neighbors for each processor and
 * determine the processor ids of each neighbor by using
 * shared edges to determine the boundary cells.  
 * This function effectively computes 
 * grid->numneighs, grid->neighs. 
 * These are both of size [numprocs][numprocs]
 * This function only shows how processor regions are
 * connected and says nothing about how cells are 
 * connected.  The criteria is for processors which share
 * the same node (i.e. cells on different processors
 * have the same node)
 *
 */
static void  BoundaryTopologyByNode(int cellbase, int cell, gridT *maingrid, 
    int numprocs, boundaryselection *list, int listsize)
{
//  printf("NODE\n");
  int nf,  cpneigh, loc, pcell, pneigh, nn, in;

  // if we don't have a good list, exit
  if(listsize < 0) {
    printf("Error!!: List for boundary node neighbors needs to be larger!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }

  // consider each node on the cells 
  for(nf=0;nf < maingrid->nfaces[cell]; nf++) {
    // get the particular node
    in = maingrid->cells[maingrid->maxfaces*cell + nf];
    // get the nodal neighbor to the particular node
    for(nn = 0; nn < maingrid->numpcneighs[in]; nn++) {
      // get nodal cell neighbor
      cpneigh = maingrid->pcneighs[in][nn];
      // if the neighbor to the cell's node (cpneigh) is not a ghost 
      if(cpneigh != -1) {
        // get the cell processors (compare with base)
        pcell = maingrid->part[cellbase];
        pneigh = maingrid->part[cpneigh];
        // if the cell proc (pcell) is not the same as it's neighbor's proc (pneigh)
        if(pcell != pneigh) {
          // loc will only show up if it has already been added to list, 
          // meaning pneigh has been recorded as a neighboor of pcell if it has, do 
          // nothing, otherwise add it to the list
          loc = IsMember(pneigh,maingrid->neighs[pcell],numprocs);
          if(loc<0) {
            // add to list where numneighs is effectively an incrementer for the end
            // of the list
            maingrid->neighs[pcell][maingrid->numneighs[pcell]]=pneigh;
            maingrid->numneighs[pcell]++;
          }
        }
        /* now process other items on the list, if any */

        // if there are more items to be considered
        if(listsize > 0) {

          // select which function to call off the first item in the list
          switch(list[1]) {
            case EDGE:
              // original selection by edge
              BoundaryTopologyByEdge(cellbase, cpneigh, maingrid, numprocs, 
                  &list[1], listsize-1);
              break;
            case NODE:
              // new selection based on nodal neighbors (more inclusive)
              BoundaryTopologyByNode(cellbase, cpneigh, maingrid, numprocs, 
                  &list[1], listsize-1);
              break;
            default:
              printf("Error in BoundaryTopologyByNode!!\n");
              break;
          }
        }
      }
    }
  }
//  printf("Terminal at node %d to go\n", listsize);
}

/*
 * Function: TransferData
 * Usage: TransferData(maingrid,localgrid,myproc);
 * -----------------------------------------------
 *
 * This function transfers cell data from the main grid onto the local 
 * grid and creates the global pointer that allows the identification of
 * the global cell index given a local index.  Because the main grid 
 * is stored on every processor, no communcation is needed. The computational
 * cells are printed first, followed by the boundary cells.
 *
 * If only edge neighbors are considered:
 * The local grid will have new ghost cells due to the cut edges.
 * There is one ghost cell for each cut edge.  
 * If nodal neighbors are considered:
 * The local grid will have new ghost cells due to cut edges where
 * all nodal neighbors (or halo neighbors) to the cut cell node
 * is considered so that there is no direct way to easily determine
 * the total number of ghost cells per cut edge.
 *
 * The local cells
 * array points to indices in the main xp,yp arrays.  But the local
 * neigh pointer points to local cell indices.
 * lcptr points from a cell in the main grid to one in the local grid.
 * mnptr points from a cell in the local grid to one in the main grid
 * and is stored on every local grid.
 *
 */
// may need to make the local cells array pointer to indicies in the 
// local array which is converted via localgrid->localtoglobalpoints 
// will need to make cells point to indicies in the local xp,yp arrays, 
// as these are loaded in ReadData and stored with node data, so the global
// index is no longer needed and a local index will suffice.  Since it is only
// used in nodal calculations it should be easily changed.
static void TransferData(gridT *maingrid, gridT **localgrid, int myproc)
{
  // allocate memory for all referential pointers and localgrid memory
  int i, j, k, n, nc, nf, ne, ng, flag, mgptr, *lcptr, *leptr, bctype, iface, grad1, grad2, enei;
  unsigned short *flagged = 
    (unsigned short *)SunMalloc(maingrid->Ne*sizeof(unsigned short),"TransferData");

  // pointers from maingrid to localgrid for cells and edges
  lcptr = (int *)SunMalloc(maingrid->Nc*sizeof(int),"TransferData");
  leptr = (int *)SunMalloc(maingrid->Ne*sizeof(int),"TransferData");

  // this GetNumCells relies upon topology and setting the number of 
  // cells for each proc
  (*localgrid)->Nc = GetNumCells(maingrid,myproc);
  (*localgrid)->maxfaces = maingrid->maxfaces;   // added part 
  // pointer to main grid for the number of cells on the local grid 
  // and weightings allocated
  (*localgrid)->mnptr = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),
      "TransferData");
  (*localgrid)->vwgt = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),
      "TransferData");

  (*localgrid)->cells = (int *)SunMalloc((*localgrid)->maxfaces*(*localgrid)->Nc*sizeof(int),
      "TransferData");

  (*localgrid)->nfaces = (int *)SunMalloc((*localgrid)->Nc*sizeof(int),"TransferData"); //added part
  (*localgrid)->xv = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"TransferData");
  (*localgrid)->yv = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"TransferData");
  (*localgrid)->xp = (REAL *)SunMalloc(maingrid->Np*sizeof(REAL),"TransferData");
  (*localgrid)->yp = (REAL *)SunMalloc(maingrid->Np*sizeof(REAL),"TransferData");
  (*localgrid)->dv = (REAL *)SunMalloc((*localgrid)->Nc*sizeof(REAL),"TransferData");

  // transfer point data to local grid
  for(j=0;j<maingrid->Np;j++) {
    (*localgrid)->xp[j] = maingrid->xp[j];
    (*localgrid)->yp[j] = maingrid->yp[j];
  }

  // this neigh critierion will need to be adapted 
  // to account for nodal neighbors, not edge neighbors
  (*localgrid)->neigh = (int *)SunMalloc((*localgrid)->maxfaces*(*localgrid)->Nc*sizeof(int),
      "TransferData");
  (*localgrid)->normal = (int *)SunMalloc((*localgrid)->maxfaces*(*localgrid)->Nc*sizeof(int),
      "TransferData");
      
  // initialize lcptr and leptr pointers to ghost (non-existent) criterion
  for(j=0;j<maingrid->Nc;j++) 
    lcptr[j]=-1;
  for(j=0;j<maingrid->Ne;j++) 
    leptr[j]=-1;

  // populate maingrid->lcptr, localgrid->mnptr, transver xv, yv, dv, vwgt, cells
  // this function utilizes the test to see if a cell is a boundary cell 
  // needed to determine if it is to be included in the list
  // currently this utilizes edge values via IsBoundaryCEllByEdge
  k=0;
  // MAXBCTYPES = 4 for now, check each boundary condition type 
  // (could be more efficient, probably
  // as the only BC type that we utilize is really the type 3 (aliased to bctype-1)
  // could just do the conditional test for all the BC types 
  // (Spell out in line IsBoundaryCell...)
  for(bctype=0;bctype<MAXBCTYPES;bctype++) {
    // loop over all the main grid cells
    for(j=0;j<maingrid->Nc;j++) {
      // check to see if the cell on the main grid and processor has
      // a consistent boundary condition type (cell BC like type 3,
      // not edge like type 2 or 4)
      if(IsBoundaryCell(j,maingrid,myproc)==bctype) {

        // check all the neighbors and for non-ghost cells 
        // which have the same processor as myproc, return flag = 1
        // flag = 1 if interprocessor boundary that's not a ghost cell
        // could be replaced by IsBoundaryCell == 2 means flag = 1

        // will need to generalize this for nodal neighbors
        // 
        // if(IsNeighborLocal(cell, maingrid, processor)) for
        // (IsNeighborByEdgeLocal and IsNeighborByNodeLocal and such)
        // then set flag is 1, else 0
        // check over all the faces
        // if this cell is a non-ghost interprocessor boundary cell

        flag = IsNeighborLocal(j,maingrid,myproc, g_halolist, g_halolistsize);

        // if we have a previous hit or the cell is on the processor, 
        //add to the structure
        // to index local cells to global cells
        // if it has been flagged because one of it's neighbors is on the 
        // local processor or the cell is on the local processor
        if(flag==1 || maingrid->part[j]==myproc) {
          // the local pointer points to list of cells 
          // on local processor from global value
          lcptr[j]=k;
          // get the pointer from the local array (k) to main array (j)
          (*localgrid)->mnptr[k]=j;
          // transfer relevant data from main grid to local grid
          (*localgrid)->xv[k]=maingrid->xv[j];
          (*localgrid)->nfaces[k]=maingrid->nfaces[j]; //added part
          (*localgrid)->yv[k]=maingrid->yv[j];
          (*localgrid)->dv[k]=maingrid->dv[j];
          (*localgrid)->vwgt[k]=maingrid->vwgt[j];
          // for each face connect cell pointer with pointers to 
          // points that make up a cell (not edges)
          for(nf=0;nf<maingrid->nfaces[j];nf++)
            // pointers to the points comprising the cell
            (*localgrid)->cells[k*(*localgrid)->maxfaces+nf]=maingrid->cells[j*maingrid->maxfaces+nf];
          k++;
        }
      }
    }
  }
//  printf("k = %d localgrid->Nc = %d\n", k, (*localgrid)->Nc);

  // populate localgrid->neigh  (just transfer values)
  // for each cell on the local grid
  for(j=0;j<(*localgrid)->Nc;j++) {
    // for each face on the local cell
    for(nf=0;nf<(*localgrid)->nfaces[j];nf++) {
      // pointer to gradient neighbor on main cell
      mgptr = maingrid->neigh[(*localgrid)->mnptr[j]*(*localgrid)->maxfaces+nf];
      // if this is not a ghost cell and is on the local grid)
      if(mgptr>=0 && IsMember(mgptr,(*localgrid)->mnptr,(*localgrid)->Nc)>=0) {
        // store the local grid neighbor info to the local grid coordinate
        (*localgrid)->neigh[j*(*localgrid)->maxfaces+nf]=lcptr[mgptr];
      }
      else {
        // otherwise it is a ghost cell
        (*localgrid)->neigh[j*maingrid->maxfaces+nf]=-1;
      }
    }
  }

  // compute the number of edges for the local grid including edges 
  // on a boundary (ghost edges)

  // allocate edge data
  (*localgrid)->Ne = GetNumEdges(*localgrid);
  (*localgrid)->edges = (int *)SunMalloc(NUMEDGECOLUMNS*(*localgrid)->Ne*sizeof(int),
      "TransferData");
  (*localgrid)->mark = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),
      "TransferData");
  (*localgrid)->eptr = (int *)SunMalloc((*localgrid)->Ne*sizeof(int),
      "TransferData");
  (*localgrid)->eneigh = (int *)SunMalloc(2*((*localgrid)->maxfaces-1)*(*localgrid)->Ne*sizeof(int),
      "TransferData");
  (*localgrid)->face = (int *)SunMalloc((*localgrid)->maxfaces*(*localgrid)->Nc*sizeof(int),
      "TransferData");
  (*localgrid)->grad = (int *)SunMalloc(2*(*localgrid)->Ne*sizeof(int),
      "TransferData");
  (*localgrid)->gradf = (int *)SunMalloc(2*(*localgrid)->Ne*sizeof(int),
      "TransferData");
  // populate maingrid->leptr, localgrid->eptr, localgrid->edges 
  // initialize all edge flags to 0
  for(j=0;j<maingrid->Ne;j++) 
    flagged[j]=0;

  k=0;
  // for each cell on the local grid
  for(i=0;i<(*localgrid)->Nc;i++) {
    // on each face on this cell
    for(nf=0;nf<(*localgrid)->nfaces[i];nf++) {
      // get the global face value (global pointer to edge)
      iface = maingrid->face[(*localgrid)->mnptr[i]*maingrid->maxfaces+nf];
      // if we already haven't considered this face
      if(!flagged[iface]) {
        // mark face as considered
        flagged[iface]=1;
        // transfer information for edges (pointers to points defining edge)
        for(j=0;j<NUMEDGECOLUMNS;j++) {
          (*localgrid)->edges[k*NUMEDGECOLUMNS+j]=
            maingrid->edges[iface*NUMEDGECOLUMNS+j];
        }
        // also connect corresponding mark
        (*localgrid)->mark[k]=maingrid->mark[iface];
        //Point from global edge index to local one
        leptr[iface]=k;
        //point from local edge index to global one
        (*localgrid)->eptr[k++]=iface;
      }
    }
  }


  // populate localgrid->grad
  // for each cell on the local grid
  for(i=0;i<(*localgrid)->Nc;i++) {
    //  for each edge on the local grid
    for(n=0;n<(*localgrid)->Ne;n++) {
      for(j=0;j<2;j++) {
        // initialize all neighbors to voronoi edge to ghost
        (*localgrid)->grad[2*n+j]=-1;
        // get voronoi edge neighbors from maingrid 
        nc = maingrid->grad[2*(*localgrid)->eptr[n]+j];
        // and if they aren't ghost make the connection between the 
        // local gradient and the local cell via the global information
        if(nc != -1)
          (*localgrid)->grad[2*n+j]=lcptr[nc];
      }
    }
  }  

  // create face and normal arrays now for the local grid 
  // (previously called for global grid)
  // check whether variables are right

  CreateFaceArray((*localgrid)->grad,(*localgrid)->gradf,(*localgrid)->neigh,
      (*localgrid)->face,(*localgrid)->nfaces,(*localgrid)->maxfaces,(*localgrid)->Nc,(*localgrid)->Ne);

  ReorderCellPoints((*localgrid)->face,(*localgrid)->edges,(*localgrid)->cells,(*localgrid)->nfaces,(*localgrid)->maxfaces,(*localgrid)->Nc);
  // there is no need to use Reordergradf here because face and gradf is already in order for maingrid 

  CreateNormalArray((*localgrid)->grad,(*localgrid)->face,(*localgrid)->normal,(*localgrid)->nfaces,(*localgrid)->maxfaces,(*localgrid)->Nc);

  // create node array which will keep track of all nodal neighboors
  // first, get number of total points
  // second, get a list of the points on each processor, this stores the 
  // values for each proc in localgrid->localpoints
  if(myproc==0 && VERBOSE>2)  printf("\t\tGetProcPoints\n");
  Tic();
  GetProcPoints(*localgrid, maingrid);
  if(myproc==0 && VERBOSE>2) printf("\t\t... time used is %f\n", Toc());

  // allocate nodal neighbor array information (list of neighbors for each point)
  // note that since the number of entries for each node is variable, 
  // these are a list to Np lists

  // -- careful here as the list for numppneighs may need to be larger
  (*localgrid)->numppneighs = (int *)SunMalloc((*localgrid)->Np*sizeof(int),"TransferData");
  (*localgrid)->ppneighs = (int **)SunMalloc((*localgrid)->Np*sizeof(int*),"TransferData");
  (*localgrid)->numpeneighs = (int *)SunMalloc((*localgrid)->Np*sizeof(int),"TransferData");
  (*localgrid)->peneighs = (int **)SunMalloc((*localgrid)->Np*sizeof(int*),"TransferData");
  (*localgrid)->numpcneighs = (int *)SunMalloc((*localgrid)->Np*sizeof(int),"TransferData");
  (*localgrid)->pcneighs = (int **)SunMalloc((*localgrid)->Np*sizeof(int*),"TransferData");
 
  // this function crashed because sampling of localgrid points is wrong 
  // as nodes are not reordered
  if(myproc==0 && VERBOSE>2)  printf("\t\tCreateNodeArray\n");
  Tic();

  CreateNodeArray((*localgrid), (*localgrid)->Np,(*localgrid)->Ne, (*localgrid)->Nc, myproc);
  if(myproc==0 && VERBOSE>2) printf("\t\t... time used is %f\n", Toc());

  // populate (*localgrid)->eneigh
  CreateMomentumCV(*localgrid);

  // free temporary information connecting global grid to local grid 
  // (will varry for each processor) and temporary flags
  free(flagged);
  free(lcptr);
  free(leptr);
}

/*
 * Function: IsNeighborLocal
 * Usage: IsNeighborLocal(cell, maingrid, myproc)
 * -------------------------
 * Determine if the neighbors to a particular cell should be
 * included on the local processor.  This effectively searches
 * over all the neighbors to the cell and determines if any
 * of these cells is on the local processor, if so, this is 
 * a processor boundary cell and should flagged.
 *
 */
static int IsNeighborLocal(int cell, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize)
{
  // if we don't have a good list, exit
  if(listsize < 1) {
    printf("Error!!: List for boundary neighbors needs to be larger!-local\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }
  // select function to call first item in list
  switch(list[0]) {
    case EDGE:
      // search by edge
      return IsNeighborLocalByEdge(cell, maingrid, myproc, &list[1], listsize-1);
      break;
    case NODE:
      // search by node 
      return IsNeighborLocalByNode(cell, maingrid, myproc, &list[1], listsize-1);
      break;
    default:
      printf("Error in IsNeighborLocal\n");
      break;
  }
}


/*
 * Function: IsNeighborLocalByEdge
 * Usage: IsNeighborLocalByEdge(cell, maingrid, myproc)
 * -------------------------
 * Determine if the neighbors to a particular cell should be
 * included on the local processor.  This effectively searches
 * over all the edge neighbors to the cell and determines if any
 * of these cells is on the local processor, if so, this is 
 * a processor boundary cell and should flagged.
 *
 */
static int IsNeighborLocalByEdge(int cell, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize)
{
  // if we don't have a good list, exit
  if(listsize < 0) {
    printf("Error!!: List for boundary neighbors needs to be larger!-localedge\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }
  int nf, nei;
  // will need to generalize this for nodal neighbors
       
  // over each face of the cell
  for(nf=0;nf<maingrid->nfaces[cell];nf++) {
    // get the neighbor
    nei = maingrid->neigh[cell*maingrid->maxfaces+nf];
    // if a neighbor is not a ghost cell
    if(nei != -1) {
      // and also if the processor on the face is on the current processor
      if(maingrid->part[nei]==myproc) {
        // this cell is a neighbor to the local processor
        return 1;
      }
      /* now process other items on the list, if any */

      // if there are more items to be considered
      if(listsize > 0) {

        // select which function to call off the first item in the list
        switch(list[0]) {
          case EDGE:
            // original selection by edge
            // check to see if we have a hit, if so we can quit and return 2
            if(IsNeighborLocalByEdge(nei, maingrid, myproc, &list[1], listsize-1)==1)
              return 1;
            break;
          case NODE:
            // new selection based on nodal neighbors (more inclusive)
            // check to see if we have a hit, if so we can quit and return 2
            if(IsNeighborLocalByNode(nei, maingrid, myproc, &list[1], listsize-1)==1)
              return 1;
            break;
          default:
            printf("Error in IsNeighborLocalByEdge!!\n");
            break;
        }
      }
    }
  }
  // there was no match so this is not a neighbor to the local processor
  return 0;
}


/*
 * Function: IsNeighborLocalByNode
 * Usage: IsNeighborLocalByNode(cell, maingrid, myproc)
 * -------------------------
 * Determine if the neighbors to a particular cell should be
 * included on the local processor.  This effectively searches
 * over all the nodal neighbors to the cell and determines if any
 * of these cells is on the local processor, if so, this is 
 * a processor boundary cell and should flagged.
 *
 */
static int IsNeighborLocalByNode(int cell, gridT *maingrid, int myproc,
    boundaryselection *list, int listsize)
{
  // if we don't have a good list, exit
  if(listsize < 0) {
    printf("Error!!: List for boundary neighbors needs to be larger!-localnode\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }
  int nf, in, nei, nn;
  // will need to generalize this for nodal neighbors

  // over each face of the cell
  for(nf=0;nf<maingrid->nfaces[cell];nf++) {
    // get each node
    in = maingrid->cells[cell*maingrid->maxfaces+nf];
    // loop over all the nodal cell neighbors
    for(nn = 0; nn < maingrid->numpcneighs[in]; nn++) {
      // get the nodal cell neighbors
      nei = maingrid->pcneighs[in][nn];
      // if a neighbor is not a ghost cell
      if(nei != -1) {
        // and also if the processor on the face is on the current processor
        if(maingrid->part[nei]==myproc) {
          // this cell is a neighbor to the local processor
          return 1;
        }
        /* now process other items on the list, if any */

        // if there are more items to be considered
        if(listsize > 0) {

          // select which function to call off the first item in the list
          switch(list[0]) {
            case EDGE:
              // original selection by edge
              // check to see if we have a hit, if so we can quit and return 2
              if(IsNeighborLocalByEdge(nei, maingrid, myproc, &list[1], listsize-1)==1)
                return 1;
              break;
            case NODE:
              // new selection based on nodal neighbors (more inclusive)
              // check to see if we have a hit, if so we can quit and return 2
              if(IsNeighborLocalByNode(nei, maingrid, myproc, &list[1], listsize-1)==1)
                return 1;
              break;
            default:
              printf("Error in IsNeighborLocalByNode!!\n");
              break;
          }
        }
      }
    }
  }
  // there was no match so this is not a neighbor to the local processor
  return 0;
}

/*
 * Function: IsNeighborGlobal
 * Usage: IsNeighborGlobal(cell, maingrid, localgrid, myproc)
 * -------------------------
 * Determine if the neighbors to a particular cell should be
 * included for interprocessor communication.  
 *
 */
static inline int IsNeighborGlobal(int cell, gridT *maingrid, 
    gridT **localgrid, int myproc, boundaryselection *list, int listsize)
{
  // if we don't have a good list, exit
  if(listsize < 1) {
    printf("Error!!: List for boundary neighbors needs to be larger!-global\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }
  // select function to call first item in list
  switch(list[0]) {
    case EDGE:
      // search by edge
      return IsNeighborGlobalByEdge(cell, maingrid, localgrid,
          myproc, list, listsize-1);
      break;
    case NODE:
      // search by node 
      return IsNeighborGlobalByNode(cell, maingrid, localgrid,
          myproc, list, listsize-1);
      break;
    default:
      printf("Error in IsNeighborGlobal\n");
      break;
  }
}


/*
 * Function: IsNeighborGlobalByEdge
 * Usage: IsNeighborGlobalByEdge(cell, maingrid, localgrid, myproc)
 * -------------------------
 * Determine if the neighbors to a particular cell should be
 * included for interprocessor communication.  
 *
 */
static int IsNeighborGlobalByEdge(int cell, gridT *maingrid, 
    gridT **localgrid, int myproc, boundaryselection *list, int listsize)
{
  // if we don't have a good list, exit
  if(listsize < 0) {
    printf("Error!!: List for boundary neighbors needs to be larger!-globaledge\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }
  int nf, nei;
  // will need to generalize this for nodal neighbors
       
  // over each face of the cell
  for(nf=0;nf<(*localgrid)->nfaces[cell];nf++) {
    // get the neighbor
    nei = (*localgrid)->neigh[cell*(*localgrid)->maxfaces+nf];
    // if a neighbor is not a ghost cell
    if(nei != -1) {
      // and also if the processor on the face is on the current processor
      if(maingrid->part[(*localgrid)->mnptr[nei]]==myproc) {
        // this cell is a neighbor to the local processor
        return 1;
      }

      /* now process other items on the list, if any */

      // if there are more items to be considered
      if(listsize > 0) {

        // select which function to call off the first item in the list
        switch(list[1]) {
          case EDGE:
            // original selection by edge
            // check to see if we have a hit, if so we can quit and return 2
            if(IsNeighborGlobalByEdge(nei, maingrid, localgrid, 
                myproc, &list[1], listsize-1)==1)
              return 1;
            break;
          case NODE:
            // new selection based on nodal neighbors (more inclusive)
            // check to see if we have a hit, if so we can quit and return 2
            if(IsNeighborGlobalByNode(nei, maingrid, localgrid, 
                myproc, &list[1], listsize-1)==1)
              return 1;
            break;
          default:
            printf("Error in IsNeighborLocalByNode!!\n");
            break;
        }
      }
    }
  }
  // there was no match so this is not a neighbor to the local processor
  return 0;
}

/*
 * Function: IsNeighborGlobalByNode
 * Usage: IsNeighborGlobalByNode(cell, maingrid, localgrid, myproc)
 * -------------------------
 * Determine if the neighbors to a particular cell should be
 * included for interprocessor communication.  
 *
 */
static int IsNeighborGlobalByNode(int cell, gridT *maingrid, 
    gridT **localgrid, int myproc, boundaryselection *list, int listsize)
{
  // if we don't have a good list, exit
  if(listsize < 0) {
    printf("Error!!: List for boundary neighbors needs to be larger!-globalnode\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);  
  }
  int nf, nei, nn, in;
       
  // over each node of the cell
  for(nf=0;nf<(*localgrid)->nfaces[cell];nf++) {
    // get the node 
    in = (*localgrid)->cells[cell*(*localgrid)->maxfaces+nf];
    // loop over all the nodal neighbors
    for(nn=0; nn < (*localgrid)->numpcneighs[in]; nn++) {
      // get the neighbor
      nei = (*localgrid)->pcneighs[in][nn];
      // if a neighbor is not a ghost cell
      if(nei != -1) {
        // and also if the processor on the face is on the current processor
        if(maingrid->part[(*localgrid)->mnptr[nei]]==myproc) {
          // this cell is a neighbor to the local processor
          return 1;
        }
        /* now process other items on the list, if any */

        // if there are more items to be considered
        if(listsize > 0) {

          // select which function to call off the first item in the list
          switch(list[1]) {
            case EDGE:
              // original selection by edge
              // check to see if we have a hit, if so we can quit and return 2
              if(IsNeighborGlobalByEdge(nei, maingrid, localgrid, 
                    myproc, &list[1], listsize-1)==1)
                return 1;
              break;
            case NODE:
              // new selection based on nodal neighbors (more inclusive)
              // check to see if we have a hit, if so we can quit and return 2
              if(IsNeighborGlobalByNode(nei, maingrid, localgrid, 
                    myproc, &list[1], listsize-1)==1)
                return 1;
              break;
            default:
              printf("Error in IsNeighborLocalByNode!!\n");
              break;
          }
        }
      }
    }
  }
  // there was no match so this is not a neighbor to the local processor
  return 0;
}

/*
 * Function: Geometry
 * Usage: Geometry(maingrid, grid, myproc)
 * -------------------------
 * Compute  geometry information such as cell area (AC), 
 * edge length (df), gradient edge length (dg) between voronoi
 * cell centers [If this is a boundary edge then dg is twice the 
 * distance between the edge and the Voronoi point on the inside 
 * of the boundary.], distange between voronoi cell center and edge
 * (def), cell normals (n1 and n2), location
 * of edge centers (xe, ye) [not neceessarily coincident upon 
 * edge mid-point], and coefficients that make up
 * tangents to compute advection (xi)
 *
 */
static void Geometry(gridT *maingrid, gridT **grid, int myproc)
{
  int n, nf, npc, ne, k, j, Nc=(*grid)->Nc, Ne=(*grid)->Ne, Np=(*grid)->Np, p1, p2,grad1,grad2,enei;
  REAL xt[(*grid)->maxfaces], yt[(*grid)->maxfaces], xc, yc, den, R0, tx, ty, tmag, xdott;
  
  (*grid)->Ac = (REAL *)SunMalloc(Nc*sizeof(REAL),"Geometry");
  (*grid)->df = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->dg = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->def = (REAL *)SunMalloc((*grid)->maxfaces*Nc*sizeof(REAL),"Geometry");
  (*grid)->n1 = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->n2 = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->xe = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->ye = (REAL *)SunMalloc(Ne*sizeof(REAL),"Geometry");
  (*grid)->xi = (REAL *)SunMalloc(2*((*grid)->maxfaces-1)*Ne*sizeof(REAL),"Geometry");

  // add on the total area of all the cells neighboring a point
  (*grid)->Actotal = (REAL **)SunMalloc(Np*sizeof(REAL*),"Geometry");
  // allocate momory over all the cells 
  for(n = 0; n < Np; n++) {
    // allocate memory for each layer
    (*grid)->Actotal[n] = (REAL *)SunMalloc((*grid)->Nkp[n]*sizeof(REAL),"Geometry");
  }

  /* Compute the area of each cell.*/
  if(myproc==0 && VERBOSE>2) printf("\t\tComputing cell areas...\n");
  for(n=0;n<Nc;n++) {
    for(nf=0;nf<(*grid)->nfaces[n];nf++) {
      xt[nf]=maingrid->xp[(*grid)->cells[n*(*grid)->maxfaces+nf]];
      yt[nf]=maingrid->yp[(*grid)->cells[n*(*grid)->maxfaces+nf]];
    }
    // compute area of cell
    (*grid)->Ac[n] = GetArea(xt,yt,(*grid)->nfaces[n]);
  }

  /* Compute the summed area of all the neighboring cells to a node */
  // this is of limited utility in many computations since Actotal will be 
  // different with real world bathymetry
  if(myproc==0 && VERBOSE>2) 
    printf("\t\tComputing neighboring cell areas (Actotal)...\n");
  for(n=0 ; n < Np; n++) {
    for(k=0; k < (*grid)->Nkp[n]; k++) {
      // initialize Actotal
      (*grid)->Actotal[n][k] = 0;
      // sum over all the relevant cell neighbors to get the total area
      for(npc=0; npc < (*grid)->numpcneighs[n] ; npc++) {
        // only consider cell neighbors at this layer or deeper
        if(k < (*grid)->Nk[(*grid)->pcneighs[n][npc]]) {
            // accumulate the area for each of the neighboring cells
            (*grid)->Actotal[n][k] += (*grid)->Ac[(*grid)->pcneighs[n][npc]];
        }
      }
    }
  }
            
  
  /* Compute the length of each edge */
  if(myproc==0 && VERBOSE>2) printf("\t\tComputing edge lengths...\n");
  for(n=0;n<Ne;n++) {
    // compute the length of the edge
    (*grid)->df[n]=
      sqrt(pow(maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n]]-
	       maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n+1]],2)+
	   pow(maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n]]-
	       maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n+1]],2));
  }

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
    //corrected part to make we can calculate xe and ye when grad[2n]==1
    if((*grid)->grad[2*n]>=0){
    xdott = ((*grid)->xv[(*grid)->grad[2*n]]-maingrid->xp[p1])*tx+
      ((*grid)->yv[(*grid)->grad[2*n]]-maingrid->yp[p1])*ty;
    (*grid)->xe[n] = maingrid->xp[p1]+xdott*tx;
    (*grid)->ye[n] = maingrid->yp[p1]+xdott*ty;
    }
    else{
    xdott = ((*grid)->xv[(*grid)->grad[2*n+1]]-maingrid->xp[p1])*tx+
      ((*grid)->yv[(*grid)->grad[2*n+1]]-maingrid->yp[p1])*ty;
    (*grid)->xe[n] = maingrid->xp[p1]+xdott*tx;
    (*grid)->ye[n] = maingrid->yp[p1]+xdott*ty;
    }
    
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
  if(myproc==0 && VERBOSE>2) printf("\t\tComputing n1, n2, and dg..\n");
  // for each edge
  for(n=0;n<Ne;n++) {
    // if both cells are not ghost cells 
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
        printf("Coincident Voronoi points on edge %d (%d,%d)!\n",n,
            (*grid)->grad[2*n],(*grid)->grad[2*n+1]);
      }
    }
    else {
      xc = (*grid)->xe[n];
      yc = (*grid)->ye[n];
      // one neighboring cell is a ghost cell
      if((*grid)->grad[2*n]==-1) {
        (*grid)->n1[n] = xc-(*grid)->xv[(*grid)->grad[2*n+1]];
        (*grid)->n2[n] = yc-(*grid)->yv[(*grid)->grad[2*n+1]];
        (*grid)->dg[n] = sqrt(pow((*grid)->n1[n],2)+
            pow((*grid)->n2[n],2));
        (*grid)->n1[n] = (*grid)->n1[n]/(*grid)->dg[n];
        (*grid)->n2[n] = (*grid)->n2[n]/(*grid)->dg[n];
        (*grid)->dg[n] = 2.0*(*grid)->dg[n];
      } 
      // other cell is ghost cell 
      else {
        (*grid)->n1[n] = (*grid)->xv[(*grid)->grad[2*n]]-xc;
        (*grid)->n2[n] = (*grid)->yv[(*grid)->grad[2*n]]-yc;
        (*grid)->dg[n] = sqrt(pow((*grid)->n1[n],2)+pow((*grid)->n2[n],2));
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
    for(nf=0;nf<(*grid)->nfaces[n];nf++) {
      ne = (*grid)->face[n*(*grid)->maxfaces+nf];
      (*grid)->def[n*(*grid)->maxfaces+nf] = 
        -(((*grid)->xv[n]-maingrid->xp[(*grid)->edges[ne*NUMEDGECOLUMNS]])*(*grid)->n1[ne]+
            ((*grid)->yv[n]-maingrid->yp[(*grid)->edges[ne*NUMEDGECOLUMNS]])*(*grid)->n2[ne])*
        (*grid)->normal[n*(*grid)->maxfaces+nf];

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
  if(myproc==0 && VERBOSE>2) printf("\t\tComputing xi coefficients...\n");
  // initialixe xi to 0 over all edges and faces
  for(n=0;n<Ne;n++) {
    grad1=(*grid)->grad[2*n];
    grad2=(*grid)->grad[2*n+1];
    enei=(*grid)->nfaces[grad1]+(*grid)->nfaces[grad2]-2;
    for(nf=0;nf<enei;nf++) {
      (*grid)->xi[2*((*grid)->maxfaces-1)*n+nf]=0;
    }
    // for each edge
    for(n=0;n<Ne;n++) {
      if((*grid)->n1[n]!=0 && (*grid)->n2[n]!=0) {
        for(nf=0;nf<2;nf++) {
          j = (*grid)->eneigh[2*((*grid)->maxfaces-1)*n+2*nf];
          k = (*grid)->eneigh[2*((*grid)->maxfaces-1)*n+2*nf+1];
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
            (*grid)->xi[2*((*grid)->maxfaces-1)*n+2*nf]=
              ((*grid)->n1[k]*(*grid)->n1[n]+(*grid)->n2[k]*(*grid)->n2[n])/den;
            (*grid)->xi[2*((*grid)->maxfaces-1)*n+2*nf+1]=
              -((*grid)->n1[j]*(*grid)->n1[n]+(*grid)->n2[j]*(*grid)->n2[n])/den;
          }
        }
      }
    }
  }
}

/*
 * Function: GetCircumcircleRadius
 * Usage: GetCircumcircleRadius(xt, yt, Nf)
 * -------------------------
 * Compute the circumcircle radius for a polygon
 * with verticies (xt,yt) and Nf faces
 *
 */
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
  
/*
 * Function: GetArea
 * Usage: GetArea(xt, yt, Nf)
 * --------------------------------
 * Compute the area of a polygon with verticies
 * (xt,yt) that has Nf faces
 *
 */
REAL GetArea(REAL *xt, REAL *yt, int Nf)
{
  REAL b,r1,r2,h,l,xt2[3],yt2[3],area;
  int i;
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
    if(Nf>3) {
      area=GetArea(xt,yt,3);
      for(i=0;i<(Nf-3);i++){        
        xt2[0]=xt[(i+2)];
        xt2[1]=xt[(i+3)];
        xt2[2]=xt[0];
        yt2[0]=yt[(i+2)];
        yt2[1]=yt[(i+3)];
        yt2[2]=yt[0];
        area+=GetArea(xt2,yt2,3);
      }
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
	 &(grid->yv[nstart]),&(grid->dv[nstart]),ncount,grid->maxfaces);

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
    
    if(nc1 != -1 && nc2 != -1 && grid->nfaces[nc1]==3 && grid->nfaces[nc2]==3) {
      xv1 = grid->xv[nc1];
      xv2 = grid->xv[nc2];
      yv1 = grid->yv[nc1];
      yv2 = grid->yv[nc2];
      xc1 = 0;
      xc2 = 0;
      yc1 = 0;
      yc2 = 0;
      for(nf=0;nf<grid->nfaces[nc1];nf++) {
	xc1 += grid->xp[grid->cells[nc1*grid->maxfaces+nf]]/grid->nfaces[nc1];
	yc1 += grid->yp[grid->cells[nc1*grid->maxfaces+nf]]/grid->nfaces[nc1];
      }
      for(nf=0;nf<grid->nfaces[nc2];nf++) {
	xc2 += grid->xp[grid->cells[nc2*grid->maxfaces+nf]]/grid->nfaces[nc2];
	yc2 += grid->yp[grid->cells[nc2*grid->maxfaces+nf]]/grid->nfaces[nc2];
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
	for(nf=0;nf<grid->nfaces[nc2];nf++) {
	  xc2 += grid->xp[grid->cells[nc2*grid->maxfaces+nf]]/grid->nfaces[nc2];
	  yc2 += grid->yp[grid->cells[nc2*grid->maxfaces+nf]]/grid->nfaces[nc2];
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
	for(nf=0;nf<grid->nfaces[nc1];nf++) {
	  xc1 += grid->xp[grid->cells[nc1*grid->maxfaces+nf]]/grid->nfaces[nc1];
	  yc1 += grid->yp[grid->cells[nc1*grid->maxfaces+nf]]/grid->nfaces[nc1];
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

static int CorrectAngles(gridT *grid, int myproc) ///questions
{
  int n, nf, nc1, nc2, numcorr, i,k;
  REAL xc, yc, xv1, xv2, yv1, yv2, xc1, xc2, xc3, yc1, yc2, yc3, xv_sum, yv_sum, dg, dg0, ang1, ang2,
    dot1, dot2, dot3, mag1, mag2, mag3, cosVoronoiRatio;
  REAL VoronoiRatio=MPI_GetValue(DATAFILE,"VoronoiRatio","CorrectVoronoi",0);
  cosVoronoiRatio = cos(VoronoiRatio*PI/180.0);

  numcorr=0;
  for(n=0;n<grid->Nc;n++) {
    xv_sum=0;
    yv_sum=0;
    k=0;
    if(grid->nfaces[n]==3){
      for(i=0;i<grid->nfaces[n]-2;i++){
        xc1 = grid->xp[grid->cells[n*grid->maxfaces+i]];
        xc2 = grid->xp[grid->cells[n*grid->maxfaces+i+1]];
        xc3 = grid->xp[grid->cells[n*grid->maxfaces+i+2]];
        yc1 = grid->yp[grid->cells[n*grid->maxfaces]];
        yc2 = grid->yp[grid->cells[n*grid->maxfaces+i+1]];
        yc3 = grid->yp[grid->cells[n*grid->maxfaces+i+2]];
        mag1 = sqrt(pow(xc2-xc1,2)+pow(yc2-yc1,2));
        mag2 = sqrt(pow(xc3-xc1,2)+pow(yc3-yc1,2));
        mag3 = sqrt(pow(xc3-xc2,2)+pow(yc3-yc2,2));
        dot1 = (mag1*mag1+mag2*mag2-mag3*mag3)/(2*mag1*mag2);
        dot2 = (mag2*mag2+mag3*mag3-mag1*mag1)/(2*mag2*mag3);
        dot3 = (mag3*mag3+mag1*mag1-mag2*mag2)/(2*mag3*mag1);

        if(dot1<=cosVoronoiRatio || dot2<=cosVoronoiRatio || dot3<=cosVoronoiRatio) {
	  k++;
          xv_sum=xv_sum+(xc1+xc2+xc3)/3;
	  yv_sum=yv_sum+(yc1+yc2+yc3)/3;
        }
      }
    }
    if(k>0){
      grid->xv[n]=xv_sum/(grid->nfaces[n]-2);
      grid->yv[n]=yv_sum/(grid->nfaces[n]-2);
      numcorr++;
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
    if(nc1>=0 && nc2>=0) {
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
    // this does not appear to be used in the function below and could be removed
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
    // this does not appear to be used in the function below and could be removed
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

  // for each neighbor
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    // get the neighboring processors
    neighproc = grid->myneighs[neigh];

    // for each of the cell to send between the pair processor
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
  
  

