/*
 * File: grid.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 09/30/02
 * --------------------------------
 * This file contains grid-based functions.
 *
 * $Id: grid.c,v 1.9 2003-04-29 00:15:06 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.8  2003/04/26 14:16:37  fringer
 * Added initialization function ReturnDepth .
 *
 * Revision 1.7  2003/04/22 02:43:42  fringer
 * Changed makepointers() in grid.c so that the ghost points include
 * extra points for the ELM interpolation.  Only the number of neihbors is
 * specified.  This may not guarantee that there are at least 3 neighbors for
 * each cell.
 *
 * Also, had to change send/recvs in 2d data transfer functions because
 * sends/recvs would hang on processor with more than 2 neighbors...
 *
 * Removed code containing jflux terms because they were causing code to
 * crash with certain number of procs.
 *
 * Revision 1.6  2003/04/21 20:25:46  fringer
 * Working version before addition of ghost cells for Kriging.
 *
 * Revision 1.5  2002/11/30 13:44:45  fringer
 * Working version for the simplified Shelf bathymetry.
 *
 * Revision 1.4  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
 * Revision 1.3  2002/11/03 02:09:42  fringer
 * Moved up to field-scale and stable!!!
 *
 * Revision 1.2  2002/11/03 00:51:42  fringer
 * Removed the filename read from ReadFileNames so that it reads the
 * file names from the main datafile defined in suntans.h
 *
 * Revision 1.1  2002/11/03 00:17:18  fringer
 * Initial revision
 *
 *
 */
#include "grid.h"
#include "util.h"
#include "initialization.h"

/*
 * Private function declarations.
 */
static void InitLocalGrid(gridT **grid);
static void VertGrid(gridT *maingrid, gridT **localgrid, MPI_Comm comm);
static void GetGraph(GraphType *graph, gridT *grid, MPI_Comm comm);
void Topology(gridT **maingrid, gridT **localgrid, int myproc, int numprocs);
static int GetNumCells(gridT *grid, int proc);
static void TransferData(gridT *maingrid, gridT **localgrid, int myproc);
static int GetNumEdges(gridT *grid);
static void Geometry(gridT *maingrid, gridT **grid, int myproc);
static REAL GetArea(REAL *xt, REAL *yt, int Nf);
static void EdgeMarkers(gridT *maingrid, gridT **localgrid, int myproc);
static void ReOrder0(gridT *maingrid, gridT **localgrid, int myproc);
static void ReOrder(gridT *grid);
static int IsCellNeighborProc(int nc, gridT *maingrid, gridT *localgrid, 
			      int myproc, int neighproc);
static int IsEdgeNeighborProc(int ne, gridT *maingrid, gridT *localgrid, 
			      int myproc, int neighproc);
static void MakePointers(gridT *maingrid, gridT **localgrid, int myproc, MPI_Comm comm);
static void MakePointers0(gridT *maingrid, gridT **localgrid, int myproc);
static void ResortBoundaries(gridT *localgrid, int myproc);
static void InterpDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static void FreeGrid(gridT *grid, int numprocs);
static void OutputData(gridT *maingrid, gridT *grid, int myproc);
static void CreateFaceArray(int *grad, int *neigh, int *face, int Nc, int Ne);
static void CreateNormalArray(int *grad, int *face, int *normal, int Nc);
void CorrectVoronoi(gridT *grid);
static void CreateNearestPointers(gridT *maingrid, gridT *localgrid, int myproc);

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
  int Np, Ne, Nc;
  gridT *maingrid;

  ReadFileNames(myproc);

  if(!TRIANGULATE) {
    Np = getsize(POINTSFILE);
    Ne = getsize(EDGEFILE);
    Nc = getsize(CELLSFILE);

    // Every processor will know about data read in from
    // triangle as well as the depth.
    if(myproc==0 && VERBOSE>0) printf("Initializing Main Grid...\n");
    InitMainGrid(&maingrid,Np,Ne,Nc);

    if(myproc==0 && VERBOSE>0) printf("Reading Grid...\n");
    ReadMainGrid(maingrid);
  } else {
    if(myproc==0 && VERBOSE>0) printf("Triangulating the point set...\n");
    GetTriangulation(&maingrid,myproc);
  }

  // Set the voronoi points to be the centroids of the
  // cells if the CORRECTVORONOI flag is 1.
  if((int)MPI_GetValue(DATAFILE,"CorrectVoronoi","ReadMainGrid",0))
    CorrectVoronoi(maingrid);

  if(myproc==0 && VERBOSE>0) printf("Getting the depth for graph weights...\n");
  GetDepth(maingrid,myproc,numprocs,comm);

  if(myproc==0 && VERBOSE>0) printf("Computing Graph...\n");
  CreateCellGraph(maingrid);

  if(myproc==0 && VERBOSE>0) printf("Computing Connectivity...\n");
  Connectivity(maingrid,myproc);

  if(myproc==0 && VERBOSE>0) printf("Partitioning...\n");
  Partition(maingrid,localgrid,comm);

  //  CheckCommunicateCells(maingrid,*localgrid,myproc,comm);
  //  CheckCommunicateEdges(maingrid,*localgrid,myproc,comm);
  //  SendRecvCellData2D((*localgrid)->dv,*localgrid,myproc,comm);

  OutputData(maingrid,*localgrid,myproc);
  FreeGrid(maingrid,numprocs);
}

void Partition(gridT *maingrid, gridT **localgrid, MPI_Comm comm)
{
  int j, n, numflag=0, wgtflag=0, options[5], edgecut;
  int myproc, numprocs, proc;
  MPI_Status status;
  GraphType graph;

  InitLocalGrid(localgrid);

  MPI_Comm_size(comm, &numprocs);
  MPI_Comm_rank(comm, &myproc);

  options[0] = 0;
  wgtflag = 2;
  numflag = 0;

  GetGraph(&graph,maingrid,comm);

  (*localgrid)->part = idxmalloc(graph.nvtxs, "TestParMetis: part");

  /*
   * Partition the graph and create the part array.
   */
  if(myproc==0 && VERBOSE>2) printf("Partitioning with ParMETIS_PartKway...\n");
  ParMETIS_PartKway(graph.vtxdist,graph.xadj,graph.adjncy,graph.vwgt,NULL,
		    &wgtflag,&numflag,&numprocs,options,&edgecut,
		    (*localgrid)->part,&comm);

  if(myproc==0 && VERBOSE>2) printf("Redistributing the partition arrays...\n");
  if(myproc!=0) 
    MPI_Send((void *)(*localgrid)->part,graph.nvtxs,MPI_INT,0,1,comm); 
  else {
    for(j=0;j<graph.nvtxs;j++)
      maingrid->part[graph.vtxdist[myproc]+j]=(*localgrid)->part[j];
    for(proc=1;proc<numprocs;proc++) 
      MPI_Recv((void *)(maingrid->part+graph.vtxdist[proc]),
	       graph.vtxdist[proc+1]-graph.vtxdist[proc],MPI_INT,proc,1,comm,&status);
  }
  MPI_Bcast((void *)maingrid->part,maingrid->Nc,MPI_INT,0,comm);

  Topology(&maingrid,localgrid,myproc,numprocs);

  if(myproc==0 && VERBOSE>2) printf("Creating nearest neighbor arrays...\n");
  CreateNearestPointers(maingrid,*localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("Transferring data...\n");
  TransferData(maingrid,localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("Creating edge markers...\n");
  EdgeMarkers(maingrid,localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("Computing edge and voronoi distances and areas...\n");
  Geometry(maingrid,localgrid,myproc);

  if(myproc==0 && VERBOSE>2) printf("Vert grid...\n");
  VertGrid(maingrid,localgrid,comm);

  if(myproc==0 && VERBOSE>2) printf("Reordering...\n");
  //  ReOrder0(maingrid,localgrid,myproc);
  //  ReOrder(*localgrid);

  if(myproc==0 && VERBOSE>2) printf("Making pointers...\n");
  MakePointers(maingrid,localgrid,myproc,comm);

  // NEED THIS!
  //  ResortBoundaries(*localgrid,myproc);

  if(VERBOSE>3) {
    ReportConnectivity(*localgrid,maingrid,myproc);
    ReportPartition(maingrid,*localgrid,myproc,comm);
  }
}

static void CreateNearestPointers(gridT *maingrid, gridT *localgrid, int myproc) {
  int i, j, k, m, Ncpart, Nclocal;
  int Nnearestcells=(int)MPI_GetValue(DATAFILE,"Nnearestcells","CreateNearestPointers",0);
  int Nnearestedges=(int)MPI_GetValue(DATAFILE,"Nnearestedges","CreateNearestPointers",0);
  int *lcptr = (int *)malloc(maingrid->Nc*sizeof(int));

  unsigned short *found=(unsigned short *)
    malloc(maingrid->Nc*sizeof(unsigned short));
  
  maingrid->Nnearestcells=Nnearestcells;
  maingrid->Nnearestedges=Nnearestedges;

  Ncpart=0;
  for(i=0;i<maingrid->Nc;i++) {
    found[i]=0;
    if(maingrid->part[i]==myproc) 
      Ncpart++;
  }

  localgrid->nearestcells=(int **)malloc(Ncpart*sizeof(int *));
  localgrid->nearestedges=(int **)malloc(Ncpart*sizeof(int *));
  for(i=0;i<Ncpart;i++)
    localgrid->nearestcells[i]=(int *)malloc(Nnearestcells*sizeof(int));
  for(i=0;i<Ncpart;i++)
    localgrid->nearestedges[i]=(int *)malloc(Nnearestedges*sizeof(int));
  
  k=0;
  Nclocal=0;
  for(i=0;i<maingrid->Nc;i++) {
    if(maingrid->part[i]==myproc) {
      FindNearest(localgrid->nearestcells[k],maingrid->xv,maingrid->yv,
      		  maingrid->Nc,Nnearestcells,maingrid->xv[i],maingrid->yv[i]);
      for(j=0;j<Nnearestcells;j++)
	if(!found[localgrid->nearestcells[k][j]]) {
	  found[localgrid->nearestcells[k][j]]=1;
	  Nclocal++;
	}
      k++;
    }
  }
  localgrid->Nc = Nclocal;

  localgrid->mnptr = (int *)malloc(localgrid->Nc*sizeof(int));
  for(i=0;i<maingrid->Nc;i++)
    found[i]=0;
  
  m=0;
  for(i=0;i<maingrid->Nc;i++) 
    if(maingrid->part[i]==myproc) {
      found[i]=1;
      lcptr[i]=m;
      localgrid->mnptr[m++]=i;
    }
  k=0;
  for(i=0;i<maingrid->Nc;i++) 
    if(maingrid->part[i]==myproc) {
      for(j=0;j<Nnearestcells;j++)
	if(!found[localgrid->nearestcells[k][j]]) {
	  found[localgrid->nearestcells[k][j]]=1;
	  lcptr[localgrid->nearestcells[k][j]]=m;
	  localgrid->mnptr[m++]=localgrid->nearestcells[k][j];
	}
      k++;
    }

  for(i=0;i<Ncpart;i++) 
    for(j=0;j<Nnearestcells;j++) 
      localgrid->nearestcells[i][j]=lcptr[localgrid->nearestcells[i][j]];

  /*
  for(i=0;i<Ncpart;i++) {
    printf("Cell %d: ",i);
    for(j=0;j<Nnearestcells;j++)
      printf("%d ",localgrid->nearestcells[i][j]);
    printf("\n");
  }
  */

  //  free(localgrid->mnptr);
  free(lcptr);
  free(found);
}

void SendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm)
{
  int n, neigh, neighproc, receivesize, sendsize, noncontig, flag;
  int senstart, senend, recstart, recend;
  REAL **recv, **send;
  MPI_Status status;

  recv = (REAL **)malloc(grid->Nneighs*sizeof(REAL *));
  send = (REAL **)malloc(grid->Nneighs*sizeof(REAL *));

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    send[neigh] = (REAL *)malloc(grid->num_cells_send[neigh]*sizeof(REAL));
    recv[neigh] = (REAL *)malloc(grid->num_cells_recv[neigh]*sizeof(REAL));

    for(n=0;n<grid->num_cells_send[neigh];n++)
      send[neigh][n]=celldata[grid->cell_send[neigh][n]];

    MPI_Send((void *)(send[neigh]),grid->num_cells_send[neigh],
	     MPI_DOUBLE,neighproc,1,comm); 
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Recv((void *)(recv[neigh]),grid->num_cells_recv[neigh],
	     MPI_DOUBLE,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    for(n=0;n<grid->num_cells_recv[neigh];n++)
      celldata[grid->cell_recv[neigh][n]]=recv[neigh][n];
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    free(send[neigh]);
    free(recv[neigh]);
  }
  free(send);
  free(recv);
}

void SendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm)
{
  int k, n, nstart, neigh, neighproc, receivesize, sendsize, noncontig, flag;
  int senstart, senend, recstart, recend, *Nsend, *Nrecv;
  REAL **recv, **send;
  MPI_Status status;

  recv = (REAL **)malloc(grid->Nneighs*sizeof(REAL *));
  send = (REAL **)malloc(grid->Nneighs*sizeof(REAL *));
  Nsend = (int *)malloc(grid->Nneighs*sizeof(int));
  Nrecv = (int *)malloc(grid->Nneighs*sizeof(int));

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];

    Nsend[neigh] = 0;
    Nrecv[neigh] = 0;
    for(n=0;n<grid->num_cells_send[neigh];n++) 
      Nsend[neigh]+=grid->Nk[grid->cell_send[neigh][n]];
    for(n=0;n<grid->num_cells_recv[neigh];n++) 
      Nrecv[neigh]+=grid->Nk[grid->cell_recv[neigh][n]];

    send[neigh] = (REAL *)malloc(Nsend[neigh]*sizeof(REAL));
    recv[neigh] = (REAL *)malloc(Nrecv[neigh]*sizeof(REAL));

    nstart=0;
    for(n=0;n<grid->num_cells_send[neigh];n++) {
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
    for(n=0;n<grid->num_cells_recv[neigh];n++) {
      for(k=0;k<grid->Nk[grid->cell_recv[neigh][n]];k++) 
	celldata[grid->cell_recv[neigh][n]][k]=recv[neigh][nstart+k];
      nstart+=grid->Nk[grid->cell_recv[neigh][n]];
    }
  }

  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    free(send[neigh]);
    free(recv[neigh]);
  }
  free(send);
  free(recv);
  free(Nsend);
  free(Nrecv);
}

void SendRecvCellData0(REAL *celldata, gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm)
{
  int n, neigh, neighproc, receivesize, sendsize, noncontig, flag;
  int senstart, senend, recstart, recend;
  REAL **recv, **send;
  MPI_Status status;

  recv = (REAL **)malloc(maingrid->numneighs[myproc]*sizeof(REAL *));
  send = (REAL **)malloc(maingrid->numneighs[myproc]*sizeof(REAL *));

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];

    send[neigh] = (REAL *)malloc(localgrid->num_cells_send[neigh]*sizeof(REAL));
    recv[neigh] = (REAL *)malloc(localgrid->num_cells_recv[neigh]*sizeof(REAL));

    for(n=0;n<localgrid->num_cells_send[neigh];n++)
      send[neigh][n]=celldata[localgrid->cell_send[neigh][n]];

    MPI_Send((void *)(send[neigh]),localgrid->num_cells_send[neigh],MPI_DOUBLE,neighproc,1,comm); 
    MPI_Recv((void *)(recv[neigh]),localgrid->num_cells_recv[neigh],MPI_DOUBLE,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    for(n=0;n<localgrid->num_cells_recv[neigh];n++)
      celldata[localgrid->cell_recv[neigh][n]]=recv[neigh][n];
  }

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    free(send[neigh]);
    free(recv[neigh]);
  }
  free(send);
  free(recv);
}

void CheckCommunicateCells(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm)
{
  int n, neigh, neighproc, receivesize, sendsize, noncontig, flag;
  int senstart, senend, recstart, recend;
  int **recv, **send;
  MPI_Status status;

  recv = (int **)malloc(maingrid->numneighs[myproc]*sizeof(int *));
  send = (int **)malloc(maingrid->numneighs[myproc]*sizeof(int *));

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];
    send[neigh] = (int *)malloc(localgrid->num_cells_send[neigh]*sizeof(int));
    recv[neigh] = (int *)malloc(localgrid->num_cells_recv[neigh]*sizeof(int));
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
    free(send[neigh]);
    free(recv[neigh]);
  }
  free(send);
  free(recv);
}

void CheckCommunicateEdges(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm)
{
  int j,n, neigh, neighproc, receivesize, sendsize, noncontig, flag;
  int senstart, senend, recstart, recend;
  int **recv, **send;
  MPI_Status status;

  recv = (int **)malloc(maingrid->numneighs[myproc]*sizeof(int *));
  send = (int **)malloc(maingrid->numneighs[myproc]*sizeof(int *));

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];

    send[neigh] = (int *)malloc(localgrid->num_edges_send[neigh]*sizeof(int));
    recv[neigh] = (int *)malloc(localgrid->num_edges_recv[neigh]*sizeof(int));

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
    free(send[neigh]);
    free(recv[neigh]);
  }
  free(send);
  free(recv);
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
  *grid = (gridT *)malloc(sizeof(gridT));
  
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
  (*grid)->xp = (REAL *)malloc((*grid)->Np*sizeof(REAL));
  (*grid)->yp = (REAL *)malloc((*grid)->Np*sizeof(REAL));
  // (x,y) coordinates of voronoi points
  (*grid)->xv = (REAL *)malloc((*grid)->Nc*sizeof(REAL));
  (*grid)->yv = (REAL *)malloc((*grid)->Nc*sizeof(REAL));

  // Pointers to xp,yp coordinates that define two endpoints of faces (0<edges<Np)
  (*grid)->edges = (int *)malloc(NUMEDGECOLUMNS*(*grid)->Ne*sizeof(int));
  // Pointers to xv,yv coordinates that define two endpoints of voronoi edges (0<grad<Np)
  (*grid)->grad = (int *)malloc(2*(*grid)->Ne*sizeof(int));
  // Pointers to xp,yp coordinates of vertices that make up polygons (0<cells<Np)
  (*grid)->cells = (int *)malloc(NFACES*(*grid)->Nc*sizeof(int));
  // Pointers to neighboring cells (0<neigh<Nc)
  (*grid)->neigh = (int *)malloc(NFACES*(*grid)->Nc*sizeof(int));
  // Dot product of unique normal with outward normal
  (*grid)->normal = (int *)malloc(NFACES*(*grid)->Ne*sizeof(int));
  // Indices of voronoi edges to cells
  (*grid)->grad = (int *)malloc(2*(*grid)->Ne*sizeof(int));
  // Indices of pointers to faces of each cell
  (*grid)->face = (int *)malloc(NFACES*(*grid)->Nc*sizeof(int));
  // Indices to edges for momentum control volume
  (*grid)->eneigh = (int *)malloc(2*(NFACES-1)*(*grid)->Ne*sizeof(int));

  // Depth at Voronoi points
  (*grid)->dv = (REAL *)malloc((*grid)->Nc*sizeof(REAL));
  // Weights for partitioning cell graph
  (*grid)->vwgt = (int *)malloc((*grid)->Nc*sizeof(int));

  // Assigned cell partition
  (*grid)->part = (int *)malloc((*grid)->Nc*sizeof(int));
  // For ordering and making sure boundary edge data is contiguous.
  (*grid)->order = (int *)malloc((*grid)->Ne*sizeof(int));
  // Edge markers
  (*grid)->mark = (int *)malloc((*grid)->Ne*sizeof(int));
  // Stores the indices to the start and end points of the adjncy array
  (*grid)->xadj = (int *)malloc(((*grid)->Nc+1)*sizeof(int));
  (*grid)->vtxdist = (int *)malloc(100*sizeof(int));
}

/*
 * Function: ReadMainGrid
 * Usage: ReadmainGrid(grid);
 * --------------------------
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
void ReadMainGrid(gridT *grid)
{
  int j, n, nei, nf, numprocs;
  char str[BUFFERLENGTH];
  FILE *ifile, *fid;

  ifile = fopen(POINTSFILE,"r");
  for(n=0;n<grid->Np;n++) {
    grid->xp[n]=getfield(ifile,str);
    grid->yp[n]=getfield(ifile,str);
    getfield(ifile,str);
  }
  fclose(ifile);
  
  ifile = fopen(EDGEFILE,"r");
  for(n=0;n<grid->Ne;n++) {
    for(j=0;j<NUMEDGECOLUMNS-1;j++) 
      grid->edges[NUMEDGECOLUMNS*n+j]=(int)getfield(ifile,str);
    grid->mark[n]=(int)getfield(ifile,str);
    for(j=0;j<2;j++) 
      grid->grad[2*n+j]=(int)getfield(ifile,str);
  }
  fclose(ifile);

  ifile = fopen(CELLSFILE,"r");
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
  char str[BUFFERLENGTH];

  MPI_GetString(POINTSFILE,DATAFILE,"points","OpenFiles",myproc);
  MPI_GetString(EDGEFILE,DATAFILE,"edges","OpenFiles",myproc);
  MPI_GetString(CELLSFILE,DATAFILE,"cells","OpenFiles",myproc);
  MPI_GetString(INPUTDEPTHFILE,DATAFILE,"depth","OpenFiles",myproc);
  MPI_GetString(CELLCENTEREDFILE,DATAFILE,"celldata","OpenFiles",myproc);
  MPI_GetString(EDGECENTEREDFILE,DATAFILE,"edgedata","OpenFiles",myproc);
  MPI_GetString(VERTSPACEFILE,DATAFILE,"vertspace","OpenFiles",myproc);
  MPI_GetString(TOPOLOGYFILE,DATAFILE,"topology","OpenFiles",myproc);
}

void GetDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm)
{
  int n, maxgridweight=100, IntDepth, Nkmax;
  REAL mindepth, maxdepth;

  Nkmax = MPI_GetValue(DATAFILE,"Nkmax","GetDepth",myproc);

  maxdepth=0.0;
  mindepth=INFTY;
  IntDepth=(int)MPI_GetValue(DATAFILE,"IntDepth","GetDepth",myproc);

  if(IntDepth) 
    InterpDepth(grid,myproc,numprocs,comm);
  else {
    for(n=0;n<grid->Nc;n++) {
      grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
    }
  }
  for(n=0;n<grid->Nc;n++) {
    if(grid->dv[n]>maxdepth)
      maxdepth = grid->dv[n];
    if(grid->dv[n]<mindepth)
      mindepth = grid->dv[n];
  }
    
  if(Nkmax>1) {
    if(mindepth!=maxdepth) 
      for(n=0;n<grid->Nc;n++) 
	grid->vwgt[n] = (int)(maxgridweight*(grid->dv[n]-mindepth)/(maxdepth-mindepth));
    else
      for(n=0;n<grid->Nc;n++) 
	grid->vwgt[n] = maxgridweight;
  } else {
      for(n=0;n<grid->Nc;n++) 
	grid->vwgt[n] = Nkmax;
  }

  if(VERBOSE>3) 
    for(n=0;n<grid->Nc;n++) 
      printf("grid->vwgt[%d]=%d\n",n,grid->vwgt[n]);
}

void CreateCellGraph(gridT *grid)
{
  int n, j, Nge;

  Nge = 0;
  for(n=0;n<grid->Nc;n++) 
    for(j=0;j<NFACES;j++)
      if(grid->neigh[NFACES*n+j]!=-1) Nge++;
  grid->Nge = Nge/2;
  grid->adjncy = (int *)malloc(Nge*sizeof(int));

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
  grid->adjncy = (int *)malloc(Nge*sizeof(int));

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

static void CreateFaceArray(int *grad, int *neigh, int *face, int Nc, int Ne)
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
	if(nb==nc2)
	  face[nc1*NFACES+nf]=n;
	nb = neigh[nc2*NFACES+nf];
	if(nb==nc1)
	  face[nc2*NFACES+nf]=n;
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
  int j, n, nf, ng, ne, neigh, nc, nc1, nc2, Nge, pcnt, *faceind;

  /* Create the face array, which contains indexes to the edges of
     each cell */
  if(myproc==0 && VERBOSE>2) printf("Creating face array and normals\n");
  CreateFaceArray(grid->grad,grid->neigh,grid->face,grid->Nc,grid->Ne);
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
  int j, nf, nei, ipart[NFACES], n1, mp;

  if(mgptr >= maingrid->Nc || mgptr < 0) {
    printf("Error in IsBoundaryCell: index out of bounds!\n");
    return -1;
  }
  if(maingrid->part[mgptr]==myproc) {
    for(nf=0;nf<NFACES;nf++) {
      nei = maingrid->neigh[mgptr*NFACES+nf];
      if(nei != -1)
	if(maingrid->part[nei] != myproc)
	  return 2;
    }
    for(nf=0;nf<NFACES;nf++) {
      if(maingrid->mark[maingrid->face[mgptr*NFACES+nf]]==3)
	return 1;
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
static void OutputData(gridT *maingrid, gridT *grid, int myproc)
{
  int j, n, nf, neigh, Np=maingrid->Np, Nc=grid->Nc, Ne=grid->Ne;
  char str[BUFFERLENGTH];
  FILE *ofile;

  if(TRIANGULATE && myproc==0) {
    if(myproc==0 && VERBOSE>2) printf("Outputting Delaunay points...\n");
    sprintf(str,"%s",POINTSFILE);
    ofile = fopen(str,"w");
    for(j=0;j<Np;j++)
      fprintf(ofile,"%f %f 0\n",maingrid->xp[j],maingrid->yp[j]);
    fclose(ofile);
  }

  if(myproc==0 && VERBOSE>2) printf("Outputting cells.dat...\n");
  sprintf(str,"%s.%d",CELLSFILE,myproc);
  ofile = fopen(str,"w");
  for(j=0;j<grid->Nc;j++) {
    fprintf(ofile,"%f %f ",grid->xv[j],grid->yv[j]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->cells[j*NFACES+nf]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->neigh[j*NFACES+nf]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  if(myproc==0 && VERBOSE>2) printf("Outputting edges.dat...\n");
  sprintf(str,"%s.%d",EDGEFILE,myproc);
  ofile = fopen(str,"w");
  for(j=0;j<grid->Ne;j++) {
    for(nf=0;nf<2;nf++)
      fprintf(ofile,"%d ",grid->edges[j*NUMEDGECOLUMNS+nf]);
    fprintf(ofile,"%d ",grid->mark[j]);
    for(nf=0;nf<2;nf++)
      fprintf(ofile,"%d ",grid->grad[2*j+nf]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  if(myproc==0 && VERBOSE>2) printf("Outputting celldata.dat...\n");
  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  ofile = fopen(str,"w");
  for(n=0;n<Nc;n++) {
    fprintf(ofile,"%f %f %f %f %d ",grid->xv[n],grid->yv[n],
	    grid->Ac[n],grid->dv[n],grid->Nk[n]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->face[NFACES*n+nf]);
    for(nf=0;nf<NFACES;nf++)
      fprintf(ofile,"%d ",grid->neigh[NFACES*n+nf]);
    for(nf=0;nf<NFACES;nf++) 
      fprintf(ofile,"%d ",grid->normal[NFACES*n+nf]);
    fprintf(ofile,"%d ",grid->mnptr[n]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  if(myproc==0 && VERBOSE>2) printf("Outputting edgedata.dat...\n");
  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  ofile = fopen(str,"w");
  for(n=0;n<Ne;n++) {
    fprintf(ofile,"%f %f 0 0 %f %f %d %d %d %d ",
	    grid->df[n],grid->dg[n],grid->n1[n],
	    grid->n2[n],grid->Nke[n],grid->grad[2*n],grid->grad[2*n+1],
	    grid->mark[n]);
    for(nf=0;nf<2*(NFACES-1);nf++)
      fprintf(ofile,"%f ",grid->xi[2*(NFACES-1)*n+nf]);
    for(nf=0;nf<2*(NFACES-1);nf++)
      fprintf(ofile,"%d ",grid->eneigh[2*(NFACES-1)*n+nf]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  if(myproc==0 && VERBOSE>2) printf("Outputting topology and boundary pointers...\n");
  sprintf(str,"%s.%d",TOPOLOGYFILE,myproc);
  ofile = fopen(str,"w");
  fprintf(ofile,"%d\n",grid->Nneighs);
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
  *grid = (gridT *)malloc(sizeof(gridT));
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
  int proc, n;

  free(grid->xp);
  free(grid->yp);
  free(grid->xv);
  free(grid->yv);
  free(grid->edges);
  free(grid->cells);
  free(grid->neigh);
  free(grid->part);
  free(grid->Nk);
  free(grid->dz);
  free(grid->face);
  free(grid->grad);
  free(grid->mark);
  free(grid->normal);
  free(grid->xadj);
  free(grid->vtxdist);

  for(proc=0;proc<numprocs;proc++) 
    free(grid->neighs[proc]);
  free(grid->neighs);
  free(grid->numneighs);

  free(grid);
}

static void VertGrid(gridT *maingrid, gridT **localgrid, MPI_Comm comm)
{
  int i, j, k, n, ne, myproc, numprocs, status, vertgridcorrect;
  REAL dz0, dmin, dmax, z, dzsmall;

  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myproc);

  
  maingrid->Nkmax = MPI_GetValue(DATAFILE,"Nkmax","VertGrid",myproc);
  vertgridcorrect = MPI_GetValue(DATAFILE,"vertgridcorrect","VertGrid",myproc);
  dzsmall = MPI_GetValue(DATAFILE,"dzsmall","VertGrid",myproc);

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
    maingrid->dz = (REAL *)malloc(maingrid->Nkmax*sizeof(REAL));

    for(k=0;k<maingrid->Nkmax;k++)
      maingrid->dz[k] = dz0;
  }
  MPI_Bcast(&(maingrid->Nkmax),1,MPI_INT,0,comm);
  if(myproc!=0)
    maingrid->dz = (REAL *)malloc(maingrid->Nkmax*sizeof(REAL));
  MPI_Bcast((void *)maingrid->dz,maingrid->Nkmax,MPI_DOUBLE,0,comm);
  MPI_Bcast(&dz0,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(&dmax,1,MPI_DOUBLE,0,comm);

  (*localgrid)->Nkmax = maingrid->Nkmax;
  (*localgrid)->dz = (REAL *)malloc((*localgrid)->Nkmax*sizeof(REAL));
  for(k=0;k<(*localgrid)->Nkmax;k++)
    (*localgrid)->dz[k]=maingrid->dz[k];
  
  maingrid->Nk = (int *)malloc(maingrid->Nc*sizeof(int));
  (*localgrid)->dztop = (REAL *)malloc((*localgrid)->Nc*sizeof(REAL));
  (*localgrid)->ctop = (int *)malloc((*localgrid)->Nc*sizeof(int));
  (*localgrid)->ctopold = (int *)malloc((*localgrid)->Nc*sizeof(int));
  (*localgrid)->etop = (int *)malloc((*localgrid)->Ne*sizeof(int));
  (*localgrid)->etopold = (int *)malloc((*localgrid)->Ne*sizeof(int));
  (*localgrid)->Nk = (int *)malloc((*localgrid)->Nc*sizeof(int));
  (*localgrid)->Nke = (int *)malloc((*localgrid)->Ne*sizeof(int));
  (*localgrid)->Nkc = (int *)malloc((*localgrid)->Ne*sizeof(int));

  // Adjust depth if bottom cell thickness is too small.  This is
  // working only for constant vertical grid spacings!!!
  // This only needs to be done if an explicit method is used for
  // vertical transport.
  /*
  for(i=0;i<maingrid->Nc;i++) 
    if(fmod(maingrid->dv[i],dz0)/dz0<dzsmall &&
       fmod(maingrid->dv[i],dz0)/dz0!=0) 
      maingrid->dv[i]-=2.0*fmod(maingrid->dv[i],dz0);
  for(i=0;i<(*localgrid)->Nc;i++)
    if(fmod((*localgrid)->dv[i],dz0)/dz0<dzsmall &&
       fmod((*localgrid)->dv[i],dz0)/dz0!=0) {
      if(WARNING) {
	printf("Warning!  Proc %d Vertical spacing at bottom boundary of cell %d too small!\n",
	       myproc,i);
	printf("Adjusting depth from %f to %f.\n",(*localgrid)->dv[i],
	       (*localgrid)->dv[i]-2.0*fmod((*localgrid)->dv[i],dz0));
      }
      (*localgrid)->dv[i]-=2.0*fmod((*localgrid)->dv[i],dz0);
    }
  */
  if((*localgrid)->Nkmax>1) {
    for(i=0;i<(*localgrid)->Nc;i++) 
      if((*localgrid)->dv[i]==dmax) 
	(*localgrid)->Nk[i] = (*localgrid)->Nkmax;
      else {
	(*localgrid)->Nk[i] = ceil((*localgrid)->dv[i]/dz0);
	//if(fabs((*localgrid)->Nk[i]*dz0 - (*localgrid)->dv[i])>=0.5*dz0)
	//	  (*localgrid)->Nk[i]--;
      }
    for(i=0;i<maingrid->Nc;i++) 
      if(maingrid->dv[i]==dmax) 
	maingrid->Nk[i] = maingrid->Nkmax;
      else {
	maingrid->Nk[i] = ceil(maingrid->dv[i]/dz0);
	//	if(fabs(maingrid->Nk[i]*dz0 - maingrid->dv[i])>=0.5*dz0) 
	//	  maingrid->Nk[i]--;
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
 * Function: GetGraph
 * Usage GetGraph(graph,grid,comm);
 * --------------------------------
 * This code was adapted from the ParMetis-2.0 code in the ParMetis
 * distribution.
 *
 */
static void GetGraph(GraphType *graph, gridT *grid, MPI_Comm comm)
{
  int i, j, k, l, numprocs, myproc;
  int nvtxs, nedges, penum, snedges, snvtxs;
  idxtype *gxadj, *gadjncy, *gvwgt;  
  idxtype *vtxdist, *sxadj, *ssize, *svwgt;
  MPI_Status status;

  MPI_Comm_size(comm, &numprocs);
  MPI_Comm_rank(comm, &myproc);

  vtxdist = graph->vtxdist = idxsmalloc(numprocs+1, 0, "ReadGraph: vtxdist");

  if (myproc == 0) {
    ssize = idxsmalloc(numprocs, 0, "ReadGraph: ssize");

    nvtxs = grid->Nc;
    gxadj = grid->xadj;
    gadjncy = grid->adjncy;
    gvwgt = (idxtype *)grid->vwgt;

    /* Construct vtxdist and send it to all the processors */
    vtxdist[0] = 0;
    for (i=0,k=nvtxs; i<numprocs; i++) {
      l = k/(numprocs-i);
      vtxdist[i+1] = vtxdist[i]+l;
      k -= l;
    }
  }

  MPI_Bcast((void *)vtxdist, numprocs+1, IDX_DATATYPE, 0, comm);

  graph->gnvtxs = vtxdist[numprocs];
  graph->nvtxs = vtxdist[myproc+1]-vtxdist[myproc];
  graph->xadj = idxmalloc(graph->nvtxs+1, "ReadGraph: xadj");
  graph->vwgt = idxmalloc(graph->nvtxs, "ReadGraph: vwgt");

  if (myproc == 0) {
    for (penum=0; penum<numprocs; penum++) {
      snvtxs = vtxdist[penum+1]-vtxdist[penum];
      sxadj = idxmalloc(snvtxs+1, "ReadGraph: sxadj");
      svwgt = idxmalloc(snvtxs, "ReadGraph: svwgt");

      idxcopy(snvtxs+1, gxadj+vtxdist[penum], sxadj);
      idxcopy(snvtxs, gvwgt+vtxdist[penum], svwgt);

      if(VERBOSE>3) 
	for(k=0;k<snvtxs;k++)
	  printf("svwgt[%d]=%d\n",k,svwgt[k]);

      for (i=snvtxs; i>=0; i--)
        sxadj[i] -= sxadj[0];

      ssize[penum] = gxadj[vtxdist[penum+1]] - gxadj[vtxdist[penum]];

      if (penum == myproc) {
        idxcopy(snvtxs+1, sxadj, graph->xadj);
	idxcopy(snvtxs, svwgt, graph->vwgt);
      } else {
        MPI_Send((void *)sxadj, snvtxs+1, IDX_DATATYPE, penum, 1, comm); 
        MPI_Send((void *)svwgt, snvtxs, IDX_DATATYPE, penum, 1, comm); 
      }

      free(sxadj);
      free(svwgt);
    }
  }
  else {
    MPI_Recv((void *)graph->xadj, graph->nvtxs+1, IDX_DATATYPE, 0, 1, comm, &status);
    MPI_Recv((void *)graph->vwgt, graph->nvtxs, IDX_DATATYPE, 0, 1, comm, &status);
  }
  if(VERBOSE>3) {
    printf("Weights on each processor after MPI_Recv\n");
    for(k=0;k<graph->nvtxs;k++)
      printf("vwgt[%d]=%d\n",k,graph->vwgt[k]);
  }

  graph->nedges = graph->xadj[graph->nvtxs];
  graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: graph->adjncy");

  if (myproc == 0) {
    for (penum=0; penum<numprocs; penum++) {
      if (penum == myproc) 
        idxcopy(ssize[penum], gadjncy+gxadj[vtxdist[penum]], graph->adjncy);
      else
        MPI_Send((void *)(gadjncy+gxadj[vtxdist[penum]]), ssize[penum], IDX_DATATYPE, penum, 1, comm); 
    }

    free(ssize);
  }
  else 
    MPI_Recv((void *)graph->adjncy, graph->nedges, IDX_DATATYPE, 0, 1, comm, &status);

  if (myproc == 0) 
    GKfree(&gxadj, &gadjncy, LTERM);

  MALLOC_CHECK(NULL);
}


static void ResortBoundaries(gridT *localgrid, int myproc)
{
  int neigh, n, *tmp;

  for(neigh=0;neigh<localgrid->Nneighs;neigh++) {
    tmp = (int *)malloc(localgrid->num_cells_send[neigh]*sizeof(int));
    for(n=0;n<localgrid->num_cells_send[neigh];n++)
      tmp[n]=localgrid->mnptr[localgrid->cell_send[neigh][n]];
    Sort(localgrid->cell_send[neigh],tmp,localgrid->num_cells_send[neigh]);
    free(tmp);

    tmp = (int *)malloc(localgrid->num_cells_recv[neigh]*sizeof(int));
    for(n=0;n<localgrid->num_cells_recv[neigh];n++)
      tmp[n]=localgrid->mnptr[localgrid->cell_recv[neigh][n]];
    Sort(localgrid->cell_recv[neigh],tmp,localgrid->num_cells_recv[neigh]);
    free(tmp);

    tmp = (int *)malloc(localgrid->num_edges_send[neigh]*sizeof(int));
    for(n=0;n<localgrid->num_edges_send[neigh];n++)
      tmp[n]=localgrid->eptr[localgrid->edge_send[neigh][n]];
    Sort(localgrid->edge_send[neigh],tmp,localgrid->num_edges_send[neigh]);
    free(tmp);

    tmp = (int *)malloc(localgrid->num_edges_recv[neigh]*sizeof(int));
    for(n=0;n<localgrid->num_edges_recv[neigh];n++)
      tmp[n]=localgrid->eptr[localgrid->edge_recv[neigh][n]];
    Sort(localgrid->edge_recv[neigh],tmp,localgrid->num_edges_recv[neigh]);
    free(tmp);
  }
}

static void MakePointers(gridT *maingrid, gridT **localgrid, int myproc, MPI_Comm comm)
{
  int i, n, nf, ne, neigh, neighproc, nc, j, k, mark, bctype, count;
  int **cell_send, **cell_recv, **edge_send, **edge_recv;
  int *num_cells_send, *num_cells_recv, *num_edges_send, *num_edges_recv;
  int *cellp, *edgep, *celldist, *edgedist, *lcptr, *leptr;
  int kcellsend, kcellrecv, kedgesend, kedgerecv;
  unsigned short *flagged;
  MPI_Status status;

  cellp = (int *)malloc((*localgrid)->Nc*sizeof(int));
  edgep = (int *)malloc((*localgrid)->Ne*sizeof(int));
  celldist = (int *)malloc((MAXBCTYPES-1)*sizeof(int));
  edgedist = (int *)malloc((MAXMARKS-1)*sizeof(int));
  lcptr = (int *)malloc(maingrid->Nc*sizeof(int));
  leptr = (int *)malloc(maingrid->Ne*sizeof(int));
  flagged = (unsigned short *)malloc((*localgrid)->Ne*sizeof(int));

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

  cell_send = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  cell_recv = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  edge_send = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  edge_recv = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  num_cells_send = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_cells_recv = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_edges_send = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_edges_recv = (int *)malloc((*localgrid)->Nneighs*sizeof(int));

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    num_cells_send[neigh]=0;
    num_cells_recv[neigh]=0;
    num_edges_send[neigh]=0;
    num_edges_recv[neigh]=0;
  }

  // Set up the pointers for receiving first and place the global
  // indices in the recv pointer array.  
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    for(n=0;n<(*localgrid)->Nc;n++) 
      if(maingrid->part[(*localgrid)->mnptr[n]]==neighproc)
	num_cells_recv[neigh]++;
    cell_recv[neigh]=(int *)malloc(num_cells_recv[neigh]*sizeof(int));

    k=0;
    for(n=0;n<(*localgrid)->Nc;n++) 
      if(maingrid->part[(*localgrid)->mnptr[n]]==neighproc)
	cell_recv[neigh][k++]=(*localgrid)->mnptr[n];
  }

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
    cell_send[neigh]=(int *)malloc(num_cells_send[neigh]*sizeof(int));
    MPI_Recv(cell_send[neigh],num_cells_send[neigh],MPI_INT,neighproc,1,comm,&status);
  }

  // Now set the indices to point to the local grid.
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    for(j=0;j<num_cells_send[neigh];j++) 
      cell_send[neigh][j]=lcptr[cell_send[neigh][j]];
    for(j=0;j<num_cells_recv[neigh];j++)
      cell_recv[neigh][j]=lcptr[cell_recv[neigh][j]];
  }

  // Now do the edges
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
    edge_send[neigh]=(int *)malloc(num_edges_send[neigh]*sizeof(int));
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
    edge_recv[neigh] = (int *)malloc(num_edges_recv[neigh]*sizeof(int));
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

static void MakePointers0(gridT *maingrid, gridT **localgrid, int myproc)
{
  int n, nf, ne, neigh, neighproc, nc, j, k, mark, bctype, count;
  int **cell_send, **cell_recv, **edge_send, **edge_recv;
  int *num_cells_send, *num_cells_recv, *num_edges_send, *num_edges_recv;
  int *cellp, *edgep, *celldist, *edgedist;
  int kcellsend, kcellrecv, kedgesend, kedgerecv;

  cellp = (int *)malloc((*localgrid)->Nc*sizeof(int));
  edgep = (int *)malloc((*localgrid)->Ne*sizeof(int));
  celldist = (int *)malloc((MAXBCTYPES-1)*sizeof(int));
  edgedist = (int *)malloc((MAXMARKS-1)*sizeof(int));

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

  cell_send = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  cell_recv = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  edge_send = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  edge_recv = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  num_cells_send = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_cells_recv = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_edges_send = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_edges_recv = (int *)malloc((*localgrid)->Nneighs*sizeof(int));

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    num_cells_send[neigh]=0;
    num_cells_recv[neigh]=0;
    num_edges_send[neigh]=0;
    num_edges_recv[neigh]=0;
  }
    
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    for(n=0;n<(*localgrid)->Nc;n++) {
      if(IsCellNeighborProc(n,maingrid,(*localgrid),myproc,neighproc)) {
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2) 
	  num_cells_send[neigh]++;
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==3)
	  num_cells_recv[neigh]++;
      }
    }
    for(n=0;n<(*localgrid)->Ne;n++)
      if((*localgrid)->mark[n]==5) 
	for(j=0;j<2;j++) {
	  nc = (*localgrid)->grad[2*n+j];
	  if(maingrid->part[(*localgrid)->mnptr[nc]]==neighproc) {
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4))
		num_edges_recv[neigh]++;
	    }
	    if(j==0)
	      nc=(*localgrid)->grad[2*n+1];
	    else
	      nc=(*localgrid)->grad[2*n];
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4))
		num_edges_send[neigh]++;
	    }
	  }
	}
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    cell_send[neigh] = (int *)malloc(num_cells_send[neigh]*sizeof(int));
    cell_recv[neigh] = (int *)malloc(num_cells_recv[neigh]*sizeof(int));
    edge_send[neigh] = (int *)malloc(num_edges_send[neigh]*sizeof(int));
    edge_recv[neigh] = (int *)malloc(num_edges_recv[neigh]*sizeof(int));
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    kcellsend=0;
    kcellrecv=0;
    neighproc=(*localgrid)->myneighs[neigh];
    for(n=0;n<(*localgrid)->Nc;n++) {
      if(IsCellNeighborProc(n,maingrid,(*localgrid),myproc,neighproc)) {
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2)
	  cell_send[neigh][kcellsend++]=n;
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==3)
	  cell_recv[neigh][kcellrecv++]=n;
      }
    }
    kedgesend=0;
    kedgerecv=0;
    for(n=0;n<(*localgrid)->Ne;n++)
      if((*localgrid)->mark[n]==5) 
	for(j=0;j<2;j++) {
	  nc = (*localgrid)->grad[2*n+j];
	  if(maingrid->part[(*localgrid)->mnptr[nc]]==neighproc) {
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4) &&
		 IsMember(ne,edge_recv[neigh],kedgerecv)==-1)
		edge_recv[neigh][kedgerecv++]=ne;
	    }
	    if(j==0)
	      nc=(*localgrid)->grad[2*n+1];
	    else
	      nc=(*localgrid)->grad[2*n];
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4) &&
		 IsMember(ne,edge_send[neigh],kedgesend)==-1)
		edge_send[neigh][kedgesend++]=ne;
	    }
	  }
	}
    num_edges_send[neigh]=kedgesend;
    num_edges_recv[neigh]=kedgerecv;
    edge_send[neigh]=ReSize(edge_send[neigh],num_edges_send[neigh]);
    edge_recv[neigh]=ReSize(edge_recv[neigh],num_edges_recv[neigh]);
  }

  /*
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    printf("%d <-> %d, Cells to Send: %d, Receive: %d\n",
	   myproc,(*localgrid)->myneighs[neigh],
	   num_cells_send[neigh],num_cells_recv[neigh]);
    printf("%d <-> %d, Edges to Send: %d, Receive: %d\n",
	   myproc,(*localgrid)->myneighs[neigh],
	   num_edges_send[neigh],num_edges_recv[neigh]);
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    printf("%d -> %d (Cells): ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_cells_send[neigh];n++)
      printf("%d ",cell_send[neigh][n]);
    printf("\n");
    printf("%d <- %d (Cells) ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_cells_recv[neigh];n++)
      printf("%d ",cell_recv[neigh][n]);
    printf("\n");
    printf("%d -> %d (Edges): ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_edges_send[neigh];n++)
      printf("%d ",edge_send[neigh][n]);
    printf("\n");
    printf("%d <- %d (Edges): ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_edges_recv[neigh];n++)
      printf("%d ",edge_recv[neigh][n]);
    printf("\n");
  }
  */

  (*localgrid)->cell_send=cell_send;
  (*localgrid)->cell_recv=cell_recv;
  (*localgrid)->edge_send=edge_send;
  (*localgrid)->edge_recv=edge_recv;
  (*localgrid)->num_cells_send=num_cells_send;
  (*localgrid)->num_cells_recv=num_cells_recv;
  (*localgrid)->num_edges_send=num_edges_send;
  (*localgrid)->num_edges_recv=num_edges_recv;
}

static void ReOrder(gridT *grid) 
{
  int n, nf, numflag, options[8], Nc = grid->Nc, Ne=grid->Ne;
  int *corder, *corderp, *eorder, *eorderp;
  REAL *tmp;

  corder = (int *)malloc(Nc*sizeof(int));
  corderp = (int *)malloc(Nc*sizeof(int));
  eorder = (int *)malloc(Ne*sizeof(int));
  eorderp = (int *)malloc(Ne*sizeof(int));
  tmp = (REAL *)malloc(2*(NFACES-1)*Ne*sizeof(REAL));

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

  grid->xadj = (int *)malloc((Nc+1)*sizeof(int));
  CreateCellGraph(grid);
  METIS_NodeND(&Nc,grid->xadj,grid->adjncy,&numflag,options,corder,corderp);
  free(grid->xadj);
  free(grid->adjncy);

  /*
  grid->xadj = (int *)malloc((Ne+1)*sizeof(int));
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

static void ReOrder0(gridT *maingrid, gridT **localgrid, int myproc)
{
  int *eorder, *eorderp, *corder, *corderp, *eflag, *cflag, *bctypeind, *marktypeind;
  int neigh, neighproc, j, k, n, nf, bctype, marktype, Nc=(*localgrid)->Nc, Ne=(*localgrid)->Ne;
  int nstart, nend, Nc_comp, Ne_comp;
  REAL *tmp, temp;

  corder=(int *)malloc(Nc*sizeof(int));
  eorder=(int *)malloc(Ne*sizeof(int));
  corderp=(int *)malloc(Nc*sizeof(int));
  eorderp=(int *)malloc(Ne*sizeof(int));
  tmp = (REAL *)malloc(2*(NFACES-1)*Ne*sizeof(REAL));

  // Place all of the computational cells before the boundary cells
  k=0;
  for(n=0;n<Nc;n++)
    corder[n]=0;
  for(n=0;n<Nc;n++) 
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)<MAXBCTYPES-1) 
      corder[k++]=n;
  (*localgrid)->Nc_comp = k;
  for(n=0;n<Nc;n++) 
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==MAXBCTYPES-1) 
      corder[k++]=n;

  for(n=0;n<Nc;n++) 
    corderp[corder[n]]=n;

  // Place all of the computational edges before the boundary edges
  k=0;
  for(n=0;n<Ne;n++)
    eorder[n]=0;
  for(n=0;n<Ne;n++)
    if((*localgrid)->mark[n]==0 || (*localgrid)->mark[n]==5)
      eorder[k++]=n;
  (*localgrid)->Ne_comp = k;
  for(n=0;n<Ne;n++)
    if((*localgrid)->mark[n]!=0 && (*localgrid)->mark[n]!=5)
      eorder[k++]=n;

  /*  Sort(&(eorder[(*localgrid)->Ne_comp]),
       &((*localgrid)->eptr[(*localgrid)->Ne_comp]),
       (*localgrid)->Ne-(*localgrid)->Ne_comp);*/
  //  Sort(&eorder,&((*localgrid)->eptr[0]),Ne);

  for(n=0;n<Ne;n++)
    eorderp[eorder[n]]=n;

  // Reorder the data corresponding to cells with the corder array
  ReOrderRealArray((*localgrid)->Ac,corder,tmp,(*localgrid)->Nc,1);
  ReOrderRealArray((*localgrid)->xv,corder,tmp,(*localgrid)->Nc,1);
  ReOrderRealArray((*localgrid)->yv,corder,tmp,(*localgrid)->Nc,1);
  ReOrderRealArray((*localgrid)->dv,corder,tmp,(*localgrid)->Nc,1);
  ReOrderIntArray((*localgrid)->cells,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->face,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->normal,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->neigh,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->vwgt,corder,(int *)tmp,(*localgrid)->Nc,1);
  ReOrderIntArray((*localgrid)->mnptr,corder,(int *)tmp,(*localgrid)->Nc,1);
  ReOrderIntArray((*localgrid)->Nk,corder,(int *)tmp,(*localgrid)->Nc,1);

  // Reorder the data corresponding to edges with the eorder array
  ReOrderRealArray((*localgrid)->df,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->dg,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->n1,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->n2,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->xi,eorder,tmp,(*localgrid)->Ne,2*(NFACES-1));
  ReOrderIntArray((*localgrid)->grad,eorder,(int *)tmp,(*localgrid)->Ne,2);
  ReOrderIntArray((*localgrid)->eneigh,eorder,(int *)tmp,(*localgrid)->Ne,2*(NFACES-1));
  ReOrderIntArray((*localgrid)->edges,eorder,(int *)tmp,(*localgrid)->Ne,NUMEDGECOLUMNS);
  ReOrderIntArray((*localgrid)->mark,eorder,(int *)tmp,(*localgrid)->Ne,1);
  ReOrderIntArray((*localgrid)->eptr,eorder,(int *)tmp,(*localgrid)->Ne,1);
  ReOrderIntArray((*localgrid)->Nke,eorder,(int *)tmp,(*localgrid)->Ne,1);

  // Now adjust the pointers to point to the new locations
  // Face and eneigh arrays point to edges
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<NFACES;nf++) 
      tmp[n*NFACES+nf]=(*localgrid)->face[n*NFACES+nf];
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<NFACES;nf++) 
      (*localgrid)->face[n*NFACES+nf]=eorderp[(int)(tmp[n*NFACES+nf])];

  for(n=0;n<Ne;n++)
    for(nf=0;nf<2*(NFACES-1);nf++)
      tmp[n*2*(NFACES-1)+nf]=(*localgrid)->eneigh[n*2*(NFACES-1)+nf];
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2*(NFACES-1);nf++)
      if((int)(tmp[n*2*(NFACES-1)+nf])!=-1)
	(*localgrid)->eneigh[n*2*(NFACES-1)+nf]=eorderp[(int)(tmp[n*2*(NFACES-1)+nf])];
      else
	(*localgrid)->eneigh[n*2*(NFACES-1)+nf]=-1;

  // Grad and neigh arrays point to cells
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2;nf++)
      tmp[2*n+nf]=(*localgrid)->grad[2*n+nf];
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2;nf++)
      if((int)(tmp[2*n+nf])!=-1)
	(*localgrid)->grad[2*n+nf]=corderp[(int)(tmp[2*n+nf])];
      else
	(*localgrid)->grad[2*n+nf]=-1;

  for(n=0;n<Nc;n++)
    for(nf=0;nf<NFACES;nf++)
      tmp[n*NFACES+nf]=(*localgrid)->neigh[n*NFACES+nf];
  for(n=0;n<Nc;n++)
    for(nf=0;nf<NFACES;nf++)
      if((int)(tmp[n*NFACES+nf])!=-1)
	(*localgrid)->neigh[n*NFACES+nf]=corderp[(int)(tmp[n*NFACES+nf])];
      else
  	(*localgrid)->neigh[n*NFACES+nf]=-1;

  free(corder);
  free(eorder);
  free(corderp);
  free(eorderp);
  free(tmp);
}

static void EdgeMarkers(gridT *maingrid, gridT **localgrid, int myproc)
{
  int j, n, ne, nc1, nc2, ne1, ne2, nf, nc, neigh, np1, np2, flag;

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

  unsigned short *flagged = (unsigned short *)malloc(grid->Nc*sizeof(unsigned short));

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
  int nf, j, k, proc, neigh1, neigh2, loc;

  (*maingrid)->numneighs=(int *)malloc(numprocs*sizeof(int));
  (*maingrid)->neighs=(int **)malloc(numprocs*sizeof(int *));

  // Initialize the arrays
  for(proc=0;proc<numprocs;proc++) {
    (*maingrid)->numneighs[proc]=0;
    (*maingrid)->neighs[proc]=(int *)malloc(numprocs*sizeof(int));    
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
  (*localgrid)->myneighs=(int *)malloc((*localgrid)->Nneighs*sizeof(int));

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
  int i, j, k, n, nc, nc1, nc2, nf, ne, ng, flag, mgptr, *lcptr, *leptr, bctype, iface;
  unsigned short *flagged = 
    (unsigned short *)malloc(maingrid->Ne*sizeof(unsigned short));

  lcptr = (int *)malloc(maingrid->Nc*sizeof(int));
  leptr = (int *)malloc(maingrid->Ne*sizeof(int));

  // Removed since set in CreateNearestPointers
  //    (*localgrid)->Nc = GetNumCells(maingrid,myproc);
  //    (*localgrid)->mnptr = (int *)malloc((*localgrid)->Nc*sizeof(int));
  (*localgrid)->vwgt = (int *)malloc((*localgrid)->Nc*sizeof(int));
  (*localgrid)->cells = (int *)malloc(NFACES*(*localgrid)->Nc*sizeof(int));
  (*localgrid)->xv = (REAL *)malloc((*localgrid)->Nc*sizeof(REAL));
  (*localgrid)->yv = (REAL *)malloc((*localgrid)->Nc*sizeof(REAL));
  (*localgrid)->dv = (REAL *)malloc((*localgrid)->Nc*sizeof(REAL));
  (*localgrid)->neigh = (int *)malloc(NFACES*(*localgrid)->Nc*sizeof(int));
  (*localgrid)->normal = (int *)malloc(NFACES*(*localgrid)->Nc*sizeof(int));

  for(j=0;j<maingrid->Nc;j++) 
    lcptr[j]=-1;
  for(j=0;j<maingrid->Ne;j++) 
    leptr[j]=-1;
  

  // Removed since set in CreateNearestPointers
  /*
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
  */

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

  for(j=0;j<(*localgrid)->Nc;j++) 
    for(nf=0;nf<NFACES;nf++) {
      mgptr = maingrid->neigh[(*localgrid)->mnptr[j]*NFACES+nf];
      if(mgptr>=0 && IsMember(mgptr,(*localgrid)->mnptr,(*localgrid)->Nc)>=0)
	(*localgrid)->neigh[j*NFACES+nf]=lcptr[mgptr];
      else
	(*localgrid)->neigh[j*NFACES+nf]=-1;
    }

  (*localgrid)->Ne = GetNumEdges(*localgrid);
  (*localgrid)->edges = (int *)malloc(NUMEDGECOLUMNS*(*localgrid)->Ne*sizeof(int));
  (*localgrid)->mark = (int *)malloc((*localgrid)->Ne*sizeof(int));
  (*localgrid)->eptr = (int *)malloc((*localgrid)->Ne*sizeof(int));
  (*localgrid)->eneigh = (int *)malloc(2*(NFACES-1)*(*localgrid)->Ne*sizeof(int));
  (*localgrid)->face = (int *)malloc(NFACES*(*localgrid)->Nc*sizeof(int));
  (*localgrid)->grad = (int *)malloc(2*(*localgrid)->Ne*sizeof(int));

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

  CreateFaceArray((*localgrid)->grad,(*localgrid)->neigh,(*localgrid)->face,
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
  int n, nf, k, j, Nc=(*grid)->Nc, Ne=(*grid)->Ne;
  REAL xt[NFACES], yt[NFACES], xc, yc, den;
  FILE *ofile;
  
  (*grid)->Ac = (REAL *)malloc(Nc*sizeof(REAL));
  (*grid)->df = (REAL *)malloc(Ne*sizeof(REAL));
  (*grid)->dg = (REAL *)malloc(Ne*sizeof(REAL));
  (*grid)->n1 = (REAL *)malloc(Ne*sizeof(REAL));
  (*grid)->n2 = (REAL *)malloc(Ne*sizeof(REAL));
  (*grid)->xi = (REAL *)malloc(2*(NFACES-1)*Ne*sizeof(REAL));

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
      xc = 0.5*(maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n]]+
		maingrid->xp[(*grid)->edges[NUMEDGECOLUMNS*n+1]]);
      yc = 0.5*(maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n]]+
		maingrid->yp[(*grid)->edges[NUMEDGECOLUMNS*n+1]]);
      if((*grid)->grad[2*n]==-1) {
        (*grid)->n1[n] = xc-(*grid)->xv[(*grid)->grad[2*n+1]];
        (*grid)->n2[n] = yc-(*grid)->yv[(*grid)->grad[2*n+1]];
	(*grid)->dg[n] = sqrt(pow((*grid)->n1[n],2)+
			      pow((*grid)->n2[n],2));
	(*grid)->n1[n] = (*grid)->n1[n]/(*grid)->dg[n];
	(*grid)->n2[n] = (*grid)->n2[n]/(*grid)->dg[n];
      } else {
        (*grid)->n1[n] = (*grid)->xv[(*grid)->grad[2*n]]-xc;
        (*grid)->n2[n] = (*grid)->yv[(*grid)->grad[2*n]]-yc;
	(*grid)->dg[n] = sqrt(pow((*grid)->n1[n],2)+
			      pow((*grid)->n2[n],2));
	(*grid)->n1[n] = (*grid)->n1[n]/(*grid)->dg[n];
	(*grid)->n2[n] = (*grid)->n2[n]/(*grid)->dg[n];
      }
    }
  }

  /* Now compute the coefficients that make up the tangents to compute advection */
  if(myproc==0 && VERBOSE>2) printf("Computing xi coefficients...\n");
  for(n=0;n<Ne;n++) 
    if((*grid)->n1[n]==0 || (*grid)->n2[n]==0)
      for(nf=0;nf<2*(NFACES-1);nf++) 
	(*grid)->xi[2*(NFACES-1)*n+nf]=0;
    else {
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

REAL GetArea(REAL *xt, REAL *yt, int Nf)
{
  int i;
  REAL n,b,r1,r2,h,l,a1,xt2[NFACES],yt2[NFACES],area;
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
}

static void InterpDepth(gridT *grid, int myproc, int numprocs, MPI_Comm comm)
{
  int n, Nd, proc, nstart, ncount, scaledepth;
  REAL *xd, *yd, *d, scaledepthfactor;
  char str[BUFFERLENGTH];
  FILE *ifile;
  MPI_Status status;

  scaledepth=(int)MPI_GetValue(DATAFILE,"scaledepth","InterpDepth",myproc);
  scaledepthfactor=MPI_GetValue(DATAFILE,"scaledepthfactor","InterpDepth",myproc);

  Nd = getsize(INPUTDEPTHFILE);
  xd = (REAL *)malloc(Nd*sizeof(REAL));
  yd = (REAL *)malloc(Nd*sizeof(REAL));
  d = (REAL *)malloc(Nd*sizeof(REAL));

  ifile = fopen(INPUTDEPTHFILE,"r");
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

  free(xd);
  free(yd);
  free(d);
}

void CorrectVoronoi(gridT *grid)
{
  int n, nf, nc1, nc2;
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
	}
      }
    }
  }
}
