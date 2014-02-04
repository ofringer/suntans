/*
 * File: sendrecv.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Functions for handling send/recv of unstructured-grid boundary data.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "sendrecv.h"
#include "timer.h"
#include "memory.h"

// Private functions (no longer used)
static void SendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm);
static void SendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm);
static void SendRecvWData(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm);
static void SendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, MPI_Comm comm);

/************************************************************************/
/*                                                                      */
/*               Public functions (used to be in grid.c)                */
/*                                                                      */
/************************************************************************/

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
void ISendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm)
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
void ISendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm)
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
void ISendRecvWData(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm)
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
void ISendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, MPI_Comm comm)
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
 * Function: CheckCommunicateCells
 * Usage: CheckCommunicateCells(maingrid,localgrid,myproc,comm);
 * -------------------------------------------------------------
 * Checks proper send/recv of interprocessor cells.
 *
 */
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

/*
 * Function: CheckCommunicateEdges
 * Usage: CheckCommunicateEdges(maingrid,localgrid,myproc,comm);
 * -------------------------------------------------------------
 * Checks proper send/recv of interprocessor edges.
 *
 */
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

/*************************************************************************/
/*                                                                       */
/* Old send/recv functions. No longer used.                              */
/*                                                                       */
/*************************************************************************/

/*
 * Function: SendRecvCellData2D
 * Usage: SendRecvCellData2D(grid->h,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the cell data back and forth between
 * processors.  
 *
 */
static void SendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm)
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
static void SendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm)
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
static void SendRecvWData(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm)
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
static void SendRecvEdgeData3D(REAL **edgedata, gridT *grid, int myproc, MPI_Comm comm)
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

