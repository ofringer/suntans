/*
 * File: report.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 10/03/02
 * --------------------------------
 * This file contains functions that print out data into
 * files as well as printing for debugging.
 *
 * $Id: report.c,v 1.2 2003-03-28 11:36:28 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2002/11/03 00:22:22  fringer
 * Initial revision
 *
 *
 */
#include "report.h"

void ParseFlags(int argc, char *argv[], int myproc)
{
  int i, j, done=0;
  GRID=0;
  SOLVE=0;
  VERBOSE=0;
  WARNING=0;
  ASCII=0;

  if(argc>1) {
    for(i=1;i<argc;i++) {
      if(argv[i][0]=='-') {
	if(argv[i][1]=='g' && strlen(argv[i])==2)
	    GRID=1;
	else if(argv[i][1]=='s' && strlen(argv[i])==2)
	  SOLVE=1;
	else if(argv[i][1]=='w' && strlen(argv[i])==2)
	  WARNING=1;
	else if(argv[i][1]=='a' && strlen(argv[i])==2)
	  ASCII=1;
	else if(argv[i][1]=='v') 
	  for(j=1;j<strlen(argv[i]);j++)
	    if(argv[i][j]=='v')
	      VERBOSE++;
	    else
	      done=1;
	else
	    done=1;
      } else 
	done=1;
    }
  }

  if(done) {
    if(myproc==0) Usage(argv[0]);
    MPI_Finalize();
    exit(0);
  }
}

void Usage(char *str)
{
  printf("Usage: %s {-g,-s,-v[vv]}\n",str);
}

/*
 * Function: ReportPartition
 * Usage: ReportPartition(maingrid,localgrid,myproc,comm);
 * -------------------------------------------------------
 * Print out the partitioning results obtained from
 * partitioning the cell-centered graph.
 *
 */
void ReportPartition(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm)
{
  int j, ncells=0, totwgt=0, wgt=0;

  totwgt=0;
  for(j=0;j<maingrid->Nc;j++)
    totwgt+=maingrid->vwgt[j];

  for(j=0;j<maingrid->Nc;j++) 
    if(maingrid->part[j] == myproc) {
      wgt+=maingrid->vwgt[j];
      ncells++;
    }

  if(myproc==0) {
    printf("\n");
    printf("Partitioning Results\n");
    printf("Main grid: %d Cells, %d Edges\n",maingrid->Nc,maingrid->Ne);
    printf("------------------------------------------------------\n");
    printf("Proc\tComp-Cells\tGhost-cells\tEdges\tWeight\n");
  }

  MPI_Barrier(comm);
  printf("%d\t%d\t\t%d\t\t%d\t%.3f\n",myproc,ncells,localgrid->Nc-ncells,
	 localgrid->Ne,(REAL)wgt/(REAL)totwgt);
}

/*
 * Function: ReportConnectivity
 * Usage: ReportConnectivity(grid,maingrid,myproc);
 * ------------------------------------------------
 * Prints out the connectivity on the current processor,
 * including the local indices as well as their corresponding
 * global indices.
 *
 */
void ReportConnectivity(gridT *grid, gridT *maingrid, int myproc)
{
  int j, nf, mgptr;

  /*
   * Check to make sure grid->mnptr exists, since this
   * must only be used for local grids.
   *
   */
  if(grid->mnptr) {
    printf("Connectivity on Processor %d (Global indices in parens)\n",myproc);
    for(j=0;j<grid->Nc;j++) {
      printf("Cell %d (BC type %d) (%d): ",j,
	     IsBoundaryCell(grid->mnptr[j],maingrid,myproc),grid->mnptr[j]);
      for(nf=0;nf<NFACES;nf++) {
	mgptr = grid->neigh[j*NFACES+nf];
	if(mgptr!=-1)
	  printf("%d (%d) ",grid->neigh[j*NFACES+nf],
		 grid->mnptr[mgptr]);
	else
	  printf("%d (-1) ",grid->neigh[j*NFACES+nf]);
      }
      printf("\n");
    }
  } else {
    printf("Error!  Function ReportConnectivity only for local grids!\n");
  }
}

/*
 * Function: ReportBoundaryDistribution
 * Usage: ReportBoundaryDistribution(maingrid,localgrid,myproc);
 * -------------------------------------------------------------
 * Print out the topology of the grid, including the edgebnddist
 * and cellbnddist arrays.
 *
 */
void ReportBoundaryDistributions(gridT *maingrid, gridT *localgrid, int myproc)
{
  int marktype, bctype, neigh;

  printf("\nCell distribution on proc %d ( neighbors ",myproc);
  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++)
    printf("%d ",maingrid->neighs[myproc][neigh]);
  printf(")\n");
  for(bctype=0;bctype<MAXBCTYPES-2;bctype++) 
    printf("%d (%d,%d) ",bctype,
	   localgrid->cellbnddist[bctype],
	   localgrid->cellbnddist[bctype+1]);
  printf("\n");
  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) 
    printf("Boundary with proc %d:  BC2: (%d,%d), BC3 (%d,%d)\n",
	   maingrid->neighs[myproc][neigh],
	   localgrid->cellbnddist[(MAXBCTYPES-2)+2*neigh],
	   localgrid->cellbnddist[(MAXBCTYPES-2)+2*neigh+1],
	   localgrid->cellbnddist[(MAXBCTYPES-2)+2*neigh+1],
	   localgrid->cellbnddist[(MAXBCTYPES-2)+2*neigh+2]);

  printf("\nEdge distribution on proc %d ( neighbors ",myproc);
  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++)
    printf("%d ",maingrid->neighs[myproc][neigh]);
  printf(")\n");
  for(marktype=0;marktype<MAXMARKS-2;marktype++) 
    printf("%d (%d,%d) ",marktype,
	   localgrid->edgebnddist[marktype],
	   localgrid->edgebnddist[marktype+1]);
  printf("\n");
  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) 
    printf("Boundary with proc %d:  BC5: (%d,%d),  BC6: (%d,%d)\n",
	   maingrid->neighs[myproc][neigh],
	   localgrid->edgebnddist[(MAXMARKS-2)+2*neigh],
	   localgrid->edgebnddist[(MAXMARKS-2)+2*neigh+1],
	   localgrid->edgebnddist[(MAXMARKS-2)+2*neigh+1],
	   localgrid->edgebnddist[(MAXMARKS-2)+2*neigh+2]);
}
  



  

