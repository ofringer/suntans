/*
 * File: gridio.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Functions for reading/writing grid data.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "gridio.h"
#include "memory.h"

// Private Variables
#define COLUMNS_IN_TRIANGLE_CELLS_FILE 8  // Number of columns in original cells.dat before hybrid version of code
#define COLUMNS_IN_TRIANGLE_EDGES_FILE 5  // Number of columns in original edges.dat before netcdf version of code that includes edge_id

// Private Function declarations
static void ReadPointsData(char *filename, gridT *grid, int myproc);
static void WritePointsData(char *filename, gridT *grid, int myproc);
static void ReadEdgesData(char *filename, gridT *grid, int myproc);
static void WriteEdgesData(char *filename, gridT *grid, int myproc, const char *withedge_id);
static void ReadEdgeCenteredData(char *filename, gridT *grid, int myproc);
static void WriteEdgeCenteredData(char *filename, gridT *grid, int myproc);
static void ReadCellsData(char *filename, gridT *grid, int myproc);
static void WriteCellsData(char *filename, gridT *grid, int myproc, const char *withfaces);
static void ReadCellCenteredData(char *filename, gridT *grid, int myproc);
static void WriteCellCenteredData(char *filename, gridT *grid, int myproc);
static void ReadNodalData(char *filename, gridT *grid, int myproc);
static void WriteNodalData(char *filename, gridT *grid, gridT *maingrid, int myproc);
static void CheckTopologyFile(char *filename, int myproc, int numprocs);
static void ReadTopologyData(char *filename, gridT *grid, int myproc);
static void WriteTopologyData(char *filename, gridT *grid, int myproc, int numprocs);
static void CheckVertSpaceFile(char *filename, int myproc, int numprocs);
static void ReadVertSpaceData(char *filename, gridT *grid, int myproc);
static void WriteVertSpaceData(char *filename, gridT *grid, int myproc);

/************************************************************************/
/*                                                                      */
/*                          Public functions.                           */
/*                                                                      */
/************************************************************************/

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
  char str[BUFFERLENGTH], str2[BUFFERLENGTH];
  FILE *ifile;

  ReadGridFileNames(myproc);

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
  CheckTopologyFile(str,myproc,numprocs);

  if(VERBOSE>2) printf("Reading %s...\n",str);
  sprintf(str,"%s.%d",TOPOLOGYFILE,myproc);
  ReadTopologyData(str,*grid,myproc);
  
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
  (*grid)->maxfaces=(int)MPI_GetValue(DATAFILE,"maxFaces","ReadGrid",myproc);
  (*grid)->neigh = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->face = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->normal = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(int),"ReadGrid");
  (*grid)->def = (REAL *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->cells = (int *)SunMalloc((*grid)->maxfaces*(*grid)->Nc*sizeof(REAL),"ReadGrid");
  (*grid)->mnptr = (int *)SunMalloc((*grid)->Nc*sizeof(int),"ReadGrid");//MR

  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Reading %s...\n",str);
  ReadCellCenteredData(str,*grid,myproc);
  
  /*
   * Now read in edge-centered data.dat
   *
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
  (*grid)->edge_id = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid"); 
  (*grid)->edges = (int *)SunMalloc((*grid)->Ne*NUMEDGECOLUMNS*sizeof(int),"ReadGrid");
  (*grid)->eptr = (int *)SunMalloc((*grid)->Ne*sizeof(int),"ReadGrid");//MR

  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Reading %s...\n",str);
  ReadEdgeCenteredData(str,*grid,myproc);

  /* 
   * Now read in node data
   *
   */
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

  sprintf(str,"%s.%d",NODEFILE,myproc);
  if(myproc==0 && VERBOSE>2) printf("Reading %s...\n",str);
  ReadNodalData(str,*grid,myproc);

  /*
   * Now read in vertical grid spacing...
   *
   */
  CheckVertSpaceFile(VERTSPACEFILE,myproc,numprocs);
  
  (*grid)->Nkmax = MPI_GetValue(DATAFILE,"Nkmax","VertGrid",myproc);
  (*grid)->dz = (REAL *)SunMalloc((*grid)->Nkmax*sizeof(REAL),"ReadGrid");

  if(myproc==0 && VERBOSE>2) printf("Reading %s...\n",VERTSPACEFILE);
  ReadVertSpaceData(VERTSPACEFILE,*grid,myproc);

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

/*
 * Function: ReadGridFileNames
 * Usage: ReadGridFileNames(myproc);
 * ---------------------------------
 * Reads the names of the files containing the grid data.
 * The file pointers are globally defined as static character arrays
 * in suntans.h.
 * 
 */
void ReadGridFileNames(int myproc)
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

/*
 * Function: ReadDepth
 * Usage: ReadDepth(maingrid,myproc);
 * ----------------------------------
 * Read in depth from depth file specified in suntans.dat, followed by "-voro", and
 * puts it into maingrid->dv.
 * 
 */
void ReadDepth(gridT *grid, int myproc) {
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
 * **MR** EDGEFILE now has 6 columns. Last column contains the edge_id for boundary segments. **MR**
 *
 * This function will populate the following:
 *
 *   grid->xp
 *   grid->yp
 *   grid->edges
 *   grid->mark
 *   grid->grad
 *   grid->edge_id
 *   grid->xv
 *   grid->yv
 *   grid->nfaces
 *   grid->cells
 *   grid->neigh
 *
 */
void ReadMainGrid(gridT *grid, int myproc)
{
  if(VERBOSE>2 && myproc==0) printf("Reading %s for main grid...\n",POINTSFILE);
  ReadPointsData(POINTSFILE,grid,myproc);

  if(VERBOSE>2 && myproc==0) printf("Reading %s for main grid...\n",EDGEFILE);
  ReadEdgesData(EDGEFILE,grid,myproc);

  if(VERBOSE>2 && myproc==0) printf("Reading %s for main grid...\n",CELLSFILE);
  ReadCellsData(CELLSFILE,grid,myproc);
}

/*
 * Function: OutputGridData
 * Usage: OutputGridData(maingrid,localgrid,myproc,numprocs);
 * ----------------------------------------------------------
 * Outputs the required grid data.
 *
 */
void OutputGridData(gridT *maingrid, gridT *grid, int myproc, int numprocs)
{
  int j, n, nf, neigh, Np=maingrid->Np, Nc=grid->Nc, Ne=grid->Ne;
  char str[BUFFERLENGTH];
  FILE *ofile;

  // Assume that the triangulation output uses the same as the original format, i.e.
  // no first column and no edge_id
  if(TRIANGULATE && myproc==0) {
    if(VERBOSE>2) printf("Outputting %s...\n",POINTSFILE);
    WritePointsData(POINTSFILE,maingrid,myproc);

    if(VERBOSE>2) printf("Outputting %s...\n",EDGEFILE);
    WriteEdgesData(EDGEFILE,maingrid,myproc,"no edge_id");

    if(VERBOSE>2) printf("Outputting %s...\n",CELLSFILE);
    WriteCellsData(CELLSFILE,maingrid,myproc,"no nfaces");
  }
  
  // On each processor, cells.dat.* always contains the number of faces in the first column,
  // even if cells.dat does not.  
  sprintf(str,"%s.%d",CELLSFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str); 
  WriteCellsData(str,grid,myproc,"with nfaces");

  // On each processor, edges.dat.* always contains the edge_id in the last column,
  // even if edges.dat does not.
  sprintf(str,"%s.%d",EDGEFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);
  WriteEdgesData(str,grid,myproc,"with edge_id");

  // celldata.dat.* always has number of faces in first column
  sprintf(str,"%s.%d",CELLCENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);
  WriteCellCenteredData(str,grid,myproc);

  sprintf(str,"%s.%d",EDGECENTEREDFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);
  WriteEdgeCenteredData(str,grid,myproc);

  sprintf(str,"%s.%d",NODEFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);
  WriteNodalData(str,grid,maingrid,myproc);

  sprintf(str,"%s.%d",TOPOLOGYFILE,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str);
  WriteTopologyData(str,grid,myproc,numprocs);

  if(myproc==0) {
    if(VERBOSE>2) printf("Outputting %s...\n",VERTSPACEFILE);
    WriteVertSpaceData(VERTSPACEFILE,grid,myproc);
  }
}

/************************************************************************/
/*                                                                      */
/*                          Private functions.                          */
/*                                                                      */
/************************************************************************/

/*
 * Function: CheckTopologyFile
 * Usage: CheckTopologyFile(filename,myproc,numprocs);
 * -------------------------------------------------
 * Check to make sure the topology file specified in the character array filename
 * is written for the right number of processors.
 *
 */
static void CheckTopologyFile(char *filename, int myproc, int numprocs) {
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadGrid",myproc);

  if(numprocs!=((int)getfield(ifile,str))) {
    if(myproc==0) {
      printf("Error! Topology file(s) %s.*\n",TOPOLOGYFILE);
      printf("is/are not written for %d processors!\n",numprocs);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  fclose(ifile);
}

/*
 * Function: ReadTopologyData
 * Usage: ReadTopologyData(filename,grid,myproc);
 * ----------------------------------------------
 * Read in the topology data on the given processor.
 *
 */
static void ReadTopologyData(char *filename, gridT *grid, int myproc) {
  int n, neigh, numprocs;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadGrid",myproc);

  numprocs=(int)getfield(ifile,str);
  grid->Nneighs=(int)getfield(ifile,str);

  grid->myneighs=(int *)SunMalloc(grid->Nneighs*sizeof(int),"ReadTopologyData");
  grid->num_cells_send=(int *)SunMalloc(grid->Nneighs*sizeof(int),"ReadTopologyData");
  grid->num_cells_recv=(int *)SunMalloc(grid->Nneighs*sizeof(int),"ReadTopologyData");
  grid->num_edges_send=(int *)SunMalloc(grid->Nneighs*sizeof(int),"ReadTopologyData");
  grid->num_edges_recv=(int *)SunMalloc(grid->Nneighs*sizeof(int),"ReadTopologyData");
  grid->cell_send=(int **)SunMalloc(grid->Nneighs*sizeof(int *),"ReadTopologyData");
  grid->cell_recv=(int **)SunMalloc(grid->Nneighs*sizeof(int *),"ReadTopologyData");
  grid->edge_send=(int **)SunMalloc(grid->Nneighs*sizeof(int *),"ReadTopologyData");
  grid->edge_recv=(int **)SunMalloc(grid->Nneighs*sizeof(int *),"ReadTopologyData");
  grid->celldist = (int *)SunMalloc((MAXBCTYPES-1)*sizeof(int),"ReadTopologyData");
  grid->edgedist = (int *)SunMalloc((MAXMARKS-1)*sizeof(int),"ReadTopologyData");
  grid->cellp = (int *)SunMalloc(grid->Nc*sizeof(int),"ReadTopologyData");
  grid->edgep = (int *)SunMalloc(grid->Ne*sizeof(int),"ReadTopologyData");

  for(neigh=0;neigh<grid->Nneighs;neigh++) 
    grid->myneighs[neigh]=(int)getfield(ifile,str);

  for(neigh=0;neigh<grid->Nneighs;neigh++) {

    grid->num_cells_send[neigh]=(int)getfield(ifile,str);
    grid->num_cells_recv[neigh]=(int)getfield(ifile,str);
    grid->num_edges_send[neigh]=(int)getfield(ifile,str);
    grid->num_edges_recv[neigh]=(int)getfield(ifile,str);

    grid->cell_send[neigh]=(int *)SunMalloc(grid->num_cells_send[neigh]*sizeof(int),"ReadTopologyData");
    grid->cell_recv[neigh]=(int *)SunMalloc(grid->num_cells_recv[neigh]*sizeof(int),"ReadTopologyData");
    grid->edge_send[neigh]=(int *)SunMalloc(grid->num_edges_send[neigh]*sizeof(int),"ReadTopologyData");
    grid->edge_recv[neigh]=(int *)SunMalloc(grid->num_edges_recv[neigh]*sizeof(int),"ReadTopologyData");

    for(n=0;n<grid->num_cells_send[neigh];n++)
      grid->cell_send[neigh][n]=(int)getfield(ifile,str);
    for(n=0;n<grid->num_cells_recv[neigh];n++)
      grid->cell_recv[neigh][n]=(int)getfield(ifile,str);
    for(n=0;n<grid->num_edges_send[neigh];n++)
      grid->edge_send[neigh][n]=(int)getfield(ifile,str);
    for(n=0;n<grid->num_edges_recv[neigh];n++)
      grid->edge_recv[neigh][n]=(int)getfield(ifile,str);
  }

  for(n=0;n<MAXBCTYPES-1;n++) 
    grid->celldist[n]=(int)getfield(ifile,str);
  for(n=0;n<MAXMARKS-1;n++) 
    grid->edgedist[n]=(int)getfield(ifile,str);
  for(n=0;n<grid->Nc;n++) 
    grid->cellp[n]=(int)getfield(ifile,str);
  for(n=0;n<grid->Ne;n++) 
    grid->edgep[n]=(int)getfield(ifile,str);

  fclose(ifile);
}

/*
 * Function: WriteTopologyData
 * Usage: WriteTopologyData(filename,grid,myproc);
 * ----------------------------------------------
 * Write the topology data on the given processor.
 *
 */
static void WriteTopologyData(char *filename, gridT *grid, int myproc, int numprocs) {
  int n, neigh;
  FILE *ofile = MPI_FOpen(filename,"w","OutputData",myproc);

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
}

/*
 * Function: CheckVertSpaceFile
 * Usage: CheckVertSpaceFile(filename,myproc,numprocs);
 * ----------------------------------------------------
 * Check to make sure the vertical space file has the right number of
 * entries.
 *
 */
static void CheckVertSpaceFile(char *filename, int myproc, int numprocs) {
  int Nkmax_specified, Nkmax_filesize;
  Nkmax_specified = MPI_GetValue(DATAFILE,"Nkmax","VertGrid",myproc);
  Nkmax_filesize = MPI_GetSize(filename,"ReadGrid",myproc);

  if(Nkmax_specified!=Nkmax_filesize) {
    printf("Error in reading in grid data!\n");
    printf("Length of %s: %d is not consistent with what is in %s: %d\n",
	   VERTSPACEFILE,Nkmax_filesize,DATAFILE,Nkmax_specified);
    MPI_Finalize();
    exit(1);
  }
}

static void ReadVertSpaceData(char *filename, gridT *grid, int myproc) {
  int n;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadGrid",myproc);

  for(n=0;n<grid->Nkmax;n++)
    grid->dz[n]=getfield(ifile,str);
  fclose(ifile);
}

static void WriteVertSpaceData(char *filename, gridT *grid, int myproc) {
  int n;
  FILE *ofile = MPI_FOpen(filename,"w","OutputData",myproc);

  for(n=0;n<grid->Nkmax;n++)
    fprintf(ofile,"%e\n",grid->dz[n]);
  fclose(ofile);
}

static void ReadPointsData(char *filename, gridT *grid, int myproc) {
  int n;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadMainGrid",myproc);

  for(n=0;n<grid->Np;n++) {
    grid->xp[n]=getfield(ifile,str);
    grid->yp[n]=getfield(ifile,str);
    getfield(ifile,str);
  }
  fclose(ifile);
}

static void WritePointsData(char *filename, gridT *grid, int myproc) {
  FILE *ofile = MPI_FOpen(POINTSFILE,"w","OutputData",myproc);

  int j, Np=grid->Np;
  for(j=0;j<Np;j++)
    fprintf(ofile,"%e %e 0\n",grid->xp[j],grid->yp[j]);
}

static void ReadEdgesData(char *filename, gridT *grid, int myproc) {
  int j, n, numColumns;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadMainGrid",myproc);

  numColumns=getNumColumns(filename);

  for(n=0;n<grid->Ne;n++) {
    for(j=0;j<NUMEDGECOLUMNS-1;j++) 
      grid->edges[NUMEDGECOLUMNS*n+j]=(int)getfield(ifile,str);
    grid->mark[n]=(int)getfield(ifile,str);
    for(j=0;j<2;j++) 
      grid->grad[2*n+j]=(int)getfield(ifile,str);
    if(numColumns>COLUMNS_IN_TRIANGLE_EDGES_FILE) {
      grid->edge_id[n]=(int)getfield(ifile,str);
    } else {
      grid->edge_id[n]=0;
    }
  }
  fclose(ifile);
}

static void WriteEdgesData(char *filename, gridT *grid, int myproc, const char *withedge_id) {
  int j, nf;
  FILE *ofile = MPI_FOpen(filename,"w","OutputData",myproc);

  for(j=0;j<grid->Ne;j++) {
    for(nf=0;nf<2;nf++)
      fprintf(ofile,"%d ",grid->edges[j*NUMEDGECOLUMNS+nf]);
    fprintf(ofile,"%d ",grid->mark[j]);
    for(nf=0;nf<2;nf++)
      fprintf(ofile,"%d ",grid->grad[2*j+nf]);
    if(!strcmp(withedge_id,"with edge_id"))
       fprintf(ofile,"%d ",grid->edge_id[j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);
}
  
static void ReadCellsData(char *filename, gridT *grid, int myproc) {
  int n, nei, nf, nfaces, numColumns;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadMainGrid",myproc);

  numColumns=getNumColumns(filename);

  for(n=0;n<grid->Nc;n++) {
 
    // Set the number of faces in each cell as follows:
    // 1) If there are more than eight columns in cells.dat, set nfaces[n] to that number only if it does 
    //    not exceed the number specified in suntans.dat
    // 2) Otherwise set it to the default value of 3
    if(numColumns>COLUMNS_IN_TRIANGLE_CELLS_FILE) {
      nfaces=getfield(ifile,str);
      if(nfaces>grid->maxfaces) {
	printf("Error!!: Number of faces %d in first column of cells file %s is larger than maximum of %d!\n",
	       nfaces,CELLSFILE,grid->maxfaces);
	MPI_Finalize();
	exit(EXIT_FAILURE);
      } else {
	grid->nfaces[n]=nfaces;
      }
    } else {
      grid->nfaces[n]=DEFAULT_NFACES;
    }
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

static void WriteCellsData(char *filename, gridT *grid, int myproc, const char *withfaces) {
  int j, nf;
  FILE *ofile = MPI_FOpen(filename,"w","OutputData",myproc);

  for(j=0;j<grid->Nc;j++) {
    if(!strcmp(withfaces,"with nfaces"))
      fprintf(ofile,"%d ",grid->nfaces[j]);
    fprintf(ofile,"%e %e ",grid->xv[j],grid->yv[j]);
    for(nf=0;nf<grid->nfaces[j];nf++)
      fprintf(ofile,"%d ",grid->cells[j*grid->maxfaces+nf]);
    for(nf=0;nf<grid->nfaces[j];nf++)
      fprintf(ofile,"%d ",grid->neigh[j*grid->maxfaces+nf]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);
}

static void ReadCellCenteredData(char *filename, gridT *grid, int myproc) {
  int n, nf;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadGrid",myproc);

  for(n=0;n<grid->Nc;n++) {
    grid->nfaces[n]=(int)getfield(ifile,str);
    if(grid->nfaces[n]>grid->maxfaces) {
      printf("Error!!: Number of faces %d in first column of cells file %s is larger than maximum of %d!\n",
 	     grid->nfaces[n],filename,grid->maxfaces);
      MPI_Finalize();
      exit(EXIT_FAILURE);    
    }    
    grid->xv[n]=getfield(ifile,str);
    grid->yv[n]=getfield(ifile,str);
    grid->Ac[n]=getfield(ifile,str);
    grid->dv[n]=getfield(ifile,str);
    grid->Nk[n]=(int)getfield(ifile,str);
    for(nf=0;nf<grid->nfaces[n];nf++)
      grid->face[grid->maxfaces*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<grid->nfaces[n];nf++)
      grid->neigh[grid->maxfaces*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<grid->nfaces[n];nf++)
      grid->normal[grid->maxfaces*n+nf]=(int)getfield(ifile,str);
    for(nf=0;nf<grid->nfaces[n];nf++)
      grid->def[grid->maxfaces*n+nf]=getfield(ifile,str);
    for(nf=0;nf<grid->nfaces[n];nf++)
      grid->cells[grid->maxfaces*n+nf]=(int)getfield(ifile,str);
    grid->mnptr[n]=(int)getfield(ifile,str);
  }
  fclose(ifile);
}

static void WriteCellCenteredData(char *filename, gridT *grid, int myproc) {
  int n, nf;
  FILE *ofile = MPI_FOpen(filename,"w","OutputData",myproc);

  for(n=0;n<grid->Nc;n++) {
    fprintf(ofile,"%d %e %e %e %e %d ",grid->nfaces[n],grid->xv[n],grid->yv[n],grid->Ac[n],grid->dv[n],grid->Nk[n]);

    for(nf=0;nf<grid->nfaces[n];nf++)
      fprintf(ofile,"%d ",grid->face[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++)
      fprintf(ofile,"%d ",grid->neigh[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++) 
      fprintf(ofile,"%d ",grid->normal[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++) 
      fprintf(ofile,"%e ",grid->def[grid->maxfaces*n+nf]);
    for(nf=0;nf<grid->nfaces[n];nf++) 
      fprintf(ofile,"%d ",grid->cells[grid->maxfaces*n+nf]);
    fprintf(ofile,"%d ",grid->mnptr[n]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);
}  

static void ReadEdgeCenteredData(char *filename, gridT *grid, int myproc) {
  int n;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadGrid",myproc);

  for(n=0;n<grid->Ne;n++) {
      grid->df[n] = getfield(ifile,str);
      grid->dg[n] = getfield(ifile,str);
      grid->n1[n] = getfield(ifile,str);
      grid->n2[n] = getfield(ifile,str);
      grid->xe[n] = getfield(ifile,str);
      grid->ye[n] = getfield(ifile,str);
      grid->Nke[n] = (int)getfield(ifile,str);
      grid->Nkc[n] = (int)getfield(ifile,str);
      grid->grad[2*n] = (int)getfield(ifile,str);
      grid->grad[2*n+1] = (int)getfield(ifile,str);
      grid->gradf[2*n] = (int)getfield(ifile,str);
      grid->gradf[2*n+1] = (int)getfield(ifile,str);
      grid->mark[n] = (int)getfield(ifile,str);
      grid->edges[NUMEDGECOLUMNS*n] = (int)getfield(ifile,str);
      grid->edges[NUMEDGECOLUMNS*n+1] = (int)getfield(ifile,str);
      grid->edge_id[n] = (int)getfield(ifile,str);//MR
      grid->eptr[n] = (int)getfield(ifile,str);//MR
  }
  fclose(ifile);
}

static void WriteEdgeCenteredData(char *filename, gridT *grid, int myproc) {
  int n, Ne=grid->Ne;
  FILE *ofile = MPI_FOpen(filename,"w","OutputData",myproc);

  for(n=0;n<Ne;n++) {
    fprintf(ofile,"%e %e %e %e %e %e %d %d %d %d %d %d %d %d %d %d %d\n",
	    grid->df[n],grid->dg[n],grid->n1[n],grid->n2[n],grid->xe[n],grid->ye[n],
	    grid->Nke[n],grid->Nkc[n],grid->grad[2*n],grid->grad[2*n+1],
	    grid->gradf[2*n],grid->gradf[2*n+1],grid->mark[n], 
	    grid->edges[n*NUMEDGECOLUMNS], grid->edges[n*NUMEDGECOLUMNS+1],grid->edge_id[n],grid->eptr[n]);
  }
  fclose(ofile);
}

static void ReadNodalData(char *filename, gridT *grid, int myproc) {
  int n, np;
  char str[BUFFERLENGTH];
  FILE *ifile = MPI_FOpen(filename,"r","ReadGrid",myproc);

  for(n=0; n < grid->Np; n++) {
    // get the pointer to the global value
    grid->localtoglobalpoints[n] = (int)getfield(ifile,str);

    // get the coordinates of the point
    grid->xp[n] = getfield(ifile,str);
    grid->yp[n] = getfield(ifile,str);

    // now get the nodal neighbors
    // now read in the number of nodal neighbors
    grid->numppneighs[n] = (int)getfield(ifile,str);

    // now allocate the memory
    grid->ppneighs[n] = (int*)SunMalloc(grid->numppneighs[n]*sizeof(int),"ReadGrid");

    // now loop over all these and get the values for the nodal nodal neighbors
    for(np=0; np < grid->numppneighs[n]; np++) {
      grid->ppneighs[n][np] = (int)getfield(ifile,str);
    }

    // now get the number of edge neighbors
    // now read in the number of edge neighbors
    grid->numpeneighs[n] = (int)getfield(ifile,str);
    // now allocate the memory
    grid->peneighs[n] = (int*)SunMalloc(grid->numpeneighs[n]*sizeof(int),"ReadGrid");
    // now loop over all these and get the values for the nodal edge neighbors
    for(np=0; np < grid->numpeneighs[n]; np++) {
      grid->peneighs[n][np] = (int)getfield(ifile,str);
    }

    // now get the number of cell neighbors
    // now read in the number of cell neighbors
    grid->numpcneighs[n] = (int)getfield(ifile,str);
    // now allocate the memory
    grid->pcneighs[n] = (int*)SunMalloc(grid->numpcneighs[n]*sizeof(int),"ReadGrid");
    // now loop over all these and get the values for the nodal cell neighbors
    for(np=0; np < grid->numpcneighs[n]; np++) {
      grid->pcneighs[n][np] = (int)getfield(ifile,str);
    }

    // now get the number of vertical layers Nkp
    grid->Nkp[n] = (int)getfield(ifile,str);
    // allocate memory for each layer
    grid->Actotal[n] = (REAL *)SunMalloc(grid->Nkp[n]*sizeof(REAL),"ReadGrid");
    // now get the total area of shared cells at this point  (Actotal)
    for(np=0; np < grid->Nkp[n]; np++) {
      grid->Actotal[n][np] = (REAL)getfield(ifile,str);
    }
  }
  fclose(ifile);
}

static void WriteNodalData(char *filename, gridT *grid, gridT *maingrid, int myproc) {
  int n, neigh;
  FILE *ofile = MPI_FOpen(filename,"w","OutputData",myproc);

  for(n=0;n<maingrid->Np;n++) { // Np needs to be local, global right now
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
}
