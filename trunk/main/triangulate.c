/*
 * File: triangulate.c
 * --------------------
 * Uses triangle libraries to create a triangulation from a specified file.
 *
 * $Id: triangulate.c,v 1.1 2003-04-26 14:20:18 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#include "suntans.h"
#include "mympi.h"
#include "grid.h"
#include "fileio.h"
#include "triangle.h"

#define TRIANGLEFORMAT 0

void GetTriangulation(gridT **grid);
void GetPoints(struct triangulateio *in, REAL *minarea);
void InitializeTriangle(struct triangulateio *mid, struct triangulateio *vorout);

/*
 * Function: GetTriangulation
 * Usage: GetTriangulation(grid);
 * ------------------------------
 * Creates a triangulation from a PLSG file using triangle libraries.
 *
 */
void GetTriangulation(gridT **grid) {
  int n, j, nf, Np, Ne, Nc;
  struct triangulateio in, out, vorout;
  REAL minarea;
  char str[BUFFERLENGTH];

  GetPoints(&in,&minarea);
  InitializeTriangle(&out,&vorout);

  //  triangulate("pzAevnqa10", &in, &out, &vorout);
  if(minarea>0)
    sprintf(str,"Qpzqevna%.5f",minarea);
  else
    sprintf(str,"Qpzqevn");

  // This is the triangulation performed by the triangle.c program.
  triangulate(str, &in, &out, &vorout);

  Np = out.numberofpoints;
  Ne = out.numberofedges;
  Nc = out.numberoftriangles;

  if(VERBOSE>0) printf("Created a triangulation with %d Cells, %d Edges, %d Delaunay points...\n",
		       Nc, Ne, Np);

  InitMainGrid(grid,Np,Ne,Nc);

  for(n=0;n<(*grid)->Np;n++) {
    (*grid)->xp[n]=out.pointlist[2*n];
    (*grid)->yp[n]=out.pointlist[2*n+1];
  }
  
  for(n=0;n<(*grid)->Ne;n++) {
    for(j=0;j<2;j++) {
      (*grid)->edges[NUMEDGECOLUMNS*n+j]=out.edgelist[2*n+j];
      (*grid)->grad[2*n+j]=vorout.edgelist[2*n+j];
    }
    (*grid)->mark[n]=out.edgemarkerlist[n];
  }

  for(n=0;n<(*grid)->Nc;n++) {
    (*grid)->xv[n] = vorout.pointlist[2*n];
    (*grid)->yv[n] = vorout.pointlist[2*n+1];
    for(nf=0;nf<NFACES;nf++) {
      (*grid)->cells[n*NFACES+nf]=out.trianglelist[n*NFACES+nf];
      (*grid)->neigh[n*NFACES+nf]=out.neighborlist[n*NFACES+nf];
    }
  }

  free(in.pointlist);
  free(in.pointmarkerlist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(in.holelist);
}

void GetPoints(struct triangulateio *in, REAL *minarea)
{
  int i, dummy, dimensions;
  char str[BUFFERLENGTH], c;
  REAL num;
  MPI_GetString(PSLGFILE,DATAFILE,"pslg","GetPoints",0);
  FILE *infile = fopen(PSLGFILE,"r");
  if(infile==NULL) {
    printf("Error in getting points for triangulation.  PSLG file %s does not exist!\n",
	   PSLGFILE);
    //    EndMPI();
  }

  in->numberofpoints = (int)getfield(infile,str);

  if(TRIANGLEFORMAT)
    for(i=0;i<3;i++) getfield(infile,str);
  in->pointlist = (REAL *) malloc(in->numberofpoints*2*sizeof(REAL));
  in->pointmarkerlist = (int *) malloc(in->numberofpoints * sizeof(int));
  for(i=0;i<in->numberofpoints;i++) {
    if(TRIANGLEFORMAT) getfield(infile,str);
    in->pointlist[2*i]= (REAL)getfield(infile,str);
    in->pointlist[2*i+1] = (REAL)getfield(infile,str);
    in->pointmarkerlist[i]= (int)getfield(infile,str);
  }

  in->numberofsegments = (int)getfield(infile,str);
  if(TRIANGLEFORMAT) getfield(infile,str);
  in->segmentlist = (int *) malloc(in->numberofsegments*2*sizeof(int));
  in->segmentmarkerlist = (int *) malloc(in->numberofsegments*sizeof(int));
  for(i=0;i<in->numberofsegments;i++) {
    if(TRIANGLEFORMAT) getfield(infile,str);
    in->segmentlist[2*i]= (int)getfield(infile,str);
    in->segmentlist[2*i+1] = (int)getfield(infile,str);
    in->segmentmarkerlist[i] = (int)getfield(infile,str);
  }

  in->numberofpointattributes = 0;
  in->numberofregions = 0;
  in->numberofholes = (int)getfield(infile,str);
  in->holelist=(REAL *)malloc(in->numberofholes*sizeof(REAL));
  for(i=0;i<in->numberofholes;i++) {
    in->holelist[2*i]=getfield(infile,str);
    in->holelist[2*i+1]=getfield(infile,str);
  }

  *minarea=getfield(infile,str);
  fclose(infile);
}

void InitializeTriangle(struct triangulateio *mid, struct triangulateio *vorout)
{
  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */

  mid->pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  mid->pointattributelist = (REAL *) NULL;
  mid->pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  mid->trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  mid->triangleattributelist = (REAL *) NULL;
  mid->neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  /* Needed only if segments are output (-p or -c) and -P not used: */
  mid->segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  mid->segmentmarkerlist = (int *) NULL;
  mid->edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  mid->edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

  vorout->pointlist = (REAL *) NULL;        /* Needed only if -v switch used. */
  /* Needed only if -v switch used and number of attributes is not zero: */
  vorout->pointattributelist = (REAL *) NULL;
  vorout->edgelist = (int *) NULL;          /* Needed only if -v switch used. */
  vorout->normlist = (REAL *) NULL;         /* Needed only if -v switch used. */
}

