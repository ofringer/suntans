/*
 * File: triangulate.c
 * --------------------
 * Uses triangle libraries to create a triangulation from a specified file.
 *
 * $Id: triangulate.c,v 1.4 2003-05-02 23:07:39 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2003/04/29 16:39:07  fringer
 * Added MPI_FOPen in place of fopen.
 *
 * Revision 1.2  2003/04/29 00:16:27  fringer
 * Fixed VERBOSE>0 line to include myproc==0, which required a redefinition
 * of the function prototype for triangulate.
 *
 * Revision 1.1  2003/04/26 14:20:18  fringer
 * Initial revision
 *
 *
 */
#include "suntans.h"
#include "mympi.h"
#include "grid.h"
#include "fileio.h"
#include "triangle.h"

#define TRIANGLEFORMAT 0

int GetTriangulation(gridT **grid, int myproc);
void GetPoints(struct triangulateio *in, REAL *minarea, int myproc);
void InitializeTriangle(struct triangulateio *mid, struct triangulateio *vorout);

/*
 * Function: GetTriangulation
 * Usage: GetTriangulation(grid);
 * ------------------------------
 * Creates a triangulation from a PLSG file using triangle libraries.
 *
 */
int GetTriangulation(gridT **grid, int myproc) {
  int n, j, nf, Np, Ne, Nc;
  struct triangulateio in, out, vorout;
  REAL minarea;
  char str[BUFFERLENGTH];

  GetPoints(&in,&minarea,myproc);
  InitializeTriangle(&out,&vorout);
  
  // Options for the triangulation:
  // Q quiet
  // z C-style numbering (start at 0 instead of 1)
  // e edge list.
  // v voronoi list.
  // n neighbor list.
  // c triangulate interior of convex hull in pslg
  // a maximum triangle area
  // q (from README triangle):
  //    Adds points to the mesh to ensure that no angles smaller than 20 degrees occur.
  if(in.numberofsegments==0)
    sprintf(str,"Qznevc"); 
  else
    sprintf(str,"Qzqnev");
  if(minarea>0)
    sprintf(str,"%sa%.5f",str,minarea);
  else
    sprintf(str,"%s",str);

  // This is the triangulation performed by the triangle.c program.
  triangulate(str, &in, &out, &vorout);

  Np = out.numberofpoints;
  Ne = out.numberofedges;
  Nc = out.numberoftriangles;

  if(VERBOSE>0 && myproc==0 && Nc>0) 
    printf("Created a triangulation with %d Cells, %d Edges, %d Delaunay points...\n",
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

  if(Nc==0)
    return 0;
  return 1;
}

void GetPoints(struct triangulateio *in, REAL *minarea, int myproc)
{
  int i, dummy, dimensions;
  char str[BUFFERLENGTH], c;
  REAL num;
  MPI_GetString(PSLGFILE,DATAFILE,"pslg","GetPoints",0);
  FILE *infile = MPI_FOpen(PSLGFILE,"r","GetPoints",myproc);

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

