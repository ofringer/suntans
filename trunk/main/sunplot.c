 /*
 * File: sunplot.c
 * ---------------
 * Description: Plot the output of suntans.
 *
 * Oliver Fringer
 * EFML Stanford University
 *
 * $Id: sunplot.c,v 1.11 2003-04-15 07:19:58 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.10  2003/04/09 17:46:08  fringer
 * Added point values using mouse click.
 *
 * Revision 1.9  2003/04/08 23:32:22  fringer
 * Added vertical velocity surface option.
 *
 * Revision 1.8  2003/04/08 17:18:44  fringer
 * Added loop for movie (the M button) and also fixed vector plots.
 *
 * Revision 1.7  2003/04/07 15:26:18  fringer
 * Added array of buttons.
 *
 * Revision 1.6  2003/04/06 23:48:51  fringer
 * Added zoom functionality and also fixed axes positioning for axis image option.
 *
 * Revision 1.5  2003/04/01 16:18:08  fringer
 * Added button functionality...
 *
 * Revision 1.4  2003/03/31 20:41:05  fringer
 * Still cleaning up...
 *
 * Revision 1.3  2003/03/31 20:23:21  fringer
 * Cleaning up...
 *
 * Revision 1.2  2003/03/31 19:57:36  fringer
 * Added log and id lines.
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include "math.h"
#include "suntans.h"
#include "fileio.h"

#define WIDTH 500
#define HEIGHT 500
#define PI 3.141592654
#define EMPTY 999999
#define INFTY 1e20
#define CMAPFILE "/home/fringer/research/SUNTANS/data/jet.cmap"
#define NUMCOLORS 66
//#define DEFAULT_FONT "-adobe-helvetica-medium-o-normal--20-140-100-100-p-98-iso8859-9"
#define DEFAULT_FONT "9x15"
#define WINLEFT .2
#define WINTOP .2
#define WINWIDTH .6
#define WINHEIGHT .7
#define BUTTONHEIGHT 20
#define ZOOMFACTOR 2.0
#define MINZOOMRATIO 1/100.0
#define MAXZOOMRATIO 100.0
#define NUMBUTTONS 21
#define POINTSIZE 2

typedef enum {
  in, out, box, none
} zoomT;

typedef enum {
  false, true
} bool;

typedef enum {
  left, center, right
} hjustifyT;

typedef enum {
  bottom, middle, top
} vjustifyT;

typedef enum {
  left_button=1, middle_button=2, right_button=3
} buttonnumT;

typedef enum {
  prevwin, gowin, nextwin,
  kupwin, kdownwin,
  prevprocwin, allprocswin, nextprocwin,
  saltwin, fswin,
  uwin, vwin, wwin, vecwin,
  depthwin, nonewin,
  edgewin, voronoiwin, delaunaywin,
  zoomwin, quitwin
} buttonName;

typedef enum {
  oneproc, allprocs
} plotProcT;

typedef enum {
  noslice, value, vertp, slice
} sliceT;

typedef enum {
  noplottype, freesurface, depth, h_d, salinity, w_velocity
} plottypeT;

typedef struct {
  char *string;
  char *mapstring;
  float l,b,w,h;
  Window butwin;
} myButtonT;

typedef struct {
  float *xc;
  float *yc;
  float **xv;
  float **yv;
  int **cells;
  int **edges;
  float **depth;

  float ***s;
  float ***u;
  float ***v;
  float ***wf;
  float ***w;
  float **h;
  float **h_d;

  float *currptr;

  int *Ne, *Nc, Np, Nkmax, nsteps, numprocs;
  int timestep, klevel;
  float umagmax;
} dataT;

void GetUMagMax(dataT *data, int klevel, int numprocs);
void DrawEdgeLines(float *xc, float *yc, int *cells, plottypeT plottype, int N);
void DrawDelaunayPoints(float *xc, float *yc, plottypeT plottype, int Np);
float *GetScalarPointer(dataT *data, plottypeT plottype, int klevel, int proc);
dataT *NewData(void);
void FreeGrid(void);
void FindNearest(dataT *data, int x, int y, int *iloc, int *procloc, int procnum, int numprocs);
void InitializeGraphics(void);
void MyDraw(dataT *data, plottypeT plottype, int procnum, int numprocs, int iloc, int procloc);
void LoopDraw(dataT *data, plottypeT plottype, int procnum, int numprocs);
void ReadGrid(int proc);
void ReadScalar(float *scal, int k, int Nk, int np, char *filename);
void ReadVelocity(float *u, float *v, int kp, int Nk, int np , char *filename);
float Min(float *x, int N);
float Max(float *x, int N);
void AxisImage(float *axes, float *data);
void Fill(XPoint *vertices, int N, int cindex, int edges);
void CAxis(dataT *data, plottypeT plottype, int klevel, int procnum, int numprocs);
void ReadColorMap(char *str);
void UnSurf(float *xc, float *yc, int *cells, float *data, int N);
void SetDataLimits(dataT *data);
void SetAxesPosition(void);
void Clf(void);
void Cla(void);
void Text(Window window, float x, float y, int boxwidth, int boxheight,
	  char *str, int fontsize, int color, 
	  hjustifyT hjustify, vjustifyT vjustify);
void DrawControls(dataT *data, int procnum, int numprocs);
Window NewButton(Window parent, char *name, int x, int y, 
		 int buttonwidth, int buttonheight, bool motion, int bordercolor);
void DrawButton(Window button, char *str);
void MapWindows(void);
void RedrawWindows(void);
void DrawZoomBox(void);
void DrawHeader(Window leftwin,Window rightwin,char *str);
void SetUpButtons(void);
void DrawVoronoiPoints(float *xv, float *yv, int Nc);
void FillCircle(int xp, int yp, int r, int ic, Window mywin);
void UnQuiver(float *xc, float *yc, int *edges, float *u, float *v, float umagmax, int Ne);
void DrawArrow(int xp, int yp, int ue, int ve, Window mywin, int ic);
void Rotate(XPoint *points, int N, int ue, int ve, int mag);
void ParseCommandLine(int N, char *str[], int *numprocs);
void ShowMessage(void);
void ReadData(dataT *data, int nstep, int numprocs);
void FreeData(dataT *data, int numprocs);
void CloseGraphics(void);

/*
 * Linux users will need to add -ldl to the Makefile to compile 
 * this example.
 *
 */
Display *dis;
myButtonT controlButtons[NUMBUTTONS];
Window win, axeswin, messagewin, controlswin;
XEvent report;
GC gc, fontgc;
XColor color;
Screen *screen;
Pixmap pix, zpix, controlspix;
Colormap colormap;
KeySym lookup;
XWindowAttributes windowAttributes;
XFontStruct *fontStruct;

int width=WIDTH, height=HEIGHT, newwidth, newheight, 
  Np0, Nc0, Ne0, n=1, k=0, keysym,
  xstart, ystart, xend, yend, lastgridread=0;
//int *cells, *edges;
//float caxis[2], *xc, *yc, *depth, *xp, axesPosition[4], dataLimits[4], buttonAxesPosition[4],
//  zoomratio, *xv, *yv, vlengthfactor=1.0;
float caxis[2], axesPosition[4], dataLimits[4], buttonAxesPosition[4],
  zoomratio, vlengthfactor=1.0;
int axisType, oldaxisType, white, black, red, blue, green, yellow, colors[NUMCOLORS];
bool edgelines, setdatalimits, pressed,   voronoipoints, delaunaypoints, vectorplot, go, goprocs,
  vertprofile, fromprofile, gridread;
char str[BUFFERLENGTH], message[BUFFERLENGTH];
zoomT zoom;
plotProcT procplottype;
sliceT sliceType;
plottypeT plottype = freesurface;

int main(int argc, char *argv[]) {
  int procnum=0, numprocs;
  dataT *data = NewData();
  bool redraw, quit=false;
  buttonnumT mousebutton;
  setdatalimits = false;
  vectorplot = false;
  zoomratio = 1;
  procplottype = allprocs;
  vertprofile = false;
  gridread = false;
  
  ParseCommandLine(argc,argv,&numprocs);  

  ReadData(data,-1,numprocs);

  InitializeGraphics();

  ReadColorMap(CMAPFILE);

  XSelectInput(dis, win, ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask );

  k=data->Nkmax/2-1;
  
  axisType='i';
  edgelines=true;
  voronoipoints=false;
  delaunaypoints=false;

  SetUpButtons();

  MapWindows();
  LoopDraw(data,plottype,procnum,numprocs);

  XMaskEvent(dis, ExposureMask, &report);
  pressed=false;
  go=false;
  goprocs=true;
  sprintf(message,"Messages are displayed here.");

  while(true) {
    fromprofile=false;
    redraw=false;
    zoom=none;
    if(go==true) {
      if(n<data->nsteps) {
	sprintf(message,"Stepping through...");
	redraw=true;
	n++;
      } else {
	redraw=false;
	go=false;
	sprintf(message,"Done!");
      }
    } else {
    XNextEvent(dis, &report);
    switch  (report.type) {
    case MotionNotify:
      xend=report.xbutton.x;
      yend=report.xbutton.y;
      XCopyArea(dis,pix,axeswin,gc,0,0,
		axesPosition[2]*width,
		axesPosition[3]*height,0,0);
      XCopyArea(dis,controlspix,controlswin,gc,0,0,
		buttonAxesPosition[2]*width,
		buttonAxesPosition[3]*height,0,0);
      DrawZoomBox();
      break;
    case ButtonPress:
      if(report.xany.window==axeswin) {
	xstart=report.xbutton.x;
	ystart=report.xbutton.y;
	pressed=true;
      }
      break;
    case ButtonRelease:
      mousebutton = report.xbutton.button;
      if(report.xany.window==controlButtons[prevwin].butwin) {
	if(mousebutton==left_button) {
	  if(n>1) { redraw = true; n--; }
	  else { redraw = false ; sprintf(message,"At n=1!"); }
	} else if(mousebutton==right_button)
	  if(n!=1) { redraw = true ; n=1; }
      } else if(report.xany.window==controlButtons[gowin].butwin) {
	if(mousebutton==left_button) {
	  go=true;
	}
      } else if(report.xany.window==controlButtons[nextwin].butwin) {
	if(mousebutton==left_button) {
	  if(n<data->nsteps) { redraw = true; n++; }
	  else { redraw = false ; sprintf(message,"At n=nsteps!"); }
	} else if(mousebutton==right_button) 
	  if(n!=data->nsteps) { redraw = true; n=data->nsteps; }
      } else if(report.xany.window==controlButtons[kdownwin].butwin) {
	if(mousebutton==left_button) {
	  if(k>0) { redraw = true; k--; }
	  else { redraw = false ; sprintf(message,"At k=1!"); }
	} else if(mousebutton==right_button)
	  if(k!=0) { redraw = true ; k=0; }
      } else if(report.xany.window==controlButtons[kupwin].butwin) {
	if(mousebutton==left_button) {
	  if(k<data->Nkmax-1) { redraw = true; k++; }
	  else { redraw = false ; sprintf(message,"At k=Nkmax!"); }
	} else if(mousebutton==right_button) 
	  if(k!=data->Nkmax-1) { redraw = true; k=data->Nkmax-1; }
      } else if(report.xany.window==controlButtons[prevprocwin].butwin && numprocs>1) {
	if(mousebutton==left_button) {
	  if(procplottype==allprocs) {
	    procnum=0;
	    procplottype=oneproc;
	    redraw=true;
	  } else if(procnum>0) { 
	    redraw = true; 
	    procnum--; 
	    procplottype=oneproc; 
	  } else { 
	    redraw = false ; 
	    sprintf(message,"At procnum==0!"); 
	  }
	} else if(mousebutton==right_button) 
	  if(procplottype==allprocs) {
	    procnum=0;
	    procplottype=oneproc;
	    redraw=true;
	  } else if(procnum!=0) { 
	    redraw = true; 
	    procnum=0; 
	    procplottype=oneproc; 
	  }
      } else if(report.xany.window==controlButtons[allprocswin].butwin
		&& mousebutton==left_button && numprocs>1) {
	sprintf(message,"Plotting all procs...");
	procplottype=allprocs;
	redraw=true;
      } else if(report.xany.window==controlButtons[nextprocwin].butwin && numprocs>1) {
	if(mousebutton==left_button) {
	  if(procplottype==allprocs) {
	    procnum=numprocs-1;
	    procplottype=oneproc;
	    redraw=true;
	  } else if(procnum<numprocs-1) { 
	    redraw = true;
	    procnum++;
	    procplottype=oneproc;
	  } else { 
	    redraw = false ; 
	    sprintf(message,"At procnum==numprocs!"); 
	  }
	} else if(mousebutton==right_button) 
	  if(procplottype==allprocs) {
	    procnum=numprocs-1;
	    procplottype=oneproc;
	    redraw=true;
	  } else if(procnum!=numprocs) { 
	    redraw = true; 
	    procnum=numprocs-1; 
	    procplottype=oneproc; 
	  }
      } else if(report.xany.window==controlButtons[saltwin].butwin && mousebutton==left_button) {
	if(plottype!=salinity) {
	  plottype=salinity;
	  sprintf(message,"Salinity selected...");
	  redraw=true;
	} else 
	  sprintf(message,"Salinity is already being displayed...");
      } else if(report.xany.window==controlButtons[fswin].butwin && mousebutton==left_button) {
	if(plottype!=freesurface) {
	  plottype=freesurface;
	  sprintf(message,"Free-surface selected...");
	  redraw=true;
	} else
	  sprintf(message,"Free-surface is already being displayed...");
      } else if(report.xany.window==controlButtons[vecwin].butwin) {
	if(vectorplot==false) {
	  vectorplot=true;
	  sprintf(message,"Vectors on...");
	  vlengthfactor=1;
	} else {
	  if(mousebutton==left_button) {
	    vectorplot=false;
	    sprintf(message,"Vectors off...");
	  } else if(mousebutton==middle_button) {
	    vlengthfactor/=2;
	    sprintf(message,"Halving the vector lengths...");
	  } else if(mousebutton==right_button) {
	    vlengthfactor*=2;
	    sprintf(message,"Doubling the vector lengths...");
	  }	
	}    
	redraw=true;
      } else if(report.xany.window==controlButtons[wwin].butwin && mousebutton==left_button) {
	if(plottype!=w_velocity) {
	  plottype=w_velocity;
	  sprintf(message,"Vertical velocity selected...");
	  redraw=true;
	} else
	  sprintf(message,"Vertical velocity is already being displayed...");
      } else if(report.xany.window==controlButtons[depthwin].butwin && mousebutton==left_button) {
	if(plottype!=depth) {
	  plottype=depth;
	  sprintf(message,"Depth selected...");
	  redraw=true;
	} else
	  sprintf(message,"Depth is already being displayed...");
      } else if(report.xany.window==controlButtons[depthwin].butwin && mousebutton==right_button) {
	if(plottype!=h_d) {
	  plottype=h_d;
	  sprintf(message,"Water height (h+d) selected...");
	  redraw=true;
	} else
	  sprintf(message,"Water height (h+d) is already being displayed...");
      } else if(report.xany.window==controlButtons[edgewin].butwin && mousebutton==left_button) {
	if(edgelines==false) {
	  sprintf(message,"Drawing edge lines");
	  edgelines=true;
	} else {
	  sprintf(message,"Removing edge lines");
	  edgelines=false;
	}
	redraw=true;
      } else if(report.xany.window==controlButtons[voronoiwin].butwin 
		&& mousebutton==left_button) {
	if(voronoipoints==false) {
	  sprintf(message,"Drawing Voronoi points");
	  voronoipoints=true;
	} else {
	  sprintf(message,"Removing Voronoi points");
	  voronoipoints=false;
	}
	redraw=true;
      } else if(report.xany.window==controlButtons[delaunaywin].butwin 
		&& mousebutton==left_button) {
	if(delaunaypoints==false) {
	  sprintf(message,"Drawing Delaunay points");
	  delaunaypoints=true;
	} else {
	  sprintf(message,"Removing Delaunay points");
	  delaunaypoints=false;
	}
	redraw=true;
      } else if(report.xany.window==controlButtons[nonewin].butwin && mousebutton==left_button) {
	if(plottype!=noplottype) {
	  plottype=noplottype;
	  sprintf(message,"Removing surface plot...");
	  redraw=true;
	}
      } else if(report.xany.window==controlButtons[zoomwin].butwin && mousebutton==left_button) {
	if(vertprofile==true) {
	  sprintf(message,"Zooming on...");
	  vertprofile=false;
	}
	else {
	  sprintf(message,"Zooming off...");
	  vertprofile=true;
	  sliceType=noslice;
	}
      } else if(report.xany.window==controlButtons[quitwin].butwin && mousebutton==left_button) {
	quit=true;
      } else if(report.xany.window==axeswin) {
	xend=report.xbutton.x;
	yend=report.xbutton.y;
	if(vertprofile==false) {
	  if(report.xbutton.button==left_button) {
	    if(xend==xstart && yend==ystart) 
	      zoom=in;
	    else 
	      zoom=box;
	    redraw=true;
	  } else if(report.xbutton.button==right_button) {
	    zoom=out;
	    redraw=true;
	  } else if(report.xbutton.button==middle_button) {
	    zoom=none;
	    zoomratio=1;
	    redraw=true;
	    setdatalimits=false;
	  }
	} else {
	  if(report.xbutton.button==left_button) {
	    sliceType=value;
	  } else if(report.xbutton.button==middle_button) {
	    sliceType=vertp;
	  } else if(report.xbutton.button==right_button) {
	    sliceType=slice;
	  }
	  redraw=true;
	  fromprofile=true;
	}
      }
      pressed=false;
      break;
    case Expose:   
      XGetWindowAttributes(dis,win,&windowAttributes);
      newwidth=windowAttributes.width;
      newheight=windowAttributes.height;
      if(newwidth!=width || newheight!=height) {
	width=newwidth;
	height=newheight;
	XResizeWindow(dis,win,width,height);
	XFreePixmap(dis,pix);
	pix = XCreatePixmap(dis,axeswin,axesPosition[2]*width,
			    axesPosition[3]*height,DefaultDepthOfScreen(screen));
      }
      RedrawWindows();
      LoopDraw(data,plottype,procnum,numprocs);
      break;
    case KeyPress:
      switch(keysym=XLookupKeysym(&report.xkey, 0)) {
      case XK_q:
	quit=true;
	break;
      case XK_k: case XK_Up:
	if(k<data->Nkmax-1) { redraw = true; k++; }
	  else { redraw = false; sprintf(message,"At k=Nkmax!"); }
	break;
      case XK_j: case XK_Down:
	if(k>0) { redraw = true; k--; }
	else { redraw = false ; sprintf(message,"At k=1!"); }
	break;
      case XK_p : case XK_Left:
	if(n>1) { redraw = true; n--; }
	else { redraw = false ; sprintf(message,"At n=1!"); }
	break;
      case XK_n: case XK_Right:
	if(n<data->nsteps) { redraw = true; n++; }
	else { redraw = false ; sprintf(message,"At n=nsteps!"); }
	break;
      case XK_s:
	if(plottype!=salinity) {
	  plottype=salinity;
	  sprintf(message,"Salinity selected...");
	  redraw=true;
	}
	break;
      case XK_h:
	if(plottype!=freesurface) {
	  plottype=freesurface;
	  sprintf(message,"Free-surface selected...");
	  redraw=true;
	}
	break;
      case XK_d:
	if(plottype!=depth) {
	  plottype=depth;
	  sprintf(message,"Depth selected...");
	  redraw=true;
	}
	break;
      case XK_a:
	if(axisType=='n') {
	  sprintf(message,"Changing plottype to image");
	  axisType='i';
	} else {
	  sprintf(message,"Changing plottype to normal");
	  axisType='n';
	}
	redraw=true;
	break;
      case XK_e:
	if(edgelines==0) {
	  sprintf(message,"Drawing edge lines");
	  edgelines=1;
	} else {
	  sprintf(message,"Removing edge lines");
	  edgelines=0;
	}
	redraw=true;
	break;
      case XK_0: case XK_1: case XK_2: case XK_3: case XK_4:
      case XK_5: case XK_6: case XK_7: case XK_8: case XK_9:
	str[0]=keysym;
	str[1]='\0';
	if(procnum!=atoi(str)) {
	  procnum=atoi(str);
	  redraw=true;
	}
	break;
      }
    }
    }
    if(quit)
      break;
    if(redraw) 
      LoopDraw(data,plottype,procnum,numprocs);
    ShowMessage();
  }
  FreeData(data,numprocs);
  CloseGraphics();
  return 0;
}

void ShowMessage(void) {
  int font_height;
  font_height=fontStruct->ascent+fontStruct->descent;

  XSetForeground(dis,gc,black);
  XFillRectangle(dis,messagewin,gc,
		 0,0,width*axesPosition[2],height*buttonAxesPosition[1]);

  XSetForeground(dis,fontgc,white);
  XDrawString(dis,messagewin,fontgc,0,font_height,
	      message,strlen(message));  
}

/*
void ReadGrid(int proc) {
  int i, ind;
  FILE *tfile,*dfile,*ifile = fopen("/home/fringer/research/SUNTANS/data/jet.cmap","r"), *efile,
    *pfile = fopen("/home/fringer/research/SUNTANS/data/points.dat","r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/cells.dat.%d",proc);
  tfile = fopen(str,"r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/celldata.dat.%d",proc);
  dfile = fopen(str,"r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/edges.dat.%d",proc);
  efile = fopen(str,"r");

  FreeGrid();

  Np=0;
  while(fgets(str,256,pfile)) Np++;
  fclose(pfile);
  pfile = fopen("/home/fringer/research/SUNTANS/data/points.dat","r"),

  Nc=0;
  while(fgets(str,256,tfile)) Nc++;
  fclose(tfile);
  sprintf(str,"/home/fringer/research/SUNTANS/data/cells.dat.%d",proc);
  tfile = fopen(str,"r"),

  Ne=0;
  while(fgets(str,256,efile)) Ne++;
  fclose(efile);
  sprintf(str,"/home/fringer/research/SUNTANS/data/edges.dat.%d",proc);
  efile = fopen(str,"r");

  xc = (float *)malloc(Np*sizeof(float));
  yc = (float *)malloc(Np*sizeof(float));
  xv = (float *)malloc(Nc*sizeof(float));
  yv = (float *)malloc(Nc*sizeof(float));
  cells = (int *)malloc(3*Nc*sizeof(int));
  edges = (int *)malloc(2*Ne*sizeof(int));
  depth = (float *)malloc(Nc*sizeof(float));

  for(i=0;i<Np;i++) {
    fscanf(pfile,"%f %f %d",&xc[i],&yc[i],&ind);
  }
  fclose(pfile);

  for(i=0;i<Ne;i++) 
    fscanf(efile,"%d %d %d %d %d",&edges[2*i],&edges[2*i+1],&ind,&ind,&ind);

  for(i=0;i<Nc;i++) {
    fscanf(dfile,"%f %f %f %f %d %d %d %d %d %d %d %d %d %d %d",
	   &xv[i],&yv[i],&xp,&depth[i],&ind,&ind,&ind,&ind,
	   &ind,&ind,&ind,&ind,&ind,&ind,&ind);
    fscanf(tfile,"%f %f %d %d %d %d %d %d",&xp,&xp,&cells[3*i],&cells[3*i+1],&cells[3*i+2],
	   &ind,&ind,&ind);
  }
  fclose(tfile);
  fclose(dfile);
  fclose(efile);

  lastgridread=proc;
  gridread=true;
}

void FreeGrid(void) {
  if(gridread) {
    free(xc);
    free(yc);
    free(xv);
    free(yv);
    free(cells);
    free(edges);
    free(depth);
  }
}
*/
/*
void ReadVelocity(float *u, float *v, int kp, int Nk, int np , char *filename) {
  int i, currptr, count;
  FILE *ufile = fopen(filename,"r");
  double *dum = (double *)malloc(Ne*sizeof(double));
  float dummy;
  
  currptr=fseek(ufile,(3*(np-1)*Ne*Nk + 3*Ne*(k-1))*sizeof(double),0);
  count = fread(dum,sizeof(double),Ne,ufile);
  for(i=0;i<Ne;i++)
    u[i]=dum[i];
  count=fread(dum,sizeof(double),Ne,ufile);
  for(i=0;i<Ne;i++)
    v[i]=dum[i];

  fclose(ufile);
  free(dum);
}

void ReadScalar(float *scal, int kp, int Nk, int np, char *filename) {
  int i, currptr;
  FILE *sfile = fopen(filename,"r");
  double *dum = (double *)malloc(Nc*sizeof(double));
  float dummy;

  for(i=0;i<Nc;i++)
    scal[i]=EMPTY;

  currptr=fseek(sfile,(Nc*(kp-1)+(np-1)*Nc*Nk)*sizeof(double),SEEK_CUR);
  fread(dum,sizeof(double),Nc,sfile);

  for(i=0;i<Nc;i++)
    scal[i]=dum[i];
  fclose(sfile);

  free(dum);
}
*/

void LoopDraw(dataT *data, plottypeT plottype, int procnum, int numprocs) {
  float umagmax;
  int iloc, procloc, proc;

  ReadData(data,n,numprocs);

  if(vertprofile && sliceType==value) {
    FindNearest(data,xend,yend,&iloc,&procloc,procnum,numprocs);
  } 

  if(plottype!=noplottype)
    CAxis(data,plottype,k,procnum,numprocs);

  if(vectorplot)
    GetUMagMax(data,k,numprocs);

  if(procplottype==allprocs) 
    for(proc=0;proc<numprocs;proc++) {
      if(proc==0 && !fromprofile) {
	SetDataLimits(data);
	SetAxesPosition();
	Cla();
      }
      MyDraw(data,plottype,proc,numprocs,iloc,procloc);
    }
  else {
    if(!fromprofile) {
      SetDataLimits(data);
      SetAxesPosition();
      Cla();
    }
    MyDraw(data,plottype,procnum,numprocs,iloc,procloc);
  }

  if(vectorplot) 
    if(procplottype==allprocs) 
      for(proc=0;proc<numprocs;proc++)
	UnQuiver(data->xc,data->yc,
		 data->edges[proc],data->u[proc][k],data->v[proc][k],
		 data->umagmax,data->Ne[proc]);
    else 
      UnQuiver(data->xc,data->yc,
	       data->edges[procnum],data->u[procnum][k],data->v[procnum][k],
	       data->umagmax,data->Ne[procnum]);
  
  XFlush(dis);
  XCopyArea(dis,pix,axeswin,gc,0,0,
	    axesPosition[2]*width,
	    axesPosition[3]*height,0,0);
  XCopyArea(dis,controlspix,controlswin,gc,0,0,
	    buttonAxesPosition[2]*width,
	    buttonAxesPosition[3]*height,0,0);
}

void MyDraw(dataT *data, plottypeT plottype, int procnum, int numprocs, int iloc, int procloc)
{
  int i;
  char tmpstr[BUFFERLENGTH];
  float *scal = GetScalarPointer(data,plottype,k,procnum);

  if(vertprofile==true && sliceType==value && procnum==procloc)
    if(plottype==freesurface || 
       plottype==depth || 
       plottype==h_d)
      sprintf(message,"Proc: %d, value(i=%d)=%f (x,y)=(%.2f,%.2f)",
	    procnum,iloc,scal[iloc],data->xv[procnum][iloc],data->yv[procnum][iloc]);
    else
      sprintf(message,"Proc: %d, value(i=%d,k=%d)=%f (x,y)=(%.2f,%.2f)",
	      procnum,iloc,k,scal[iloc],data->xv[procnum][iloc],
	      data->yv[procnum][iloc]);

  if(!fromprofile) {
    if(plottype!=noplottype)
      UnSurf(data->xc,data->yc,data->cells[procnum],scal,data->Nc[procnum]);
  
    if(delaunaypoints)
      DrawDelaunayPoints(data->xc,data->yc,plottype,data->Np);

    if(voronoipoints) 
      DrawVoronoiPoints(data->xv[procnum],data->yv[procnum],data->Nc[procnum]);
  
    if(edgelines)
      DrawEdgeLines(data->xc,data->yc,data->cells[procnum],plottype,data->Nc[procnum]);

    DrawControls(data,procnum,numprocs);
  }     
}

void UnQuiver(float *xc, float *yc, int *edges, float *u, float *v, float umagmax, int Ne) {
  int j, ic;
  float xe, ye, umag, l, lmax, n1, n2;
  int xp, yp, ue, ve, vlength;

  if(plottype==noplottype)
    ic = white;
  else
    ic = black;

  lmax = 0;
  for(j=0;j<Ne;j++) {
    l=sqrt(pow(xc[edges[2*j]]-xc[edges[2*j+1]],2)+
	   pow(yc[edges[2*j]]-yc[edges[2*j+1]],2));
    if(l>lmax) lmax=l;
  }
  vlength = vlengthfactor*(int)(lmax/(dataLimits[1]-dataLimits[0])*
		  axesPosition[2]*width);

  for(j=0;j<Ne;j++) {
    if(u[j]!=EMPTY) {
      xe = 0.5*(xc[edges[2*j]]+xc[edges[2*j+1]]);
      ye = 0.5*(yc[edges[2*j]]+yc[edges[2*j+1]]);
      
      xp = axesPosition[2]*width*(xe-dataLimits[0])/
	(dataLimits[1]-dataLimits[0]);
      yp = 	axesPosition[3]*height*(1-(ye-dataLimits[2])/
					(dataLimits[3]-dataLimits[2]));
      
      umag = sqrt(u[j]*u[j]+v[j]*v[j]);
      n1 = 0;
      if(u[j]!=0)
	n1 = umag/u[j];
      n2 = 0;
      if(v[j]!=0)
	n2 = umag/v[j];
      
      ue = vlength*u[j]/umagmax;
      ve = vlength*v[j]/umagmax;
      
      DrawArrow(xp,yp,ue,ve,pix,ic);
    }
  }
  
}

void DrawArrow(int xp, int yp, int ue, int ve, Window mywin, int ic) {
  int i, mag;
  float r=0.5, angle=30;
  XPoint points[4];

  mag = (int)sqrt(ue*ue+ve*ve);

  points[0].x = r*mag;
  points[1].x = 0;
  points[2].x = 0;
  points[3].x = r*mag;
  points[0].y = 0;
  points[1].y = -r*mag*sin(angle*M_PI/180);
  points[2].y = r*mag*sin(angle*M_PI/180);
  points[3].y = 0;
  Rotate(points,4,ue,-ve,mag);
  
  for(i=0;i<4;i++) {
    points[i].x+=(xp+(1-r)*ue);
    points[i].y+=(yp-(1-r)*ve);
  }
  
  XSetForeground(dis,gc,ic);
  XDrawLine(dis,mywin,gc,xp,yp,xp+ue,yp-ve);
  XFillPolygon(dis,mywin,gc,points,4,Convex,CoordModeOrigin);
}

void Rotate(XPoint *points, int N, int ue, int ve, int mag) {
  int i, xr, yr;
  float s = (float)ve/(float)mag, c = (float)ue/(float)mag;

  for(i=0;i<N;i++) {
    xr = points[i].x*c + points[i].y*s;
    yr = points[i].x*s - points[i].y*c;
    points[i].x = xr;
    points[i].y = yr;
  }
}    
  
void DrawDelaunayPoints(float *xc, float *yc, plottypeT plottype, int Np) {
  int i, ci, xp, yp;

  if(plottype != noplottype) 
    ci = black;
  else
    ci = white;

  for(i=0;i<Np;i++) {
    xp = axesPosition[2]*width*(xc[i]-dataLimits[0])/
      (dataLimits[1]-dataLimits[0]);
    yp = 	axesPosition[3]*height*(1-(yc[i]-dataLimits[2])/
					(dataLimits[3]-dataLimits[2]));
    FillCircle(xp,yp,POINTSIZE,ci,pix);
  }
}
  
void DrawVoronoiPoints(float *xv, float *yv, int Nc) {
  int xp, yp, i, ic;

  for(i=0;i<Nc;i++) {
    xp = axesPosition[2]*width*(xv[i]-dataLimits[0])/
      (dataLimits[1]-dataLimits[0]);
    yp = 	axesPosition[3]*height*(1-(yv[i]-dataLimits[2])/
					(dataLimits[3]-dataLimits[2]));
    FillCircle(xp,yp,POINTSIZE,red,pix);
  }
}

void FillCircle(int xp, int yp, int r, int ic, Window mywin) {
  XSetForeground(dis,gc,ic);
  XFillArc(dis,mywin,gc,xp-r/2,yp-r/2,r,r,0,360*64);
}
  
float Min(float *x, int N) {
  int i;
  float xmin = INFTY;
  for(i=0;i<N;i++)
    if(x[i]<xmin) xmin=x[i];

  return xmin;
}

float Max(float *x, int N) {
  int i;
  float xmax = -INFTY;
  for(i=0;i<N;i++)
    if(x[i]>xmax) xmax=x[i];

  return xmax;
}

void AxisImage(float *axes, float *data) {
  float h, w, dx, dy, densx, densy;

  dx = data[1]-data[0];
  dy = data[3]-data[2];
  
  densx = width*axes[2]/dx;
  densy = height*axes[3]/dy;

  if(densx>densy) {
    w = axes[2]*densy/densx;
    axes[0] = axes[0]+0.5*(axes[2]-w);
    axes[2] = w;
  } else {
    h = axes[3]*densx/densy;
    axes[1] = axes[1]+0.5*(axes[3]-h);
    axes[3] = h;
  }
}

void Fill(XPoint *vertices, int N, int cindex, int edges) {
  int i, ic;
  if(plottype!=noplottype) {
    XSetForeground(dis,gc,colors[cindex]);
    XFillPolygon(dis,pix,gc,vertices,3,Convex,CoordModeOrigin);
  }
}

float *GetScalarPointer(dataT *data, plottypeT plottype, int klevel, int proc) {
  float *scal;

  switch(plottype) {
  case freesurface:
    return data->h[proc];
    break;
  case depth:
    return data->depth[proc];
    break;
  case salinity:
    return data->s[proc][klevel];
    break;
  case w_velocity:
    return data->w[proc][klevel];
    break;
  case h_d:
    return data->h_d[proc];
    break;
  }
  return NULL;
}
  
void CAxis(dataT *data, plottypeT plottype, int klevel, int procnum, int numprocs) {
  int i, is, ie, ni;
  float *scal;

  caxis[0] = EMPTY;
  caxis[1] = -EMPTY;

  if(procplottype==allprocs) {
    is=0;
    ie=numprocs;
  } else {
    is=procnum;
    ie=procnum+1;
  }

  for(i=is;i<ie;i++) {
    for(ni=0;ni<data->Nc[i];ni++) {
      scal = GetScalarPointer(data,plottype,klevel,i);
      if(scal[ni]<=caxis[0] && scal[ni]!=0 && scal[ni]!=EMPTY) caxis[0]=scal[ni];
      if(scal[ni]>=caxis[1] && scal[ni]!=0 && scal[ni]!=EMPTY) caxis[1]=scal[ni];
    }
  }
}

void ReadColorMap(char *str) {
  int i;
  FILE *ifile = fopen(str,"r");
  float r, g, b;

  for(i=0;i<NUMCOLORS-2;i++) {
    fscanf(ifile,"%f %f %f\n",&r, &g, &b);
    color.red = r * 0xffff;
    color.green = g * 0xffff;
    color.blue = b * 0xffff;
    XAllocColor(dis,colormap,&color);  
    colors[i] = color.pixel;
  }
  colors[NUMCOLORS-2]=white;
  colors[NUMCOLORS-1]=black;
  fclose(ifile);
}

void DrawEdgeLines(float *xc, float *yc, int *cells, plottypeT plottype, int N) {
  int i, j, ind;
  XPoint *vertices = (XPoint *)malloc(4*sizeof(XPoint));

  for(i=0;i<N;i++) {
    for(j=0;j<3;j++) {
      vertices[j].x = 
	axesPosition[2]*width*(xc[cells[3*i+j]]-dataLimits[0])/
	(dataLimits[1]-dataLimits[0]);
      vertices[j].y = 
	axesPosition[3]*height*(1-(yc[cells[3*i+j]]-dataLimits[2])/
	(dataLimits[3]-dataLimits[2]));
    }

    vertices[3].x = vertices[0].x;
    vertices[3].y = vertices[0].y;

    if(plottype!=noplottype)
      XSetForeground(dis,gc,black);
    else
      XSetForeground(dis,gc,white);
    XDrawLines(dis,pix,gc,vertices,4,0);
  }
  free(vertices);
}

void UnSurf(float *xc, float *yc, int *cells, float *data, int N) {
  int i, j, ind;
  XPoint *vertices = (XPoint *)malloc(4*sizeof(XPoint));

  for(i=0;i<N;i++) {
    for(j=0;j<3;j++) {
      vertices[j].x = 
	axesPosition[2]*width*(xc[cells[3*i+j]]-dataLimits[0])/
	(dataLimits[1]-dataLimits[0]);
      vertices[j].y = 
	axesPosition[3]*height*(1-(yc[cells[3*i+j]]-dataLimits[2])/
	(dataLimits[3]-dataLimits[2]));
    }

    vertices[3].x = vertices[0].x;
    vertices[3].y = vertices[0].y;

    ind = (data[i]-caxis[0])/(caxis[1]-caxis[0])*(NUMCOLORS-3);
    if(data[i]==EMPTY || (plottype=='D' && data[i]==0))
      ind = NUMCOLORS-1;
    if(ind==NUMCOLORS-2)
      printf("index = %d\n",i);

    Fill(vertices,3,ind,edgelines);
  }
  free(vertices);
}

void InitializeGraphics(void) {
  int screen_number;
  XGCValues fontvalues;
  XColor temp1, temp2;

  dis = XOpenDisplay(NULL);
  screen = ScreenOfDisplay(dis,0);
  screen_number = XScreenNumberOfScreen(screen);

  black = BlackPixel(dis,screen_number);
  white = WhitePixel(dis,screen_number);

  width=XDisplayWidth(dis,screen_number);
  height=XDisplayHeight(dis,screen_number);
  
  win = XCreateSimpleWindow(dis, RootWindow(dis, 0), 
			    WINLEFT*width,WINTOP*height,
			    WINWIDTH*width,WINHEIGHT*height, 
			    0,black,black);
  width = WINWIDTH*width;
  height = WINHEIGHT*height;

  //  pix = XCreatePixmap(dis,win,width,height,DefaultDepthOfScreen(screen));

  XMapWindow(dis, win);

  colormap = DefaultColormap(dis,screen_number);
  gc = XCreateGC(dis, win,0,0);

  XAllocNamedColor(dis,colormap,"red",&temp1,&temp2);
  red = temp2.pixel;
  XAllocNamedColor(dis,colormap,"blue",&temp1,&temp2);
  blue = temp2.pixel;
  XAllocNamedColor(dis,colormap,"green",&temp1,&temp2);
  green = temp2.pixel;
  XAllocNamedColor(dis,colormap,"yellow",&temp1,&temp2);
  yellow = temp2.pixel;
  
  fontStruct = XLoadQueryFont(dis,DEFAULT_FONT);
  if(!fontStruct) {
    printf("Font \"%s\" does not exist!\n",DEFAULT_FONT);
    exit(0);
  } else {
    fontvalues.background = black;
    fontvalues.font = fontStruct->fid;
    fontvalues.fill_style = FillSolid;
    fontvalues.line_width = 1;
    fontgc = XCreateGC(dis,win, GCForeground | GCBackground |
		       GCFont | GCLineWidth | GCFillStyle, &fontvalues);
  }

  axesPosition[0] = 0.25;
  axesPosition[1] = 0.05;
  axesPosition[2] = 0.7;
  axesPosition[3] = 0.9;

  buttonAxesPosition[0]=.1*axesPosition[0];
  buttonAxesPosition[1]=axesPosition[1];
  buttonAxesPosition[2]=.8*axesPosition[0];
  buttonAxesPosition[3]=axesPosition[3];
}

void SetDataLimits(dataT *data) {
  float dx, dy, Xc, Yc, x1, y1, x2, y2, a1, a2, xmin, xmax, ymin, ymax;
  dx = dataLimits[1]-dataLimits[0];
  dy = dataLimits[3]-dataLimits[2];

  if(!setdatalimits) {
    dataLimits[0] = Min(data->xc,data->Np);
    dataLimits[1] = Max(data->xc,data->Np);
    dataLimits[2] = Min(data->yc,data->Np);
    dataLimits[3] = Max(data->yc,data->Np);
    setdatalimits=true;
  } else {
    switch(zoom) {
    case in:
      if(zoomratio>=MINZOOMRATIO) {
	Xc = dataLimits[0]+(float)xstart/((float)axesPosition[2]*width)*dx;
	Yc = dataLimits[3]-(float)ystart/((float)axesPosition[3]*height)*dy;
	dataLimits[0]=Xc-dx/ZOOMFACTOR/2;
	dataLimits[1]=Xc+dx/ZOOMFACTOR/2;
	dataLimits[2]=Yc-dy/ZOOMFACTOR/2;
	dataLimits[3]=Yc+dy/ZOOMFACTOR/2;
	zoomratio/=2;
      } else 
	printf("Cannot zoom further.  Beyond zoom limits!\n");
      break;
    case out:
      if(zoomratio<=MAXZOOMRATIO) {
	Xc = dataLimits[0]+(float)xstart/((float)axesPosition[2]*width)*dx;
	Yc = dataLimits[3]-(float)ystart/((float)axesPosition[3]*height)*dy;
	dataLimits[0]=Xc-dx*ZOOMFACTOR/2;
	dataLimits[1]=Xc+dx*ZOOMFACTOR/2;
	dataLimits[2]=Yc-dy*ZOOMFACTOR/2;
	dataLimits[3]=Yc+dy*ZOOMFACTOR/2;
	zoomratio*=2;
      } else 
	printf("Cannot zoom further.  Beyond zoom limits!\n");
      break;
    case box:
      a1 = dx*dy;
      x1 = dataLimits[0]+(float)xstart/((float)axesPosition[2]*width)*dx;
      x2 = dataLimits[0]+(float)xend/((float)axesPosition[2]*width)*dx;
      y1 = dataLimits[2]+(1-(float)ystart/((float)axesPosition[3]*height))*dy;
      y2 = dataLimits[2]+(1-(float)yend/((float)axesPosition[3]*height))*dy;
      a2 = abs(x2-x1)*abs(y2-y1);
      if(x2>x1) {
	xmin = x1;
	xmax = x2;
      }	else {
	xmin = x2;
	xmax = x1;
      }

      if(y2>y1) {
	ymin = y1;
	ymax = y2;
      }	else {
	ymin = y2;
	ymax = y1;
      }
      zoomratio*=sqrt(a2/a1);

      if(zoomratio>=MINZOOMRATIO) {
	dataLimits[0]=xmin;
	dataLimits[1]=xmax;
	dataLimits[2]=ymin;
	dataLimits[3]=ymax;
      } else 
	printf("Too fine a zooming area.  Please try again!\n");
      break;
    }
  }
}

void SetAxesPosition(void) {
  axesPosition[0] = 0.25;
  axesPosition[1] = 0.05;
  axesPosition[2] = 0.7;
  axesPosition[3] = 0.9;

  if(axisType=='i') 
    AxisImage(axesPosition,dataLimits);

  XMoveResizeWindow(dis,axeswin,
		    axesPosition[0]*width,
		    axesPosition[1]*height,
		    axesPosition[2]*width,
		    axesPosition[3]*height);
  XFreePixmap(dis,pix);
  pix = XCreatePixmap(dis,axeswin,axesPosition[2]*width,
		      axesPosition[3]*height,DefaultDepthOfScreen(screen));
}

void Cla(void) {
  XSetForeground(dis,gc,black);
  XFillRectangle(dis,pix,gc,0,0,width*axesPosition[2],height*axesPosition[3]);
}

void Clf(void) {
  XSetForeground(dis,gc,black);
  XFillRectangle(dis,win,gc,0,0,width,height);
}

void Text(Window window, float x, float y, int boxwidth, int boxheight,
	  char *str, int fontsize, int color, 
	  hjustifyT hjustify, vjustifyT vjustify) {

  char **fonts;
  XFontStruct **fontsReturn;
  int i, xf, yf, font_width, font_height, numfonts;

  font_width=XTextWidth(fontStruct,str,strlen(str));
  font_height=fontStruct->ascent+fontStruct->descent;

  switch(hjustify) {
  case left:
    xf=x*boxwidth;
    break;
  case center:
    xf=x*boxwidth-font_width/2;
    break;
  case right:
    xf=x*boxwidth-font_width;
    break;
  default:
    printf("Error!  Unknown horizontal justification type!\n");
    break;
  }

  switch(vjustify) {
  case bottom:
    yf=(1-y)*boxheight;
    break;
  case middle:
    yf=y*boxheight+font_height/2;
    break;
  case top:
    yf=y*boxheight-font_height;
    break;
  default:
    printf("Error!  Unknown vertical justification type!\n");
    break;
  }

  XSetForeground(dis,fontgc,color);
  XDrawString(dis,window,fontgc,xf,yf,str,strlen(str));  
  /*
  fonts = XListFontsWithInfo(dis,"*cour*",200,&numfonts,fontsReturn);
  for(i=0;i<numfonts;i++)
    printf("Found %s\n",fonts[i]);
  exit(0);
  */
}

void RedrawWindows(void) {
  int buttonnum;

  XMoveResizeWindow(dis,controlswin,
		    buttonAxesPosition[0]*width,
		    buttonAxesPosition[1]*height,
		    buttonAxesPosition[2]*width,
		    buttonAxesPosition[3]*height);
  XFreePixmap(dis,controlspix);
  controlspix = XCreatePixmap(dis,controlswin,buttonAxesPosition[2]*width,
			      buttonAxesPosition[3]*height,DefaultDepthOfScreen(screen));

  XMoveResizeWindow(dis,axeswin,
		      axesPosition[0]*width,
		      axesPosition[1]*height,
		      axesPosition[2]*width,
		      axesPosition[3]*height);

  XMoveResizeWindow(dis,messagewin,
		    axesPosition[0]*width,
		    (1.2*buttonAxesPosition[1]+buttonAxesPosition[3])*height,
		    axesPosition[2]*width,
		    buttonAxesPosition[1]*height);

  /*
  XMoveResizeWindow(dis,prevwin,
		    (0*buttonAxesPosition[0]+.05*buttonAxesPosition[2])*width,
		    (0*buttonAxesPosition[1]+.12*buttonAxesPosition[2])*height,
		    buttonAxesPosition[2]*.4*width,
		    BUTTONHEIGHT);
  */
  for(buttonnum=0;buttonnum<NUMBUTTONS;buttonnum++)
    XMoveResizeWindow(dis,controlButtons[buttonnum].butwin,
		      controlButtons[buttonnum].l*buttonAxesPosition[2]*width,
		      controlButtons[buttonnum].b*buttonAxesPosition[2]*height,
		      controlButtons[buttonnum].w*buttonAxesPosition[2]*width,
		      controlButtons[buttonnum].h);
  /*
  XMoveResizeWindow(dis,nextwin,
		    (0*buttonAxesPosition[0]+.55*buttonAxesPosition[2])*width,
		    (0*buttonAxesPosition[1]+.12*buttonAxesPosition[2])*height,
		    buttonAxesPosition[2]*.4*width,
		    BUTTONHEIGHT);

  XMoveResizeWindow(dis,kdownwin,
		    (0*buttonAxesPosition[0]+.05*buttonAxesPosition[2])*width,
		    (0*buttonAxesPosition[1]+.37*buttonAxesPosition[2])*height,
		    buttonAxesPosition[2]*.4*width,
		    BUTTONHEIGHT);
  XMoveResizeWindow(dis,kupwin,
		    (0*buttonAxesPosition[0]+.55*buttonAxesPosition[2])*width,
		    (0*buttonAxesPosition[1]+.37*buttonAxesPosition[2])*height,
		    buttonAxesPosition[2]*.4*width,
		    BUTTONHEIGHT);
  */
}
  
void MapWindows(void) {
  int buttonnum; 

  controlswin = NewButton(win,"controls",
			  buttonAxesPosition[0]*width,
			  buttonAxesPosition[1]*height,
			  buttonAxesPosition[2]*width,
			  buttonAxesPosition[3]*height,false,white);
  controlspix = XCreatePixmap(dis,controlswin,
			      buttonAxesPosition[2]*width,
			      buttonAxesPosition[3]*height,
			      DefaultDepthOfScreen(screen));
  XMapWindow(dis,controlswin);

  axeswin = NewButton(win,"axes",
		      axesPosition[0]*width,
		      axesPosition[1]*height,
		      axesPosition[2]*width,
		      axesPosition[3]*height,true,black);
  pix = XCreatePixmap(dis,axeswin,
		      axesPosition[2]*width,
		      axesPosition[3]*height,DefaultDepthOfScreen(screen));
  XMapWindow(dis,axeswin);

  messagewin = NewButton(win,"message",
			 axesPosition[0]*width,
			 (1.2*buttonAxesPosition[1]+buttonAxesPosition[3])*height,
			 axesPosition[2]*width,
			 buttonAxesPosition[1]*height,false,black);
  XMapWindow(dis,messagewin);

  for(buttonnum=0;buttonnum<NUMBUTTONS;buttonnum++) {
    controlButtons[buttonnum].butwin=NewButton(controlswin,
			      controlButtons[buttonnum].mapstring,
			      controlButtons[buttonnum].l*buttonAxesPosition[2]*width,
			      controlButtons[buttonnum].b*buttonAxesPosition[2]*height,
			      controlButtons[buttonnum].w*buttonAxesPosition[2]*width,
			      controlButtons[buttonnum].h,false,black);
    XMapWindow(dis,controlButtons[buttonnum].butwin);
  }
}

void DrawControls(dataT *data, int procnum, int numprocs) {
  int buttonnum;
  XPoint *vertices = (XPoint *)malloc(5*sizeof(XPoint));
  
  XSetForeground(dis,gc,black);
  XFillRectangle(dis,controlspix,gc,
		 0,0,
		 buttonAxesPosition[2]*width,
		 buttonAxesPosition[3]*height);

  XSetForeground(dis,gc,white);
  XDrawRectangle(dis,controlspix,gc,0,0,
		 buttonAxesPosition[2]*width,
		 buttonAxesPosition[3]*height);

  for(buttonnum=0;buttonnum<NUMBUTTONS;buttonnum++) 
    DrawButton(controlButtons[buttonnum].butwin,controlButtons[buttonnum].string);

  sprintf(str,"Step: %d of %d",n,data->nsteps);
  DrawHeader(controlButtons[prevwin].butwin,controlButtons[nextwin].butwin,str);

  sprintf(str,"Level: %d of %d",k+1,data->Nkmax);
  DrawHeader(controlButtons[kdownwin].butwin,controlButtons[kupwin].butwin,str);

  if(procplottype==allprocs || numprocs==1)
    sprintf(str,"Showing All %d Procs", numprocs);
  else
    sprintf(str,"Processor: %d of %d",procnum+1,numprocs);

  DrawHeader(controlButtons[prevprocwin].butwin,controlButtons[nextprocwin].butwin,str);
}

Window NewButton(Window parent, char *name, int x, int y, 
		 int buttonwidth, int buttonheight, bool motion, int bordercolor) {
  char c;

  XSetWindowAttributes attr;
  XSizeHints hints;
  Window button;

  attr.background_pixel = black;
  attr.border_pixel = bordercolor;
  attr.backing_store = NotUseful;
  if(motion) 
    attr.event_mask = ButtonReleaseMask | ButtonPressMask | Button1MotionMask ;
  else
    attr.event_mask = ButtonReleaseMask | ButtonPressMask;
  attr.bit_gravity = SouthWestGravity;
  attr.win_gravity = SouthWestGravity;
  attr.save_under = False;
  button = XCreateWindow(dis,parent, x, y, buttonwidth, buttonheight,
                         2, 0, InputOutput, CopyFromParent,
                         CWBackPixel | CWBorderPixel | CWEventMask |
                         CWBitGravity | CWWinGravity | CWBackingStore |
                         CWSaveUnder, &attr);
  hints.width = buttonwidth;
  hints.height = buttonheight - 4;
  hints.min_width = 0;
  hints.min_height = buttonheight - 4;
  hints.max_width = buttonwidth;
  hints.max_height = buttonheight - 4;
  hints.width_inc = 1;
  hints.height_inc = 1;
  hints.flags = PMinSize | PMaxSize | PSize | PResizeInc;
  XSetStandardProperties(dis,button,name,"SUNTANS",None,(char **)NULL,0,&hints);

  return button;
}

void DrawButton(Window button, char *str) {
  int x, y, w, h, d, b, font_width, font_height, border=3;
  Window root;

  font_width=XTextWidth(fontStruct,str,strlen(str));
  font_height=fontStruct->ascent+fontStruct->descent;

  XGetGeometry(dis,button,&root, &x, &y, &w, &h, &d, &b);

  XSetForeground(dis,gc,white);  
  XFillRectangle(dis,button, gc, 0, 0,w,h);

  XSetForeground(dis,gc,black);  
  XDrawRectangle(dis,button, gc,border,border,w-2*border,h-2*border);
  XDrawLine(dis,button,gc,0,0,border,border);
  XDrawLine(dis,button,gc,w,0,w-border,border);
  XDrawLine(dis,button,gc,0,h,border,h-border);
  XDrawLine(dis,button,gc,w,h,w-border,h-border);

  XSetForeground(dis,fontgc,black);  
  XDrawString(dis,button,fontgc,(int)(w/2-font_width/2),(int)font_height,str,strlen(str));  
}

void DrawHeader(Window leftwin,Window rightwin,char *str) {
  int x1, y1, w1, h1, d1, b1, font_width, font_height,
    x2, y2, w2, h2, d2, b2, headsep=5;

  font_width=XTextWidth(fontStruct,str,strlen(str));
  font_height=fontStruct->ascent+fontStruct->descent;

  XGetGeometry(dis,leftwin,&RootWindow(dis,0), &x1, &y1, &w1, &h1, &d1, &b1);
  XGetGeometry(dis,rightwin,&RootWindow(dis,0), &x2, &y2, &w2, &h2, &d2, &b2);

  XSetForeground(dis,fontgc,white);  
  XDrawString(dis,controlspix,
	      fontgc,(int)((x1+x2+w2)/2)-font_width/2,y1-headsep,str,strlen(str));  
}

void DrawZoomBox(void) {
  int x1, y1, boxwidth, boxheight;
  boxwidth=abs(xend-xstart);
  boxheight=abs(yend-ystart);
  if(xend>xstart)
    x1 = xstart;
  else
    x1 = xend;
  if(yend>ystart)
    y1 = ystart;
  else
    y1 = yend;
  XSetForeground(dis,gc,white);
  XDrawRectangle(dis,axeswin,gc,
		 x1,y1,boxwidth,boxheight);
}

void SetUpButtons(void) {
  int buttonnum;
  float dist=0.25;

  controlButtons[prevwin].string="<--";
  controlButtons[prevwin].mapstring="prev";
  controlButtons[prevwin].l=0.05;
  controlButtons[prevwin].b=0.12;
  controlButtons[prevwin].w=0.35;
  controlButtons[prevwin].h=(float)BUTTONHEIGHT;

  controlButtons[gowin].string="M";
  controlButtons[gowin].mapstring="go";
  controlButtons[gowin].l=0.45;
  controlButtons[gowin].b=0.12;
  controlButtons[gowin].w=0.1;
  controlButtons[gowin].h=(float)BUTTONHEIGHT;

  controlButtons[nextwin].string="-->";
  controlButtons[nextwin].mapstring="next";
  controlButtons[nextwin].l=0.6;
  controlButtons[nextwin].b=0.12;
  controlButtons[nextwin].w=0.35;
  controlButtons[nextwin].h=(float)BUTTONHEIGHT;

  controlButtons[kdownwin].string="<--";
  controlButtons[kdownwin].mapstring="kdownwin";
  controlButtons[kdownwin].l=0.05;
  controlButtons[kdownwin].b=controlButtons[nextwin].b+dist;
  controlButtons[kdownwin].w=0.4;
  controlButtons[kdownwin].h=(float)BUTTONHEIGHT;

  controlButtons[kupwin].string="-->";
  controlButtons[kupwin].mapstring="kupwin";
  controlButtons[kupwin].l=0.55;
  controlButtons[kupwin].b=controlButtons[nextwin].b+dist;
  controlButtons[kupwin].w=0.4;
  controlButtons[kupwin].h=(float)BUTTONHEIGHT;

  controlButtons[prevprocwin].string="<--";
  controlButtons[prevprocwin].mapstring="prevprocwin";
  controlButtons[prevprocwin].l=0.05;
  controlButtons[prevprocwin].b=controlButtons[nextwin].b+2*dist;
  controlButtons[prevprocwin].w=0.25;
  controlButtons[prevprocwin].h=(float)BUTTONHEIGHT;

  controlButtons[allprocswin].string="ALL";
  controlButtons[allprocswin].mapstring="allprocswin";
  controlButtons[allprocswin].l=0.375;
  controlButtons[allprocswin].b=controlButtons[nextwin].b+2*dist;
  controlButtons[allprocswin].w=0.25;
  controlButtons[allprocswin].h=(float)BUTTONHEIGHT;

  controlButtons[nextprocwin].string="-->";
  controlButtons[nextprocwin].mapstring="nextprocwin";
  controlButtons[nextprocwin].l=0.7;
  controlButtons[nextprocwin].b=controlButtons[nextwin].b+2*dist;
  controlButtons[nextprocwin].w=0.25;
  controlButtons[nextprocwin].h=(float)BUTTONHEIGHT;

  controlButtons[saltwin].string="Salinity";
  controlButtons[saltwin].mapstring="saltwin";
  controlButtons[saltwin].l=0.05;
  controlButtons[saltwin].b=controlButtons[nextwin].b+3*dist;
  controlButtons[saltwin].w=0.4;
  controlButtons[saltwin].h=(float)BUTTONHEIGHT;

  controlButtons[fswin].string="Surface";
  controlButtons[fswin].mapstring="saltwin";
  controlButtons[fswin].l=0.55;
  controlButtons[fswin].b=controlButtons[nextwin].b+3*dist;
  controlButtons[fswin].w=0.4;
  controlButtons[fswin].h=(float)BUTTONHEIGHT;

  controlButtons[uwin].string="U";
  controlButtons[uwin].mapstring="uwin";
  controlButtons[uwin].l=0.05;
  controlButtons[uwin].b=controlButtons[nextwin].b+4*dist;
  controlButtons[uwin].w=0.1;
  controlButtons[uwin].h=(float)BUTTONHEIGHT;

  controlButtons[vwin].string="V";
  controlButtons[vwin].mapstring="vwin";
  controlButtons[vwin].l=0.2;
  controlButtons[vwin].b=controlButtons[nextwin].b+4*dist;
  controlButtons[vwin].w=0.1;
  controlButtons[vwin].h=(float)BUTTONHEIGHT;

  controlButtons[wwin].string="W";
  controlButtons[wwin].mapstring="wwin";
  controlButtons[wwin].l=0.35;
  controlButtons[wwin].b=controlButtons[nextwin].b+4*dist;
  controlButtons[wwin].w=0.1;
  controlButtons[wwin].h=(float)BUTTONHEIGHT;

  controlButtons[vecwin].string="Vectors";
  controlButtons[vecwin].mapstring="vecwin";
  controlButtons[vecwin].l=0.55;
  controlButtons[vecwin].b=controlButtons[nextwin].b+4*dist;
  controlButtons[vecwin].w=0.4;
  controlButtons[vecwin].h=(float)BUTTONHEIGHT;

  controlButtons[depthwin].string="Depth";
  controlButtons[depthwin].mapstring="depthwin";
  controlButtons[depthwin].l=0.05;
  controlButtons[depthwin].b=controlButtons[nextwin].b+5*dist;
  controlButtons[depthwin].w=0.4;
  controlButtons[depthwin].h=(float)BUTTONHEIGHT;

  controlButtons[nonewin].string="None";
  controlButtons[nonewin].mapstring="nonewin";
  controlButtons[nonewin].l=0.55;
  controlButtons[nonewin].b=controlButtons[nextwin].b+5*dist;
  controlButtons[nonewin].w=0.4;
  controlButtons[nonewin].h=(float)BUTTONHEIGHT;

  controlButtons[edgewin].string="Edges";
  controlButtons[edgewin].mapstring="edgewin";
  controlButtons[edgewin].l=0.05;
  controlButtons[edgewin].b=controlButtons[nextwin].b+6*dist;
  controlButtons[edgewin].w=0.25;
  controlButtons[edgewin].h=(float)BUTTONHEIGHT;

  controlButtons[voronoiwin].string="Voro";
  controlButtons[voronoiwin].mapstring="voronoiwin";
  controlButtons[voronoiwin].l=0.375;
  controlButtons[voronoiwin].b=controlButtons[nextwin].b+6*dist;
  controlButtons[voronoiwin].w=0.25;
  controlButtons[voronoiwin].h=(float)BUTTONHEIGHT;

  controlButtons[delaunaywin].string="Dela";
  controlButtons[delaunaywin].mapstring="delaunaywin";
  controlButtons[delaunaywin].l=0.7;
  controlButtons[delaunaywin].b=controlButtons[nextwin].b+6*dist;
  controlButtons[delaunaywin].w=0.25;
  controlButtons[delaunaywin].h=(float)BUTTONHEIGHT;

  controlButtons[zoomwin].string="ZOOM";
  controlButtons[zoomwin].mapstring="zoomwin";
  controlButtons[zoomwin].l=0.05;
  controlButtons[zoomwin].b=controlButtons[nextwin].b+7*dist;
  controlButtons[zoomwin].w=0.9;
  controlButtons[zoomwin].h=(float)BUTTONHEIGHT;

  controlButtons[quitwin].string="QUIT";
  controlButtons[quitwin].mapstring="quitwin";
  controlButtons[quitwin].l=0.05;
  controlButtons[quitwin].b=controlButtons[nextwin].b+8*dist;
  controlButtons[quitwin].w=0.9;
  controlButtons[quitwin].h=(float)BUTTONHEIGHT;
}

void ParseCommandLine(int N, char *str[], int *numprocs) {
  int i;

  switch(N) {
  case 1:
    *numprocs=1;
    break;
  case 3:
    if(!strcmp(str[1],"-np")) {
      for(i=0;i<strlen(str[2]);i++)
	if(!isdigit(str[2][i]))
	  break;
      *numprocs=atoi(str[2]);
    }
    break;
  default:
    printf("Usage: %s [-np 2]\n",str[0]);
    break;
  }
}

void FindNearest(dataT *data, int x, int y, int *iloc, int *procloc, int procnum, int numprocs) {
  int proc, i, procstart, procend;
  float mindist, dist, xval, yval;

  xval = dataLimits[0]+
    (float)x*(dataLimits[1]-dataLimits[0])/(axesPosition[2]*(float)width);
  yval = dataLimits[2]+
    (float)(height*axesPosition[3]-y)*
    (dataLimits[3]-dataLimits[2])/(axesPosition[3]*(float)height);

  if(procplottype==allprocs) {
    procstart=0;
    procend=numprocs;
  } else {
    procstart=procnum;
    procend=procnum+1;
  }

  for(proc=procstart;proc<procend;proc++) {
    for(i=0;i<data->Nc[proc];i++) {
      dist = sqrt(pow(xval-data->xv[proc][i],2)+pow(yval-data->yv[proc][i],2));
      if(i==0 && proc==procstart) mindist=dist;
      if(dist<mindist) {
	mindist=dist;
	*iloc=i;
	*procloc=proc;
      }
    }
  }
}

dataT *NewData(void) {
  dataT *data = (dataT *)malloc(sizeof(data));
  data->timestep=-1;
  data->klevel=-1;
  return data;
}

void CloseGraphics(void) {
  XFreePixmap(dis,pix);  
  XFreePixmap(dis,controlspix);  
  XDestroyWindow(dis,win);
  XCloseDisplay(dis);
}

void FreeData(dataT *data, int numprocs) {
  int proc, j;
    for(proc=0;proc<numprocs;proc++) {
      free(data->cells[proc]);
      free(data->edges[proc]);
      free(data->xv[proc]);
      free(data->yv[proc]);
      free(data->depth[proc]);
      free(data->h[proc]);
      free(data->h_d[proc]);

      for(j=0;j<data->Nkmax;j++) {
	free(data->s[proc][j]);
	free(data->u[proc][j]);
	free(data->v[proc][j]);
	free(data->wf[proc][j]);
	free(data->w[proc][j]);
      }

      free(data->s[proc]);
      free(data->u[proc]);
      free(data->v[proc]);
      free(data->wf[proc]);
      free(data->w[proc]);
    }
    free(data->cells);
    free(data->edges);
    free(data->xv);
    free(data->yv);
    free(data->depth);
    free(data->h);
    free(data->h_d);
    free(data->s);
    free(data->u);
    free(data->v);
    free(data->wf);
    free(data->w);
    free(data->Ne);
    free(data->Nc);
    free(data->xc);
    free(data->yc);
}

void ReadData(dataT *data, int nstep, int numprocs) {
  int i, j, proc, status, ind, ind0;
  float xind;
  double *dummy;
  char string[BUFFERLENGTH];
  FILE *fid;

  GetString(POINTSFILE,DATAFILE,"points",&status);
  GetString(EDGEFILE,DATAFILE,"edges",&status);
  GetString(CELLSFILE,DATAFILE,"cells",&status);
  GetString(CELLCENTEREDFILE,DATAFILE,"celldata",&status);

  if(nstep==-1) {
    data->Np=getsize(POINTSFILE);

    data->cells = (int **)malloc(numprocs*sizeof(int *));
    data->edges = (int **)malloc(numprocs*sizeof(int *));

    data->xv = (float **)malloc(numprocs*sizeof(float *));
    data->yv = (float **)malloc(numprocs*sizeof(float *));

    data->depth = (float **)malloc(numprocs*sizeof(float *));
    data->h = (float **)malloc(numprocs*sizeof(float *));
    data->h_d = (float **)malloc(numprocs*sizeof(float *));
    data->s = (float ***)malloc(numprocs*sizeof(float **));
    data->u = (float ***)malloc(numprocs*sizeof(float **));
    data->v = (float ***)malloc(numprocs*sizeof(float **));
    data->wf = (float ***)malloc(numprocs*sizeof(float **));
    data->w = (float ***)malloc(numprocs*sizeof(float **));

    data->Ne = (int *)malloc(numprocs*sizeof(int));
    data->Nc = (int *)malloc(numprocs*sizeof(int));
    
    data->timestep=-1;

    data->Nkmax=(int)GetValue("suntans.dat","Nkmax",&status);
    data->nsteps=(int)GetValue(DATAFILE,"nsteps",&status);
    data->numprocs=numprocs;

    for(proc=0;proc<numprocs;proc++) {
      sprintf(string,"%s.%d",EDGEFILE,proc);
      data->Ne[proc]=getsize(string);
      sprintf(string,"%s.%d",CELLSFILE,proc);
      data->Nc[proc]=getsize(string);
    }
    
    data->xc = (float *)malloc(data->Np*sizeof(float *));
    data->yc = (float *)malloc(data->Np*sizeof(float *));

    for(proc=0;proc<numprocs;proc++) {

      data->cells[proc]=(int *)malloc(3*data->Nc[proc]*sizeof(int));
      data->edges[proc]=(int *)malloc(2*data->Ne[proc]*sizeof(int));

      data->xv[proc]=(float *)malloc(data->Nc[proc]*sizeof(float));
      data->yv[proc]=(float *)malloc(data->Nc[proc]*sizeof(float));

      data->depth[proc]=(float *)malloc(data->Nc[proc]*sizeof(float));
      data->h[proc]=(float *)malloc(data->Nc[proc]*sizeof(float));
      data->h_d[proc]=(float *)malloc(data->Nc[proc]*sizeof(float));
      data->s[proc]=(float **)malloc(data->Nkmax*sizeof(float *));
      data->u[proc]=(float **)malloc(data->Nkmax*sizeof(float *));
      data->v[proc]=(float **)malloc(data->Nkmax*sizeof(float *));
      data->wf[proc]=(float **)malloc(data->Nkmax*sizeof(float *));
      data->w[proc]=(float **)malloc(data->Nkmax*sizeof(float *));

      for(i=0;i<data->Nkmax;i++) {
	data->s[proc][i]=(float *)malloc(data->Nc[proc]*sizeof(float));
	data->u[proc][i]=(float *)malloc(data->Ne[proc]*sizeof(float));
	data->v[proc][i]=(float *)malloc(data->Ne[proc]*sizeof(float));
	data->wf[proc][i]=(float *)malloc(data->Ne[proc]*sizeof(float));
	data->w[proc][i]=(float *)malloc(data->Nc[proc]*sizeof(float));
	}
    }

    fid = fopen(POINTSFILE,"r");
    for(i=0;i<data->Np;i++) 
      fscanf(fid,"%f %f %d",&(data->xc[i]),&(data->yc[i]),&ind);
    fclose(fid);  

    for(proc=0;proc<numprocs;proc++) {
      sprintf(string,"%s.%d",CELLSFILE,proc);
      fid = fopen(string,"r");
      for(i=0;i<data->Nc[proc];i++) 
	fscanf(fid,"%f %f %d %d %d %d %d %d",&xind,&xind,
	       &(data->cells[proc][3*i]),
	       &(data->cells[proc][3*i+1]),
	       &(data->cells[proc][3*i+2]),
	       &ind,&ind,&ind);
      fclose(fid);

      sprintf(string,"%s.%d",EDGEFILE,proc);
      fid = fopen(string,"r");
      for(i=0;i<data->Ne[proc];i++) {
	fscanf(fid,"%d %d %d %d %d",
	       &(data->edges[proc][2*i]),
	       &(data->edges[proc][2*i+1]),
	       &ind,&ind,&ind);
      }
      fclose(fid);
      
      sprintf(string,"%s.%d",CELLCENTEREDFILE,proc);
      fid = fopen(string,"r");
      for(i=0;i<data->Nc[proc];i++) {
	fscanf(fid,"%f %f %f %f %d %d %d %d %d %d %d %d %d %d %d",
	       &(data->xv[proc][i]),
	       &(data->yv[proc][i]),
	       &xind,
	       &(data->depth[proc][i]),
	       &ind,&ind,&ind,&ind,&ind,&ind,&ind,&ind,&ind,&ind,&ind);
      }
      fclose(fid);
    }
  } else if(data->timestep != nstep) {
    data->timestep=nstep;

    for(proc=0;proc<numprocs;proc++) {
      dummy=(double *)malloc(data->Ne[proc]*sizeof(double));

      sprintf(string,"/home/fringer/research/SUNTANS/data/s.dat.%d",proc);
      fid = fopen(string,"r");
      fseek(fid,(nstep-1)*data->Nc[proc]*data->Nkmax*sizeof(double),0);
      for(i=0;i<data->Nkmax;i++) {
	fread(dummy,sizeof(double),data->Nc[proc],fid);      
	for(j=0;j<data->Nc[proc];j++) 
	  data->s[proc][i][j]=dummy[j];
      }
      fclose(fid);

      sprintf(string,"/home/fringer/research/SUNTANS/data/u.dat.%d",proc);
      fid = fopen(string,"r");
      fseek(fid,3*(nstep-1)*data->Ne[proc]*data->Nkmax*sizeof(double),0);
      for(i=0;i<data->Nkmax;i++) {
	fread(dummy,sizeof(double),data->Ne[proc],fid);      
	for(j=0;j<data->Nc[proc];j++) 
	  data->u[proc][i][j]=dummy[j];
	fread(dummy,sizeof(double),data->Ne[proc],fid);      
	for(j=0;j<data->Nc[proc];j++) 
	  data->v[proc][i][j]=dummy[j];
	fread(dummy,sizeof(double),data->Ne[proc],fid);      
	for(j=0;j<data->Nc[proc];j++) 
	  data->wf[proc][i][j]=dummy[j];
      }
      fclose(fid);

      sprintf(string,"/home/fringer/research/SUNTANS/data/w.dat.%d",proc);
      fid = fopen(string,"r");
      fseek(fid,(nstep-1)*data->Nc[proc]*data->Nkmax*sizeof(double),0);
      for(i=0;i<data->Nkmax;i++) {
	fread(dummy,sizeof(double),data->Nc[proc],fid);      
	for(j=0;j<data->Nc[proc];j++) 
	  data->w[proc][i][j]=dummy[j];
      }
      fclose(fid);

      sprintf(string,"/home/fringer/research/SUNTANS/data/fs.dat.%d",proc);
      fid = fopen(string,"r");
      fseek(fid,(nstep-1)*data->Nc[proc]*sizeof(double),0);
      fread(dummy,sizeof(double),data->Nc[proc],fid);      
      for(i=0;i<data->Nc[proc];i++)
	data->h[proc][i]=dummy[i];
      fclose(fid);

      for(i=0;i<data->Nc[proc];i++)
	data->h_d[proc][i]=data->h[proc][i]+data->depth[proc][i];

      free(dummy);
    }
  }
}
    
void GetUMagMax(dataT *data, int klevel, int numprocs) {
  int j, proc;
  float umag, umagmax=0, ud, vd;
  
  if(data->klevel!=klevel) {
    data->klevel=klevel;
    for(proc=0;proc<numprocs;proc++) 
      for(j=0;j<data->Ne[proc];j++) {
	ud = data->u[proc][klevel][j];
	vd = data->v[proc][klevel][j];
	umag=sqrt(pow(ud,2)+pow(vd,2));
	if(umag>umagmax) umagmax=umag;
      }
    data->umagmax=umagmax;
  }
}


    
    

    
