/*
 * File: sunplot.c
 * ---------------
 * Description: Plot the output of suntans.
 *
 * Oliver Fringer
 * EFML Stanford University
 *
 * $Id: sunplot.c,v 1.8 2003-04-08 17:18:44 fringer Exp $
 * $Log: not supported by cvs2svn $
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

#define WIDTH 500
#define HEIGHT 500
#define PI 3.141592654
#define EMPTY 99999
#define INFTY 1e20
#define BUFFER 256
#define CMAPFILE "jet.cmap"
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
#define NUMBUTTONS 20
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
  quitwin
} buttonName;

typedef enum {
  oneproc, allprocs
} plotProcT;

typedef struct {
  char *string;
  char *mapstring;
  float l,b,w,h;
  Window butwin;
} myButtonT;

void InitializeGraphics(void);
void MyDraw(char plottype, int procnum);
void ReadGrid(int proc);
void ReadScalar(float *scal, int k, int Nk, int np, char *filename);
void ReadVelocity(float *u, float *v, int kp, int Nk, int np , char *filename);
float Min(float *x, int N);
float Max(float *x, int N);
void AxisImage(float *axes, float *data);
void Fill(XPoint *vertices, int N, int cindex, int edges);
void CAxis(float *caxis, float *data, int Nc);
void ReadColorMap(char *str);
void UnSurf(float *xc, float *yc, int *cells, float *data, int N);
void SetDataLimits(void);
void SetAxesPosition(void);
void Clf(void);
void Cla(void);
void Text(Window window, float x, float y, int boxwidth, int boxheight,
	  char *str, int fontsize, int color, 
	  hjustifyT hjustify, vjustifyT vjustify);
void DrawControls(void);
Window NewButton(Window parent, char *name, int x, int y, 
		 int buttonwidth, int buttonheight, bool motion);
void DrawButton(Window button, char *str);
void MapWindows(void);
void RedrawWindows(void);
void DrawZoomBox(void);
void DrawHeader(Window leftwin,Window rightwin,char *str);
void SetUpButtons(void);
void DrawVoronoiPoints(float *xv, float *yv, int Nc);
void FillCircle(int xp, int yp, int r, int ic, Window mywin);
void UnQuiver(float *xc, float *yc, int *edges, float *u, float *v, int Ne);
void DrawArrow(int xp, int yp, int ue, int ve, Window mywin, int ic);
void Rotate(XPoint *points, int N, int ue, int ve, int mag);

/*
 * Linux users will need to add -ldl to the Makefile to compile 
 * this example.
 *
 */
Display *dis;
myButtonT controlButtons[NUMBUTTONS];
Window win, axeswin, controlswin;
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
  plottype='h',Np, Nc, Ne, n=1, nsteps=21, k=1, Nkmax=20, keysym,
  xstart, ystart, xend, yend;
int *cells, *edges;
float caxis[2], *xc, *yc, *depth, *xp, axesPosition[4], dataLimits[4], buttonAxesPosition[4],
  zoomratio, *xv, *yv;
int axisType, oldaxisType, white, black, red, blue, green, yellow, 
  procnum=0, numprocs=1,colors[NUMCOLORS];
bool edgelines, setdatalimits, pressed,   voronoipoints, delaunaypoints, vectorplot, go;
char str[BUFFER];
zoomT zoom;
plotProcT procplottype;

int main() {
  bool redraw;
  buttonnumT mousebutton;
  setdatalimits = false;
  vectorplot = true;
  zoomratio = 1;
  procplottype = allprocs; // or allprocs

  InitializeGraphics();

  ReadColorMap(CMAPFILE);

  XSelectInput(dis, win, ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask );

  k=Nkmax/2;
  
  axisType='i';
  edgelines=true;
  voronoipoints=false;
  delaunaypoints=false;
  ReadGrid(procnum);

  SetUpButtons();

  MapWindows();
  MyDraw(plottype,procnum);

  XMaskEvent(dis, ExposureMask, &report);
  pressed=false;
  go=false;
  while(true) {
    redraw=false;
    zoom=none;
    if(go==true) {
      if(n<nsteps) {
	redraw=true;
	n++;
      } else {
	redraw=false;
	go=false;
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
	  else { redraw = false ; printf("At n=1!\n"); }
	} else if(mousebutton==right_button)
	  if(n!=1) { redraw = true ; n=1; }
      } else if(report.xany.window==controlButtons[gowin].butwin) {
	if(mousebutton==left_button) {
	  go=true;
	}
      } else if(report.xany.window==controlButtons[nextwin].butwin) {
	if(mousebutton==left_button) {
	  if(n<nsteps) { redraw = true; n++; }
	  else { redraw = false ; printf("At n=nsteps!\n"); }
	} else if(mousebutton==right_button) 
	  if(n!=nsteps) { redraw = true; n=nsteps; }
      } else if(report.xany.window==controlButtons[kdownwin].butwin) {
	if(mousebutton==left_button) {
	  if(k>1) { redraw = true; k--; }
	  else { redraw = false ; printf("At k=1!\n"); }
	} else if(mousebutton==right_button)
	  if(k!=1) { redraw = true ; k=1; }
      } else if(report.xany.window==controlButtons[kupwin].butwin) {
	if(mousebutton==left_button) {
	  if(k<Nkmax) { redraw = true; k++; }
	  else { redraw = false ; printf("At k=Nkmax!\n"); }
	} else if(mousebutton==right_button) 
	  if(k!=Nkmax) { redraw = true; k=Nkmax; }
      } else if(report.xany.window==controlButtons[prevprocwin].butwin && numprocs>1) {
	if(mousebutton==left_button) {
	  if(procnum>0) { redraw = true; procnum--; procplottype=oneproc; }
	  else { redraw = false ; printf("At procnum==0!\n"); }
	} else if(mousebutton==right_button) 
	  if(procnum!=0) { redraw = true; procnum=0; procplottype=oneproc; }
      } else if(report.xany.window==controlButtons[allprocswin].butwin
		&& mousebutton==left_button && numprocs>1) {
	if(procplottype!=allprocs) redraw=true; 
	procplottype=allprocs; 
      } else if(report.xany.window==controlButtons[nextprocwin].butwin && numprocs>1) {
	if(mousebutton==left_button) {
	  if(procnum<numprocs-1) { redraw = true; procnum++; procplottype=oneproc; }
	  else { redraw = false ; printf("At procnum==numprocs!\n"); }
	} else if(mousebutton==right_button) 
	  if(procnum!=numprocs) { redraw = true; procnum=numprocs-1; procplottype=oneproc; }
      } else if(report.xany.window==controlButtons[saltwin].butwin && mousebutton==left_button) {
	if(plottype!='s') {
	  plottype='s';
	  printf("Salinity selected...\n");
	  redraw=true;
	}
      } else if(report.xany.window==controlButtons[fswin].butwin && mousebutton==left_button) {
	if(plottype!='h') {
	  plottype='h';
	  printf("Free-surface selected...\n");
	  redraw=true;
	}
      } else if(report.xany.window==controlButtons[vecwin].butwin && mousebutton==left_button) {
	if(vectorplot==false) {
	  vectorplot=true;
	  printf("Vectors on...\n");
	} else {
	  vectorplot=false;
	  printf("Vectors off...\n");
	}
	redraw=true;
      } else if(report.xany.window==controlButtons[depthwin].butwin && mousebutton==left_button) {
	if(plottype!='d') {
	  plottype='d';
	  printf("Depth selected...\n");
	  redraw=true;
	}
      } else if(report.xany.window==controlButtons[depthwin].butwin && mousebutton==right_button) {
	if(plottype!='D') {
	  plottype='D';
	  printf("h+d selected...\n");
	  redraw=true;
	}
      } else if(report.xany.window==controlButtons[edgewin].butwin && mousebutton==left_button) {
	if(edgelines==false) {
	  printf("Drawing edge lines\n");
	  edgelines=true;
	} else {
	  printf("Removing edge lines\n");
	  edgelines=false;
	}
	redraw=true;
      } else if(report.xany.window==controlButtons[voronoiwin].butwin 
		&& mousebutton==left_button) {
	if(voronoipoints==false) {
	  printf("Drawing Voronoi points\n");
	  voronoipoints=true;
	} else {
	  printf("Removing Voronoi points\n");
	  voronoipoints=false;
	}
	redraw=true;
      } else if(report.xany.window==controlButtons[delaunaywin].butwin 
		&& mousebutton==left_button) {
	if(delaunaypoints==false) {
	  printf("Drawing Delaunay points\n");
	  delaunaypoints=true;
	} else {
	  printf("Removing Delaunay points\n");
	  delaunaypoints=false;
	}
	redraw=true;
      } else if(report.xany.window==controlButtons[nonewin].butwin && mousebutton==left_button) {
	if(plottype!='n') {
	  plottype='n';
	  printf("Removing surface plot...\n");
	  redraw=true;
	}
      } else if(report.xany.window==controlButtons[quitwin].butwin && mousebutton==left_button) {
	exit(0);
      } else if(report.xany.window==axeswin) {
	xend=report.xbutton.x;
	yend=report.xbutton.y;
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
      MyDraw(plottype,procnum);
      break;
    case KeyPress:
      switch(keysym=XLookupKeysym(&report.xkey, 0)) {
      case XK_q:
	exit(0);
	break;
      case XK_k: case XK_Up:
	if(k<Nkmax) { redraw = true; k++; }
	  else { redraw = false; printf("At k=Nkmax!\n"); }
	break;
      case XK_j: case XK_Down:
	if(k>1) { redraw = true; k--; }
	else { redraw = false ; printf("At k=1!\n"); }
	break;
      case XK_p : case XK_Left:
	if(n>1) { redraw = true; n--; }
	else { redraw = false ; printf("At n=1!\n"); }
	break;
      case XK_n: case XK_Right:
	if(n<nsteps) { redraw = true; n++; }
	else { redraw = false ; printf("At n=nsteps!\n"); }
	break;
      case XK_s:
	if(plottype!='s') {
	  plottype='s';
	  printf("Salinity selected...\n");
	  redraw=true;
	}
	break;
      case XK_h:
	if(plottype!='h') {
	  plottype='h';
	  printf("Free-surface selected...\n");
	  redraw=true;
	}
	break;
      case XK_d:
	if(plottype!='d') {
	  plottype='d';
	  printf("Depth selected...\n");
	  redraw=true;
	}
	break;
      case XK_a:
	if(axisType=='n') {
	  printf("Changing plottype to image\n");
	  axisType='i';
	} else {
	  printf("Changing plottype to normal\n");
	  axisType='n';
	}
	redraw=true;
	break;
      case XK_e:
	if(edgelines==0) {
	  printf("Drawing edge lines\n");
	  edgelines=1;
	} else {
	  printf("Removing edge lines\n");
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
	  ReadGrid(procnum);
	  redraw=true;
	}
	break;
      }
    }
    }
    if(redraw) 
      MyDraw(plottype,procnum);
  }
  return 0;
}

void ReadGrid(int proc) {
  int i, ind;
  FILE *tfile,*dfile,*ifile = fopen("jet.cmap","r"), *efile,
    *pfile = fopen("/home/fringer/research/SUNTANS/data/points.dat","r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/cells.dat.%d",proc);
  tfile = fopen(str,"r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/celldata.dat.%d",proc);
  dfile = fopen(str,"r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/edges.dat.%d",proc);
  efile = fopen(str,"r");

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
  printf("Ne = %d\n",Ne);

  xc = (float *)malloc(Np*sizeof(float));
  yc = (float *)malloc(Np*sizeof(float));
  xv = (float *)malloc(Nc*sizeof(float));
  yv = (float *)malloc(Nc*sizeof(float));
  cells = (int *)malloc(3*Nc*sizeof(int));
  edges = (int *)malloc(2*Ne*sizeof(int));
  depth = (float *)malloc(Nc*sizeof(float));

  for(i=0;i<Np;i++) {
    fscanf(pfile,"%f %f %d\n",&xc[i],&yc[i],&ind);
  }
  fclose(pfile);

  for(i=0;i<Ne;i++) 
    fscanf(efile,"%d %d %d %d %d\n",&edges[2*i],&edges[2*i+1],&ind,&ind,&ind);

  for(i=0;i<Nc;i++) {
    fscanf(dfile,"%f %f %f %f %d %d %d %d %d %d %d %d %d %d %d\n",
	   &xv[i],&yv[i],&xp,&depth[i],&ind,&ind,&ind,&ind,
	   &ind,&ind,&ind,&ind,&ind,&ind,&ind);
    fscanf(tfile,"%f %f %d %d %d %d %d %d\n",&xp,&xp,&cells[3*i],&cells[3*i+1],&cells[3*i+2],
	   &ind,&ind,&ind);
  }
  fclose(tfile);
  fclose(dfile);
  fclose(efile);
}

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

void MyDraw(char plottype, int procnum)
{
  int i;
  float *scal, *u, *v;
  FILE *sfile, *hfile, *ufile;

  sprintf(str,"/home/fringer/research/SUNTANS/data/s.dat.%d",procnum);
  sfile = fopen(str,"r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/h.dat.%d",procnum);
  hfile = fopen(str,"r");
  sprintf(str,"/home/fringer/research/SUNTANS/data/u.dat.%d",procnum);
  ufile = fopen(str,"r");

  scal = (float *)malloc(Nc*sizeof(float));
  u = (float *)malloc(Ne*sizeof(float));
  v = (float *)malloc(Ne*sizeof(float));

  if(vectorplot)
    ReadVelocity(u,v,k,Nkmax,n,"/home/fringer/research/SUNTANS/data/u.dat.0");

  switch(plottype) {
  case 's':
    ReadScalar(scal,k,Nkmax,n,"/home/fringer/research/SUNTANS/data/s.dat.0");
    break;
  case 'h':
    ReadScalar(scal,1,1,n,"/home/fringer/research/SUNTANS/data/fs.dat.0");
    break;
  case 'd':
    for(i=0;i<Nc;i++)
      scal[i]=depth[i];
    break;
  case 'D':
    ReadScalar(scal,1,1,n,"/home/fringer/research/SUNTANS/data/fs.dat.0");
    for(i=0;i<Nc;i++)
      scal[i]+=depth[i];
    break;
  }
  
  SetDataLimits();
  
  SetAxesPosition();
  
  CAxis(caxis,scal,Nc);
  
  UnSurf(xc,yc,cells,scal,Nc);
  
  if(voronoipoints) 
    DrawVoronoiPoints(xv,yv,Nc);

  if(vectorplot)
    UnQuiver(xc,yc,edges,u,v,Ne);

  DrawControls();
  XFlush(dis);

  XCopyArea(dis,pix,axeswin,gc,0,0,
	    axesPosition[2]*width,
	    axesPosition[3]*height,0,0);
  XCopyArea(dis,controlspix,controlswin,gc,0,0,
	    buttonAxesPosition[2]*width,
	    buttonAxesPosition[3]*height,0,0);

  free(scal);
  free(u);
  free(v);
}

void UnQuiver(float *xc, float *yc, int *edges, float *u, float *v, int Ne) {
  int j, ic;
  float xe, ye, umag, umagmax, l, lmax, n1, n2;
  int xp, yp, ue, ve, vlength=100;

  if(plottype=='n')
    ic = white;
  else
    ic = black;

  lmax = 0;
  for(j=0;j<Ne;j++) {
    l=sqrt(pow(xc[edges[2*j]]-xc[edges[2*j+1]],2)+
	   pow(yc[edges[2*j]]-yc[edges[2*j+1]],2));
    if(l>lmax) lmax=l;
  }
  vlength = (int)(lmax/(dataLimits[1]-dataLimits[0])*
		  axesPosition[2]*width);

  umagmax=0;
  for(j=0;j<Ne;j++) {
    umag=sqrt(u[j]*u[j]+v[j]*v[j]);
    if(umag>umagmax) umagmax=umag;
  }

  for(j=0;j<Ne;j++) {
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
  if(plottype!='n') {
    XSetForeground(dis,gc,colors[cindex]);
    XFillPolygon(dis,pix,gc,vertices,3,Convex,CoordModeOrigin);
  }
  if(edges) {
    if(plottype=='n')
      XSetForeground(dis,gc,white);
    else
	XSetForeground(dis,gc,black);
    XDrawLines(dis,pix,gc,vertices,4,0);
  }
  if(delaunaypoints) {
    if(plottype=='n')
      ic = white;
    else
      ic = black;
    for(i=0;i<N-1;i++)
      FillCircle(vertices[i].x,vertices[i].y,POINTSIZE,ic,pix);
  }
}

void CAxis(float *caxis, float *data, int N) {
  int i;

  caxis[0] = EMPTY;
  caxis[1] = -EMPTY;

  for(i=0;i<N;i++) {
    if(data[i]<=caxis[0] && data[i]!=0) caxis[0]=data[i];
    if(data[i]>=caxis[1] && data[i]!=0) caxis[1]=data[i];
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

void UnSurf(float *xc, float *yc, int *cells, float *data, int N) {
  int i, j, ind;
  XPoint *vertices = (XPoint *)malloc(4*sizeof(XPoint));

  Clf();

  for(i=0;i<N;i++) {
    for(j=0;j<3;j++) {
      vertices[j].x = 
	axesPosition[2]*width*(xc[cells[3*i+j]]-dataLimits[0])/
	(dataLimits[1]-dataLimits[0]);
      vertices[j].y = 
	axesPosition[3]*height*(1-(yc[cells[3*i+j]]-dataLimits[2])/
	(dataLimits[3]-dataLimits[2]));
      /*
      vertices[j].x = axesPosition[0]*width+
	axesPosition[2]*width*(xc[cells[3*i+j]]-dataLimits[0])/
	(dataLimits[1]-dataLimits[0]);
      vertices[j].y = axesPosition[1]*height+
	axesPosition[3]*height*(1-(yc[cells[3*i+j]]-dataLimits[2])/
	(dataLimits[3]-dataLimits[2]));
      */
    }

    vertices[3].x = vertices[0].x;
    vertices[3].y = vertices[0].y;

    ind = (data[i]-caxis[0])/(caxis[1]-caxis[0])*(NUMCOLORS-2);
    if(data[i]==EMPTY || (plottype=='D' && data[i]==0))
      ind = NUMCOLORS-1;

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
  printf("Screen is %d by %d pixels.\n",width,height);
  
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

void SetDataLimits(void) {
  float dx, dy, Xc, Yc, x1, y1, x2, y2, a1, a2, xmin, xmax, ymin, ymax;
  dx = dataLimits[1]-dataLimits[0];
  dy = dataLimits[3]-dataLimits[2];

  if(!setdatalimits) {
    dataLimits[0] = Min(xc,Np);
    dataLimits[1] = Max(xc,Np);
    dataLimits[2] = Min(yc,Np);
    dataLimits[3] = Max(yc,Np);
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

void Clf(void) {
  XSetForeground(dis,gc,black);
  XFillRectangle(dis,pix,gc,0,0,width*axesPosition[2],height*axesPosition[3]);
}

void Cla(void) {
  XSetForeground(dis,gc,black);
  XFillRectangle(dis,pix,gc,
		 0,0,
		 axesPosition[2]*width,
		 axesPosition[3]*height);
}

void Text(Window window, float x, float y, int boxwidth, int boxheight,
	  char *str, int fontsize, int color, 
	  hjustifyT hjustify, vjustifyT vjustify) {

  char **fonts;
  XFontStruct ** fontsReturn;
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
			  buttonAxesPosition[3]*height,false);
  controlspix = XCreatePixmap(dis,controlswin,
			      buttonAxesPosition[2]*width,
			      buttonAxesPosition[3]*height,
			      DefaultDepthOfScreen(screen));
  XMapWindow(dis,controlswin);

  axeswin = NewButton(win,"axes",
		      axesPosition[0]*width,
		      axesPosition[1]*height,
		      axesPosition[2]*width,
		      axesPosition[3]*height,true);
  pix = XCreatePixmap(dis,axeswin,
		      axesPosition[2]*width,
		      axesPosition[3]*height,DefaultDepthOfScreen(screen));
  XMapWindow(dis,axeswin);
  /*
  prevwin = NewButton(controlswin,"prev",
		      (0*buttonAxesPosition[0]+.05*buttonAxesPosition[2])*width,
		      (0*buttonAxesPosition[1]+.12*buttonAxesPosition[2])*height,
		      buttonAxesPosition[2]*.4*width,
		      BUTTONHEIGHT,false);
  XMapWindow(dis,prevwin);
  */
  for(buttonnum=0;buttonnum<NUMBUTTONS;buttonnum++) {
    controlButtons[buttonnum].butwin=NewButton(controlswin,
			      controlButtons[buttonnum].mapstring,
			      controlButtons[buttonnum].l*buttonAxesPosition[2]*width,
			      controlButtons[buttonnum].b*buttonAxesPosition[2]*height,
			      controlButtons[buttonnum].w*buttonAxesPosition[2]*width,
			      controlButtons[buttonnum].h,false);
    XMapWindow(dis,controlButtons[buttonnum].butwin);
  }
  /*
  nextwin = NewButton(controlswin,"next",
		      (0*buttonAxesPosition[0]+.55*buttonAxesPosition[2])*width,
		      (0*buttonAxesPosition[1]+.12*buttonAxesPosition[2])*height,
		      buttonAxesPosition[2]*.4*width,
		      BUTTONHEIGHT,false);
  XMapWindow(dis,nextwin);

  kdownwin = NewButton(controlswin,"kup",
		      (0*buttonAxesPosition[0]+.05*buttonAxesPosition[2])*width,
		      (0*buttonAxesPosition[1]+.37*buttonAxesPosition[2])*height,
		      buttonAxesPosition[2]*.4*width,
		      BUTTONHEIGHT,false);
  XMapWindow(dis,kdownwin);

  kupwin = NewButton(controlswin,"kdown",
		      (0*buttonAxesPosition[0]+.55*buttonAxesPosition[2])*width,
		      (0*buttonAxesPosition[1]+.37*buttonAxesPosition[2])*height,
		      buttonAxesPosition[2]*.4*width,
		      BUTTONHEIGHT,false);
  XMapWindow(dis,kupwin);
  */
}

void DrawControls(void) {
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

  sprintf(str,"Step: %d of %d",n,nsteps);
  DrawHeader(controlButtons[prevwin].butwin,controlButtons[nextwin].butwin,str);

  sprintf(str,"Level: %d of %d",k,Nkmax);
  DrawHeader(controlButtons[kdownwin].butwin,controlButtons[kupwin].butwin,str);

  if(procplottype==allprocs || numprocs==1)
    sprintf(str,"Showing All %d Procs", numprocs);
  else
    sprintf(str,"Processor: %d of %d",procnum+1,numprocs);

  DrawHeader(controlButtons[prevprocwin].butwin,controlButtons[nextprocwin].butwin,str);
}

Window NewButton(Window parent, char *name, int x, int y, 
		 int buttonwidth, int buttonheight, bool motion) {
  char c;

  XSetWindowAttributes attr;
  XSizeHints hints;
  Window button;

  attr.background_pixel = black;
  attr.border_pixel = white;
  attr.backing_store = NotUseful;
  if(motion) 
    attr.event_mask = ExposureMask | ButtonReleaseMask | ButtonPressMask | Button1MotionMask ;
  else
    attr.event_mask = ExposureMask | ButtonReleaseMask | ButtonPressMask;
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
  int x, y, w, h, d, b, font_width, font_height;
  Window root;

  font_width=XTextWidth(fontStruct,str,strlen(str));
  font_height=fontStruct->ascent+fontStruct->descent;

  XGetGeometry(dis,button,&root, &x, &y, &w, &h, &d, &b);

  XSetForeground(dis,gc,white);  
  XFillRectangle(dis,button, gc, 0, 0,w,h);
  
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

  controlButtons[quitwin].string="QUIT";
  controlButtons[quitwin].mapstring="quitwin";
  controlButtons[quitwin].l=0.05;
  controlButtons[quitwin].b=controlButtons[nextwin].b+7*dist;
  controlButtons[quitwin].w=0.9;
  controlButtons[quitwin].h=(float)BUTTONHEIGHT;
}
