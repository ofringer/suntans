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

void mydraw(char plottype, int procnum);
void ReadGrid(int proc);
void ReadScalar(float *scal, int k, int Nk, int np, char *filename);
float Min(float *x, int N);
float Max(float *x, int N);
void AxisImage(float *axes, float *data);
void Fill(XPoint *vertices, int N, int cindex, int edges);

/*Linux users will need to add -ldl to the Makefile to compile 
 *this example.
 */
Display *dis;
Window win;
XEvent report;
GC gc;
XColor color;
Screen *screen;
Pixmap pix;
Colormap colormap;

int plottype='s',Np, Nc, n=1, nsteps=21, k=1, Nkmax=70;
int *cells;
float cmap[64][3], *xc, *yc, *depth, *xp;
int axisType, oldaxisType, edgelines, white, black, procnum=0;
char str[BUFFER];

int main() {
  int i, ind, redraw;
  FILE *ifile = fopen("jet.cmap","r");
  float cmap[64][3], caxis[2], val;
  for(i=0;i<64;i++) {
    fscanf(ifile,"%f %f %f\n",&cmap[i][0],&cmap[i][1],&cmap[i][2]);
  }
  fclose(ifile);
  caxis[0] = 0;
  caxis[1] = 1;

  dis = XOpenDisplay(NULL);
  win = XCreateSimpleWindow(dis, RootWindow(dis, 0), 
			    1, 1, WIDTH, HEIGHT, 0, BlackPixel (dis, 0), BlackPixel(dis, 0));
  screen = ScreenOfDisplay(dis,0);
  pix = XCreatePixmap(dis,win,WIDTH,HEIGHT,DefaultDepthOfScreen(screen));

  XMapWindow(dis, win);

  colormap = DefaultColormap(dis, 0);
  gc = XCreateGC(dis, win, 0, 0);
  //  black = BlackPixel(dis,screen);
  white = WhitePixel(dis,0);

  XSelectInput(dis, win, ExposureMask | KeyPressMask | ButtonPressMask);

  k=Nkmax/2;
  
  axisType='i';
  edgelines=1;
  ReadGrid(procnum);
  //  mydraw(plottype,procnum);
  //  XCopyArea(dis,pix,win,gc,0,0,WIDTH,HEIGHT,0,0);

  while (1)  {
    redraw=0;
    XNextEvent(dis, &report);
    switch  (report.type) {
    case Expose:   
      printf("Exposed...\n");
      mydraw(plottype,procnum);
      XCopyArea(dis,win,pix,gc,0,0,WIDTH,HEIGHT,0,0);
      break;
    case KeyPress:
      switch(XLookupKeysym(&report.xkey, 0)) {
      case XK_q:
	exit(0);
	break;
      case XK_k: case XK_Up:
	if(k<Nkmax) { redraw = 1; k++; }
	  else { redraw = 0; printf("At k=Nkmax!\n"); }
	break;
      case XK_j: case XK_Down:
	if(k>1) { redraw = 1; k--; }
	else { redraw = 0 ; printf("At k=1!\n"); }
	break;
      case XK_p : case XK_Left:
	if(n>1) { redraw = 1; n--; }
	else { redraw = 0 ; printf("At n=1!\n"); }
	break;
      case XK_n: case XK_Right:
	if(n<nsteps) { redraw = 1; n++; }
	else { redraw = 0 ; printf("At n=nsteps!\n"); }
	break;
      case XK_s:
	if(plottype!='s') {
	  plottype='s';
	  printf("Salinity selected...\n");
	  redraw=1;
	}
	break;
      case XK_h:
	if(plottype!='h') {
	  plottype='h';
	  printf("Free-surface selected...\n");
	  redraw=1;
	}
	break;
      case XK_d:
	if(plottype!='d') {
	  plottype='d';
	  printf("Depth selected...\n");
	  redraw=1;
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
	redraw=1;
	break;
      case XK_e:
	if(edgelines==0) {
	  printf("Drawing edge lines\n");
	  edgelines=1;
	} else {
	  printf("Removing edge lines\n");
	  edgelines=0;
	}
	redraw=1;
	break;
      default:
	if(isdigit(XLookupKeysym(&report.xkey, 0))) {
	  str[0]=XLookupKeysym(&report.xkey, 0);
	  str[1]='\0';
	  if(procnum!=atoi(str))
	    procnum=atoi(str);
	    ReadGrid(procnum);
	    redraw=1;
	  }
	break;
      }
      if(redraw) {
	mydraw(plottype,procnum);
	XCopyArea(dis,win,pix,gc,0,0,WIDTH,HEIGHT,0,0);
      }
      break;
    }
  }
  return 0;
}

void ReadGrid(int proc) {
  int i, ind;
  char str1[256];
  FILE *tfile,*dfile,*ifile = fopen("jet.cmap","r"), 
    *pfile = fopen("/home/fringer/research/SUNTANS/data/points.dat","r");
  sprintf(str1,"/home/fringer/research/SUNTANS/data/cells.dat.%d",proc);
  tfile = fopen(str1,"r");
  sprintf(str1,"/home/fringer/research/SUNTANS/data/celldata.dat.%d",proc);
  dfile = fopen(str1,"r");

  Np=0;
  while(fgets(str1,256,pfile)) Np++;
  fclose(pfile);
  pfile = fopen("/home/fringer/research/SUNTANS/data/points.dat","r"),

  Nc=0;
  while(fgets(str1,256,tfile)) Nc++;
  fclose(tfile);
  sprintf(str1,"/home/fringer/research/SUNTANS/data/cells.dat.%d",proc);
  tfile = fopen(str1,"r"),

  xc = (float *)malloc(Np*sizeof(float));
  yc = (float *)malloc(Np*sizeof(float));
  cells = (int *)malloc(3*Nc*sizeof(int));
  depth = (float *)malloc(Nc*sizeof(float));

  for(i=0;i<64;i++) {
    fscanf(ifile,"%f %f %f\n",&cmap[i][0],&cmap[i][1],&cmap[i][2]);
  }
  fclose(ifile);
  for(i=0;i<Np;i++) {
    fscanf(pfile,"%f %f %d\n",&xc[i],&yc[i],&ind);
  }
  fclose(pfile);

  for(i=0;i<Nc;i++) {
    fscanf(dfile,"%f %f %f %f %d %d %d %d %d %d %d %d %d %d %d\n",
	   &xp,&xp,&xp,&depth[i],&ind,&ind,&ind,&ind,
	   &ind,&ind,&ind,&ind,&ind,&ind,&ind);
    fscanf(tfile,"%f %f %d %d %d %d %d %d\n",&xp,&xp,&cells[3*i],&cells[3*i+1],&cells[3*i+2],
	   &ind,&ind,&ind);
  }
  fclose(tfile);
  fclose(dfile);
}

void ReadScalar(float *scal, int k, int Nk, int np, char *filename) {
  int i, kp, currptr;
  FILE *sfile = fopen(filename,"r");
  double *dum = (double *)malloc(Nc*sizeof(double));
  float dummy;

  for(i=0;i<Nc;i++)
    scal[i]=EMPTY;

  if(Nk==1) {
    /*
    for(np=1;np<n;np++)
      for(i=0;i<Nc;i++)
	fscanf(sfile,"%f\n",&dummy);
    */
    fseek(sfile,(np-1)*Nc*sizeof(float),SEEK_SET);
    /*
  for(i=0;i<Nc;i++) 
    fscanf(sfile,"%f\n",&scal[i]);
    */
    fread(dum,sizeof(double),Nc,sfile);
  } else {
    /*
    for(np=1;np<n;np++)
      for(i=0;i<Nc;i++)
	for(kp=0;kp<Nk;kp++)
	  fscanf(sfile,"%f\n",&dummy);
    */
    //    printf("np = %d\n",np);
    currptr=fseek(sfile,(np-1)*Nc*Nk*sizeof(double),SEEK_CUR);
    //    printf("Current pointer = %d\n",ftell(sfile));
    /*
    for(i=0;i<Nc;i++) {
      for(kp=0;kp<k;kp++)
	fscanf(sfile,"%f\n",&dummy);
      fscanf(sfile,"%f\n",&scal[i]);
      for(kp=k+1;kp<Nk;kp++)
	fscanf(sfile,"%f\n",&dummy);
    }
    */
    //    printf("Nk = %d, Nc= %d, k=%d\n",Nk,Nc,k);
    currptr=fseek(sfile,Nc*(k-1)*sizeof(double),SEEK_CUR);
    //    printf("Current pointer = %d\n",ftell(sfile));
    fread(dum,sizeof(double),Nc,sfile);
    for(i=0;i<Nc;i++)
      scal[i]=dum[i];
  }
  fclose(sfile);
}

void mydraw(char plottype, int procnum)
{
  int x1, y1, x2, y2, N=3, mode=1, x, y, width, height, font_height, Nx, Ny, kp, np;
  XPoint *vertices = (XPoint *)malloc(4*sizeof(XPoint));
  char  str1[256];
  XFontStruct *font_info;
  int i, ind, j;
  FILE *sfile, *hfile;
  float caxis[2], val, dx, dy, xp, yp, xmax, xmin, ymax, ymin, zmax, zmin,
    xp1, yp1, *salt, axesPosition[4], dataLimits[4];

  sprintf(str1,"/home/fringer/research/SUNTANS/data/s.dat.%d",procnum);
  sfile = fopen(str1,"r");
  sprintf(str1,"/home/fringer/research/SUNTANS/data/h.dat.%d",procnum);
  hfile = fopen(str1,"r");

  salt = (float *)malloc(Nc*sizeof(float));

  switch(plottype) {
  case 's':
    ReadScalar(salt,k,Nkmax,n,"/home/fringer/research/SUNTANS/data/s.dat.0");
    break;
  case 'h':
    ReadScalar(salt,k,1,n,"/home/fringer/research/SUNTANS/data/fs.dat.0");
    break;
  case 'd':
    for(i=0;i<Nc;i++)
      salt[i]=depth[i];
    break;
  }

  caxis[0] = EMPTY;
  caxis[1] = -EMPTY;

  strcpy(str1,"RECTANGLE"); 

  vertices[0].x = 50;
  vertices[0].y = 50;
  vertices[1].x = 200;
  vertices[1].y = 50;
  vertices[2].x = 75;
  vertices[2].y = 100;

  color.flags = DoRed | DoGreen | DoBlue;

  dataLimits[0] = Min(xc,Np);
  dataLimits[1] = Max(xc,Np);
  dataLimits[2] = Min(yc,Np);
  dataLimits[3] = Max(yc,Np);
  //  printf("Xmin = %f, Xmax = %f, Ymin = %f, Ymax = %f\n",
  //	 dataLimits[0],dataLimits[1],dataLimits[2],dataLimits[3]);

  axesPosition[0] = 0.1;
  axesPosition[1] = 0.1;
  axesPosition[2] = 0.8;
  axesPosition[3] = 0.8;

  if(axisType=='i')
    AxisImage(axesPosition,dataLimits);

  Nx = 500;
  Ny = 500;
  dx = (float)(xmax-xmin)/Ny;
  dy = (float)(ymax-ymin)/Nx;
  zmax = 0;
  zmin = 0;
  for(i=0;i<Nc;i++) {
    //salt[i]=depth[i];
    if(salt[i]<=caxis[0] && salt[i]!=0) caxis[0]=salt[i];
    if(salt[i]>=caxis[1] && salt[i]!=0) caxis[1]=salt[i];
  }

  //  printf("%f %f %f\n",caxis[0],caxis[1],caxis[0]-caxis[1]);

  color.red = 0;
  color.green = 0;
  color.blue = 0;
  XAllocColor(dis, colormap, &color);
  XSetForeground(dis,gc,color.pixel);
  XFillRectangle(dis,win,gc,0,0,WIDTH,HEIGHT);
  for(i=0;i<Nc;i++) {
    for(j=0;j<3;j++) {
      vertices[j].x = axesPosition[0]*WIDTH+
	axesPosition[2]*WIDTH*(xc[cells[3*i+j]]-dataLimits[0])/
	(dataLimits[1]-dataLimits[0]);
      vertices[j].y = axesPosition[1]*HEIGHT+
	axesPosition[3]*HEIGHT*(1-(yc[cells[3*i+j]]-dataLimits[2])/
	(dataLimits[3]-dataLimits[2]));
    }

    vertices[3].x = vertices[0].x;
    vertices[3].y = vertices[0].y;

    ind = (salt[i]-caxis[0])/(caxis[1]-caxis[0])*64;
    if(salt[i]==0)
      ind = 0;
    /*
    color.red = cmap[ind][0] * 0xffff;
    color.green = cmap[ind][1] * 0xffff;
    color.blue = cmap[ind][2] * 0xffff;
    XAllocColor(dis, colormap, &color);
    XSetForeground(dis,gc,color.pixel);
    */
    Fill(vertices,3,ind,edgelines);
    //XFillPolygon(dis,win,gc,vertices,3,Convex,CoordModeOrigin);
    //XDrawLines(dis,win,gc,vertices,4,0);
  }
  
  sprintf(str1,"Time step: %d, K-level: %d",n,k); 
  XSetForeground(dis,gc,white);
  XDrawString(dis,win,gc,WIDTH/4,HEIGHT-10,str1,strlen(str1));

  XFlush(dis);
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
  float height, width, dx, dy;
  height = axes[3];
  width = axes[2];
  dx = data[1]-data[0];
  dy = data[3]-data[2];

  if(dx>dy) {
    height=width*dy/dx;
    axes[1] = axes[1]+0.5*(axes[3]-height);
    axes[3] = height;
  } else {
    width=height*dx/dy;
    axes[0] = axes[0]+0.5*(axes[2]-width);
    axes[2] = width;
  }
}

void Fill(XPoint *vertices, int N,int cindex, int edges) {
    color.red = cmap[cindex][0] * 0xffff;
    color.green = cmap[cindex][1] * 0xffff;
    color.blue = cmap[cindex][2] * 0xffff;
    XAllocColor(dis, colormap,&color);
    XSetForeground(dis,gc,color.pixel);
    XFillPolygon(dis,win,gc,vertices,3,Convex,CoordModeOrigin);

    if(edges) {
      color.red = 0;
      color.green = 0;
      color.blue = 0;
      XAllocColor(dis, colormap,&color);
      XSetForeground(dis,gc,color.pixel);
      XDrawLines(dis,win,gc,vertices,4,0);
    }
}

