/*
 * File: phys.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 10/21/02
 * --------------------------------
 * This file contains physically-based functions.
 *
 * $Id: phys.c,v 1.17 2003-06-10 02:30:23 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.16  2003/05/12 00:06:23  fringer
 * Added ComputeTangentialVelocity(gridT *grid, physT *phys), which
 * is needed for advection of horizontal velocity (used to be computed
 * for outputting data only).
 *
 * Added ComputeTraceBack, which computes the xd,yd, and zd points of
 * the Lagrangian traceback.  Right now it is a simple xd=x-u*dt.
 *
 * Added FindNearestEdge, which finds the nearest edge as a member of
 * the nearestedges array.  This reduces the search time required by
 * a factor of Ne/Nnearestedges, since the entire grid does not have to
 * be searched.  This was not the origingal intention of the nearestedges
 * arrays but turned out to work well!
 *
 * Added Quadratic, which performs a quadratic interpolation in the
 * vertical.  As of now, vertical interpolation is performed with
 * this function to the level of the traceback, and then Kriging
 * interpolation with linear drift is performed with the Kriging
 * function.
 *
 * Added Kriging.  This function computes the value of the horizontal
 * normal to a face given the normal components at the other edges
 * in the nearestedges array.
 *
 * Added InterpolateEdge, which computes the traceback, determines i
 * f it is within the domain, and then interpolates with Kriging.
 *
 * FindNearestPlane was added but was not used since the interpolation is
 * not done using the cells around the traceback point but from the
 * departure point.
 *
 * To do: precompute the kriging matrix and its inverse, and more
 * accurately compute the vertical interpolation.  Quadratic is very
 * inaccurate because it does not take into account different vertical
 * grid spacings at the fs and bottom.
 *
 * Revision 1.15  2003/05/05 01:27:01  fringer
 * Added AdvectHorizontalVelocity which only sets phys->utmp[i]=0.
 * Changed SendRecvCellData2D to ISendRecvCellData2D (non-blocking send/recv)
 *
 * Revision 1.14  2003/05/01 00:29:05  fringer
 * Changed SendRecvCellData2D/3D so that only the first cells are
 * transferred in the CG solver while all of them are transferred only
 * once when data needs to be output.
 *
 * Revision 1.13  2003/04/29 16:39:29  fringer
 * Added MPI_FOPen in place of fopen.
 *
 * Revision 1.12  2003/04/29 00:12:10  fringer
 * Added ReturnSalinity function and removed salinity initialization from this file.
 * Added BGSalinityFID file pointer to store background salinity in a file.
 * Fixed baroclinic term so that horizontal density gradients do not
 * induce a flow.
 * Removed no-grad boundary condition from ScalarsImp function becuase it
 * was causing unphysical flow at the boundaries.  For now the salinity
 * values in cells neighboring the boundary are not updated.
 *
 * Revision 1.11  2003/04/26 14:16:22  fringer
 * Added initialization function ReturnFreeSurface
 *
 * Revision 1.10  2003/04/22 02:42:32  fringer
 * Changed makepointers() in grid.c so that the ghost points include
 * extra points for the ELM interpolation.  Only the number of neihbors is
 * specified.  This may not guarantee that there are at least 3 neighbors for
 * each cell.
 *
 * Revision 1.9  2003/04/21 20:26:08  fringer
 * Working version before addition of ghost cells for kriging.
 *
 * Revision 1.8  2003/04/08 23:32:46  fringer
 * Added binary output functionality so far for fs,u,v,w.  Use the -a option to enforce ASCII output.
 *
 * Revision 1.7  2002/12/01 10:39:39  fringer
 * Added initial density vector s0 as well as Flux terms to compute internal
 * wave energy fluxes at specified faces.
 *
 * Revision 1.6  2002/11/30 13:44:26  fringer
 * Working version for the simplified Shelf bathymetry.
 *
 * Revision 1.5  2002/11/05 01:31:17  fringer
 * Added baroclinic term
 *
 * Revision 1.4  2002/11/03 20:36:33  fringer
 * Added parameters to check for mass and volume conservation
 *
 * Revision 1.3  2002/11/03 20:25:08  fringer
 * Cleaned up the udpatedz routine.
 *
 * Revision 1.2  2002/11/03 02:09:56  fringer
 * Moved up to field-scale and stable!!
 *
 * Revision 1.1  2002/11/03 00:18:49  fringer
 * Initial revision
 *
 *
 */
#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "util.h"
#include "initialization.h"
#include "memory.h"

/*
 * Private Function declarations.
 *
 */
static void UpdateDZ(gridT *grid, physT *phys, int option);
static void BarotropicPredictor(gridT *grid, physT *phys, 
				propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void CGSolve(gridT *grid, physT *phys, propT *prop, 
		    int myproc, int numprocs, MPI_Comm comm);
static void UpdateU(gridT *grid, physT *phys, propT *prop);
static void HydroW(gridT *grid, physT *phys);
static void WtoVerticalFace(gridT *grid, physT *phys);
static void UpdateScalars(gridT *grid, physT *phys, propT *prop);
static void UpdateScalarsImp(gridT *grid, physT *phys, propT *prop);
static void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs,
			  MPI_Comm comm);
static void Check(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void OutputData(gridT *grid, physT *phys, propT *prop,
		int myproc, int numprocs, MPI_Comm comm);
static void ReadProperties(propT **prop, int myproc);
static void OpenFiles(propT *prop, int myproc);
static void Progress(propT *prop, int myproc);
static void EddyViscosity(gridT *grid, physT *phys, propT *prop);
static void ReadProperties(propT **prop, int myproc);
static void OpenFiles(propT *prop, int myproc);
static void AdvectHorizontalVelocity(gridT *grid, physT *phys, propT *prop,
				     int myproc, int numprocs);
static void ComputeTangentialVelocity(REAL **u, REAL **ut, gridT *grid);
static void ComputeTraceBack(REAL x0, REAL y0, REAL z0, 
			     REAL *xd, REAL *yd, REAL *zd, 
			     REAL un, REAL ut, 
			     REAL W, REAL n1, REAL n2, REAL dt);
static int FindNearestEdge(int j0, REAL xd, REAL yd, REAL zd, gridT *grid, physT *phys);
static void InterpolateEdge(REAL *ui, REAL *vi, REAL *wi,
			    REAL xd, REAL yd, REAL zd,  
			    int j0, int k0, REAL z0,  
			    int jnear, int knear, 
			    gridT *grid, physT *phys, propT *prop, REAL Cab);
static REAL Quadratic(REAL f1, REAL f2, REAL f3, REAL r);
static REAL Bilinear(REAL f1, REAL f2, REAL f3, REAL r);

static void FindNearestPlane(REAL zd, REAL *dz, REAL *znear, int *knear, int Nkmax);
static void InitializeKriging(gridT *grid, int pow);

void AllocatePhysicalVariables(gridT *grid, physT **phys)
{
  int flag=0, ci, i, j, k, Nc=grid->Nc, Ne=grid->Ne;

  *phys = (physT *)SunMalloc(sizeof(physT),"AllocatePhysicalVariables");

  (*phys)->u = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->D = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->utmp = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->ut = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn = (REAL ***)SunMalloc(3*sizeof(REAL **),"AllocatePhysicalVariables");
  for(ci=0;ci<3;ci++)
    (*phys)->Cn[ci] = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wf = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");

  for(j=0;j<Ne;j++) {
    if(grid->Nkc[j]<grid->Nke[j]) {
      printf("Error!  Nkc(=%d)<Nke(=%d) at edge %d\n",grid->Nkc[j],grid->Nke[j],j);
      flag = 1;
    }
    (*phys)->u[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->ut[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    for(ci=0;ci<3;ci++)
      (*phys)->Cn[ci][j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wf[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
  }
  if(flag) {
    MPI_Finalize();
    exit(0);
  }

  (*phys)->h = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->hold = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  
  (*phys)->w = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wtmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->q = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s0 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->nu_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->tau_T = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->tau_B = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdT = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdB = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  
  for(i=0;i<Nc;i++) {
    (*phys)->w[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->q[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s0[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->nu_tv[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  (*phys)->ap = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->am = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->bp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->bm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->a = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->b = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->c = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->d = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
}

void FreePhysicalVariables(gridT *grid, physT *phys)
{
  int i, j, k, ci, Nc=grid->Nc, Ne=grid->Ne;
  
  for(j=0;j<Ne;j++) {
    free(phys->u[j]);
    free(phys->utmp[j]);
    free(phys->ut[j]);
    for(ci=0;ci<3;ci++)
      free(phys->Cn[ci][j]);
    free(phys->wf[j]);
  }

  for(i=0;i<Nc;i++) {
    free(phys->w[i]);
    free(phys->wtmp[i]);
    free(phys->q[i]);
    free(phys->s[i]);
    free(phys->s0[i]);
    free(phys->stmp[i]);
    free(phys->nu_tv[i]);
  }

  free(phys->h);
  free(phys->htmp);
  free(phys->w);
  free(phys->wtmp);
  free(phys->wf);
  free(phys->q);
  free(phys->s);
  free(phys->s0);
  free(phys->stmp);
  free(phys->nu_tv);
  free(phys->tau_T);
  free(phys->tau_B);
  free(phys->CdT);
  free(phys->CdB);
  free(phys->u);
  free(phys->D);
  free(phys->utmp);
  free(phys->ut);
  for(ci=0;ci<3;ci++)
    free(phys->Cn[ci]);
  free(phys->Cn);

  free(phys->ap);
  free(phys->am);
  free(phys->bp);
  free(phys->bm);
  free(phys->a);
  free(phys->b);
  free(phys->c);
  free(phys->d);

  free(phys);
}
    
void InitializePhysicalVariables(gridT *grid, physT *phys)
{
  int i, j, jptr, k, n, nc1, nc2, ne, nf, Nc=grid->Nc, Ne=grid->Ne;
  REAL r, u, v, xc, yc, hf, hfmax, z;

  for(i=0;i<Nc;i++) {
    phys->h[i]=ReturnFreeSurface(grid->xv[i],grid->yv[i],grid->dv[i]);
    if(phys->h[i]<-grid->dv[i])
      phys->h[i]=-grid->dv[i];
  }

  UpdateDZ(grid,phys,1);

  for(i=0;i<Nc;i++) {
    phys->w[i][grid->Nk[i]]=0;
    for(k=0;k<grid->Nk[i];k++) {
      phys->w[i][k]=0;
      phys->q[i][k]=0;
      phys->s[i][k]=0;
      phys->s0[i][k]=0;
    }
  }

  for(i=0;i<Nc;i++) {
    z = 0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //z-=grid->dzz[i][k]/2;
      z-=grid->dz[k]/2;
      phys->s[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z);
      phys->s0[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z);
      //      z-=grid->dzz[i][k]/2;
      z-=grid->dz[k]/2;
    }
  }

  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nkc[j];k++) 
      phys->u[j][k]=0;
    for(k=0;k<grid->Nke[j];k++)
      phys->wf[j][k]=0;
  }

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    xc=0.5*(grid->xv[nc1]+grid->xv[nc2]);
    yc=0.5*(grid->yv[nc1]+grid->yv[nc2]);
    u = -(yc-7.5)/7.5;
    v = (xc-7.5)/7.5;
    z = 0;
    for(k=0;k<grid->Nkc[j];k++) {
      z-=0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
      phys->u[j][k]=ReturnHorizontalVelocity(xc,yc,grid->n1[j],grid->n2[j],z);
      z-=0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
    }
  }

  /*
   * Determine minimum and maximum scalar values.
   */
  phys->smin=phys->s[0][0];
  phys->smax=phys->s[0][0];
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) {
      if(phys->s[i][k]<phys->smin) phys->smin=phys->s[i][k];      
      if(phys->s[i][k]>phys->smax) phys->smax=phys->s[i][k];      
    }
}

void InitializeVerticalGrid(gridT **grid)
{
  int i, k, Nc=(*grid)->Nc;

  (*grid)->dzz = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->dzzold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");

  for(i=0;i<Nc;i++) {
    (*grid)->dzz[i]=(REAL *)SunMalloc(((*grid)->Nk[i])*sizeof(REAL),"InitializeVerticalGrid");
    (*grid)->dzzold[i]=(REAL *)SunMalloc(((*grid)->Nk[i])*sizeof(REAL),"InitializeVerticalGrid");
    for(k=0;k<(*grid)->Nk[i];k++) {
      (*grid)->dzz[i][k]=(*grid)->dz[k];  
      (*grid)->dzzold[i][k]=(*grid)->dz[k];  
    }
  }
}

static void UpdateDZ(gridT *grid, physT *phys, int option)
{
  int i, j, k, ne1, ne2, Nc=grid->Nc, Ne=grid->Ne, flag;
  REAL z;

  if(!option) {
    for(j=0;j<Ne;j++)
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->Nk[i];k++)
	grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }

  for(i=0;i<Nc;i++) {
    z = 0;
    for(k=0;k<grid->Nk[i];k++)
      z-=grid->dz[k];
    grid->dzz[i][grid->Nk[i]-1]=grid->dz[grid->Nk[i]-1]+grid->dv[i]+z;
  }
  
  //  for(k=grid->ctop[8];k<grid->Nk[8];k++)
  //    printf("dzz[8][%d]=%f\n",k,grid->dzz[8][k]);

  if(grid->Nkmax>1) {
    for(i=0;i<Nc;i++) {
      z = 0;
      flag = 0;

      for(k=0;k<grid->Nk[i];k++) {
	z-=grid->dz[k];
	if(phys->h[i]>=z) 
	  if(!flag) {
	    if(k==grid->Nk[i]-1) {
	      grid->dzz[i][k]=phys->h[i]+grid->dv[i];
	      grid->ctop[i]=k;
	    } else if(phys->h[i]==z) {
	      grid->dzz[i][k]=0;
	      grid->ctop[i]=k+1;
	    } else {
	      grid->dzz[i][k]=phys->h[i]-z;
	      grid->ctop[i]=k;
	    }
	    flag=1;
	  } else {
	    if(k==grid->Nk[i]-1) 
	      grid->dzz[i][k]=grid->dz[k]+grid->dv[i]+z;
	    else 
	      if(z<-grid->dv[i])
		grid->dzz[i][k]=0;
	      else 
		grid->dzz[i][k]=grid->dz[k];
	  } 
	else 
	  grid->dzz[i][k]=0;
      }
    }
  } else 
    for(i=0;i<Nc;i++) 
      grid->dzz[i][0]=grid->dv[i]+phys->h[i];


  //  for(k=grid->ctop[8];k<grid->Nk[8];k++)
  //    printf("dzz[8][%d]=%f\n",k,grid->dzz[8][k]);

  /*
  for(i=0;i<Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      printf("dzz[%d][%d]=%f, h+d=%f\n",i,k,grid->dzz[i][k],phys->h[i]+grid->dv[i]);
  */
  /*
  for(i=0;i<Nc;i++) {
    z=0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      z+=grid->dzz[i][k];
    if(grid->Nk[i]*grid->dz[0]<grid->dv[i])
      printf("Nk[%d]=%d,Nkdz=%f,dv=%f,dv+h=%f,sum=%f\n",i,grid->Nk[i],grid->Nk[i]*grid->dz[0],
	     grid->dv[i],grid->dv[i]+phys->h[i],z);	  
  }
  */
  /*
  // Correction to make sure dz doesn't get too small.
  for(i=0;i<Nc;i++) {
    if(grid->dzz[i][grid->ctop[i]]/grid->dz[grid->ctop[i]]<0.1) 
      if(grid->ctop[i]<grid->Nk[i]-1) {
	grid->ctop[i]++;
	grid->dzz[i][grid->ctop[i]]+=grid->dzz[i][grid->ctop[i]-1];
	grid->dzz[i][grid->ctop[i]-1]=0;
      }
  }
  */

  for(j=0;j<grid->Ne;j++) {
    ne1 = grid->grad[2*j];
    ne2 = grid->grad[2*j+1];
    if(ne1 == -1)
      grid->etop[j]=grid->ctop[ne2];
    else if(ne2 == -1)
      grid->etop[j]=grid->ctop[ne1];
    else if(grid->ctop[ne1]<grid->ctop[ne2])
      grid->etop[j]=grid->ctop[ne1];
    else
      grid->etop[j]=grid->ctop[ne2];
  }
  if(option) {
    for(j=0;j<Ne;j++) 
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->Nk[i];k++)
	grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }
}

void Solve(gridT *grid, physT *phys, int myproc, int numprocs, MPI_Comm comm)
{
  int i, j, k, n;
  extern int TotSpace;
  propT *prop;
  //  FILE *fid;// = fopen("../../data/hleft.dat","w");

  ReadProperties(&prop,myproc);
  OpenFiles(prop,myproc);

  prop->n=0;
  ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);

  if(prop->nonlinear) {
    if(VERBOSE>1) printf("Initializing the covariance matrix for Kriging...\n");
    InitializeKriging(grid,prop->kriging_cov);
  }

  if(VERBOSE>1) printf("Processor %d,  Total memory: %d bytes\n",myproc,TotSpace);

  for(n=1;n<=prop->nsteps;n++) {
    prop->n = n;
    prop->rtime = (n-1)*prop->dt;
    if(prop->nsteps>0) {

      EddyViscosity(grid,phys,prop);

      AdvectHorizontalVelocity(grid,phys,prop,myproc,numprocs);
      BarotropicPredictor(grid,phys,prop,myproc,numprocs,comm);
      ISendRecvCellData2D(phys->h,grid,myproc,comm,first);

      HydroW(grid,phys);

      UpdateScalarsImp(grid,phys,prop);
      SendRecvCellData3D(phys->s,grid,myproc,comm,first);

      // We need w on the vertical faces as well as the tangential
      // components of velocity on the faces in order to compute the
      // tracebacks and interpolation for the advection of the horizontal velocity
      // at the faces.
      WtoVerticalFace(grid,phys);
      ComputeTangentialVelocity(phys->u,phys->ut,grid);
    }

    Progress(prop,myproc);
    OutputData(grid,phys,prop,myproc,numprocs,comm);
    Check(grid,phys,prop,myproc,numprocs,comm);
  }
  //  free(prop);
}

/*
 * Function: AdvectHorizontalVelocity
 * Usage: AdvectHorizontalVelocity(grid,phys,prop,myproc,numprocs);
 * ----------------------------------------------------------------
 * This function computes the tracebacks from each computational
 * edge and stores the interpolated values in the array
 * phys->utmp.
 *
 */
static void AdvectHorizontalVelocity(gridT *grid, physT *phys, propT *prop,
				     int myproc, int numprocs) {
  int j, jptr, k, nc1, nc2, jnear, knear, n, numiters=1;
  REAL x0, y0, z0, xd, yd, zd, zd0, h0, dv0, *ustar, *vstar, *wstar, Cab;

  ustar = phys->a;
  vstar = phys->b;
  wstar = phys->c;

  // Set utmp to zero for all Nke.
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nke[j];k++)
      phys->utmp[j][k]=0;
  }

  // Set utmp=u if this is a linear calculation.
  if(!prop->nonlinear) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->utmp[j][k]=phys->u[j][k];
    }
  } else {    
    // Interpolate the u.  Since etop[j] may contain newly-wetted
    // edges this is okay since the velocity field is interpolated
    // from the previous time step.
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      
      if(nc1==-1 || nc2==-1) { 
	printf("Error!  Trying to compute advection on a boundary edge!\n");
	exit(1);
      }
      
      // Assume face center is half-way between voronoi points.
      h0 = 0.5*(phys->h[nc1]+phys->h[nc2]);
      dv0 = 0.5*(grid->dv[nc1]+grid->dv[nc2]);
      x0 = 0.5*(grid->xv[nc1]+grid->xv[nc2]);
      y0 = 0.5*(grid->yv[nc1]+grid->yv[nc2]);
      
      for(k=0;k<grid->etop[j];k++)
	phys->utmp[j][k]=0;
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	ustar[k]=phys->u[j][k]*grid->n1[j]-phys->ut[j][k]*grid->n2[j];
	vstar[k]=phys->u[j][k]*grid->n2[j]+phys->ut[j][k]*grid->n1[j];
	wstar[k]=phys->wf[j][k];
      }

      for(n=0;n<numiters;n++) {
	z0 = h0;
	if(n==numiters-1)
	  Cab=0;
	else
	  Cab=0.5;

	for(k=grid->etop[j];k<grid->Nke[j];k++) {
	  z0-=0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
	  
	  xd=x0-prop->dt*ustar[k];
	  yd=y0-prop->dt*vstar[k];
	  zd=z0-prop->dt*wstar[k];

	  if(FindNearestEdge(j,xd,yd,zd,grid,phys)>=0) {
	    FindNearestPlane(zd,grid->dz,&zd0,&knear,grid->Nkmax);
	    InterpolateEdge(&(ustar[k]),&(vstar[k]),&(wstar[k]),
			    xd,yd,zd,j,k,zd0,j,knear,grid,phys,prop,Cab);
	  } else {
	    ustar[k]=0;
	    vstar[k]=0;
	    wstar[k]=0;
	  } 

	  z0-=0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
	}
      }

      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	phys->utmp[j][k]=ustar[k]*grid->n1[j]+vstar[k]*grid->n2[j];
    }
  }

  // Hold on to the current velocity field to create V* using Adams-Bashforth.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
      
    for(k=0;k<grid->etop[j];k++) {
      phys->Cn[0][j][k]=0;
      phys->Cn[1][j][k]=0;
      phys->Cn[2][j][k]=0;
    }
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->Cn[0][j][k]=phys->u[j][k];
      phys->Cn[1][j][k]=phys->ut[j][k];
      phys->Cn[2][j][k]=phys->wf[j][k];
    }
  }
}

static void InterpolateEdge(REAL *ui, REAL *vi, REAL *wi,
			    REAL xd, REAL yd, REAL zd,  
			    int j0, int k0, REAL z0,  
			    int jnear, int knear, 
			    gridT *grid, physT *phys, propT *prop, REAL Cab) {
  int ci, n, j, k, je, nk, status, nc1, nc2, sign;
  REAL U, V, un1, un2, dz,  *x, *y, **uf, us[3];

  uf = (REAL **)SunMalloc(3*sizeof(REAL *),"InterpolateEdge");
  for(ci=0;ci<3;ci++)
    uf[ci]=(REAL *)SunMalloc(grid->Nnearestedges*sizeof(REAL),"InterpolateEdge");

  x = (REAL *)SunMalloc(grid->Nnearestedges*sizeof(REAL),"InterpolateEdge");
  y = (REAL *)SunMalloc(grid->Nnearestedges*sizeof(REAL),"InterpolateEdge");

  if(zd<z0) {
    nk=1;
    sign=-1;
  } else {
    nk=-1;
    sign=1;
  }

  if(knear==grid->etop[j0])
    if(zd<=z0) 
      nk=1;
    else
      nk=0;
  else if(knear==grid->Nke[j0]-1)
    if(zd>=z0) 
      nk=-1;
    else
      nk=0;

  if(grid->etop[j0]==grid->Nke[j0]-1)
    nk=0;

  for(ci=0;ci<3;ci++) 
    for(n=0;n<grid->Nnearestedges;n++) {
      je = grid->nearestedges[j0][n];
      
      nc1=grid->grad[2*je];
      nc2=grid->grad[2*je+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;
      
      dz = 0.25*(grid->dzz[nc1][knear]+grid->dzz[nc2][knear]+
		 grid->dzz[nc1][knear+nk]+grid->dzz[nc2][knear+nk]);
      
      if(ci==0) {
	un1=(1+Cab)*(phys->u[je][knear]*grid->n1[je]-phys->ut[je][knear]*grid->n2[je])-
	  Cab*(phys->Cn[0][je][knear]*grid->n1[je]-phys->Cn[1][je][knear]*grid->n2[je]);
	un2=(1+Cab)*(phys->u[je][knear+nk]*grid->n1[je]-phys->ut[je][knear+nk]*grid->n2[je])-
	  Cab*(phys->Cn[0][je][knear+nk]*grid->n1[je]-phys->Cn[1][je][knear+nk]*grid->n2[je]);
      } else if(ci==1) {
	un1=(1+Cab)*(phys->u[je][knear]*grid->n2[je]+phys->ut[je][knear]*grid->n1[je])-
	  Cab*(phys->Cn[0][je][knear]*grid->n2[je]+phys->Cn[1][je][knear]*grid->n1[je]);
	un2=(1+Cab)*(phys->u[je][knear+nk]*grid->n2[je]+phys->ut[je][knear+nk]*grid->n1[je])-
	  Cab*(phys->Cn[0][je][knear+nk]*grid->n2[je]+phys->Cn[1][je][knear+nk]*grid->n1[je]);
      } else {
	un1=(1+Cab)*phys->wf[je][knear]-Cab*phys->Cn[2][je][knear];
	un2=(1+Cab)*phys->wf[je][knear+nk]-Cab*phys->Cn[2][je][knear+nk];
      }
      
      // If dz=0 this means that the interpolation is being performed in
      // a dry zone, so set uf[n]=un1 if this is so.
      if(dz==0)
	uf[ci][n] = un1;
      else 
	uf[ci][n] = un1 + sign*(zd-z0)*(un2-un1)/dz;
    }
  
  for(n=0;n<grid->Nnearestedges;n++) {
    je = grid->nearestedges[j0][n];
    x[n] = grid->xe[je];
    y[n] = grid->ye[je];
  }

  kriging(xd,yd,us,x,y,uf,3,prop->kriging_cov,grid->Nnearestedges,grid->Kriging[j0]);

  *ui=us[0];
  *vi=us[1];
  *wi=us[2];

  free(x);
  free(y);
  free(uf[0]);
  free(uf[1]);
  free(uf[2]);
  free(uf);
}
 
static REAL Bilinear(REAL f1, REAL f2, REAL f3, REAL r) {
  return f2 + 0.5*(r*(f1-f2) + r*(f3-f2));
}
   
static REAL Quadratic(REAL f1, REAL f2, REAL f3, REAL r) {
  return f1 + 0.5*r*(f3-f2) + 0.5*pow(r,2)*(f3-2*f2+f1);
}

static void ComputeTraceBack(REAL x0, REAL y0, REAL z0, 
			     REAL *xd, REAL *yd, REAL *zd, 
			     REAL un, REAL ut, 
			     REAL W, REAL n1, REAL n2, REAL dt) {

  REAL U, V;

  // U is the x component of velocity at the face
  // V is the y component of velocity at the face.
  // U^2 + V^2 = un^2 + ut^2
  U = un*n1-ut*n2;
  V = un*n2+ut*n1;

  *xd=x0-dt*U;
  *yd=y0-dt*V;
  *zd=z0-dt*W;
}
  
static int FindNearestEdge(int j0, REAL xd, REAL yd, REAL zd, gridT *grid, physT *phys) {

  int j, je, jnearest, departurecell;
  REAL dot, dist, dist0;

  dist0 = INFTY;

  for(j=0;j<grid->Nnearestedges;j++) {
    je = grid->nearestedges[j0][j];
    if(je==-1) { 
      printf("Error! jnearest[%d][%d]=-1\n",j0,j,je);
      exit(0);
    }
    dist = pow(grid->xe[je]-xd,2)+pow(grid->ye[je]-yd,2);
    if(dist<dist0) {
      jnearest=je;
      dist0=dist;
    }
  }

  // This is the dot product of the normal vector with the position vector
  // of the departure point relative to the center of the edge.  Since
  // the normal vector always points towards grad[2*jnearest] then if the dot product
  // is positive then this is its cell.
  dot = 
    (xd-grid->xe[jnearest])*grid->n1[jnearest]+
    (yd-grid->ye[jnearest])*grid->n2[jnearest];

  if(dot>0)
    departurecell=grid->grad[2*jnearest];
  else
    departurecell=grid->grad[2*jnearest+1];
  
  // If the grad pointer is a ghost cell then the departure point is outside
  // the boundary at the top-most cell.
  if(departurecell==-1) 
    return -1;
  
  // If it is within the boundary at the top-most cell, now make sure it
  // is within the upper and lower bounds of the domain.
  if(zd>=phys->h[departurecell] || zd<=-grid->dv[departurecell]) 
    return -1;

  // Otherwise it is a valid edge so return it!
  return jnearest;
}

static void FindNearestPlane(REAL zd, REAL *dz, REAL *znear, int *knear, int Nkmax) {
  int k;
  REAL ztop=0, zbot=-dz[0];
  
  *knear=0;
  *znear=0.5*(zbot+ztop);
  for(k=0;k<Nkmax;k++) {
    if(zd<ztop && zd>zbot) {
      *knear=k;
      *znear=0.5*(zbot+ztop);
      break;
    }
    if(k<Nkmax-1) {
      ztop=zbot;
      zbot-=dz[k+1];
    }
  }
}

static void EddyViscosity(gridT *grid, physT *phys, propT *prop)
{
  int i, j, jptr, k;
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    phys->tau_T[j]=grid->n1[j]*prop->tau_T;
    phys->tau_B[j]=0;
    phys->CdT[j]=prop->CdT;
    phys->CdB[j]=prop->CdB;
  }
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      phys->nu_tv[i][k]=0;
  }
}

static void UpdateU(gridT *grid, physT *phys, propT *prop)
{
  int j, k, jptr;
  REAL dt = prop->dt;
  REAL omega = prop->omega;
  REAL rtime = prop->rtime;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    for(k=0;k<grid->etop[j];k++)
      phys->u[j][k]=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->u[j][k]-=GRAV*dt*(phys->h[grid->grad[2*j]]-phys->h[grid->grad[2*j+1]])/grid->dg[j];
    }
  }
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    for(k=0;k<grid->Nke[j];k++) 
      phys->u[j][k]=prop->flux*cos(prop->omega*prop->rtime);

  }

  HydroW(grid,phys);
}    

/*
 * Function: BarotropicPredictor 
 * Usage: BarotropicPredictor(grid,phys,prop,myproc,numprocs,comm);
 * ----------------------------------------------------------------
 * Update the free surface.
 *
 */
static void BarotropicPredictor(gridT *grid, physT *phys, 
				propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, j, jptr, ne, nf, normal, nc1, nc2, k, upwind, k0;
  REAL sum, dt=prop->dt, theta=prop->theta, thetaFS=prop->thetaFS, fluxheight;
  REAL *a, *b, *c, *d, *e1, **E;

  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;
  e1 = phys->ap;
  E = phys->ut;

  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=0;
    //    for(k=grid->etop[j];k<grid->Nke[j];k++)
    //      phys->utmp[j][k]=0;
  }

  // First create U**.  This is where the semi-Lagrangian formulation is added
  // This is commented out since the ELM scheme is interpolated in 
  // AdvectHorizontalVelocity.
  /*
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    for(k=0;k<grid->etop[j];k++)
      phys->utmp[j][k]=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]=phys->u[j][k];
  }
  */

  // phys->utmp contains the interpolated velocity field which is computed
  // in AdvectHorizontalVelocity.

  // Add on the baroclinic term
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      for(k0=grid->etop[j];k0<k;k0++)
	phys->utmp[j][k]-=0.5*GRAV*prop->beta*dt*(phys->s[nc1][k0]-phys->s[nc2][k0])*
	  (grid->dzz[nc1][k0]+grid->dzz[nc2][k0])/grid->dg[j];
    }
    /*
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->utmp[j][k]-=0.5*dt*GRAV*prop->beta*
	(phys->s[nc1][k]*grid->dzz[nc1][k]-
	 phys->s[nc2][k]*grid->dzz[nc2][k])/grid->dg[j];
      for(k0=grid->etop[j];k0<k;k0++)
	phys->utmp[j][k]-=dt*GRAV*prop->beta*(phys->s[nc1][k0]*grid->dzz[nc1][k0]-
					      phys->s[nc2][k0]*grid->dzz[nc2][k0])/grid->dg[j];
    }
    */
  }	

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // Add the explicit part of the free-surface to U**.
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]-=GRAV*(1-thetaFS)*dt*(phys->h[nc1]-phys->h[nc2])/grid->dg[j];

    // Add the shear stress from the top cell
    phys->utmp[j][grid->etop[j]]+=dt*phys->tau_T[j];

    /*
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      printf("utmp[%d][%d]=%f\n",j,k,phys->utmp[j][k]);
    */

    // Create the tridiagonal entries and formulate U***
    if(!(grid->dzz[nc1][grid->etop[j]]==0 && grid->dzz[nc2][grid->etop[j]]==0)) {

      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	a[k]=0;
	b[k]=0;
	c[k]=0;
	d[k]=0;
      }

      // Vertical eddy-viscosity interpolated to faces since it is stored
      // at cell-centers.
      for(k=grid->etop[j]+1;k<grid->Nke[j];k++) 
	c[k]=0.25*(phys->nu_tv[nc1][k-1]+phys->nu_tv[nc2][k-1]+
		   phys->nu_tv[nc1][k]+phys->nu_tv[nc2][k]);
      
      // Coefficients for the viscous terms.  Face heights are taken as
      // the average of the face heights on either side of the face (not upwinded).
      for(k=grid->etop[j]+1;k<grid->Nke[j];k++) 
	a[k]=2.0*(prop->nu+c[k])/(0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k])*
				  (grid->dzz[nc1][k-1]+grid->dzz[nc2][k-1]+
				   grid->dzz[nc1][k]+grid->dzz[nc2][k]));
      
      for(k=grid->etop[j];k<grid->Nke[j]-1;k++) {
	b[k]=2.0*(prop->nu+c[k+1])/(0.25*(grid->dzz[nc1][k]+grid->dzz[nc2][k])*
				    (grid->dzz[nc1][k]+grid->dzz[nc2][k]+
				     grid->dzz[nc1][k+1]+grid->dzz[nc2][k+1]));
      }
      /*
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	printf("a[%d]=%f,b[%d]=%f\n",k,a[k],k,b[k]);
      */

      if(grid->Nke[j]-grid->etop[j]>1) {

	// Explicit part of the viscous term
	for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++)
	  phys->utmp[j][k]+=dt*(1-theta)*(a[k]*phys->u[j][k-1]-
					  (a[k]+b[k])*phys->u[j][k]+
					  b[k]*phys->u[j][k+1]);
	
	// Top cell
	phys->utmp[j][grid->etop[j]]+=dt*(1-theta)*(-(b[grid->etop[j]]+2.0*phys->CdT[j]*
						      fabs(phys->u[j][grid->etop[j]])/
						      (grid->dzz[nc1][grid->etop[j]]+
						       grid->dzz[nc2][grid->etop[j]]))*
						    phys->u[j][grid->etop[j]]
						    +b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);
	
	// Bottom cell
	phys->utmp[j][grid->Nke[j]-1]+=dt*(1-theta)*(a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2]-
						     (a[grid->Nke[j]-1]+2.0*phys->CdB[j]*
						      fabs(phys->u[j][grid->Nke[j]-1])/
						      (grid->dzz[nc1][grid->Nke[j]-1]+
						       grid->dzz[nc2][grid->Nke[j]-1]))*
						     phys->u[j][grid->Nke[j]-1]);
      
      } else 
	phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdB[j]+phys->CdT[j])/
	  (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
	  fabs(phys->u[j][grid->etop[j]])*phys->u[j][grid->etop[j]];

      // Now set up the coefficients for the tridiagonal inversion for the
      // implicit part.  These are given from the arrays above in the discrete operator
      // d^2U/dz^2 = -theta dt a_k U_{k-1} + (1+theta dt (a_k+b_k)) U_k - theta dt b_k U_{k+1}

      // Right hand side U** is given by d[k] here.
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	e1[k]=1.0;
	d[k]=phys->utmp[j][k];
      }
      
      if(grid->Nke[j]-grid->etop[j]>1) {
	// Top cells
	c[grid->etop[j]]=-theta*dt*b[grid->etop[j]];
	b[grid->etop[j]]=1.0+theta*dt*(b[grid->etop[j]]+
				       2.0*phys->CdT[j]*fabs(phys->u[j][grid->etop[j]])/
				       (grid->dzz[nc1][grid->etop[j]]+
					grid->dzz[nc2][grid->etop[j]]));
	a[grid->etop[j]]=0;
	
	// Bottom cell
	c[grid->Nke[j]-1]=0;
	b[grid->Nke[j]-1]=1.0+theta*dt*(a[grid->Nke[j]-1]+
					2.0*phys->CdB[j]*fabs(phys->u[j][grid->Nke[j]-1])/
					(grid->dzz[nc1][grid->Nke[j]-1]+
					 grid->dzz[nc2][grid->Nke[j]-1]));
	a[grid->Nke[j]-1]=-theta*dt*a[grid->Nke[j]-1];
      
	// Interior cells
	for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
	  c[k]=-theta*dt*b[k];
	  b[k]=1.0+theta*dt*(a[k]+b[k]);
	  a[k]=-theta*dt*a[k];
	}
      } else {
	b[grid->etop[j]]=1.0+2.0*theta*dt*fabs(phys->u[j][grid->etop[j]])/
	  (grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]])*
	  (phys->CdB[j]+phys->CdT[j]);
      }	  

      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	if(grid->dzz[nc1][k]==0 && grid->dzz[nc2][k]==0) {
	  printf("Exiting because dzz[%d][%d]=%f or dzz[%d][%d]=%f\n",
		 nc1,k,grid->dzz[nc1][k],nc2,k,grid->dzz[nc2][k]);
	  exit(0);
	}
	if(a[k]!=a[k]) printf("a[%d] problems, dzz[%d][%d]=%f\n",k,j,k,grid->dzz[j][k]);
	if(b[k]!=b[k] || b[k]==0) printf("b[%d] problems\n",k);
	if(c[k]!=c[k]) printf("c[%d] problems\n",k);
      }

      /*
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	if(b[k]!=1 || c[k]!=0 || a[k]!=0 || d[k]!=0) 
	  printf("j=%d, a[%d]=%f,b[%d]=%f,c[%d]=%f,d[%d]=%f\n",j,k,a[k],k,b[k],k,c[k],k,d[k]);
      */

      // Now utmp will have U*** in it, which is given by A^{-1}U**, and E will have
      // A^{-1}e1, where e1 = [1,1,1,1,1,...,1]^T 
      if(grid->Nke[j]-grid->etop[j]>1) {
	//	for(k=grid->etop[j];k<grid->Nke[j];k++) 
	//	  printf("Before: e1[%d]=%f E[%d][%d]=%f\n",k,e1[k],j,k,E[j][k]);
	TriSolve(&(a[grid->etop[j]]),&(b[grid->etop[j]]),&(c[grid->etop[j]]),
		 &(d[grid->etop[j]]),&(phys->utmp[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);
	TriSolve(&(a[grid->etop[j]]),&(b[grid->etop[j]]),&(c[grid->etop[j]]),
		 &(e1[grid->etop[j]]),&(E[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);	
	//	for(k=grid->etop[j];k<grid->Nke[j];k++) 
	//	  printf("After: E[%d][%d]=%f\n",j,k,E[j][k]);
      } else {
	phys->utmp[j][grid->etop[j]]/=b[grid->etop[j]];
	E[j][grid->etop[j]]=1.0/b[grid->etop[j]];
      }

      // Now vertically integrate E to create the vertically integrated flux-face
      // values that comprise the coefficients of the free-surface solver.  This
      // will create the D vector, where D=DZ^T E (which should be given by the
      // depth when there is no viscosity.
      phys->D[j]=0;
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	fluxheight=UpWind(phys->u[j][k],grid->dzz[nc1][k],grid->dzz[nc2][k]);

	phys->D[j]+=fluxheight*E[j][k];
      }
    }
  }

  for(j=0;j<grid->Ne;j++) 
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      if(phys->utmp[j][k]!=phys->utmp[j][k]) {
	printf("Error at %d %d (U***=nan)\n",j,k);
	exit(1);
      }

  // So far we have U*** and D.  Now we need to create h* in htmp.   This
  // will comprise the source term for the free-surface solver.
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];


    sum = 0;
    for(nf=0;nf<NFACES;nf++) {
      
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      for(k=grid->etop[ne];k<grid->Nke[ne];k++) {
	fluxheight=UpWind(phys->u[ne][k],grid->dzz[nc1][k],grid->dzz[nc2][k]);

	sum+=((1-thetaFS)*phys->u[ne][k]+thetaFS*phys->utmp[ne][k])*
	  fluxheight*grid->df[ne]*normal;
      }
    }
    phys->htmp[i]=phys->h[i]-dt/grid->Ac[i]*sum;
  }

  // Set the free-surface values at the boundary.  These are set for h and not for
  // htmp since the boundary values for h will be used in the cg solver.
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    phys->h[i] = prop->amp*sin(prop->omega*prop->rtime);
  }

  // Now we have the required components for the CG solver for the free-surface:
  //
  // h^{n+1} - g*(theta*dt)^2/Ac * Sum_{faces} D_{face} dh^{n+1}/dn df N = htmp
  //
  // L(h) = b
  //
  // L(h) = h + 1/Ac * Sum_{faces} D_{face} dh^{n+1}/dn N
  // b = htmp
  //
  // As the initial guess let h^{n+1} = h^n, so just leave it as it is to
  // begin the solver.
  CGSolve(grid,phys,prop,myproc,numprocs,comm);

  // Add on the implicit barotropic term to obtain the hydrostatic horizontal velocity field.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=grid->etop[j];k<grid->Nke[j];k++)
      phys->u[j][k]=phys->utmp[j][k]-GRAV*thetaFS*dt*E[j][k]*
	(phys->h[nc1]-phys->h[nc2])/grid->dg[j];

    if(grid->etop[j]==grid->Nke[j]-1 && grid->dzz[nc1][grid->etop[j]]==0 &&
       grid->dzz[nc2][grid->etop[j]]==0) 
      phys->u[j][grid->etop[j]]=0;
  }

  UpdateDZ(grid,phys,0);

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(grid->etop[j]>grid->etopold[j]) 
      for(k=0;k<grid->etop[j];k++)
	phys->u[j][k]=0;
    else 
      for(k=grid->etop[j];k<grid->etopold[j];k++)
	phys->u[j][k]=phys->utmp[j][k]-GRAV*thetaFS*dt*
	  (phys->h[nc1]-phys->h[nc2])/grid->dg[j];
  }

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      if(IsNan(phys->u[j][k]))
	printf("Nk[%d]=%d,Nk[%d]=%d,h[%d]=%f, h[%d]=%f, u[%d][%d]=NaN, etop: %d, etopold: %d, Nke: %d\n",
	       grid->grad[2*j],grid->Nk[grid->grad[2*j]],
	       grid->grad[2*j+1],grid->Nk[grid->grad[2*j+1]],
	       grid->grad[2*j],phys->h[grid->grad[2*j]],
	       grid->grad[2*j+1],phys->h[grid->grad[2*j+1]],
	       j,k,grid->etop[j],grid->etopold[j],grid->Nke[j]);
  }
						
  // Set the flux values at boundary cells if specified

  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k]=prop->flux*cos(prop->omega*prop->rtime);
  }
}

static void CGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, j, iptr, nf, ne, nc1, nc2, n, niters, *N, Nc=grid->Nc, *counts, count;
  REAL *h, *hold, *D, *hsrc, myresid, resid, residold, val, tmp, **M, relax, myNsrc, Nsrc, coef;
  //  FILE *fid=fopen("../../data/cg.m","w");

  h = phys->h;
  hold = phys->hold;
  D = phys->D;
  hsrc = phys->htmp;
  N = grid->normal;

  tmp = GRAV*pow(prop->thetaFS*prop->dt,2);

  /*
  M = (REAL **)SunMalloc(Nc*sizeof(REAL *),"CGSolve");
  for(i=0;i<Nc;i++) {
    M[i] = (REAL *)SunMalloc(Nc*sizeof(REAL),"CGSolve");
    for(j=0;j<Nc;j++)
      M[i][j]=0;
  }
  */
  /*
  for(i=0;i<Nc;i++)
    printf("%d: Edges: %d %d %d, Neighs: %d %d %d\n",
	   i,grid->face[3*i],grid->face[3*i+1],grid->face[3*i+2],
	   grid->neigh[3*i],grid->neigh[3*i+1],grid->neigh[3*i+2]);
  */
  /*
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    M[i][i]=1;
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      if(grid->neigh[i*NFACES+nf]!=-1) {
	M[grid->neigh[i*NFACES+nf]][i]=-tmp*D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	M[i][i]+=M[i][grid->neigh[i*NFACES+nf]];
      }
    }
  }
  */

  ISendRecvCellData2D(h,grid,myproc,comm,first);

  relax = prop->relax;
  niters = prop->maxiters;
  resid=0;
  myresid=0;

  myNsrc=0;
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    myNsrc+=pow(hsrc[i],2);
  }
  MPI_Reduce(&myNsrc,&(Nsrc),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&Nsrc,1,MPI_DOUBLE,0,comm);
  Nsrc=sqrt(Nsrc);

  /*
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    printf("%f %f\n",grid->xv[i],hsrc[i]);
  }
  */

  for(n=0;n<niters;n++) {

    for(i=0;i<grid->Nc;i++) {
      hold[i] = h[i];
    }

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      h[i] = hsrc[i];

      coef=1;
      for(nf=0;nf<NFACES;nf++) 
	if(grid->neigh[i*NFACES+nf]!=-1) {
	  ne = grid->face[i*NFACES+nf];

	  coef+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	  h[i]+=relax*tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
	    phys->h[grid->neigh[i*NFACES+nf]]/grid->Ac[i];
	}
      //      if(coef<=0) { printf("COEF Problems\n"); exit(0); }
      if(coef!=coef && numprocs==1) { printf("COEF is NAN\n"); exit(0); }
      h[i]/=coef;
    }

    residold=resid;
    myresid=0;
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      hold[i] = hsrc[i];

      coef=1;
      for(nf=0;nf<NFACES;nf++) 
	if(grid->neigh[i*NFACES+nf]!=-1) {
	  ne = grid->face[i*NFACES+nf];
	  coef+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	  hold[i]+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
	    phys->h[grid->neigh[i*NFACES+nf]]/grid->Ac[i];
	}
      myresid+=pow(hold[i]/coef-h[i],2);
    }
    MPI_Reduce(&myresid,&(resid),1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Bcast(&resid,1,MPI_DOUBLE,0,comm);
    resid=sqrt(resid);

    //    printf("Proc %d: %e\n",myproc,resid);

    ISendRecvCellData2D(h,grid,myproc,comm,first);
    MPI_Barrier(comm);

    if(fabs(resid)<prop->epsilon)
      break;
  }
  if(n==niters && myproc==0) 
    printf("Warning... Iteration not converging after %d steps! RES=%e\n",n,resid);
  
  /*
  fprintf(fid,"M=[");
  for(i=0;i<Nc;i++) {
    for(j=0;j<Nc;j++)
      fprintf(fid,"%f ",M[i][j]);
    fprintf(fid,"\n");
  }
  fprintf(fid,"];\n");
  fprintf(fid,"b=[");
  for(i=0;i<Nc;i++)
    fprintf(fid,"%f ",hsrc[i]);
  fprintf(fid,"]';\n");
  fprintf(fid,"x0=[");
  for(i=0;i<Nc;i++)
    fprintf(fid,"%f ",h[i]);
  fprintf(fid,"]';\n");
  fprintf(fid,"xv=[");
  for(i=0;i<Nc;i++)
    fprintf(fid,"%f ",grid->xv[i]);
  fprintf(fid,"]';");
  fclose(fid);

  for(i=0;i<Nc;i++) 
    free(M[i]);
  free(M);
  */

  for(i=0;i<grid->Nc;i++)
    if(h[i]!=h[i]) 
      printf("NaN h[%d] in cgsolve!\n",i);

}


static void UpdateScalarsImp(gridT *grid, physT *phys, propT *prop)
{
  int i, j, jptr, iptr, k, nf, kp, km, kstart, kend, ktop;
  int Nc=grid->Nc, Ne=grid->Ne, normal, nc1, nc2, ne;
  REAL wup, wum, wlp, wlm, df, Ac, dt=prop->dt, tmp;
  REAL *a, *b, *c, *d, *ap, *am, dznew, dzold, mass, rat;

  ap = phys->ap;
  am = phys->am;
  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;

  for(i=0;i<Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->s[i][k];

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    Ac = grid->Ac[i];
    
    if(grid->ctop[i]>=grid->ctopold[i]) {
      ktop=grid->ctop[i];
      dznew=grid->dzz[i][ktop];
    } else {
      ktop=grid->ctopold[i];
      dznew=0;
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++)
	dznew+=grid->dzz[i][k];      
    }

    for(k=0;k<grid->Nk[i]+1;k++) {
      ap[k] = 0.5*(phys->w[i][k]+fabs(phys->w[i][k]));
      am[k] = 0.5*(phys->w[i][k]-fabs(phys->w[i][k]));
    }
    for(k=ktop+1;k<grid->Nk[i];k++) {
      a[k-ktop]=dt*am[k];
      b[k-ktop]=grid->dzz[i][k]+dt*(ap[k]-am[k+1]);
      c[k-ktop]=-dt*ap[k+1];
    }

    a[0]=0;
    b[0]=dznew-dt*am[ktop+1];
    c[0]=-dt*ap[ktop+1];
    b[(grid->Nk[i]-1)-ktop]+=c[(grid->Nk[i]-1)-ktop];

    for(k=ktop+1;k<grid->Nk[i];k++) 
      d[k-ktop]=grid->dzzold[i][k]*phys->stmp[i][k];

    d[0]=0;
    if(grid->ctopold[i]<=grid->ctop[i])
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
	d[0]+=grid->dzzold[i][k]*phys->stmp[i][k];
    else
      d[0]=grid->dzzold[i][ktop]*phys->stmp[i][ktop];

    for(nf=0;nf<NFACES;nf++) {
    
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      for(k=0;k<grid->Nkc[ne];k++) 
	ap[k] = 0.5*(phys->u[ne][k]+fabs(phys->u[ne][k]))*phys->stmp[nc2][k]*grid->dzzold[nc2][k]
	  +0.5*(phys->u[ne][k]-fabs(phys->u[ne][k]))*phys->stmp[nc1][k]*grid->dzzold[nc1][k];

      for(k=ktop+1;k<grid->Nk[i];k++) 
      	d[k-ktop]-=dt*ap[k]*df*normal/Ac;

      for(k=0;k<=ktop;k++)
	  d[0]-=dt*ap[k]*df*normal/Ac;
    }

    if(grid->Nk[i]-ktop>1) 
      TriSolve(a,b,c,d,&(phys->s[i][ktop]),grid->Nk[i]-ktop);
    else if(b[0]!=0)
      phys->s[i][ktop]=d[0]/b[0];

    if(dznew!=0) {
      mass = phys->s[i][ktop]*dznew;
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++)
	phys->s[i][k]=mass/dznew;
    }

    for(k=0;k<grid->ctop[i];k++)
      phys->s[i][k]=0;
  }

  /*
   * Need to modify the boundary cells here.  For now they
   * are left unmodified.
   *
   */
  /*
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      tmp=phys->s[i][k];
      phys->s[i][k]=0;
      km=0;
      for(nf=0;nf<NFACES;nf++) {
	if(grid->neigh[i*NFACES+nf]!=-1) 
	  if(k<=grid->Nk[grid->neigh[i*NFACES+nf]] &&
	     k>=grid->ctop[grid->neigh[i*NFACES+nf]]) {
	    phys->s[i][k]+=phys->s[grid->neigh[i*NFACES+nf]][k];
	    km++;
	  }
      }
      if(km>0)
	phys->s[i][k]/=km;
      else
	phys->s[i][k]=tmp;
    }
   }
  */
}

static void UpdateScalars(gridT *grid, physT *phys, propT *prop)
{
  int i, j, jptr, iptr, k, nf, kp, km, kstart, kend;
  int Nc=grid->Nc, Ne=grid->Ne, normal, nc1, nc2, ne;
  REAL wup, wum, wlp, wlm, df, Ac, dt=prop->dt;
  REAL *ap, *am, dznew, dzold, mass, rat;

  ap = phys->ap;
  am = phys->am;

  for(i=0;i<Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->s[i][k];

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->Nk[i];k++) 
      phys->s[i][k]*=(grid->dzzold[i][k]*grid->Ac[i]);
  }

  for(j=0;j<grid->Ne;j++) {
    
    df = grid->df[j];
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    for(k=0;k<grid->Nkc[j];k++) {
      ap[k] = 0.5*(phys->u[j][k]+fabs(phys->u[j][k]));
      am[k] = 0.5*(phys->u[j][k]-fabs(phys->u[j][k]));
    }      

    for(k=0;k<grid->Nkc[j];k++)
      phys->utmp[j][k]=(ap[k]*phys->stmp[nc2][k]*grid->dzzold[nc2][k]+
			am[k]*phys->stmp[nc1][k]*grid->dzzold[nc1][k])*df;
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    Ac = grid->Ac[i];

    for(nf=0;nf<NFACES;nf++) {

      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];

      for(k=0;k<grid->Nk[i];k++) 
      	phys->s[i][k]-=dt*phys->utmp[ne][k]*normal;
    }

    for(k=0;k<grid->Nk[i]+1;k++) {
      ap[k] = 0.5*(phys->w[i][k]+fabs(phys->w[i][k]));
      am[k] = 0.5*(phys->w[i][k]-fabs(phys->w[i][k]));
    }

    for(k=1;k<grid->Nk[i]-1;k++) {
      phys->s[i][k]-=dt*Ac*(ap[k]*phys->stmp[i][k]+
			    am[k]*phys->stmp[i][k-1]-
			    ap[k+1]*phys->stmp[i][k+1]-
			    am[k+1]*phys->stmp[i][k]);
    }
    if(grid->Nk[i]>1) {
      // Top cell only has flux through lower face.
      phys->s[i][0]+=dt*Ac*(ap[1]*phys->stmp[i][1]+
			    am[1]*phys->stmp[i][0]);
      // Bottom cell only has flux through upper face.
      k = grid->Nk[i]-1;
      phys->s[i][grid->Nk[i]-1]-=dt*Ac*(ap[grid->Nk[i]-1]*phys->stmp[i][grid->Nk[i]-1]+
					am[grid->Nk[i]-1]*phys->stmp[i][grid->Nk[i]-2]);
    }
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    mass=0;
    for(k=0;k<=grid->ctop[i];k++) {
      mass+=phys->s[i][k];
      phys->s[i][k]=0;
    }
    //    if(grid->dzz[i][grid->ctop[i]]!=0 && grid->dzz[i][grid->ctop[i]]!=0) 
    if(grid->dzz[i][grid->ctop[i]]!=0) 
      phys->s[i][grid->ctop[i]]=mass/(grid->Ac[i]*grid->dzz[i][grid->ctop[i]]);
    else
      phys->s[i][grid->ctop[i]]=0;

    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      phys->s[i][k]/=(grid->Ac[i]*grid->dzz[i][k]);
  }
}

static void UpdateScalars0(gridT *grid, physT *phys, propT *prop)
{
  int i, j, jptr, iptr, k, nf, kp, km, kstart, kend;
  int Nc=grid->Nc, Ne=grid->Ne, normal, nc1, nc2, ne;
  REAL wup, wum, wlp, wlm, df, Ac, dt=prop->dt;
  REAL *ap, *am, dznew, dzold;

  ap = phys->ap;
  am = phys->am;

  for(i=0;i<Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) {
      phys->stmp[i][k]=phys->s[i][k];
      phys->wtmp[i][k]=0;
    }

  /*
  for(i=0;i<Nc;i++)
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      phys->s[i][k]*=grid->dzzold[i][k]/grid->dzz[i][k];
  */

  for(j=0;j<grid->Ne;j++) {
    
    df = grid->df[j];
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    for(k=0;k<grid->Nke[j];k++) {
      ap[k] = 0.5*(phys->u[j][k]+fabs(phys->u[j][k]));
      am[k] = 0.5*(phys->u[j][k]-fabs(phys->u[j][k]));
    }      

    for(k=0;k<grid->Nke[j];k++)
      phys->utmp[j][k]=(ap[k]*phys->stmp[nc2][k]*grid->dzzold[nc2][k]+
			am[k]*phys->stmp[nc1][k]*grid->dzzold[nc1][k])*df;
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    Ac = grid->Ac[i];

    for(nf=0;nf<NFACES;nf++) {

      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];

      for(k=0;k<grid->ctopold[i];k++)
     	phys->s[i][grid->ctopold[i]]-=dt*phys->utmp[ne][k]*normal/
	  (Ac*grid->dzzold[i][grid->ctopold[i]]);
      for(k=grid->ctopold[i];k<grid->Nk[i];k++)
	phys->s[i][k]-=dt*phys->utmp[ne][k]*normal/
	  (Ac*grid->dzzold[i][k]);
    }

    for(k=0;k<grid->Nk[i]+1;k++) {
      ap[k] = 0.5*(phys->w[i][k]+fabs(phys->w[i][k]));
      am[k] = 0.5*(phys->w[i][k]-fabs(phys->w[i][k]));
    }

    for(k=grid->ctopold[i]+1;k<grid->Nk[i]-1;k++) {
      phys->s[i][k]-=dt*(ap[k]*phys->stmp[i][k]+
			    am[k]*phys->stmp[i][k-1]-
			    ap[k+1]*phys->stmp[i][k+1]-
			    am[k+1]*phys->stmp[i][k])/grid->dzzold[i][k];
    }
    if(grid->Nk[i]-grid->ctopold[i]>1) {
      // Top cell only has flux through lower face.
      k = grid->ctop[i];
      phys->s[i][k]+=dt*(ap[k+1]*phys->stmp[i][k+1]+
			 am[k+1]*phys->stmp[i][k])/grid->dzzold[i][k];
      // Bottom cell only has flux through upper face.
      k = grid->Nk[i]-1;
      phys->s[i][k]-=dt*(ap[k]*phys->stmp[i][k]+
			 am[k]*phys->stmp[i][k-1])/grid->dzzold[i][k];
    }
  }
  for(i=0;i<Nc;i++) {
    if(grid->ctop[i]>=grid->ctopold[i]) {
      // Emptying cells
      dzold = 0;
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++) {
	dzold+=phys->s[i][k]*grid->dzzold[i][k];
	phys->s[i][k]=0;
      }
      phys->s[i][grid->ctop[i]]=dzold/grid->dzz[i][grid->ctop[i]];
    } else {
      // Filling cells
      dznew = 0;
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++)
	dznew+=grid->dzz[i][k];
      //      for(k=grid->ctop[i];k<=grid->ctopold[i];k++)
      //	phys->s[i][k]=phys->s[i][grid->ctopold[i]]*
      //	  grid->dzzold[i][grid->ctopold[i]]/dznew;
    }
  }
  /*   
  for(i=0;i<Nc;i++) {
    if(grid->ctop[i]>=grid->ctopold[i]) {
      dzold = 0;
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
	dzold+=grid->dzzold[i][k];
      phys->s[i][grid->ctop[i]]=phys->s[i][grid->ctop[i]]*
      	(dzold/grid->dzz[i][grid->ctop[i]]);
    } else {
      dznew = 0;
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++)
	dznew+=grid->dzz[i][k];
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++)
	phys->s[i][k]=phys->s[i][grid->ctopold[i]]*
	  grid->dzzold[i][grid->ctopold[i]]/dznew;
    }
  }
  */
} 

static void HydroW0(gridT *grid, physT *phys)
{
  int i, k, nf, nc1, nc2, iptr, ne, normal;
  REAL *ap, *am, df, Ac;

  ap = phys->ap;
  am = phys->am;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    Ac = grid->Ac[i];

    for(k=0;k<grid->Nk[i]+1;k++) {
      phys->wtmp[i][k]=0;
      phys->w[i][k] = 0;
    }

    for(nf=0;nf<NFACES;nf++) {

      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;
      df = grid->df[ne];

      for(k=0;k<grid->Nk[i];k++) {
	ap[k]=0.5*(phys->u[i][k]+fabs(phys->u[i][k]));
	am[k]=0.5*(phys->u[i][k]-fabs(phys->u[i][k]));
      }

      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) 
	phys->wtmp[i][k]+=(ap[k]*grid->dzz[nc2][k]+
			   am[k]*grid->dzz[nc1][k])*df*normal/Ac;
    }

    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      phys->w[i][k]=phys->w[i][k+1]+phys->wtmp[i][k];
    }
  }
}

static void HydroW(gridT *grid, physT *phys)
{
  int i, k, nf, iptr, ne, nc1, nc2;
  REAL ap, am;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    for(k=0;k<grid->ctop[i];k++)
      phys->w[i][k] = 0;
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      phys->w[i][k] = phys->w[i][k+1];
      for(nf=0;nf<NFACES;nf++) {
	ne = grid->face[i*NFACES+nf];
	nc1 = grid->grad[2*ne];
	nc2 = grid->grad[2*ne+1];
	if(nc1==-1) nc1=nc2;
	if(nc2==-1) nc2=nc1;

	ap = 0.5*(phys->u[ne][k]+fabs(phys->u[ne][k]));
	am = 0.5*(phys->u[ne][k]-fabs(phys->u[ne][k]));
	phys->w[i][k]-=(ap*grid->dzz[nc2][k]+am*grid->dzz[nc1][k])*
	  grid->df[ne]*grid->normal[i*NFACES+nf]/grid->Ac[i];
      }
    }
  }
}

static void WtoVerticalFace(gridT *grid, physT *phys)
{
  int j, k, nf, nc1, nc2;
  for(j=0;j<grid->Ne;j++) {
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1 != -1 && nc2 != -1)
      for(k=0;k<grid->Nke[j];k++)
	phys->wf[j][k]=0.25*(phys->w[nc1][k]+
			     phys->w[nc1][k+1]+
			     phys->w[nc2][k]+
			     phys->w[nc2][k+1]);
    else if(nc1 == -1)
      for(k=0;k<grid->Nke[j];k++)
	phys->wf[j][k]=0.5*(phys->w[nc2][k]+
			    phys->w[nc2][k+1]);
    else
      for(k=0;k<grid->Nke[j];k++)
	phys->wf[j][k]=0.5*(phys->w[nc1][k]+
			    phys->w[nc1][k+1]);
  }
}
     
static void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs,
			  MPI_Comm comm)
{
  int i, nc1, nc2, jflux, iptr, j, k, k0, Nc=grid->Nc, Ne=grid->Ne;
  REAL mass, volume, volh, height, Ep, u_barotropic, r_baro;

  if(myproc==0) phys->mass=0;
  if(myproc==0) phys->volume=0;
  if(myproc==0) phys->Ep=0;

  mass=0;
  volume=0;
  volh=0;
  Ep=0;

  phys->Eflux1 = 0;
  phys->Eflux2 = 0;
  phys->Eflux3 = 0;
  phys->Eflux4 = 0;

  //  jflux = 25;
  /*
  jflux = 0;
  u_barotropic = 0;
  nc1 = grid->grad[2*jflux];
  nc2 = grid->grad[2*jflux+1];

  for(k=grid->etop[jflux];k<grid->Nke[jflux];k++)
    u_barotropic+=0.5*phys->u[jflux][k]*(grid->dzz[nc1][k]+
					 grid->dzz[nc2][k]);
  u_barotropic/=(0.5*(grid->dv[nc1]+grid->dv[nc2]));

  for(k=grid->etop[jflux];k<grid->Nke[jflux];k++) {
    r_baro = 0;
    for(k0=grid->etop[jflux];k0<k;k0++)
      r_baro+=0.5*GRAV*prop->beta*((phys->s[nc1][k0]-phys->s0[nc1][k0])*grid->dzz[nc1][k0]+
				   (phys->s[nc2][k0]-phys->s0[nc2][k0])*grid->dzz[nc2][k0]);
    phys->Eflux1+=0.5*r_baro*(phys->u[jflux][k]-u_barotropic)*(grid->dzz[nc1][k]+
							    grid->dzz[nc2][k]);
    phys->Eflux3+=0.5*r_baro*phys->u[jflux][k]*(grid->dzz[nc1][k]+
						grid->dzz[nc2][k]);
  }
  phys->Eflux1*=grid->df[jflux];
  phys->Eflux3*=grid->df[jflux];
  
  //  jflux = 225;
  jflux = 0;
  u_barotropic = 0;
  nc1 = grid->grad[2*jflux];
  nc2 = grid->grad[2*jflux+1];

  for(k=grid->etop[jflux];k<grid->Nke[jflux];k++)
    u_barotropic+=0.5*phys->u[jflux][k]*(grid->dzz[nc1][k]+
					 grid->dzz[nc2][k]);
  u_barotropic/=(0.5*(grid->dv[nc1]+grid->dv[nc2]));

  for(k=grid->etop[jflux];k<grid->Nke[jflux];k++) {
    r_baro = 0;
    for(k0=grid->etop[jflux];k0<k;k0++)
      r_baro+=0.5*GRAV*prop->beta*((phys->s[nc1][k0]-phys->s0[nc1][k0])*grid->dzz[nc1][k0]+
				   (phys->s[nc2][k0]-phys->s0[nc2][k0])*grid->dzz[nc2][k0]);
    phys->Eflux2+=0.5*r_baro*(phys->u[jflux][k]-u_barotropic)*(grid->dzz[nc1][k]+
							    grid->dzz[nc2][k]);
    phys->Eflux4+=0.5*r_baro*phys->u[jflux][k]*(grid->dzz[nc1][k]+
						grid->dzz[nc2][k]);
  }
  phys->Eflux2*=grid->df[jflux];
  phys->Eflux4*=grid->df[jflux];
  */
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    height = 0;
    volh+=grid->Ac[i]*(grid->dv[i]+phys->h[i]);
    Ep+=0.5*grid->Ac[i]*(pow(phys->h[i],2)-pow(grid->dv[i],2));
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      height += grid->dzz[i][k];
      volume+=grid->Ac[i]*grid->dzz[i][k];
      //      Ep+=0.5*grid->Ac[i]*pow(grid->dzz[i][k],2);
      mass+=phys->s[i][k]*grid->Ac[i]*grid->dzz[i][k];
    }
    Ep/=volume;
  
    /*
    printf("d+h=%f, sum=%f, dv=%f, h=%f, Nkdz=%f, dz0=%f\n",
	   (grid->dv[i]+phys->h[i]),height,grid->dv[i],phys->h[i],grid->Nk[i]*grid->dz[0],
	   grid->dz[0]);
    */
  }

  MPI_Reduce(&mass,&(phys->mass),1,MPI_DOUBLE,MPI_SUM,0,comm);
  //MPI_Reduce(&volh,&(phys->volume),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&volume,&(phys->volume),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&Ep,&(phys->Ep),1,MPI_DOUBLE,MPI_SUM,0,comm);

  if(myproc==0) {
    if(prop->n==0) {
      phys->volume0 = phys->volume;
      phys->mass0 = phys->mass;
      phys->Ep0 = phys->Ep;
    } else {
      if(fabs((phys->volume-phys->volume0)/phys->volume0)>CONSERVED && prop->volcheck)
	printf("Warning! Not volume conservative! V(0)=%e, V(t)=%e\n",
	       phys->volume0,phys->volume);
      if(fabs((phys->mass-phys->mass0)/phys->volume0)>CONSERVED && prop->masscheck) 
	printf("Warning! Not mass conservative! M(0)=%e, M(t)=%e\n",
	       phys->mass0,phys->mass);
    }
  }
}
    
static void Check(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, k, ic, kc, Nc=grid->Nc, Ne=grid->Ne;
  int uflag=1, sflag=1, hflag=1, myalldone, alldone;
  REAL C, Cmax;

  for(i=0;i<Nc;i++) 
    if(phys->h[i]!=phys->h[i]) {
	hflag=0;
	break;
    }

  for(i=0;i<Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      if(phys->s[i][k]!=phys->s[i][k]) {
	sflag=0;
	break;
      }
    if(!sflag)
      break;
  }

  for(i=0;i<Ne;i++) {
    for(k=0;k<grid->Nke[i];k++)
      if(phys->u[i][k]!=phys->u[i][k]) {
	uflag=0;
	break;
      }
    if(!uflag)
      break;
  }

  Cmax=0;
  for(i=0;i<Ne;i++) 
    for(k=0;k<grid->Nke[i];k++) {
      C = fabs(phys->u[i][k])*prop->dt/grid->dg[i];
      if(C>Cmax) {
	ic = i;
	kc = k;
	Cmax = C;
      }
    }

  /*
  if(Cmax>prop->Cmax) {
    printf("Warning! U(%d,%d)=%f, dx(%d)=%f, Time step changed from %f to ",
	   ic,kc,phys->u[ic][kc],ic,grid->dg[ic],prop->dt);
    prop->dt = prop->dt*prop->Cmax/Cmax;
    Cmax = 1.0;
    printf("%f\n",prop->dt);
  }
  */

  myalldone=0;
  if(!uflag || !sflag || !hflag || Cmax>prop->Cmax) {
    printf("Time step %d: Processor %d, Run is blowing up!\n",prop->n,myproc);
    
    if(Cmax>prop->Cmax)
      printf("Courant number problems at (%d,%d), Umax=%f, dx=%f Cmax=%.2f > %.2f\n",
	     ic,kc,phys->u[ic][kc],grid->dg[ic],Cmax,prop->Cmax);
    else
      printf("Courant number is okay: Cmax=%.2f < %.2f\n",Cmax,prop->Cmax);
    if(!uflag)
      printf("U is divergent.\n");
    else
      printf("U is okay.\n");
    if(!sflag)
      printf("Scalar is divergent.\n");
    else
      printf("Scalar is okay.\n");
    if(!hflag)
      printf("Free-surface is divergent.\n");
    else
      printf("Free-surface is okay.\n");
    
    myalldone=1;
  }

  MPI_Reduce(&myalldone,&alldone,1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&alldone,1,MPI_INT,0,comm);

  if(alldone) {
    MPI_Finalize();
    exit(0);
  }
}

static void ComputeTangentialVelocity(REAL **u, REAL **ut, gridT *grid) {
  
  int eneigh, j, k, ne;
  // ut will store the tangential component of velocity on each face
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nke[j];k++) 
      ut[j][k]=0;
    for(ne=0;ne<2*(NFACES-1);ne++) {
      eneigh=grid->eneigh[2*(NFACES-1)*j+ne];
      if(eneigh!=-1)
	for(k=0;k<grid->Nke[j];k++)
	  ut[j][k] += grid->xi[2*(NFACES-1)*j+ne]*
	    u[grid->eneigh[2*(NFACES-1)*j+ne]][k];
    }
    if(grid->grad[2*j]!=-1 && grid->grad[2*j+1]!=-1)
      for(k=0;k<grid->Nke[j];k++)
	ut[j][k]*=0.5;
  }
  
}

static void OutputData(gridT *grid, physT *phys, propT *prop,
		int myproc, int numprocs, MPI_Comm comm)
{
  int i, j, k, ne, eneigh, nwritten;
  REAL *tmp = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"OutputData");

  if(!(prop->n%prop->ntconserve)) {
    ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);
    if(myproc==0)
      fprintf(prop->ConserveFID,"%e %e %e %e %e %e %e %e\n",prop->rtime,phys->mass,phys->volume,
	      phys->Ep-phys->Ep0,phys->Eflux1,phys->Eflux2,phys->Eflux3,phys->Eflux4);
  }

  if(!(prop->n%prop->ntout) || prop->n==1) {

    if(myproc==0 && VERBOSE>1) printf("Outputting data at step %d of %d\n",prop->n,prop->nsteps);

    // Transfer back and forth all the data to all of the ghost points.
    // This is mostly so that ghost cells do not have empty values in the data.
    // It is NOT necessarily required for the computation!
    ISendRecvCellData2D(phys->h,grid,myproc,comm,all);
    SendRecvCellData3D(phys->s,grid,myproc,comm,all);

    if(ASCII) 
      for(i=0;i<grid->Nc;i++)
	fprintf(prop->FreeSurfaceFID,"%f\n",phys->h[i]);
    else {
      nwritten=fwrite(phys->h,sizeof(REAL),grid->Nc,prop->FreeSurfaceFID);
      if(nwritten!=grid->Nc) {
	printf("Error outputting free-surface data!\n");
	exit(EXIT_WRITING);
      }
    }
    fflush(prop->FreeSurfaceFID);

    // ut stores the tangential component of velocity on the faces.
    if(ASCII) 
      for(j=0;j<grid->Ne;j++) {
	for(k=0;k<grid->Nke[j];k++) 
	  fprintf(prop->HorizontalVelocityFID,"%e %e %e\n",
		  phys->u[j][k]*grid->n1[j]-phys->ut[j][k]*grid->n2[j],
		  phys->u[j][k]*grid->n2[j]+phys->ut[j][k]*grid->n1[j],
		  phys->wf[j][k]);
	for(k=grid->Nke[j];k<grid->Nkmax;k++)
	  fprintf(prop->HorizontalVelocityFID,"0.0 0.0 0.0\n");
      }
    else 
      for(k=0;k<grid->Nkmax;k++) {
	for(j=0;j<grid->Ne;j++) {
	  if(k<grid->Nke[j]) 
	    tmp[j]=phys->u[j][k]*grid->n1[j]-0*phys->ut[j][k]*grid->n2[j];
	  else
	    tmp[j]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Ne,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Ne) {
	  printf("Error outputting Horizontal Velocity data!\n");
	  exit(EXIT_WRITING);
	}
	for(j=0;j<grid->Ne;j++) {
	  if(k<grid->Nke[j])
	    tmp[j]=phys->u[j][k]*grid->n2[j]+0*phys->ut[j][k]*grid->n1[j];
	  else
	    tmp[j]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Ne,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Ne) {
	  printf("Error outputting Horizontal Velocity data!\n");
	  exit(EXIT_WRITING);
	}
	for(j=0;j<grid->Ne;j++) {
	  if(k<grid->Nke[j])
	    tmp[j]=phys->wf[j][k];
	  else
	    tmp[j]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Ne,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Ne) {
	  printf("Error outputting Face Velocity data!\n");
	  exit(EXIT_WRITING);
	}
      }
    fflush(prop->HorizontalVelocityFID);

    if(ASCII)
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i]+1;k++)
	  fprintf(prop->VerticalVelocityFID,"%e\n",phys->w[i][k]);
	for(k=grid->Nk[i]+1;k<grid->Nkmax+1;k++)
	  fprintf(prop->VerticalVelocityFID,"0.0\n");
      }
    else {
      for(k=0;k<grid->Nkmax+1;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i]+1)
	    phys->htmp[i]=phys->w[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->VerticalVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting vertical velocity data!\n");
	  exit(EXIT_WRITING);
	}
      }
    }
    fflush(prop->VerticalVelocityFID);

    if(ASCII) {
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++)
	  fprintf(prop->SalinityFID,"%e\n",phys->s[i][k]);
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->SalinityFID,"0.0\n");
      }
      if(prop->n==1) {
	for(i=0;i<grid->Nc;i++) {
	  for(k=0;k<grid->Nk[i];k++)
	    fprintf(prop->BGSalinityFID,"%e\n",phys->s0[i][k]);
	  for(k=grid->Nk[i];k<grid->Nkmax;k++)
	    fprintf(prop->BGSalinityFID,"0.0\n");
	}
      }
    } else {
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    phys->htmp[i]=phys->s[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->SalinityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting salinity data!\n");
	  exit(EXIT_WRITING);
	}
      }
      if(prop->n==1) 
	for(k=0;k<grid->Nkmax;k++) {
	  for(i=0;i<grid->Nc;i++) {
	    if(k<grid->Nk[i])
	      phys->htmp[i]=phys->s0[i][k];
	    else
	      phys->htmp[i]=EMPTY;
	  }
	  nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->BGSalinityFID);
	  if(nwritten!=grid->Nc) {
	    printf("Error outputting background salinity data!\n");
	    exit(EXIT_WRITING);
	  }
	}
    }
    fflush(prop->SalinityFID);
    
    for(i=0;i<grid->Nc;i++) {
      for(k=0;k<grid->Nk[i];k++)
	fprintf(prop->VerticalGridFID,"%e\n",grid->dzz[i][k]);
      for(k=grid->Nk[i];k<grid->Nkmax;k++)
	fprintf(prop->VerticalGridFID,"0.0\n");
    }
  }
  fflush(prop->VerticalGridFID);

  if(prop->n==1)
    fclose(prop->BGSalinityFID);

  if(prop->n==prop->nsteps) {
    fclose(prop->FreeSurfaceFID);
    fclose(prop->HorizontalVelocityFID);
    fclose(prop->VerticalVelocityFID);
    fclose(prop->SalinityFID);
    fclose(prop->VerticalGridFID);
    if(myproc==0) fclose(prop->ConserveFID);
  }

  SunFree(tmp,grid->Ne*sizeof(REAL),"OutputData");
}

static void ReadProperties(propT **prop, int myproc)
{
  *prop = (propT *)SunMalloc(sizeof(propT),"ReadProperties");
  
  (*prop)->theta = MPI_GetValue(DATAFILE,"theta","ReadProperties",myproc);
  (*prop)->thetaAB = MPI_GetValue(DATAFILE,"thetaAB","ReadProperties",myproc);
  (*prop)->thetaFS = MPI_GetValue(DATAFILE,"thetaFS","ReadProperties",myproc);
  (*prop)->beta = MPI_GetValue(DATAFILE,"beta","ReadProperties",myproc);
  (*prop)->nu = MPI_GetValue(DATAFILE,"nu","ReadProperties",myproc);
  (*prop)->tau_T = MPI_GetValue(DATAFILE,"tau_T","ReadProperties",myproc);
  (*prop)->CdT = MPI_GetValue(DATAFILE,"CdT","ReadProperties",myproc);
  (*prop)->CdB = MPI_GetValue(DATAFILE,"CdB","ReadProperties",myproc);
  (*prop)->dt = MPI_GetValue(DATAFILE,"dt","ReadProperties",myproc);
  (*prop)->Cmax = MPI_GetValue(DATAFILE,"Cmax","ReadProperties",myproc);
  (*prop)->nsteps = (int)MPI_GetValue(DATAFILE,"nsteps","ReadProperties",myproc);
  (*prop)->ntout = (int)MPI_GetValue(DATAFILE,"ntout","ReadProperties",myproc);
  (*prop)->ntprog = (int)MPI_GetValue(DATAFILE,"ntprog","ReadProperties",myproc);
  (*prop)->ntconserve = (int)MPI_GetValue(DATAFILE,"ntconserve","ReadProperties",myproc);
  (*prop)->maxiters = (int)MPI_GetValue(DATAFILE,"maxiters","ReadProperties",myproc);
  (*prop)->epsilon = MPI_GetValue(DATAFILE,"epsilon","ReadProperties",myproc);
  (*prop)->relax = MPI_GetValue(DATAFILE,"relax","ReadProperties",myproc);
  (*prop)->amp = MPI_GetValue(DATAFILE,"amp","ReadProperties",myproc);
  (*prop)->omega = MPI_GetValue(DATAFILE,"omega","ReadProperties",myproc);
  (*prop)->flux = MPI_GetValue(DATAFILE,"flux","ReadProperties",myproc);
  (*prop)->volcheck = MPI_GetValue(DATAFILE,"volcheck","ReadProperties",myproc);
  (*prop)->masscheck = MPI_GetValue(DATAFILE,"masscheck","ReadProperties",myproc);
  (*prop)->kriging_cov = MPI_GetValue(DATAFILE,"kriging_cov","ReadProperties",myproc);
  (*prop)->nonlinear = MPI_GetValue(DATAFILE,"nonlinear","ReadProperties",myproc);
}

static void OpenFiles(propT *prop, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  MPI_GetString(filename,DATAFILE,"FreeSurfaceFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->FreeSurfaceFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetString(filename,DATAFILE,"HorizontalVelocityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->HorizontalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetString(filename,DATAFILE,"VerticalVelocityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->VerticalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetString(filename,DATAFILE,"SalinityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->SalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetString(filename,DATAFILE,"BGSalinityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->BGSalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetString(filename,DATAFILE,"VerticalGridFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->VerticalGridFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  if(myproc==0) {
    MPI_GetString(filename,DATAFILE,"ConserveFile","OpenFiles",myproc);
    sprintf(str,"%s",filename);
    prop->ConserveFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  }
}

static void Progress(propT *prop, int myproc) 
{
  int progout, prog;

  if(myproc==0 && prop->ntprog>0 && VERBOSE>0) {
    progout = (int)(prop->nsteps*(double)prop->ntprog/100);
    prog=(int)(100.0*(double)prop->n/(double)prop->nsteps);
    if(progout>0)
      if(!(prop->n%progout))
	printf("%d%% Complete.\n",prog);
  }
}

static void InitializeKriging(gridT *grid, int pow) {
  int j, jptr, n;
  REAL *x = (REAL *)SunMalloc(grid->Nnearestedges*sizeof(REAL),"InitializeKriging");
  REAL *y = (REAL *)SunMalloc(grid->Nnearestedges*sizeof(REAL),"InitializeKriging");

  grid->Kriging = (REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"InitializeKriging");

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    for(n=0;n<grid->Nnearestedges;n++) {
      x[n] = grid->xe[grid->nearestedges[j][n]];
      y[n] = grid->ye[grid->nearestedges[j][n]];
    }
    
    grid->Kriging[j]=(REAL *)inversekriging(x,y,pow,grid->Nnearestedges);
  }

  SunFree(x,grid->Nnearestedges*sizeof(REAL),"InitializeKriging");
  SunFree(y,grid->Nnearestedges*sizeof(REAL),"InitializeKriging");
}

