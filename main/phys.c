/*
 * File: phys.c
 * Author: Oliver Fringer
 * Institution: Stanford University
 * Date: 10/21/02
 * --------------------------------
 * This file contains physically-based functions.
 *
 * $Id: phys.c,v 1.50 2004-05-15 00:05:36 fringer Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.49  2004/05/14 02:26:15  fringer
 * Removed the theta-AB method, also only sum up until <k and not <=k
 * in the baroclinic term.
 *
 * Added a ramping down of the value of theta over an amount of time
 * given by thetaramptime to reduce transient oscillations.
 *
 * Revision 1.48  2004/04/24 00:36:50  fringer
 * Working out some bugs with wetting/drying...
 *
 * Revision 1.47  2004/04/22 04:54:32  fringer
 * Added wtmp2 and stmp3 to store velocity and salinity at time step
 * n in order to implement theta method for buoyancy term in
 * momentum equation.  Added StoreVariables to store u,w, and s at
 * time step n, and also added updatescalars before advecthorizontal
 * velocity to compute the scalar qty at step n+1 with an explicit
 * update (theta=0 in updatescalars).  Also added the scalold term
 * in the updatescalars function to use the old scalar value in
 * updatescalars.
 *
 * Also changed the explicit terms in updatescalars to use the
 * values at time step n instead of step n+1.
 *
 * Revision 1.46  2004/04/22 04:15:27  fringer
 * Works for wetting/drying and Eulerian advection:
 * 1) Horizontal advection computed only for edges k>=etop[j]+1
 * 2) NewCells() sets new cell velocities equal to u[etopold[j]]. This
 * is effectively a zeroth order extrapolation.
 * 3) Vertical advection only in cells k>=ctop[i]+2
 * 4) Nonhydrostatic correction only performed for columns
 * with etop[j]<grid->Nke[j]-1
 *
 * Revision 1.45  2004/04/21 02:32:49  fringer
 * Added the k0=grid->etop[j] to the baroclinic term and had the
 * summation for the integral changed so that it includes the bottom
 * cell in the summation.
 *
 * Revision 1.44  2004/04/14 13:31:22  fringer
 * Added  QCoefficients to compute coefficients for pressure solver
 * at each time step in order to speed up CG (preconditioned) solution
 *
 * Revision 1.43  2004/04/06 00:18:30  fringer
 * Added computation of flux at boundary faces (of type 2) in
 * nonlinear advection of momentum.
 *
 * Revision 1.42  2004/03/16 17:59:39  fringer
 * Set ut to 0 in the beginning of advecthorizontal velocity since
 * this was causing horizontal advection errors.
 *
 * Revision 1.41  2004/03/13 02:48:52  fringer
 * Added ability to switch between upwind and central in advection
 * scheme (nonlinear: 1, upwind ; 2, central).  Also fixed some advection
 * bugs so that fluxes are not computed outside of the grid at
 * surface faces.
 *
 * Revision 1.40  2004/03/13 00:43:00  fringer
 * Working lock-exchange case -> Nonhydrostatic, nonlinear, central differencing
 * for momentum, upwind for scalar.
 *
 * Revision 1.39  2004/03/12 06:24:41  fringer
 * Removed OpenBoundaryFluxes function and added it to its own .c file
 * boundaries.c, where all boundary values are defined.
 *
 * Revision 1.38  2004/03/12 06:14:22  fringer
 * Working version for Huntington beach case.
 *
 * Revision 1.37  2004/03/05 19:47:01  fringer
 * Vertical advection not working for estuarine case.
 *
 * Revision 1.36  2004/01/27 05:25:49  fringer
 * Working version for Monterey Bay case.  OpenBoundaryFluxes still not working.
 * Had to add send/recv routines in advecthorizvelocity so that stmp/2 terms were
 * defined at the boundary nodes.
 *
 * Revision 1.35  2003/12/11 01:42:12  fringer
 * Changed momentum advection scheme so that source terms stmp and stmp2
 * are conservative when summed up.  Also added ability to perform
 * advection using either upwind or central differencing.
 *
 * Changed operatorq so that it works with vertical walls using kmin and
 * kmax.
 *
 * Added NewCells routine which is an attempt to conserve momentum due
 * to changes in volume of the upper layer, but this is still not working
 * very well.
 *
 * Revision 1.34  2003/12/08 20:02:16  fringer
 * Changed theta to 1-theta in hydrow when adding the explicit part
 * of the nonhydrostatic pressure gradient to the open boundary
 * fluxes.
 *
 * Revision 1.33  2003/12/08 20:01:06  fringer
 * Added OpenBoundaryFluxes for BC type 2, which forces the open boundary
 * with velocity and damps out fluctuations.  Also set q=0 at open boundaries
 * so this required changes to operatorq and corrector, as well as
 * hydrow, which adds the explicit part of the pressure gradient to
 * the velocities at the open boundary.
 *
 * Revision 1.32  2003/12/03 04:05:52  fringer
 * Moved the volume correction for the horizontal velocity to
 * the Corrector function so that the post-pressure velocity field
 * is corrected.
 *
 * Revision 1.31  2003/12/03 02:07:02  fringer
 * Added output to a progress file that shows the current time step.
 * This file is specified as the ProgressFile in suntans.dat
 *
 * Also added lines for open boundary types (marker=2) and forced
 * flux types (marker=4).
 *
 * Also set the boundary condition for the cgsolveq to no gradient
 * by setting D[j]=df[j]/dg[j] and setting D[j] to 0 on faces wehre
 * there should be no q gradient.
 *
 * Revision 1.30  2003/12/02 23:33:57  fringer
 * Added error checking code in updatedz and also updated cgsolve (q and h)
 * so that the iteration stops if eps0=0.
 *
 * Revision 1.29  2003/12/02 04:38:27  fringer
 * Changed send/recv routines to be nonblocking.
 *
 * Revision 1.28  2003/12/02 02:31:44  fringer
 * Removed ability to transfer first or all in send/recv routines.
 *
 * Revision 1.27  2003/12/02 02:15:54  fringer
 * Removed all traces of kriging, including findnearestplane and
 * findnearestedges.
 *
 * Revision 1.26  2003/12/02 01:05:54  fringer
 * Working nonhydrostatic version.  FreeGrid and FreePhysicalVariables are
 * causing problems.
 *
 * Revision 1.25  2003/10/29 18:11:11  fringer
 * Changed boundary condition adjustment for free-surface flux so
 * that it computes fluxes for mark=3 edges as opposed to grad!=-1
 * and n1!=0.
 *
 * Revision 1.24  2003/10/28 04:08:07  fringer
 * Implemented theta method for scalar advection UpdataUImp.
 * Used to be theta=1.  Theta value is given by prop->thetaS.
 *
 * Revision 1.23  2003/10/27 21:44:02  fringer
 * Made explicit part of advection of scalar 2nd order with AB2.
 * Added Cn_U,Cn_R,and Cn_T to store the time n-1 terms for the AB2.
 *
 * Revision 1.22  2003/10/27 17:22:58  fringer
 * Added spong layer and changed loops in momentum advection so that
 * stmp and stmp2 are computed in the interprocessor ghost cells.
 *
 * Revision 1.21  2003/10/25 02:04:08  fringer
 * Removed Kriging terms and added baroclinic and cariolis terms to AB vector
 * Cn[][][].
 *
 * Revision 1.20  2003/10/25 00:55:09  fringer
 * Added Eulerian scheme with AB for horizontal advection.
 * Before removal of Kriging terms.
 *
 * Revision 1.19  2003/09/16 22:29:51  fringer
 * Added ComputeVelocityVector to compute velocity vector at faces and also changed outputdata so that horizontal velocity is output at cell centers.
 *
 * Revision 1.18  2003/06/10 03:28:46  fringer
 * Changed MPI_GetString to MPI_GetFile
 * Changed output of total memory to be in units of Mb
 *
 * Revision 1.17  2003/06/10 02:30:23  fringer
 * Changed ELM scheme so that kriging matrices are computed during the
 * first time step.  Inverses are therefore not computed at each time
 * step, but rather, they are stored in the 2d array grid->kriging[][]
 *
 * Also changed malloc/free to SunMalloc/SunFree.
 *
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
static void Corrector(REAL **qc, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs);
static void CGSolve(gridT *grid, physT *phys, propT *prop, 
		    int myproc, int numprocs, MPI_Comm comm);
static void CGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys, propT *prop, 
		     int myproc, int numprocs, MPI_Comm comm);
static void ConditionQ(REAL **x, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
static void GuessQ(REAL **q, REAL **wold, REAL **w, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void GSSolve(gridT *grid, physT *phys, propT *prop, 
		    int myproc, int numprocs, MPI_Comm comm);
static REAL InnerProduct(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static REAL InnerProduct3(REAL **x, REAL **y, gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static void OperatorH(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop);
static void OperatorQC(REAL **coef, REAL **fcoef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop);
static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, physT *phys, propT *prop);
static void OperatorQ(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop);
static void HydroW(REAL **w, gridT *grid, physT *phys, propT *prop);
static void WtoVerticalFace(gridT *grid, physT *phys);
static void UpdateScalars(gridT *grid, physT *phys, propT *prop, REAL **scalold, REAL **scal, REAL **Cn, REAL theta);
static void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs,
			  MPI_Comm comm);
static void Check(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void Progress(propT *prop, int myproc);
static void EddyViscosity(gridT *grid, physT *phys, propT *prop);
static void AdvectHorizontalVelocity(gridT *grid, physT *phys, propT *prop,
				     int myproc, int numprocs, MPI_Comm comm);
static void StoreVariables(gridT *grid, physT *phys);
static void NewCells(gridT *grid, physT *phys, propT *prop);
static void AdvectVerticalVelocity(gridT *grid, physT *phys, propT *prop,
				   int myproc, int numprocs);
static void ComputeTangentialVelocity(REAL **u, REAL **ut, gridT *grid);
static void AdjustBoundaryFluxes(REAL **u, REAL **uc, REAL **vc, gridT *grid);
static void ComputeNormalVelocity(REAL **u, REAL **uc, REAL **vc, gridT *grid);
static void ComputeVelocityVector(REAL **u, REAL **uc, REAL **vc, gridT *grid);
static void OutputData(gridT *grid, physT *phys, propT *prop,
		int myproc, int numprocs, int blowup, MPI_Comm comm);

void AllocatePhysicalVariables(gridT *grid, physT **phys)
{
  int flag=0, i, j, Nc=grid->Nc, Ne=grid->Ne, nf;

  *phys = (physT *)SunMalloc(sizeof(physT),"AllocatePhysicalVariables");

  (*phys)->u = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->uc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->uold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->D = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->utmp = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->utmp2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->ut = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_U = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");

  (*phys)->wf = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");

  for(j=0;j<Ne;j++) {
    if(grid->Nkc[j]<grid->Nke[j]) {
      printf("Error!  Nkc(=%d)<Nke(=%d) at edge %d\n",grid->Nkc[j],grid->Nke[j],j);
      flag = 1;
    }
    (*phys)->u[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->ut[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_U[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
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
  (*phys)->wtmp2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_W = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->q = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->qtmp = (REAL **)SunMalloc(NFACES*Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s0 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_R = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp3 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->nu_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->tau_T = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->tau_B = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdT = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdB = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  
  for(i=0;i<Nc;i++) {
    (*phys)->uc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->uold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp2[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_W[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->q[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    for(nf=0;nf<NFACES;nf++)
      (*phys)->qtmp[i*NFACES+nf] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s0[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_R[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp2[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp3[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
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
  int i, j, Nc=grid->Nc, Ne=grid->Ne, nf;
  
  for(j=0;j<Ne;j++) {
    free(phys->u[j]);
    free(phys->utmp[j]);
    free(phys->utmp2[j]);
    free(phys->ut[j]);
    free(phys->Cn_U[j]);
    free(phys->wf[j]);
  }

  for(i=0;i<Nc;i++) {
    free(phys->uc[i]);
    free(phys->vc[i]);
    free(phys->uold[i]);
    free(phys->vold[i]);
    free(phys->w[i]);
    free(phys->wtmp[i]);
    free(phys->wtmp2[i]);
    free(phys->Cn_W[i]);
    free(phys->q[i]);
    for(nf=0;nf<NFACES;nf++)
      free(phys->qtmp[i*NFACES+nf]);
    free(phys->s[i]);
    free(phys->T[i]);
    free(phys->s0[i]);
    free(phys->Cn_R[i]);
    free(phys->Cn_T[i]);
    free(phys->stmp[i]);
    free(phys->stmp2[i]);
    free(phys->stmp3[i]);
    free(phys->nu_tv[i]);
  }

  free(phys->h);
  free(phys->htmp);
  free(phys->uc);
  free(phys->vc);
  free(phys->w);
  free(phys->wtmp);
  free(phys->wtmp2);
  free(phys->Cn_W);
  free(phys->wf);
  free(phys->q);
  free(phys->qtmp);
  free(phys->s);
  free(phys->T);
  free(phys->s0);
  free(phys->Cn_R);
  free(phys->Cn_T);
  free(phys->stmp);
  free(phys->stmp2);
  free(phys->stmp3);
  free(phys->nu_tv);
  free(phys->tau_T);
  free(phys->tau_B);
  free(phys->CdT);
  free(phys->CdB);
  free(phys->u);
  free(phys->D);
  free(phys->utmp);
  free(phys->ut);
  free(phys->Cn_U);

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
    
/*
 * Function: ReadPhysicalVariables
 * Usage: ReadPhysicalVariables(grid,phys,prop,myproc);
 * ----------------------------------------------------
 * 
 * Need to take care of uold here!
 *
 */
void ReadPhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc) {

  int i, j;

  if(VERBOSE>1 && myproc==0) printf("Reading from rstore...\n");
    
  fread(&(prop->nstart),sizeof(int),1,prop->StartFID);

  fread(phys->h,sizeof(REAL),grid->Nc,prop->StartFID);
  for(j=0;j<grid->Ne;j++) 
    fread(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(j=0;j<grid->Ne;j++) 
    fread(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  for(i=0;i<grid->Nc;i++) 
    fread(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StartFID);
  fclose(prop->StartFID);

  UpdateDZ(grid,phys,0);
  ComputeVelocityVector(phys->u,phys->uc,phys->vc,grid);
}


void InitializePhysicalVariables(gridT *grid, physT *phys, propT *prop)
{
  int i, j, k, Nc=grid->Nc;
  REAL z;

  prop->nstart=0;

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
      phys->T[i][k]=0;
      phys->s0[i][k]=0;
    }
  }

  for(i=0;i<Nc;i++) {
    z = 0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //z-=grid->dzz[i][k]/2;
      z-=grid->dz[k]/2;
      phys->T[i][k]=ReturnTemperature(grid->xv[i],grid->yv[i],z,grid->dv[i]);
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
  ComputeVelocityVector(phys->u,phys->uc,phys->vc,grid);
  ComputeVelocityVector(phys->u,phys->uold,phys->vold,grid);

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

void Solve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int n;
  extern int TotSpace;
  FILE *fid = fopen("/tmp/fs.dat","w");

  prop->n=0;
  ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);

  if(VERBOSE>1) printf("Processor %d,  Total memory: %d Mb\n",myproc,(int)(TotSpace/(1024*1e3)));
  
  prop->theta0=prop->theta;

  for(n=prop->nstart+1;n<=prop->nsteps+prop->nstart;n++) {
    prop->n = n;
    prop->rtime = (n-1)*prop->dt;
    if(prop->nsteps>0) {

      if(prop->thetaramptime!=0)
	prop->theta=(1-exp(-prop->rtime/prop->thetaramptime))*prop->theta0+
	  exp(-prop->rtime/prop->thetaramptime);

      EddyViscosity(grid,phys,prop);

      StoreVariables(grid,phys);

      if(VERBOSE>2 && myproc==0) printf("Creating horizontal source term...\n");
      //      if(prop->beta) UpdateScalars(grid,phys,prop,phys->stmp3,phys->s,phys->Cn_R,0);
      AdvectHorizontalVelocity(grid,phys,prop,myproc,numprocs,comm);

      if(VERBOSE>2 && myproc==0) printf("Solving for h and hydrostatic u (predictor)...\n");
      BarotropicPredictor(grid,phys,prop,myproc,numprocs,comm);
      ISendRecvCellData2D(phys->h,grid,myproc,comm);

      if(prop->newcells) NewCells(grid,phys,prop);
      ISendRecvEdgeData3D(phys->u,grid,myproc,comm);

      if(prop->nonhydrostatic) {
	//	  WtoVerticalFace(grid,phys);
	AdvectVerticalVelocity(grid,phys,prop,myproc,numprocs);

	ComputeQSource(phys->stmp,grid,phys,prop,myproc,numprocs);
	//	GuessQ(phys->stmp2,phys->w,phys->uold,grid,phys,prop,myproc,numprocs,comm);
	CGSolveQ(phys->stmp2,phys->stmp,phys->uold,grid,phys,prop,myproc,numprocs,comm);
	Corrector(phys->stmp2,grid,phys,prop,myproc,numprocs,comm);
      }

      HydroW(phys->w,grid,phys,prop);
      ISendRecvWData(phys->w,grid,myproc,comm);

      if(prop->beta) {
	UpdateScalars(grid,phys,prop,phys->stmp3,phys->s,phys->Cn_R,prop->thetaS);
	ISendRecvCellData3D(phys->s,grid,myproc,comm);
      }
      if(prop->gamma) {
	UpdateScalars(grid,phys,prop,phys->T,phys->T,phys->Cn_T,prop->thetaS);
	ISendRecvCellData3D(phys->T,grid,myproc,comm);
      }

      // utmp2 contains the velocity field at time step n, u contains
      // it at time step n+1.  This is so that at the next time step
      // phys->uold contains velocity at time step n-1 and phys->uc contains
      // that at time step n.
      ComputeVelocityVector(phys->utmp2,phys->uold,phys->vold,grid);
      ComputeVelocityVector(phys->u,phys->uc,phys->vc,grid);
    }
    fprintf(fid,"%f %f\n",prop->rtime,phys->h[0]);
    fflush(fid);

    Progress(prop,myproc);
    OutputData(grid,phys,prop,myproc,numprocs,0,comm);
    Check(grid,phys,prop,myproc,numprocs,comm);
  }
}

static void StoreVariables(gridT *grid, physT *phys) {
  int i, j, k, iptr, jptr;

  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) {
      phys->stmp3[i][k]=phys->s[i][k];
      phys->wtmp2[i][k]=phys->w[i][k];
    }

  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=0;
    for(k=0;k<grid->Nke[j];k++)
      phys->utmp2[j][k]=phys->u[j][k];
  }
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
				     int myproc, int numprocs, MPI_Comm comm) {
  int i, nf, j, jptr, k, nc, nc1, nc2, ne, k0, kmin, kmax, wetdry_offset;
  REAL *a, *b, *c, fab, sum;

  if(prop->wetdry)
    wetdry_offset=1;
  else
    wetdry_offset=0;

  a = phys->a;
  b = phys->b;
  c = phys->c;

  if(prop->n==1) {
    fab=1;
    for(j=0;j<grid->Ne;j++)
      for(k=0;k<grid->Nke[j];k++)
	phys->Cn_U[j][k]=0;
  } else
    fab=1.5;

  // Set utmp and ut to zero.
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nke[j];k++) {
      phys->utmp[j][k]=0;
      phys->ut[j][k]=0;
    }
  }

  // Set utmp=u if this is a linear calculation.
  if(0) { //!prop->nonlinear) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];

      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->utmp[j][k]=phys->u[j][k];

      // New cells attain the velocity of the cell they emerged from.
      // For the nonlinear case the velocity field is advected into it.
      if(grid->etop[j]<grid->etopold[j])
	for(k=grid->etop[j];k<grid->etopold[j];k++)
	  phys->utmp[j][k]=phys->u[j][grid->etopold[j]];
    }
  } else {    

    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];

      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->utmp[j][k]=(1-fab)*phys->Cn_U[j][k]+0*phys->u[j][k]
	  -(1-prop->theta)*prop->dt/grid->dg[j]*(phys->q[nc1][k]-phys->q[nc2][k])
	  +(1.0-prop->dt*exp(-0.5*(grid->xv[nc1]+grid->xv[nc2])/prop->sponge_distance)/
	    prop->sponge_decay)*phys->u[j][k];

	phys->Cn_U[j][k]=0;
      }
    }

    // Add the Coriolis Terms
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	phys->Cn_U[j][k]+=prop->dt*prop->Coriolis_f*(0.5*(phys->vc[nc1][k]+phys->vc[nc2][k])*grid->n1[j]-
						     0.5*(phys->uc[nc1][k]+phys->uc[nc2][k])*grid->n2[j]);
    }
    
    // Add on the baroclinic term and the pressure gradient from the previous time step
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	k0=grid->etop[j];
	//	phys->Cn_U[j][k]-=0.5*prop->dt*GRAV*prop->beta*
	//	  (prop->theta*(phys->s[nc1][k0]+phys->s[nc2][k0])+(1-prop->theta)*(phys->stmp3[nc1][k0]+phys->stmp3[nc2][k0]))*
	//	  (grid->dzz[nc1][k0]-grid->dzz[nc2][k0])/grid->dg[j];

	for(k0=grid->etop[j];k0<k;k0++)
	  phys->Cn_U[j][k]-=0.5*GRAV*prop->beta*prop->dt*(phys->s[nc1][k0]-phys->s[nc2][k0])*
	    (grid->dzz[nc1][k0]+grid->dzz[nc2][k0])/grid->dg[j];
      }
      /*
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	for(k0=grid->etop[j];k0<=k;k0++)
	  phys->Cn_U[j][k]-=GRAV*prop->beta*prop->dt*phys->s[nc1][k0]*grid->dzz[nc1][k0]/grid->dg[j];
	for(k0=grid->etop[j];k0<=k;k0++)
	  phys->Cn_U[j][k]+=GRAV*prop->beta*prop->dt*phys->s[nc2][k0]*grid->dzz[nc2][k0]/grid->dg[j];
      }
      */
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

    // Compute Eulerian advection of momentum (nonlinear!=0)
    if(prop->nonlinear) {

      // U-fluxes at boundary cells
      for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
	j = grid->edgep[jptr];
	
	for(k=grid->etop[j];k<grid->Nke[j];k++)
	  phys->ut[j][k]=phys->u[j][k]*grid->dzz[grid->grad[2*j]][k]*grid->n1[j];
      }

      // Compute the u-component fluxes at the faces
      if(prop->nonlinear==1)  // Upwind
	for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	  j = grid->edgep[jptr];
	  
	  nc1 = grid->grad[2*j];
	  nc2 = grid->grad[2*j+1];
	  
	  if(grid->ctop[nc1]>grid->ctop[nc2])
	    kmin = grid->ctop[nc1];
	  else
	    kmin = grid->ctop[nc2];

	  for(k=0;k<kmin;k++)
	    phys->ut[j][k]=0;

	  for(k=kmin;k<grid->Nke[j];k++) {
	    if(phys->u[j][k]>0)
	      nc=nc2;
	    else
	      nc=nc1;
	    phys->ut[j][k]=phys->uc[nc][k]*grid->dzz[nc][k];
	  }
	}
      else // Central
	for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	  j = grid->edgep[jptr];
	  
	  nc1 = grid->grad[2*j];
	  nc2 = grid->grad[2*j+1];
	  
	  if(grid->ctop[nc1]>grid->ctop[nc2])
	    kmin = grid->ctop[nc1];
	  else
	    kmin = grid->ctop[nc2];

	  for(k=0;k<kmin;k++)
	    phys->ut[j][k]=0;

	  for(k=kmin;k<grid->Nke[j];k++) 
	    phys->ut[j][k]=0.5*(phys->uc[nc1][k]*grid->dzz[nc1][k]+
				phys->uc[nc2][k]*grid->dzz[nc2][k]);
	}
      
      // Now compute the cell-centered source terms and put them into stmp and also
      // add on diffusion of u.
      for(i=0;i<grid->Nc;i++) {
	
	for(k=0;k<grid->Nk[i];k++) 
	  phys->stmp[i][k]=0;
	
	for(nf=0;nf<NFACES;nf++) {
	  
	  ne = grid->face[i*NFACES+nf];
	  nc = grid->neigh[i*NFACES+nf];

	  if(grid->Nk[nc]<grid->Nk[i])
	    kmax = grid->Nk[nc];
	  else
	    kmax = grid->Nk[i];

	  if(grid->ctop[nc]>grid->ctop[i])
	    kmin = grid->ctop[nc];
	  else
	    kmin = grid->ctop[i];

	  for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	    phys->stmp[i][k]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	      (grid->Ac[i]*grid->dzz[i][k]);
	  
	  if(nc!=-1) {
	    for(k=kmin;k<kmax;k++)
	      phys->stmp[i][k]-=prop->nu_H*(phys->uc[nc][k]-phys->uc[i][k])*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	    for(k=kmax;k<grid->Nkc[ne];k++)
	      phys->stmp[i][k]+=2*prop->nu_H*phys->uc[i][k]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	  }

	  // Top cell is filled with momentum from neighboring cells
	  for(k=grid->etop[ne];k<=grid->ctop[i];k++) 
	    phys->stmp[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	      (grid->Ac[i]*grid->dzz[i][grid->ctop[i]]);
	}
      }

      // V-fluxes at boundary cells
      for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
	j = grid->edgep[jptr];
	
	for(k=grid->etop[j];k<grid->Nke[j];k++)
	  phys->ut[j][k]=phys->u[j][k]*grid->dzz[grid->grad[2*j]][k]*grid->n2[j];
      }
      
      // Compute the v-component fluxes at the faces
      if(prop->nonlinear==1)  // Upwind
	for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	  j = grid->edgep[jptr];
	  
	  nc1 = grid->grad[2*j];
	  nc2 = grid->grad[2*j+1];
	  
	  if(grid->ctop[nc1]>grid->ctop[nc2])
	    kmin = grid->ctop[nc1];
	  else
	    kmin = grid->ctop[nc2];

	  for(k=0;k<kmin;k++)
	    phys->ut[j][k]=0;
	  for(k=kmin;k<grid->Nke[j];k++) {
	    if(phys->u[j][k]>0)
	      nc=nc2;
	    else
	      nc=nc1;
	    phys->ut[j][k]=phys->vc[nc][k]*grid->dzz[nc][k];
	  }
	}
      else // Central
	for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	  j = grid->edgep[jptr];
	  
	  nc1 = grid->grad[2*j];
	  nc2 = grid->grad[2*j+1];
	  
	  if(grid->ctop[nc1]>grid->ctop[nc2])
	    kmin = grid->ctop[nc1];
	  else
	    kmin = grid->ctop[nc2];

	  for(k=0;k<kmin;k++)
	    phys->ut[j][k]=0;

	  for(k=kmin;k<grid->Nke[j];k++) 
	    phys->ut[j][k]=0.5*(phys->vc[nc1][k]*grid->dzz[nc1][k]+
				phys->vc[nc2][k]*grid->dzz[nc2][k]);
	}

      // Now compute the cell-centered source terms and put them into stmp and also
      // add on diffusion of v.
      for(i=0;i<grid->Nc;i++) {
	
	for(k=0;k<grid->Nk[i];k++) 
	  phys->stmp2[i][k]=0;
	
	for(nf=0;nf<NFACES;nf++) {
	  
	  ne = grid->face[i*NFACES+nf];
	  nc = grid->neigh[i*NFACES+nf];

	  if(grid->Nk[nc]<grid->Nk[i])
	    kmax = grid->Nk[nc];
	  else
	    kmax = grid->Nk[i];

	  if(grid->ctop[nc]>grid->ctop[i])
	    kmin = grid->ctop[nc];
	  else
	    kmin = grid->ctop[i];
	  
	  for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	    phys->stmp2[i][k]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	      (grid->Ac[i]*grid->dzz[i][k]);
	  
	  if(nc!=-1) {
	    for(k=kmin;k<kmax;k++)
	      phys->stmp2[i][k]-=prop->nu_H*(phys->vc[nc][k]-phys->vc[i][k])*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	    for(k=kmax;k<grid->Nkc[ne];k++)
	      phys->stmp2[i][k]+=2*prop->nu_H*phys->vc[i][k]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
	  }

	  // Top cell is filled with momentum from neighboring cells
	  for(k=grid->etop[ne];k<=grid->ctop[i];k++) 
	    phys->stmp2[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	      (grid->Ac[i]*grid->dzz[i][grid->ctop[i]]);
	}
      }
      
      // Now do vertical advection
      for(i=0;i<grid->Nc;i++) {
	
	if(prop->nonlinear==1)  // Upwind
	  for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	    a[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
			(phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
	    b[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
			(phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
	  }
	else  // Central
	  for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	    a[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
				   grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1]));
	    b[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
				   grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1]));
	  }

	for(k=grid->ctop[i]+2;k<grid->Nk[i]-1;k++) {
	  phys->stmp[i][k]+=(a[k]-a[k+1])/grid->dzz[i][k];
	  phys->stmp2[i][k]+=(b[k]-b[k+1])/grid->dzz[i][k];
	}
	
	if(grid->ctop[i]!=grid->Nk[i]-1) {
	  if(grid->dzz[i][grid->ctop[i]]/grid->dzz[i][grid->ctop[i]+1]<.1) {
	    phys->stmp[i][grid->ctop[i]+1]-=a[grid->ctop[i]+2]/
	      (grid->dzz[i][grid->ctop[i]]+grid->dzz[i][grid->ctop[i]+1]);
	    phys->stmp2[i][grid->ctop[i]+1]-=b[grid->ctop[i]+2]/
	      (grid->dzz[i][grid->ctop[i]]+grid->dzz[i][grid->ctop[i]+1]);
	    
	    phys->stmp[i][grid->ctop[i]]-=a[grid->ctop[i]+2]/
	      (grid->dzz[i][grid->ctop[i]]+grid->dzz[i][grid->ctop[i]+1]);
	    phys->stmp2[i][grid->ctop[i]]-=b[grid->ctop[i]+2]/
	      (grid->dzz[i][grid->ctop[i]]+grid->dzz[i][grid->ctop[i]+1]);
	  } else {
	    phys->stmp[i][grid->ctop[i]+1]+=(a[grid->ctop[i]+1]-a[grid->ctop[i]+2])/grid->dzz[i][grid->ctop[i]+1];
	    phys->stmp2[i][grid->ctop[i]+1]+=(b[grid->ctop[i]+1]-b[grid->ctop[i]+2])/grid->dzz[i][grid->ctop[i]+1];
	    
	    phys->stmp[i][grid->ctop[i]]-=a[grid->ctop[i]+1]/grid->dzz[i][grid->ctop[i]];
	    phys->stmp2[i][grid->ctop[i]]-=b[grid->ctop[i]+1]/grid->dzz[i][grid->ctop[i]];
	  }
	  
	  // Bottom - advection only comes in through top of cell.
	  phys->stmp[i][grid->Nk[i]-1]+=a[grid->Nk[i]-1]/grid->dzz[i][grid->Nk[i]-1];
	  phys->stmp2[i][grid->Nk[i]-1]+=b[grid->Nk[i]-1]/grid->dzz[i][grid->Nk[i]-1];
	}
      }
      // Check to make sure integrated fluxes are 0 for conservation
      sum=0;
      for(i=0;i<grid->Nc;i++) {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  sum+=grid->Ac[i]*phys->stmp[i][k]*grid->dzz[i][k];
      }
      if(fabs(sum)>1e-5 && WARNING) printf("Warning, not U-momentum conservative!\n");
      
      sum=0;
      for(i=0;i<grid->Nc;i++) {
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  sum+=grid->Ac[i]*phys->stmp2[i][k]*grid->dzz[i][k];
      }
      if(fabs(sum)>1e-5 && WARNING) printf("Warning, not V-momentum conservative!\n");
      
      // Send/recv stmp and stmp2 to account for advective fluxes in ghost cells at
      // interproc boundaries.
      ISendRecvCellData3D(phys->stmp,grid,myproc,comm);
      ISendRecvCellData3D(phys->stmp2,grid,myproc,comm);
      
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr]; 

	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	if(nc1==-1) nc1=nc2;
	if(nc2==-1) nc2=nc1;

	if(grid->ctop[nc1]>grid->ctop[nc2])
	  k0=grid->ctop[nc1]+wetdry_offset;
	else
	  k0=grid->ctop[nc2]+wetdry_offset;

	for(k=k0;k<grid->Nk[nc1];k++) 
	  phys->Cn_U[j][k]-=grid->def[nc1*NFACES+grid->gradf[2*j]]/grid->dg[j]
	    *prop->dt*(phys->stmp[nc1][k]*grid->n1[j]+phys->stmp2[nc1][k]*grid->n2[j]);
	
	for(k=k0;k<grid->Nk[nc2];k++) 
	  phys->Cn_U[j][k]-=grid->def[nc2*NFACES+grid->gradf[2*j+1]]/grid->dg[j]
	    *prop->dt*(phys->stmp[nc2][k]*grid->n1[j]+phys->stmp2[nc2][k]*grid->n2[j]);
      }
    }

    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->utmp[j][k]+=fab*phys->Cn_U[j][k];
    }
  }
}

static void NewCells(gridT *grid, physT *phys, propT *prop) {
 
  int j, jptr, k, nc1, nc2;
  REAL dz;
 
  // Correct the velocity resulting from changes in volume
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
 
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
 
    if(grid->etop[j]<=grid->etopold[j]) {
      dz = 0;
      for(k=grid->etop[j];k<=grid->etopold[j];k++)
        dz+=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
 
      for(k=grid->etop[j];k<=grid->etopold[j];k++)
        phys->u[j][k]=phys->u[j][grid->etopold[j]]/dz*
          (0.5*(grid->dzzold[nc1][k]+grid->dzzold[nc2][k]));
    }
  }
}

static void AdvectVerticalVelocity(gridT *grid, physT *phys, propT *prop,
				     int myproc, int numprocs) {
  int i, iptr, j, jptr, k, ne, nf, nc, nc1, nc2, kmin, wetdry_offset;
  REAL fab, sum;

  if(prop->wetdry)
    wetdry_offset=1;
  else
    wetdry_offset=0;

  if(prop->n==1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
	phys->Cn_W[i][k]=0;
  } else
    fab=1.5;

  // Add on the nonhydrostatic pressure gradient from the previous time
  // step.
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->wtmp[i][k]=phys->w[i][k]+(1-fab)*phys->Cn_W[i][k];
      phys->Cn_W[i][k]=0;
    }
    
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      phys->wtmp[i][k]-=2.0*(1-prop->theta)*prop->dt/(grid->dzz[i][k-1]+grid->dzz[i][k])*
	(phys->q[i][k-1]-phys->q[i][k]);
    phys->wtmp[i][grid->ctop[i]]+=2.0*(1-prop->theta)*prop->dt/grid->dzz[i][grid->ctop[i]]*
      phys->q[i][grid->ctop[i]];
  }

  // Compute Eulerian advection (nonlinear!=0)
  if(prop->nonlinear) {
    // Compute the w-component fluxes at the faces

    if(prop->nonlinear==1) // Upwind
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	
	for(k=kmin;k<grid->Nke[j];k++) {
	  if(phys->u[j][k]>0)
	    nc = nc2;
	  else
	    nc = nc1;

	  phys->ut[j][k]=0.5*(phys->w[nc][k]+phys->w[nc][k+1])*grid->dzz[nc][k];
	}
      }
    else // Central 
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	j = grid->edgep[jptr];
	
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(grid->ctop[nc1]>grid->ctop[nc2])
	  kmin = grid->ctop[nc1];
	else
	  kmin = grid->ctop[nc2];
	
	for(k=0;k<kmin;k++)
	  phys->ut[j][k]=0;
	for(k=kmin;k<grid->Nke[j];k++) 
	  phys->ut[j][k]=0.5*(0.5*(phys->w[nc1][k]+phys->w[nc1][k+1])*grid->dzz[nc1][k]+
			      0.5*(phys->w[nc2][k]+phys->w[nc2][k+1])*grid->dzz[nc2][k]);
      }
    
    for(i=0;i<grid->Nc;i++) {
      
      for(k=0;k<grid->Nk[i];k++) 
	phys->stmp[i][k]=0;
      
      for(nf=0;nf<NFACES;nf++) {
	
	ne = grid->face[i*NFACES+nf];
	
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	  phys->stmp[i][k]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][k]);
	
	// Top cell is filled with momentum from neighboring cells
	for(k=grid->etop[ne];k<=grid->ctop[i];k++) 
	  phys->stmp[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*NFACES+nf]/
	    (grid->Ac[i]*grid->dzz[i][grid->ctop[i]]);
      }
      
      // Vertical advection
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	phys->stmp[i][k]+=(pow(phys->w[i][k],2)-pow(phys->w[i][k+1],2))/grid->dzz[i][k];
      }
      // Top cell
      phys->stmp[i][grid->ctop[i]]-=pow(phys->w[i][grid->ctop[i]+1],2)/grid->dzz[i][grid->ctop[i]];
    }
    
    // Check to make sure integrated fluxes are 0 for conservation
    sum=0;
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	sum+=grid->Ac[i]*phys->stmp[i][k]*grid->dzz[i][k];
    }
    if(fabs(sum)>1e-5 && WARNING) printf("Warning, not W-momentum conservative!\n");
    
    //Now use the cell-centered advection terms to update the advection at the faces
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr]; 
      
      for(k=grid->ctop[i]+1+wetdry_offset;k<grid->Nk[i];k++) 
	phys->Cn_W[i][k]-=prop->dt*(grid->dzz[i][k-1]*phys->stmp[i][k-1]+grid->dzz[i][k]*phys->stmp[i][k])/
	  (grid->dzz[i][k-1]+grid->dzz[i][k]);
      
      // Top flux advection consists only of top cell
      if(!prop->wetdry) {
	k=grid->ctop[i];
	phys->Cn_W[i][k]-=prop->dt*phys->stmp[i][k];
      }
    }
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->w[i][k]=phys->wtmp[i][k]+fab*phys->Cn_W[i][k];
  }
}

static void Corrector(REAL **qc, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr, j, jptr, k;

  // Correct the horizontal velocity only if this is not a boundary point.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 
    if(phys->D[j]!=0 && grid->etop[j]<grid->Nke[j]-1)
      for(k=grid->etop[j];k<grid->Nke[j];k++)
	phys->u[j][k]-=prop->theta*prop->dt/grid->dg[j]*
	  (qc[grid->grad[2*j]][k]-qc[grid->grad[2*j+1]][k]);
  }
  
  // Now at the open boundaries
  /*
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr]; 

    for(k=grid->etop[j];k<grid->Nke[j];k++)
      phys->u[j][k]-=2.0*prop->thetaB*prop->dt/grid->dg[j]*
      	qc[grid->grad[2*j]][k];
  }
  */

  // Update the pressure since qc is a pressure correction
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    if(grid->ctop[i]<grid->Nk[i]-1)
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	phys->q[i][k]+=qc[i][k];
  }

  // Send q to the boundary cells now that it has been updated with qc
  ISendRecvCellData3D(phys->q,grid,myproc,comm);

  /*
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      grid->df[ne]=phys->D[ne];
    }
  }
  */
}

static void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs) {
  
  int i, iptr, j, jptr, k, nf, ne, nc1, nc2;
  REAL *ap=phys->a, *am=phys->b;

  /*
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr]; 

    for(k=grid->etop[j];k<grid->Nke[j];k++)
      phys->u[j][k]*=prop->theta/prop->thetaB;
  }
  */

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k] = grid->Ac[i]*(phys->w[i][k]-phys->w[i][k+1]);

    for(nf=0;nf<NFACES;nf++) {

      ne = grid->face[i*NFACES+nf];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	ap[k] = 0.5*(phys->u[ne][k]+fabs(phys->u[ne][k]));
	am[k] = 0.5*(phys->u[ne][k]-fabs(phys->u[ne][k]));
      }

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	src[i][k]+=(grid->dzz[nc2][k]*ap[k]+grid->dzz[nc1][k]*am[k])*
	  grid->normal[i*NFACES+nf]*grid->df[ne];
    }

    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k]/=(prop->theta*prop->dt);
  }

  /*
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr]; 

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k]/=(prop->theta/prop->thetaB);
  }
  */

  // D[j] is used in OperatorQ, and it must be zero to ensure no gradient
  // at the hydrostatic faces.
  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=grid->df[j]/grid->dg[j];
  }

  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(nf=0;nf<NFACES;nf++)
      phys->D[grid->face[i*NFACES+nf]]=0;
  }
}

static void CGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr, k, n, niters;

  REAL **x, **r, **p, **z, mu, nu, eps, eps0;

  z = phys->uc;
  x = q;
  r = phys->vc;
  p = src;

  // Preconditioner
  if(prop->qprecond) {
    ConditionQ(c,grid,phys,prop,myproc,comm);
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	p[i][k]/=c[i][k];
	x[i][k]*=0*c[i][k];
      }
    }
  }

  niters = prop->qmaxiters;

  // Create coefficients for OperatorQ
  QCoefficients(phys->wtmp,phys->qtmp,c,grid,phys,prop);

  // Initialization for CG
  if(prop->qprecond) OperatorQC(phys->wtmp,phys->qtmp,x,z,c,grid,phys,prop);
  else OperatorQ(phys->wtmp,x,z,c,grid,phys,prop);
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      r[i][k] = p[i][k]-z[i][k];
      p[i][k] = r[i][k];
    }
  }
  eps = eps0 = InnerProduct3(p,p,grid,myproc,numprocs,comm);

  for(n=0;n<niters && eps0!=0 && !IsNan(eps0);n++) {

    ISendRecvCellData3D(p,grid,myproc,comm);
    if(prop->qprecond) OperatorQC(phys->wtmp,phys->qtmp,p,z,c,grid,phys,prop);
    else OperatorQ(phys->wtmp,p,z,c,grid,phys,prop);

    mu = 1/eps;
    nu = eps/InnerProduct3(p,z,grid,myproc,numprocs,comm);

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	x[i][k] += nu*p[i][k];
	r[i][k] -= nu*z[i][k];
      }
    }
    eps = InnerProduct3(r,r,grid,myproc,numprocs,comm);
    mu*=eps;

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	p[i][k] = r[i][k] + mu*p[i][k];
    }
    //    printf("%e\n",sqrt(eps/eps0));
    if(VERBOSE>3) printf("CGSolve Pressure Iteration: %d, resid=%e, proc=%d\n",n,sqrt(eps/eps0),myproc);
    if(sqrt(eps/eps0)<prop->qepsilon) 
      break;
  }

  if(n==niters && myproc==0 && WARNING) 
    printf("Warning... Time step %d, Pressure iteration not converging after %d steps! RES=%e\n",prop->n,n,sqrt(eps/eps0));
  else
    if(VERBOSE>2 && myproc==0) printf("Time step %d, CGSolve pressure converged after %d iterations, res=%e\n",prop->n,n,sqrt(eps/eps0));

  // Rescale the preconditioned solution 
  if(prop->qprecond) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	x[i][k]/=c[i][k];
    }
  }

  ISendRecvCellData3D(x,grid,myproc,comm);
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
  int i, iptr, j, jptr, ne, nf, nf1, normal, nc1, nc2, k;
  REAL sum, dt=prop->dt, theta=prop->theta, fluxheight;
  REAL *a, *b, *c, *d, *e1, **E;

  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;
  e1 = phys->ap;
  E = phys->ut;

  // Set D[j] = 0 
  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i]+1;k++) 
      phys->wtmp2[i][k]=phys->w[i][k];

  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=0;
    for(k=0;k<grid->Nke[j];k++)
      phys->utmp2[j][k]=phys->u[j][k];
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

  // phys->u contains the velocity specified at the open boundaries
  // It is also the velocity at time step n.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]=phys->u[j][k];
  }

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // Add the explicit part of the free-surface to U**.
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]-=GRAV*(1-theta)*dt*(phys->h[nc1]-phys->h[nc2])/grid->dg[j];

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
	printf("Error in function BarotropicPredictor at j=%d k=%d (U***=nan)\n",j,k);
	exit(1);
      }

  // So far we have U*** and D.  Now we need to create h* in htmp.   This
  // will comprise the source term for the free-surface solver.  Before we
  // do this we need to set the new velocity at the open boundary faces and
  // place them into utmp.  These need the velocity from the old time step uold
  // As well as the current vectors.
  OpenBoundaryFluxes(phys->q,phys->utmp,phys->utmp2,grid,phys,prop);
  
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

	sum+=((1-theta)*phys->u[ne][k]+theta*phys->utmp[ne][k])*
	  fluxheight*grid->df[ne]*normal;
      }
    }
    phys->htmp[i]=phys->h[i]-dt/grid->Ac[i]*sum;
  }

  // Set the free-surface values at the boundary.  These are set for h and not for
  // htmp since the boundary values for h will be used in the cg solver.
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    //phys->htmp[i] = phys->h[i]+prop->dt*prop->amp*prop->omega*cos(prop->omega*prop->rtime);
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
  if(prop->cgsolver==0)
    GSSolve(grid,phys,prop,myproc,numprocs,comm);
  else if(prop->cgsolver==1)
    CGSolve(grid,phys,prop,myproc,numprocs,comm);

  // Add on the implicit barotropic term to obtain the hydrostatic horizontal velocity field.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k]=phys->utmp[j][k]-GRAV*theta*dt*E[j][k]*
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
	phys->u[j][k]=phys->utmp[j][k]-GRAV*theta*dt*
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

  // Set the flux values at boundary cells if specified (marker=4)
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k]=prop->flux*cos(prop->omega*prop->rtime);
  }

  // Set the flux values at the open boundary (marker=2).  These
  // were set to utmp previously in OpenBoundaryFluxes.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k] = phys->utmp[j][k];
  }

  // Now set the fluxes at the free-surface boundary by assuming dw/dz = 0
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      if(grid->mark[ne]==3) {
	for(k=grid->etop[ne];k<grid->Nke[ne];k++) {
	  phys->u[ne][k] = 0;
	  sum=0;
	  for(nf1=0;nf1<NFACES;nf1++)
	    sum+=phys->u[grid->face[i*NFACES+nf1]][k]*grid->df[grid->face[i*NFACES+nf1]]*grid->normal[i*NFACES+nf1];
	  phys->u[ne][k]=-sum/grid->df[ne]/grid->normal[i*NFACES+nf];
	}
      }
    }
  }
}

static void CGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, n, niters;
  REAL *x, *r, *D, *p, *z, mu, nu, eps, eps0;

  z = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"CGSolve");
  x = phys->h;
  r = phys->hold;
  D = phys->D;
  p = phys->htmp;

  niters = prop->maxiters;

  // Initialization for CG

  // First need to adjust the source term for the cells just inside
  // the forced boundary.  This is done only if there are forced boundary
  // points, i.e. if celldist[2]!=celldist[1].
  if(grid->celldist[2]!=grid->celldist[1]) {
    for(i=0;i<grid->Nc;i++) {
      z[i] = 0;
    }
    
    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      
      z[i] = x[i];
    }
    OperatorH(z,r,grid,phys,prop);
    
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      p[i]-=r[i];
    }

    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      
      p[i]=0;
    }  
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    r[i] = p[i];
    x[i] = 0;
  }
  eps0 = eps = InnerProduct(r,r,grid,myproc,numprocs,comm);

  for(n=0;n<niters && eps0!=0;n++) {

    ISendRecvCellData2D(p,grid,myproc,comm);
    OperatorH(p,z,grid,phys,prop);

    mu = 1/eps;
    nu = eps/InnerProduct(p,z,grid,myproc,numprocs,comm);

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      x[i] += nu*p[i];
      r[i] -= nu*z[i];
    }
    eps = InnerProduct(r,r,grid,myproc,numprocs,comm);
    mu*=eps;

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      p[i] = r[i] + mu*p[i];
    }

    if(VERBOSE>3) printf("CGSolve Iteration: %d, resid=%e, proc=%d\n",n,sqrt(eps/eps0),myproc);
    if(sqrt(eps/eps0)<prop->epsilon) 
      break;
  }

  if(n==niters && myproc==0 && WARNING) 
    printf("Warning... Time step %d, Iteration not converging after %d steps! RES=%e\n",prop->n,n,sqrt(eps/eps0));
  else
    if(VERBOSE>2 && myproc==0) printf("Time step %d, CGSolve converged after %d iterations, res=%e\n",prop->n,n,sqrt(eps/eps0));

  ISendRecvCellData2D(x,grid,myproc,comm);
  SunFree(z,grid->Nc*sizeof(REAL),"CGSolve");
}

static REAL InnerProduct(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, MPI_Comm comm) {
  
  int i, iptr;
  REAL sum, mysum=0;
  
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    mysum+=x[i]*y[i];
  }
  MPI_Reduce(&mysum,&(sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&sum,1,MPI_DOUBLE,0,comm);

  return sum;
}  

static REAL InnerProduct3(REAL **x, REAL **y, gridT *grid, int myproc, int numprocs, MPI_Comm comm) {
  
  int i, k, iptr;
  REAL sum, mysum=0;
  
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      mysum+=x[i][k]*y[i][k];
  }
  MPI_Reduce(&mysum,&(sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&sum,1,MPI_DOUBLE,0,comm);

  return sum;
}  

static void OperatorH(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop) {
  
  int i, iptr, ne, nf;
  REAL tmp = GRAV*pow(prop->theta*prop->dt,2);

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      y[i] = x[i];
      for(nf=0;nf<NFACES;nf++) 
	if(grid->neigh[i*NFACES+nf]!=-1) {
	  ne = grid->face[i*NFACES+nf];

	  y[i]+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
	    (x[i]-x[grid->neigh[i*NFACES+nf]])/grid->Ac[i];
	}
  }
}

static void OperatorQC(REAL **coef, REAL **fcoef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop) {

  int i, iptr, k, ne, nf, nc, kmin, kmax;
  REAL *a = phys->a;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	y[i][k]=-x[i][k];

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {

	    ne = grid->face[i*NFACES+nf];

	    if(grid->ctop[nc]>grid->ctop[i])
	      kmin = grid->ctop[nc];
	    else
	      kmin = grid->ctop[i];
	    
	    for(k=kmin;k<grid->Nke[ne];k++) 
	      y[i][k]+=x[nc][k]*fcoef[i*NFACES+nf][k];
	  }

      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	y[i][k]+=coef[i][k]*x[i][k-1]+coef[i][k+1]*x[i][k+1];

      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      y[i][k]+=coef[i][k+1]*x[i][k+1];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      y[i][k]+=coef[i][k]*x[i][k-1];
  }
  /*
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j = grid->edgep[jptr];

      i = grid->grad[2*j];

      ne = grid->face[i*NFACES+nf];

      for(k=grid->ctop[i];k<grid->Nke[j];k++) 
	y[i][k]-=2.0*x[i][k]*grid->dzz[i][k]*phys->D[j];
  }
  */
}

static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, physT *phys, propT *prop) {

  int i, iptr, k, kmin, nf, nc, ne;

  if(prop->qprecond) 
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]]/c[i][grid->ctop[i]];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
	coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1])/(c[i][k]*c[i][k-1]);

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {

	    ne = grid->face[i*NFACES+nf];

	    if(grid->ctop[nc]>grid->ctop[i])
	      kmin = grid->ctop[nc];
	    else
	      kmin = grid->ctop[i];
	    
	    for(k=kmin;k<grid->Nke[ne];k++) 
	      fcoef[i*NFACES+nf][k]=grid->dzz[i][k]*phys->D[ne]/(c[i][k]*c[nc][k]);
	}
    }
  else
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
	coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
  }
}

static void OperatorQ(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop) {

  int i, iptr, k, ne, nf, nc, kmin, kmax;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	y[i][k]=0;

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {
	  
	    ne = grid->face[i*NFACES+nf];

	    if(grid->ctop[nc]>grid->ctop[i])
	      kmin = grid->ctop[nc];
	    else
	      kmin = grid->ctop[i];

	    for(k=kmin;k<grid->Nke[ne];k++) 
	      y[i][k]+=(x[nc][k]-x[i][k])*grid->dzz[i][k]*phys->D[ne];
	  }

      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	y[i][k]+=coef[i][k]*x[i][k-1]-(coef[i][k]+coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      y[i][k]+=(-2*coef[i][k]-coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      y[i][k]+=coef[i][k]*x[i][k-1]-coef[i][k]*x[i][k];
  }
  /*
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j = grid->edgep[jptr];

      i = grid->grad[2*j];

      ne = grid->face[i*NFACES+nf];

      for(k=grid->ctop[i];k<grid->Nke[j];k++) 
	y[i][k]-=2.0*x[i][k]*grid->dzz[i][k]*phys->D[j];
  }
  */
}

static void GuessQ(REAL **q, REAL **wold, REAL **w, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {
  
  int i, iptr, k;
  REAL qerror;

  HydroW(phys->w,grid,phys,prop);

  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
  
    q[i][grid->ctop[i]]=0;
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      q[i][k]=q[i][k-1]-(grid->dzz[i][k]+grid->dzz[i][k-1])/(2.0*prop->dt*prop->theta)*
	(wold[i][k]-w[i][k]);
    }
    
    // Adjust so that q[i][grid->ctop[i]]+q[i][grid->ctop[i]+1]=0;
    qerror = 0.5*(q[i][grid->ctop[i]]+q[i][grid->ctop[i]+1]);
    q[i][grid->ctop[i]]-=qerror;
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      q[i][k]-=qerror;

    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      q[i][k]=-q[i][k];
  }
}

static void ConditionQ(REAL **x, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {

  int i, iptr, k, ne, nf, nc, kmin, kmax, warn=0;
  REAL *a = phys->a;

  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	x[i][k]=0;

      for(nf=0;nf<NFACES;nf++) 
	if((nc=grid->neigh[i*NFACES+nf])!=-1) {

	  ne = grid->face[i*NFACES+nf];
	  
	  if(grid->Nk[nc]<grid->Nk[i])
	    kmax = grid->Nk[nc];
	  else
	    kmax = grid->Nk[i];
	  
	  if(grid->ctop[nc]>grid->ctop[i])
	    kmin = grid->ctop[nc];
	  else
	    kmin = grid->ctop[i];
	  
	  for(k=kmin;k<kmax;k++) 
	    x[i][k]+=grid->dzz[i][k]*phys->D[ne];
	}
      
      a[grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
	a[k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1]);

      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	x[i][k]+=(a[k]+a[k+1]);

      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      x[i][k]+=2*a[k]+a[k+1];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      x[i][k]+=a[k];
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	if(x[i][k]<=0) {
	  x[i][k]=1;
	  warn=1;
	}
	x[i][k]=sqrt(x[i][k]);
      }
  }
  if(WARNING && warn) printf("Warning...invalid preconditioner!\n");
  /*
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j = grid->edgep[jptr];

      i = grid->grad[2*j];

      ne = grid->face[i*NFACES+nf];

      for(k=grid->ctop[i];k<grid->Nke[j];k++) 
	x[i][k]+=2.0*grid->dzz[i][k]*phys->D[j];
  }
  */
  // Send the preconditioner to the neighboring processors.
  ISendRecvCellData3D(x,grid,myproc,comm);
}

static void GSSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, nf, ne, n, niters, *N;
  REAL *h, *hold, *D, *hsrc, myresid, resid, residold, tmp, relax, myNsrc, Nsrc, coef;
  //  FILE *fid=fopen("../../data/cg.m","w");

  h = phys->h;
  hold = phys->hold;
  D = phys->D;
  hsrc = phys->htmp;
  N = grid->normal;

  tmp = GRAV*pow(prop->theta*prop->dt,2);

  ISendRecvCellData2D(h,grid,myproc,comm);

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

    //printf("Proc %d: %e\n",myproc,resid);

    ISendRecvCellData2D(h,grid,myproc,comm);
    MPI_Barrier(comm);

    if(fabs(resid)<prop->epsilon)
      break;
  }
  if(n==niters && myproc==0 && WARNING) 
    printf("Warning... Iteration not converging after %d steps! RES=%e\n",n,resid);
  
  for(i=0;i<grid->Nc;i++)
    if(h[i]!=h[i]) 
      printf("NaN h[%d] in cgsolve!\n",i);

}

static void UpdateScalars(gridT *grid, physT *phys, propT *prop, REAL **scalold, REAL **scal, REAL **Cn, REAL theta)
{
  int i, iptr, k, nf, ktop;
  int Nc=grid->Nc, normal, nc1, nc2, ne;
  REAL df, Ac, dt=prop->dt, fab, *a, *b, *c, *d, *ap, *am, dznew, mass;

  ap = phys->ap;
  am = phys->am;
  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;

  if(prop->n==1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
	Cn[i][k]=0;
  } else
    fab=1.5;

  for(i=0;i<Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=scalold[i][k];

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

    // These are the advective components of the tridiagonal
    // at the new time step.
    for(k=0;k<grid->Nk[i]+1;k++) {
      ap[k] = 0.5*(phys->w[i][k]+fabs(phys->w[i][k]));
      am[k] = 0.5*(phys->w[i][k]-fabs(phys->w[i][k]));
    }
    for(k=ktop+1;k<grid->Nk[i];k++) {
      a[k-ktop]=theta*dt*am[k];
      b[k-ktop]=grid->dzz[i][k]+theta*dt*(ap[k]-am[k+1]);
      c[k-ktop]=-theta*dt*ap[k+1];
    }

    // Top cell
    a[0]=0;
    b[0]=dznew-theta*dt*am[ktop+1];
    c[0]=-theta*dt*ap[ktop+1];

    // Bottom cell no-flux boundary condition
    b[(grid->Nk[i]-1)-ktop]+=c[(grid->Nk[i]-1)-ktop];

    // Explicit part into source term d[]
    for(k=ktop+1;k<grid->Nk[i];k++) 
      d[k-ktop]=grid->dzzold[i][k]*phys->stmp[i][k];

    d[0]=0;
    if(grid->ctopold[i]<=grid->ctop[i])
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
	d[0]+=grid->dzzold[i][k]*phys->stmp[i][k];
    else
      d[0]=grid->dzzold[i][ktop]*phys->stmp[i][ktop];

    // These are the advective components of the tridiagonal
    // at the old time step.
    for(k=0;k<grid->Nk[i]+1;k++) {
      ap[k] = 0.5*(phys->wtmp2[i][k]+fabs(phys->wtmp2[i][k]));
      am[k] = 0.5*(phys->wtmp2[i][k]-fabs(phys->wtmp2[i][k]));
    }
    for(k=ktop+1;k<grid->Nk[i]-1;k++) 
      d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
			       (ap[k]-am[k+1])*phys->stmp[i][k]-
			       ap[k+1]*phys->stmp[i][k+1]);

    if(ktop<grid->Nk[i]-1) {
      //Flux through bottom of top cell
      k=ktop;
      d[0]-=(1-theta)*dt*(-am[k+1]*phys->stmp[i][k]-
			   ap[k+1]*phys->stmp[i][k+1]);

      // Through top of bottom cell
      k=grid->Nk[i]-1;
      d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
			       ap[k]*phys->stmp[i][k]);
    }

    /* Comment this out for AB2 */
    /*
    for(k=0;k<grid->Nk[i];k++)
      Cn[i][k]=0;
    fab=1;
    */

    // First add on the source term from the previous time step.
    for(k=0;k<grid->Nk[i]-ktop;k++)
      d[k]+=(1-fab)*Cn[i][k];

    for(k=0;k<grid->Nk[i];k++)
      Cn[i][k]=0;

    // Now create the source term for the current time step
    for(nf=0;nf<NFACES;nf++) {
    
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      for(k=0;k<grid->Nkc[ne];k++) 
	ap[k] = dt*df*normal/Ac*(0.5*(phys->utmp2[ne][k]+fabs(phys->utmp2[ne][k]))*
				 phys->stmp[nc2][k]*grid->dzzold[nc2][k]
	  +0.5*(phys->utmp2[ne][k]-fabs(phys->utmp2[ne][k]))*
				 phys->stmp[nc1][k]*grid->dzzold[nc1][k]);

      for(k=ktop+1;k<grid->Nk[i];k++) 
      	Cn[i][k-ktop]-=ap[k];

      for(k=0;k<=ktop;k++)
	Cn[i][0]-=ap[k];
    }

    // Add on the source from the current time step to the rhs.
    for(k=0;k<grid->Nk[i]-ktop;k++)
      d[k]+=fab*Cn[i][k];

    if(grid->Nk[i]-ktop>1) 
      TriSolve(a,b,c,d,&(scal[i][ktop]),grid->Nk[i]-ktop);
    else if(b[0]!=0)
      scal[i][ktop]=d[0]/b[0];

    if(dznew!=0) {
      mass = scal[i][ktop]*dznew;
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++)
	scal[i][k]=mass/dznew;
    }

    for(k=0;k<grid->ctop[i];k++)
      scal[i][k]=0;
  }

  /*
   * Need to modify the boundary cells here.  For now they
   * are left unmodified.
   *
   */
  /*
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];

    for(k=grid->ctopold[i];k<grid->Nk[i];k++) {
      tmp=scal[i][k];
      scal[i][k]=0;
      km=0;
      for(nf=0;nf<NFACES;nf++) {
	if(grid->neigh[i*NFACES+nf]!=-1) 
	  if(k<=grid->Nk[grid->neigh[i*NFACES+nf]] &&
	     k>=grid->ctop[grid->neigh[i*NFACES+nf]]) {
	    scal[i][k]+=scal[grid->neigh[i*NFACES+nf]][k];
	    km++;
	  }
      }
      if(km>0)
	scal[i][k]/=km;
      else
	scal[i][k]=tmp;
    }
  }
  */
}

static void HydroW(REAL **w, gridT *grid, physT *phys, propT *prop)
{
  int i, k, nf, iptr, ne, nc1, nc2;
  REAL ap, am, d;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->ctop[i];k++) 
      w[i][k] = 0;

    w[i][grid->Nk[i]] = 0;
    d=grid->dzz[i][grid->Nk[i]-1];
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      w[i][k] = w[i][k+1];
      for(nf=0;nf<NFACES;nf++) {
	ne = grid->face[i*NFACES+nf];
	nc1 = grid->grad[2*ne];
	nc2 = grid->grad[2*ne+1];
	if(nc1==-1) nc1=nc2;
	if(nc2==-1) nc2=nc1;

	ap = 0.5*(phys->u[ne][k]+fabs(phys->u[ne][k]));
	am = 0.5*(phys->u[ne][k]-fabs(phys->u[ne][k]));

	w[i][k]-=(ap*grid->dzz[nc2][k]+am*grid->dzz[nc1][k])*
	  grid->df[ne]*grid->normal[i*NFACES+nf]/grid->Ac[i];
      }
    }
  }
}

static void WtoVerticalFace(gridT *grid, physT *phys)
{
  int j, k, nc1, nc2;

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
  int i, iptr, k;
  REAL mass, volume, volh, height, Ep;

  if(myproc==0) phys->mass=0;
  if(myproc==0) phys->volume=0;
  if(myproc==0) phys->Ep=0;

  mass=0;
  volume=0;
  volh=0;
  Ep=0;

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
  int i, k, icu, kcu, icw, kcw, Nc=grid->Nc, Ne=grid->Ne, ih, is, ks, iu, ku;
  int uflag=1, sflag=1, hflag=1, myalldone, alldone;
  REAL C, CmaxU, CmaxW;

  icu=kcu=icw=kcw=ih=is=ks=iu=ku=0;

  for(i=0;i<Nc;i++) 
    if(phys->h[i]!=phys->h[i]) {
	hflag=0;
	ih=i;
	break;
    }

  for(i=0;i<Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      if(phys->s[i][k]!=phys->s[i][k]) {
	sflag=0;
	is=i;
	ks=k;
	break;
      }
    if(!sflag)
      break;
  }

  for(i=0;i<Ne;i++) {
    for(k=0;k<grid->Nke[i];k++)
      if(phys->u[i][k]!=phys->u[i][k]) {
	uflag=0;
	iu=i;
	ku=k;
	break;
      }
    if(!uflag)
      break;
  }

  CmaxU=0;
  for(i=0;i<Ne;i++) 
    for(k=grid->etop[i]+1;k<grid->Nke[i];k++) {
      C = fabs(phys->u[i][k])*prop->dt/grid->dg[i];
      if(C>CmaxU) {
	icu = i;
	kcu = k;
	CmaxU = C;
      }
    }

  CmaxW=0;
  for(i=0;i<Nc;i++) 
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      C = 0.5*fabs(phys->w[i][k]+phys->w[i][k+1])*prop->dt/grid->dzz[i][k];
      if(C>CmaxW) {
	icw = i;
	kcw = k;
	CmaxW = C;
      }
    }

  myalldone=0;
  if(!uflag || !sflag || !hflag || CmaxU>prop->Cmax || CmaxW>prop->Cmax) {
    printf("Time step %d: Processor %d, Run is blowing up!\n",prop->n,myproc);
    
    if(CmaxU>prop->Cmax)
      printf("Courant number problems at (%d,%d), Umax=%f, dx=%f CmaxU=%.2f > %.2f\n",
	     icu,kcu,phys->u[icu][kcu],grid->dg[icu],CmaxU,prop->Cmax);
    else if(CmaxW>prop->Cmax)
      printf("Courant number problems at (%d,%d), Wmax=%f, dz=%f CmaxW=%.2f > %.2f\n",
	     icw,kcw,0.5*(phys->w[icw][kcw]+phys->w[icw][kcw+1]),grid->dzz[icw][kcw],CmaxW,prop->Cmax);
    else
      printf("Courant number is okay: CmaxU=%.2f,CmaxW=%.2f < %.2f\n",CmaxU,CmaxW,prop->Cmax);

    if(!uflag)
      printf("U is divergent at (%d,%d)\n",iu,ku);
    else
      printf("U is okay.\n");
    if(!sflag) 
      printf("Scalar is divergent at (%d,%d).\n",is,ks);
    else
      printf("Scalar is okay.\n");
    if(!hflag)
      printf("Free-surface is divergent at (%d)\n",ih);
    else
      printf("Free-surface is okay.\n");
    
    OutputData(grid,phys,prop,myproc,numprocs,1,comm);
    myalldone=1;
  }

  MPI_Reduce(&myalldone,&alldone,1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&alldone,1,MPI_INT,0,comm);

  if(alldone) {
    MPI_Finalize();
    exit(0);
  }
}

static void ComputeVelocityVector(REAL **u, REAL **uc, REAL **vc, gridT *grid) {
  
  int k, n, ne, nf;
  for(n=0;n<grid->Nc;n++) {
    for(k=0;k<grid->Nk[n];k++) {
      uc[n][k]=0;
      vc[n][k]=0;
    }
    for(k=grid->ctop[n];k<grid->Nk[n];k++) {
      for(nf=0;nf<NFACES;nf++) {
	ne = grid->face[n*NFACES+nf];
	uc[n][k]+=u[ne][k]*grid->n1[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
	vc[n][k]+=u[ne][k]*grid->n2[ne]*grid->def[n*NFACES+nf]*grid->df[ne];
      }
      uc[n][k]/=grid->Ac[n];
      vc[n][k]/=grid->Ac[n];
      ne = grid->face[n*NFACES+1];
    }
  }
}

static void AdjustBoundaryFluxes(REAL **u, REAL **uc, REAL **vc, gridT *grid) {
  
  int nc1, nc2, j, jptr, k, kstart, kend, iin, iout;
  REAL db, z;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc1!=-1 && nc2!=-1) {
      if(grid->Nk[nc1]>grid->Nk[nc2]) {
	iin = nc1;
	iout = nc2;
	kstart = grid->Nk[nc2];
	kend = grid->Nk[nc1];
      } else {
	iin = nc2;
	iout = nc1;
	kstart = grid->Nk[nc1];
	kend = grid->Nk[nc2];
      }

      kstart++;
      if(kstart!=kend) {
	z=grid->dzz[iout][kstart]/2;
	
	for(k=kstart;k<kend;k++) {
	  db = grid->dg[j]*(1-z/fabs(grid->dv[nc1]-grid->dv[nc2]));
	  if(db==0) {
	    printf("db = 0!\n");
	    exit(0);
	  }
	  u[j][k] = (grid->n1[j]*uc[iin][k]+grid->n2[j]*vc[iin][k])*(1-0.5*grid->dg[j]/db);

	  if(k<kend-1) z+=0.5*(grid->dzz[iout][k]+grid->dzz[iout][k+1]);
	}
      }
    }
  }
}
      
static void ComputeNormalVelocity(REAL **u, REAL **uc, REAL **vc, gridT *grid) {
  
  int j, jptr, k;
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      u[j][k]=0.5*(grid->n1[j]*(uc[grid->grad[2*j]][k]+uc[grid->grad[2*j+1]][k])+
		   grid->n2[j]*(vc[grid->grad[2*j]][k]+vc[grid->grad[2*j+1]][k]));
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
		int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int i, j, k, nwritten;
  REAL *tmp = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"OutputData");

  if(!(prop->n%prop->ntconserve)) {
    ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);
    if(myproc==0)
      fprintf(prop->ConserveFID,"%e %e %e %e %e %e %e %e\n",prop->rtime,phys->mass,phys->volume,
	      phys->Ep-phys->Ep0,phys->Eflux1,phys->Eflux2,phys->Eflux3,phys->Eflux4);
  }

  if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) {

    if(myproc==0 && VERBOSE>1) printf("Outputting data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);

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
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++) 
	  fprintf(prop->HorizontalVelocityFID,"%e %e %e\n",
		  phys->uc[i][k],phys->vc[i][k],0.5*(phys->w[i][k]+phys->w[i][k+1]));
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->HorizontalVelocityFID,"0.0 0.0 0.0\n");
      }
    else 
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i]) 
	    tmp[i]=phys->uc[i][k];
	  else
	    tmp[i]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Nc,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting Horizontal Velocity data!\n");
	  exit(EXIT_WRITING);
	}
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    tmp[i]=phys->vc[i][k];
	  else
	    tmp[i]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Nc,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting Horizontal Velocity data!\n");
	  exit(EXIT_WRITING);
	}
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    tmp[i]=0.5*(phys->w[i][k]+phys->w[i][k+1]);
	  else
	    tmp[i]=0;
	}
	nwritten=fwrite(tmp,sizeof(REAL),grid->Nc,prop->HorizontalVelocityFID);
	if(nwritten!=grid->Nc) {
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
      if(prop->n==1+prop->nstart) {
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
      if(prop->n==1+prop->nstart) {
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
    }
    fflush(prop->SalinityFID);

    if(ASCII) 
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++)
	  fprintf(prop->TemperatureFID,"%e\n",phys->T[i][k]);
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->TemperatureFID,"0.0\n");
      }
    else 
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    phys->htmp[i]=phys->T[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->TemperatureFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting temperature data!\n");
	  exit(EXIT_WRITING);
	}
      }
    fflush(prop->TemperatureFID);

    if(ASCII) 
      for(i=0;i<grid->Nc;i++) {
	for(k=0;k<grid->Nk[i];k++)
	  fprintf(prop->PressureFID,"%e\n",phys->q[i][k]);
	for(k=grid->Nk[i];k<grid->Nkmax;k++)
	  fprintf(prop->PressureFID,"0.0\n");
      }
    else 
      for(k=0;k<grid->Nkmax;k++) {
	for(i=0;i<grid->Nc;i++) {
	  if(k<grid->Nk[i])
	    phys->htmp[i]=phys->q[i][k];
	  else
	    phys->htmp[i]=EMPTY;
	}
	nwritten=fwrite(phys->htmp,sizeof(REAL),grid->Nc,prop->PressureFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting pressure data!\n");
	  exit(EXIT_WRITING);
	}
      }
    fflush(prop->PressureFID);
    
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

  if(prop->n==prop->nsteps+prop->nstart) {
    fclose(prop->FreeSurfaceFID);
    fclose(prop->HorizontalVelocityFID);
    fclose(prop->VerticalVelocityFID);
    fclose(prop->SalinityFID);
    fclose(prop->VerticalGridFID);
    if(myproc==0) fclose(prop->ConserveFID);
  }

  if(prop->n==prop->nsteps+prop->nstart || blowup) {
    if(VERBOSE>1 && myproc==0) printf("Writing to rstore...\n");
    
    nwritten=fwrite(&(prop->n),sizeof(int),1,prop->StoreFID);

    fwrite(phys->h,sizeof(REAL),grid->Nc,prop->StoreFID);
    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(j=0;j<grid->Ne;j++) 
      fwrite(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    for(i=0;i<grid->Nc;i++) 
      fwrite(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StoreFID);
    fclose(prop->StoreFID);
  }

  SunFree(tmp,grid->Ne*sizeof(REAL),"OutputData");
}

void ReadProperties(propT **prop, int myproc)
{
  *prop = (propT *)SunMalloc(sizeof(propT),"ReadProperties");
  
  (*prop)->thetaramptime = MPI_GetValue(DATAFILE,"thetaramptime","ReadProperties",myproc);
  (*prop)->theta = MPI_GetValue(DATAFILE,"theta","ReadProperties",myproc);
  (*prop)->thetaS = MPI_GetValue(DATAFILE,"thetaS","ReadProperties",myproc);
  (*prop)->thetaB = MPI_GetValue(DATAFILE,"thetaB","ReadProperties",myproc);
  (*prop)->beta = MPI_GetValue(DATAFILE,"beta","ReadProperties",myproc);
  (*prop)->gamma = MPI_GetValue(DATAFILE,"gamma","ReadProperties",myproc);
  (*prop)->nu = MPI_GetValue(DATAFILE,"nu","ReadProperties",myproc);
  (*prop)->nu_H = MPI_GetValue(DATAFILE,"nu_H","ReadProperties",myproc);
  (*prop)->tau_T = MPI_GetValue(DATAFILE,"tau_T","ReadProperties",myproc);
  (*prop)->CdT = MPI_GetValue(DATAFILE,"CdT","ReadProperties",myproc);
  (*prop)->CdB = MPI_GetValue(DATAFILE,"CdB","ReadProperties",myproc);
  (*prop)->dt = MPI_GetValue(DATAFILE,"dt","ReadProperties",myproc);
  (*prop)->Cmax = MPI_GetValue(DATAFILE,"Cmax","ReadProperties",myproc);
  (*prop)->nsteps = (int)MPI_GetValue(DATAFILE,"nsteps","ReadProperties",myproc);
  (*prop)->ntout = (int)MPI_GetValue(DATAFILE,"ntout","ReadProperties",myproc);
  (*prop)->ntprog = (int)MPI_GetValue(DATAFILE,"ntprog","ReadProperties",myproc);
  (*prop)->ntconserve = (int)MPI_GetValue(DATAFILE,"ntconserve","ReadProperties",myproc);
  (*prop)->nonhydrostatic = (int)MPI_GetValue(DATAFILE,"nonhydrostatic","ReadProperties",myproc);
  (*prop)->cgsolver = (int)MPI_GetValue(DATAFILE,"cgsolver","ReadProperties",myproc);
  (*prop)->maxiters = (int)MPI_GetValue(DATAFILE,"maxiters","ReadProperties",myproc);
  (*prop)->qmaxiters = (int)MPI_GetValue(DATAFILE,"qmaxiters","ReadProperties",myproc);
  (*prop)->qprecond = (int)MPI_GetValue(DATAFILE,"qprecond","ReadProperties",myproc);
  (*prop)->epsilon = MPI_GetValue(DATAFILE,"epsilon","ReadProperties",myproc);
  (*prop)->qepsilon = MPI_GetValue(DATAFILE,"qepsilon","ReadProperties",myproc);
  (*prop)->relax = MPI_GetValue(DATAFILE,"relax","ReadProperties",myproc);
  (*prop)->amp = MPI_GetValue(DATAFILE,"amp","ReadProperties",myproc);
  (*prop)->omega = MPI_GetValue(DATAFILE,"omega","ReadProperties",myproc);
  (*prop)->timescale = MPI_GetValue(DATAFILE,"timescale","ReadProperties",myproc);
  (*prop)->flux = MPI_GetValue(DATAFILE,"flux","ReadProperties",myproc);
  (*prop)->volcheck = MPI_GetValue(DATAFILE,"volcheck","ReadProperties",myproc);
  (*prop)->masscheck = MPI_GetValue(DATAFILE,"masscheck","ReadProperties",myproc);
  (*prop)->nonlinear = MPI_GetValue(DATAFILE,"nonlinear","ReadProperties",myproc);
  (*prop)->newcells = MPI_GetValue(DATAFILE,"newcells","ReadProperties",myproc);
  (*prop)->wetdry = MPI_GetValue(DATAFILE,"wetdry","ReadProperties",myproc);
  (*prop)->Coriolis_f = MPI_GetValue(DATAFILE,"Coriolis_f","ReadProperties",myproc);
  (*prop)->sponge_distance = MPI_GetValue(DATAFILE,"sponge_distance","ReadProperties",myproc);
  (*prop)->sponge_decay = MPI_GetValue(DATAFILE,"sponge_decay","ReadProperties",myproc);
}

void OpenFiles(propT *prop, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  MPI_GetFile(filename,DATAFILE,"FreeSurfaceFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->FreeSurfaceFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"HorizontalVelocityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->HorizontalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"VerticalVelocityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->VerticalVelocityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"SalinityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->SalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"BGSalinityFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->BGSalinityFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"TemperatureFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->TemperatureFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"PressureFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->PressureFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"VerticalGridFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->VerticalGridFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"StoreFile","OpenFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  prop->StoreFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  if(RESTART) {
    MPI_GetFile(filename,DATAFILE,"StartFile","OpenFiles",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    prop->StartFID = MPI_FOpen(str,"r","OpenFiles",myproc);
  }

  if(myproc==0) {
    MPI_GetFile(filename,DATAFILE,"ConserveFile","OpenFiles",myproc);
    sprintf(str,"%s",filename);
    prop->ConserveFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  }
}

static void Progress(propT *prop, int myproc) 
{
  int progout, prog;
  char filename[BUFFERLENGTH];
  FILE *fid;
  
  MPI_GetFile(filename,DATAFILE,"ProgressFile","Progress",myproc);
  
  if(myproc==0) {
    fid = fopen(filename,"w");
    fprintf(fid,"On %d of %d, t=%.2f (%d%% Complete, %d output)",
	    prop->n,prop->nstart+prop->nsteps,prop->rtime,100*(prop->n-prop->nstart)/prop->nsteps,
	    1+(prop->n-prop->nstart)/prop->ntout);
    fclose(fid);
  }

  if(myproc==0 && prop->ntprog>0 && VERBOSE>0) {
    progout = (int)(prop->nsteps*(double)prop->ntprog/100);
    prog=(int)(100.0*(double)(prop->n-prop->nstart)/(double)prop->nsteps);
    if(progout>0)
      if(!(prop->n%progout))
	printf("%d%% Complete.\n",prog);
  }
}

