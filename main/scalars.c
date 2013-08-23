/*
 * File: scalars.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * ----------------------------------------
 * This file contains the scalar transport function.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "scalars.h"
#include "util.h"
#include "tvd.h"
#include "initialization.h"

#define SMALL_CONSISTENCY 1e-5

REAL smin_value, smax_value;

/*
 * Function: UpdateScalars
 * Usage: UpdateScalars(grid,phys,prop,wnew,scalar,Cn,kappa,kappaH,kappa_tv,theta);
 * ---------------------------------------------------------------------------
 * Update the scalar quantity stored in the array denoted by scal using the
 * theta method for vertical advection and vertical diffusion and Adams-Bashforth
 * for horizontal advection and diffusion.
 *
 * Cn must store the AB terms from time step n-1 for this scalar
 * kappa denotes the vertical scalar diffusivity
 * kappaH denotes the horizontal scalar diffusivity
 * kappa_tv denotes the vertical turbulent scalar diffusivity
 *
 */
void UpdateScalars(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **scal, REAL **boundary_scal, REAL **Cn, 
    REAL kappa, REAL kappaH, REAL **kappa_tv, REAL theta,
    REAL **src1, REAL **src2, REAL *Ftop, REAL *Fbot, int alpha_top, int alpha_bot,
    MPI_Comm comm, int myproc, int checkflag, int TVDscheme) 
{
  int i, iptr, j, jptr, ib, k, nf, ktop;
  int Nc=grid->Nc, normal, nc1, nc2, ne;
  REAL df, dg, Ac, dt=prop->dt, fab, *a, *b, *c, *d, *ap, *am, *bd, *uflux, dznew, mass, *sp, *temp;
  REAL smin, smax, div_local, div_da;
  int k1, k2, kmin, imin, kmax, imax, mincount, maxcount, allmincount, allmaxcount, flag;

  prop->TVD = TVDscheme;
  // These are used mostly debugging to turn on/off vertical and horizontal TVD.
  prop->horiTVD = 1;
  prop->vertTVD = 1;

  ap = phys->ap;
  am = phys->am;
  bd = phys->bp;
  temp = phys->bm;
  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;

  // Never use AB2
  if(1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
        Cn[i][k]=0;
  } else
    fab=1.5;

  for(i=0;i<Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=scal[i][k];

  // Add on boundary fluxes, using stmp2 as the temporary storage
  // variable
  //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      phys->stmp2[i][k]=0;
  }

  if(boundary_scal) {
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];

      ib = grid->grad[2*j];

      // Set the value of stmp2 adjacent to the boundary to the value of the boundary.
      // This will be used to add the boundary flux when stmp2 is used again below.
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++)
        phys->stmp2[ib][k]=boundary_scal[jptr-grid->edgedist[2]][k];
    }
  }

  // Compute the scalar on the vertical faces (for horiz. advection)

  if(prop->TVD && prop->horiTVD)
    HorizontalFaceScalars(grid,phys,prop,scal,boundary_scal,prop->TVD,comm,myproc); 

  //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
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
    if(!(prop->TVD && prop->vertTVD))
      for(k=0;k<grid->Nk[i]+1;k++) {
        ap[k] = 0.5*(wnew[i][k]+fabs(wnew[i][k]));
        am[k] = 0.5*(wnew[i][k]-fabs(wnew[i][k]));
      }
    else  // Compute the ap/am for TVD schemes
      GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
          wnew,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);

    for(k=ktop+1;k<grid->Nk[i];k++) {
      a[k-ktop]=theta*dt*am[k];
      b[k-ktop]=grid->dzz[i][k]+theta*dt*(ap[k]-am[k+1]);
      c[k-ktop]=-theta*dt*ap[k+1];
    }

    // Top cell advection
    a[0]=0;
    b[0]=dznew-theta*dt*am[ktop+1];
    c[0]=-theta*dt*ap[ktop+1];

    // Bottom cell no-flux boundary condition for advection
    b[(grid->Nk[i]-1)-ktop]+=c[(grid->Nk[i]-1)-ktop];

    // Implicit vertical diffusion terms
    for(k=ktop+1;k<grid->Nk[i];k++)
      bd[k]=(2.0*kappa+kappa_tv[i][k-1]+kappa_tv[i][k])/
        (grid->dzz[i][k-1]+grid->dzz[i][k]);

    for(k=ktop+1;k<grid->Nk[i]-1;k++) {
      a[k-ktop]-=theta*dt*bd[k];
      b[k-ktop]+=theta*dt*(bd[k]+bd[k+1]);
      c[k-ktop]-=theta*dt*bd[k+1];
    }
    if(src1)
      for(k=ktop;k<grid->Nk[i];k++)
        b[k-ktop]+=grid->dzz[i][k]*src1[i][k]*theta*dt;

    // Diffusive fluxes only when more than 1 layer
    if(ktop<grid->Nk[i]-1) {
      // Top cell diffusion
      b[0]+=theta*dt*(bd[ktop+1]+2*alpha_top*bd[ktop+1]);
      c[0]-=theta*dt*bd[ktop+1];

      // Bottom cell diffusion
      a[(grid->Nk[i]-1)-ktop]-=theta*dt*bd[grid->Nk[i]-1];
      b[(grid->Nk[i]-1)-ktop]+=theta*dt*(bd[grid->Nk[i]-1]+2*alpha_bot*bd[grid->Nk[i]-1]);
    }

    // Explicit part into source term d[] 
    for(k=ktop+1;k<grid->Nk[i];k++) 
      d[k-ktop]=grid->dzzold[i][k]*phys->stmp[i][k];
    if(src1)
      for(k=ktop+1;k<grid->Nk[i];k++) 
        d[k-ktop]-=src1[i][k]*(1-theta)*dt*grid->dzzold[i][k]*phys->stmp[i][k];

    d[0]=0;
    if(grid->ctopold[i]<=grid->ctop[i]) {
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
        d[0]+=grid->dzzold[i][k]*phys->stmp[i][k];
      if(src1)
        for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
          d[0]-=src1[i][k]*(1-theta)*dt*grid->dzzold[i][k]*phys->stmp[i][k];
    } else {
      d[0]=grid->dzzold[i][ktop]*phys->stmp[i][ktop];
      if(src1)
        d[0]-=src1[i][ktop]*(1-theta)*dt*grid->dzzold[i][ktop]*phys->stmp[i][k];
    }

    // These are the advective components of the tridiagonal
    // that use the new velocity
    if(!(prop->TVD && prop->vertTVD))
      for(k=0;k<grid->Nk[i]+1;k++) {
        ap[k] = 0.5*(phys->wtmp2[i][k]+fabs(phys->wtmp2[i][k]));
        am[k] = 0.5*(phys->wtmp2[i][k]-fabs(phys->wtmp2[i][k]));
      }
    else // Compute the ap/am for TVD schemes
      GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
          phys->wtmp2,grid->dzzold,phys->stmp,i,grid->Nk[i],ktop,prop->dt,prop->TVD);

    // Explicit advection and diffusion
    for(k=ktop+1;k<grid->Nk[i]-1;k++) 
      d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
          (ap[k]-am[k+1])*phys->stmp[i][k]-
          ap[k+1]*phys->stmp[i][k+1])-
        (1-theta)*dt*(bd[k]*phys->stmp[i][k-1]
            -(bd[k]+bd[k+1])*phys->stmp[i][k]
            +bd[k+1]*phys->stmp[i][k+1]);

    if(ktop<grid->Nk[i]-1) {
      //Flux through bottom of top cell
      k=ktop;
      d[0]=d[0]-(1-theta)*dt*(-am[k+1]*phys->stmp[i][k]-
          ap[k+1]*phys->stmp[i][k+1])+
        (1-theta)*dt*(-(2*alpha_top*bd[k+1]+bd[k+1])*phys->stmp[i][k]+
            bd[k+1]*phys->stmp[i][k+1]);
      if(Ftop) d[0]+=dt*(1-alpha_top+2*alpha_top*bd[k+1])*Ftop[i];

      // Through top of bottom cell
      k=grid->Nk[i]-1;
      d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
          ap[k]*phys->stmp[i][k])-
        (1-theta)*dt*(bd[k]*phys->stmp[i][k-1]-
            (bd[k]+2*alpha_bot*bd[k])*phys->stmp[i][k]);
      if(Fbot) d[k-ktop]+=dt*(-1+alpha_bot+2*alpha_bot*bd[k])*Fbot[i];
    }

    // First add on the source term from the previous time step.
    if(grid->ctop[i]<=grid->ctopold[i]) {
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++) 
        d[0]+=(1-fab)*Cn[i][grid->ctopold[i]]/(1+fabs(grid->ctop[i]-grid->ctopold[i]));
      for(k=grid->ctopold[i]+1;k<grid->Nk[i];k++) 
        d[k-grid->ctopold[i]]+=(1-fab)*Cn[i][k];
    } else {
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++) 
        d[0]+=(1-fab)*Cn[i][k];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
        d[k-grid->ctop[i]]+=(1-fab)*Cn[i][k];
    }

    for(k=0;k<grid->ctop[i];k++)
      Cn[i][k]=0;

    if(src2)
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        Cn[i][k-ktop]=dt*src2[i][k]*grid->dzzold[i][k];
    else
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        Cn[i][k]=0;

    // Now create the source term for the current time step
    for(k=0;k<grid->Nk[i];k++)
      ap[k]=0;

    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) {
        nc2=nc1;
        if(boundary_scal && (grid->mark[ne]==2 || grid->mark[ne]==3))
          sp=phys->stmp2[nc1];
        else
          sp=phys->stmp[nc1];
      } else 
        sp=phys->stmp[nc2];

      if(!(prop->TVD && prop->horiTVD)) {
        for(k=0;k<grid->Nke[ne];k++) 
          temp[k]=UpWind(phys->utmp2[ne][k],
              phys->stmp[nc1][k],
              sp[k]);
      } else {
        for(k=0;k<grid->Nke[ne];k++) 
          if(phys->utmp2[ne][k]>0)
            temp[k]=phys->SfHp[ne][k];
          else
            temp[k]=phys->SfHm[ne][k];	    
      }

      for(k=0;k<grid->Nke[ne];k++)
        ap[k] += dt*df*normal/Ac*(theta*phys->u[ne][k]+(1-theta)*phys->utmp2[ne][k])
          *temp[k]*grid->dzf[ne][k];
    }

    for(k=ktop+1;k<grid->Nk[i];k++) 
      Cn[i][k-ktop]-=ap[k];

    for(k=0;k<=ktop;k++) 
      Cn[i][0]-=ap[k];

    // Add on the source from the current time step to the rhs.
    for(k=0;k<grid->Nk[i]-ktop;k++) 
      d[k]+=fab*Cn[i][k];

    // Add on the volume correction if h was < -d
    /*
       if(grid->ctop[i]==grid->Nk[i]-1)
       d[grid->Nk[i]-ktop-1]+=phys->hcorr[i]*phys->stmp[i][grid->ctop[i]];
       */

    for(k=ktop;k<grid->Nk[i];k++)
      ap[k]=Cn[i][k-ktop];
    for(k=0;k<=ktop;k++)
      Cn[i][k]=0;
    for(k=ktop+1;k<grid->Nk[i];k++)
      Cn[i][k]=ap[k];
    for(k=grid->ctop[i];k<=ktop;k++)
      Cn[i][k]=ap[ktop]/(1+fabs(grid->ctop[i]-ktop));

    if(grid->Nk[i]-ktop>1) 
      TriSolve(a,b,c,d,&(scal[i][ktop]),grid->Nk[i]-ktop);
    else if(prop->n>1) {
      if(b[0]>0 && phys->active[i])
        scal[i][ktop]=d[0]/b[0];
      else 
        scal[i][ktop]=0;
    }

    for(k=0;k<grid->ctop[i];k++)
      scal[i][k]=0;

    for(k=grid->ctop[i];k<grid->ctopold[i];k++) 
      scal[i][k]=scal[i][ktop];
  }

  // Code to check divergence change CHECKCONSISTENCY to 1 in suntans.h
  if(CHECKCONSISTENCY && checkflag) {

    if(prop->n==1+prop->nstart) {
      smin=INFTY;
      smax=-INFTY;
      for(i=0;i<grid->Nc;i++) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(phys->stmp[i][k]>smax) { 
            smax=phys->stmp[i][k]; 
            imax=i; 
            kmax=k; 
          }
          if(phys->stmp[i][k]<smin) { 
            smin=phys->stmp[i][k]; 
            imin=i; 
            kmin=k; 
          }
        }
      }
      MPI_Reduce(&smin,&smin_value,1,MPI_DOUBLE,MPI_MIN,0,comm);
      MPI_Reduce(&smax,&smax_value,1,MPI_DOUBLE,MPI_MAX,0,comm);
      MPI_Bcast(&smin_value,1,MPI_DOUBLE,0,comm);
      MPI_Bcast(&smax_value,1,MPI_DOUBLE,0,comm);

      if(myproc==0)
        printf("Minimum scalar: %.2f, maximum: %.2f\n",smin_value,smax_value);
    }      

    //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      flag=0;
      for(nf=0;nf<NFACES;nf++) {
	if(grid->mark[grid->face[i*NFACES+nf]]==2 || 
	   grid->mark[grid->face[i*NFACES+nf]]==3) {
          flag=1;
          break;
        }
      }

      if(!flag) {
        div_da=0;

        for(k=0;k<grid->Nk[i];k++) {
          div_da+=grid->Ac[i]*(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;

          div_local=0;
	  for(nf=0;nf<NFACES;nf++) {
	    ne=grid->face[i*NFACES+nf];
            div_local+=(theta*phys->u[ne][k]+(1-theta)*phys->utmp2[ne][k])
	      *grid->dzf[ne][k]*grid->normal[i*NFACES+nf]*grid->df[ne];
          }
          div_da+=div_local;
          div_local+=grid->Ac[i]*(theta*(wnew[i][k]-wnew[i][k+1])+
              (1-theta)*(phys->wtmp2[i][k]-phys->wtmp2[i][k+1]));

          if(k>=grid->ctop[i]) {
            if(fabs(div_local)>SMALL_CONSISTENCY && grid->dzz[imin][0]>DRYCELLHEIGHT) 
              printf("Step: %d, proc: %d, locally-divergent at %d, %d, div=%e\n",
                  prop->n,myproc,i,k,div_local);
          }
        }
        if(fabs(div_da)>SMALL_CONSISTENCY && phys->h[i]+grid->dv[i]>DRYCELLHEIGHT)
          printf("Step: %d, proc: %d, Depth-Ave divergent at i=%d, div=%e\n",
              prop->n,myproc,i,div_da);
      }
    }

    mincount=0;
    maxcount=0;
    smin=INFTY;
    smax=-INFTY;
    //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      flag=0;
      for(nf=0;nf<NFACES;nf++) {
	if(grid->mark[grid->face[i*NFACES+nf]]==2 || grid->mark[grid->face[i*NFACES+nf]]==3) {
          flag=1;
          break;
        }
      }

      if(!flag) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(scal[i][k]>smax) { 
            smax=scal[i][k]; 
            imax=i; 
            kmax=k; 
          }
          if(scal[i][k]<smin) { 
            smin=scal[i][k]; 
            imin=i; 
            kmin=k; 
          }

          if(scal[i][k]>smax_value+SMALL_CONSISTENCY && grid->dzz[i][k]>DRYCELLHEIGHT)
            maxcount++;
          if(scal[i][k]<smin_value-SMALL_CONSISTENCY && grid->dzz[i][k]>DRYCELLHEIGHT)
            mincount++;
        }
      }
    }
    MPI_Reduce(&mincount,&allmincount,1,MPI_INT,MPI_SUM,0,comm);
    MPI_Reduce(&maxcount,&allmaxcount,1,MPI_INT,MPI_SUM,0,comm);

    if(mincount!=0 || maxcount!=0) 
      printf("Not CWC, step: %d, proc: %d, smin = %e at i=%d,H=%e, smax = %e at i=%d,H=%e\n",
          prop->n,myproc,
          smin,imin,phys->h[imin]+grid->dv[imin],
          smax,imax,phys->h[imax]+grid->dv[imax]);

    if(myproc==0 && (allmincount !=0 || allmaxcount !=0))
      printf("Total number of CWC violations (all procs): s<s_min: %d, s>s_max: %d\n",
          allmincount,allmaxcount);
  }
  }
