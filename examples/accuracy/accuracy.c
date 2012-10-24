/**********************************************************************
 *
 * File: accuracy.c
 * Description:  This program analyzes the output from the internal
 * seiche runs in accuracy.sh and prints out the time accuracy results and
 * compares them to the reference solution.  The error results are output
 * in the form error(t=dtmin)/error(t=dt), so that a 2nd order accurate
 * method should result in an increase in the error by a factor of
 * 64 or larger for the case when dt=dtmin*8.
 *
 * Oliver Fringer
 * Stanford University
 * 11 July 2005
 *
 ***********************************************************************/

#include<stdio.h>
#include<math.h>
#include<errno.h>

#define Nkmax 50
#define Nc 50
#define numruns 6

#define DATA "data"
#define REFDATA "rundata/errors.dat"
#define UFILE "u.dat.0"
#define WFILE "w.dat.0"
#define SFILE "s.dat.0"
#define QFILE "q.dat.0"
#define HFILE "fs.dat.0"
#define DTFILE "data/dt.txt"

#define REFERENCE 0

#define REAL double
#define BUFFERSIZE 256

int main(void) {
  int i, k, n, nsteps[numruns];
  float dt[numruns];
  REAL errorU[numruns],errorW[numruns],errorS[numruns],
    errorQ[numruns],errorQ0[numruns],errorH[numruns],
    errorU_ref[numruns],errorW_ref[numruns],errorS_ref[numruns],
    errorQ_ref[numruns],errorQ0_ref[numruns],errorH_ref[numruns];

  FILE *fid;
  char str[BUFFERSIZE];
  REAL uref[3*Nkmax][Nc], wref[Nkmax+1][Nc], sref[Nkmax][Nc],
    qref[Nkmax][Nc], q0ref[Nkmax][Nc], href[Nc];
  REAL u[3*Nkmax][Nc], w[Nkmax+1][Nc], s[Nkmax][Nc],
    q[Nkmax][Nc], q0[Nkmax][Nc], h[Nc];

  fid = fopen(DTFILE,"r");
  for(n=0;n<numruns;n++) 
    fscanf(fid,"%f %d",&dt[n],&nsteps[n]);
  fclose(fid);

  loadData(uref,numruns,UFILE,sizeof(REAL),3*Nc*Nkmax,nsteps[numruns-1]);
  loadData(wref,numruns,WFILE,sizeof(REAL),Nc*(Nkmax+1),nsteps[numruns-1]);
  loadData(sref,numruns,SFILE,sizeof(REAL),Nc*Nkmax,nsteps[numruns-1]);
  loadData(href,numruns,HFILE,sizeof(REAL),Nc,nsteps[numruns-1]);
  loadData(qref,numruns,QFILE,sizeof(REAL),Nc*Nkmax,nsteps[numruns-1]-1);
  loadData(q0ref,numruns,QFILE,sizeof(REAL),Nc*Nkmax,nsteps[numruns-1]);

  for(k=0;k<Nkmax;k++) 
    for(i=0;i<Nc;i++) 
      qref[k][i] = 1.5*q0ref[k][i]-0.5*qref[k][i];

  for(n=0;n<numruns;n++) {
    errorU[n]=0;
    errorW[n]=0;
    errorS[n]=0;
    errorQ[n]=0;
    errorQ0[n]=0;
    errorH[n]=0;
  }

  for(i=0;i<Nc;i++) 
    errorH[numruns-1]+=href[i]*href[i];
  for(k=0;k<Nkmax;k++) 
    for(i=0;i<Nc;i++) {
      errorS[numruns-1]+=sref[k][i]*sref[k][i];
      errorQ[numruns-1]+=qref[k][i]*qref[k][i];
      errorQ0[numruns-1]+=q0ref[k][i]*q0ref[k][i];
    }
  
  for(k=0;k<Nkmax;k++) 
    for(i=0;i<Nc;i++) {
      errorU[numruns-1]+=uref[3*k][i]*uref[3*k][i];
      errorW[numruns-1]+=uref[3*k+2][i]*uref[3*k+2][i];
    }

  for(n=0;n<numruns-1;n++) {
    loadData(u,n+1,UFILE,sizeof(REAL),3*Nc*Nkmax,nsteps[n]);
    loadData(w,n+1,WFILE,sizeof(REAL),Nc*(Nkmax+1),nsteps[n]);
    loadData(s,n+1,SFILE,sizeof(REAL),Nc*Nkmax,nsteps[n]);
    loadData(h,n+1,HFILE,sizeof(REAL),Nc,nsteps[n]);
    loadData(q,n+1,QFILE,sizeof(REAL),Nc*Nkmax,nsteps[n]-1);
    loadData(q0,n+1,QFILE,sizeof(REAL),Nc*Nkmax,nsteps[n]);
    
    for(k=0;k<Nkmax;k++) 
      for(i=0;i<Nc;i++) 
	q[k][i] = 1.5*q0[k][i]-0.5*q[k][i];
    
    for(i=0;i<Nc;i++) 
      errorH[n]+=(h[i]-href[i])*(h[i]-href[i]);
    for(k=0;k<Nkmax;k++) 
      for(i=0;i<Nc;i++) {
	errorS[n]+=(s[k][i]-sref[k][i])*(s[k][i]-sref[k][i]);
	errorQ[n]+=(q[k][i]-qref[k][i])*(q[k][i]-qref[k][i]);
	errorQ0[n]+=(q0[k][i]-q0ref[k][i])*(q0[k][i]-q0ref[k][i]);
      }
    
    for(k=0;k<Nkmax;k++) 
      for(i=0;i<Nc;i++) {
	errorU[n]+=(u[3*k][i]-uref[3*k][i])*(u[3*k][i]-uref[3*k][i]);
	errorW[n]+=(u[3*k+2][i]-uref[3*k+2][i])*(u[3*k+2][i]-uref[3*k+2][i]);
      }
  }

  if(REFERENCE) {
    fid = fopen(REFDATA,"w");
    fwrite(errorU,sizeof(REAL),numruns,fid);
    fwrite(errorW,sizeof(REAL),numruns,fid);
    fwrite(errorS,sizeof(REAL),numruns,fid);
    fwrite(errorQ,sizeof(REAL),numruns,fid);
    fwrite(errorQ0,sizeof(REAL),numruns,fid);
    fwrite(errorH,sizeof(REAL),numruns,fid);
    fclose(fid);
  } else {
    if((fid = fopen(REFDATA,"r"))==NULL) {
      sprintf(str,"Error opening %s:",REFDATA);
      perror(str);
      return 1;
    }
    fread(errorU_ref,sizeof(REAL),numruns,fid);
    fread(errorW_ref,sizeof(REAL),numruns,fid);
    fread(errorS_ref,sizeof(REAL),numruns,fid);
    fread(errorQ_ref,sizeof(REAL),numruns,fid);
    fread(errorQ0_ref,sizeof(REAL),numruns,fid);
    fread(errorH_ref,sizeof(REAL),numruns,fid);
  }
    
  fprintf(stdout,"Error results (Error(n)/Error(0)): \n\n");
  if(REFERENCE) 
    fprintf(stdout,"Reference results:\n");
  else
    fprintf(stdout,"Your results:\n");
  fprintf(stdout,"--------------------------------------------------\n");
  fprintf(stdout,"dt0/dt\tU\tW\tS\tQ\tQ0\th\n");
  fprintf(stdout,"--------------------------------------------------\n");

  for(n=0;n<numruns-1;n++) { 
    fprintf(stdout,"%.0f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",
	    dt[0]/dt[n],
	    sqrt(errorU[0]/errorU[n]),
	    sqrt(errorW[0]/errorW[n]),
	    sqrt(errorS[0]/errorS[n]),
	    sqrt(errorQ[0]/errorQ[n]),
	    sqrt(errorQ0[0]/errorQ0[n]),
	    sqrt(errorH[0]/errorH[n]));
  }

  if(!REFERENCE) {
    fprintf(stdout,"\nReference results:\n");
    for(n=0;n<numruns-1;n++) { 
      fprintf(stdout,"%.0f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",
	      dt[0]/dt[n],
	      sqrt(errorU_ref[0]/errorU_ref[n]),
	      sqrt(errorW_ref[0]/errorW_ref[n]),
	      sqrt(errorS_ref[0]/errorS_ref[n]),
	      sqrt(errorQ_ref[0]/errorQ_ref[n]),
	      sqrt(errorQ0_ref[0]/errorQ0_ref[n]),
	      sqrt(errorH_ref[0]/errorH_ref[n]));
    }

    fprintf(stdout,"\nDifference (relative):\n");
    for(n=0;n<numruns-1;n++) { 
      fprintf(stdout,"%.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
	      dt[0]/dt[n],
	      (1-sqrt(errorU[0]/errorU[n])/sqrt(errorU_ref[0]/errorU_ref[n])),
	      (1-sqrt(errorW[0]/errorW[n])/sqrt(errorW_ref[0]/errorW_ref[n])),
	      (1-sqrt(errorS[0]/errorS[n])/sqrt(errorS_ref[0]/errorS_ref[n])),
	      (1-sqrt(errorQ[0]/errorQ[n])/sqrt(errorQ_ref[0]/errorQ_ref[n])),
	      (1-sqrt(errorQ0[0]/errorQ0[n])/sqrt(errorQ0_ref[0]/errorQ0_ref[n])),
	      (1-sqrt(errorH[0]/errorH[n])/sqrt(errorH_ref[0]/errorH_ref[n])));
    }
  }

  return 0;
}

int loadData(void *array, int run, char *file, int size, int numels, int nstep) {
  int n;
  char str[BUFFERSIZE];
  FILE *fid;

  sprintf(str,"%s/t%d/%s",DATA,run,file);
  fid = fopen(str,"r");
  fseek(fid,(nstep-1)*size*numels,SEEK_SET);
  fread(array,size,numels,fid);
  fclose(fid);
}
