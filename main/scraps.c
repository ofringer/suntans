/*
 * $Id: scraps.c,v 1.2 2004-05-29 20:28:50 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
static void MakePointers0(gridT *maingrid, gridT **localgrid, int myproc);
static void ReOrder0(gridT *maingrid, gridT **localgrid, int myproc);
void SendRecvCellData0(REAL *celldata, gridT *maingrid, gridT *localgrid, 
		       int myproc, MPI_Comm comm);

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

static void MakePointers0(gridT *maingrid, gridT **localgrid, int myproc)
{
  int n, nf, ne, neigh, neighproc, nc, j, k, mark, bctype, count;
  int **cell_send, **cell_recv, **edge_send, **edge_recv;
  int *num_cells_send, *num_cells_recv, *num_edges_send, *num_edges_recv;
  int *cellp, *edgep, *celldist, *edgedist;
  int kcellsend, kcellrecv, kedgesend, kedgerecv;

  cellp = (int *)malloc((*localgrid)->Nc*sizeof(int));
  edgep = (int *)malloc((*localgrid)->Ne*sizeof(int));
  celldist = (int *)malloc((MAXBCTYPES-1)*sizeof(int));
  edgedist = (int *)malloc((MAXMARKS-1)*sizeof(int));

  k=0;
  celldist[0]=0;
  // Put computational cells first
  for(n=0;n<(*localgrid)->Nc;n++)
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==0 ||
       IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2)
	cellp[k++]=n;
  celldist[1]=k;
  // Followed by boundary cells
  for(n=0;n<(*localgrid)->Nc;n++)
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==1)
      cellp[k++]=n;
  celldist[2]=k;

  k=0;
  edgedist[0]=0;
  // Put computational edges first
  for(n=0;n<(*localgrid)->Ne;n++)
    if((*localgrid)->mark[n]==0 || (*localgrid)->mark[n]==5)
      edgep[k++]=n;
  edgedist[1]=k;
  // Followed by boundary cells
  for(mark=1;mark<MAXMARKS-2;mark++) {
    for(n=0;n<(*localgrid)->Ne;n++)
      if((*localgrid)->mark[n]==mark)
	edgep[k++]=n;
    edgedist[mark+1]=k;
  }

  (*localgrid)->cellp=cellp;
  (*localgrid)->edgep=edgep;
  (*localgrid)->celldist=celldist;
  (*localgrid)->edgedist=edgedist;

  /*
  for(bctype=0;bctype<MAXBCTYPES-2;bctype++) {
    printf("Cells of type %d: ",bctype);
    for(j=celldist[bctype];j<celldist[bctype+1];j++)
      printf("%d ",cellp[j]);
    printf("\n");
  }
  for(mark=0;mark<MAXMARKS-2;mark++) {
    printf("Edges of type %d: ",mark);
    for(j=edgedist[mark];j<edgedist[mark+1];j++)
      printf("%d ",edgep[j]);
    printf("\n");
  }
  */

  cell_send = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  cell_recv = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  edge_send = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  edge_recv = (int **)malloc((*localgrid)->Nneighs*sizeof(int *));
  num_cells_send = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_cells_recv = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_edges_send = (int *)malloc((*localgrid)->Nneighs*sizeof(int));
  num_edges_recv = (int *)malloc((*localgrid)->Nneighs*sizeof(int));

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    num_cells_send[neigh]=0;
    num_cells_recv[neigh]=0;
    num_edges_send[neigh]=0;
    num_edges_recv[neigh]=0;
  }
    
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    neighproc=(*localgrid)->myneighs[neigh];
    for(n=0;n<(*localgrid)->Nc;n++) {
      if(IsCellNeighborProc(n,maingrid,(*localgrid),myproc,neighproc)) {
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2) 
	  num_cells_send[neigh]++;
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==3)
	  num_cells_recv[neigh]++;
      }
    }
    for(n=0;n<(*localgrid)->Ne;n++)
      if((*localgrid)->mark[n]==5) 
	for(j=0;j<2;j++) {
	  nc = (*localgrid)->grad[2*n+j];
	  if(maingrid->part[(*localgrid)->mnptr[nc]]==neighproc) {
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4))
		num_edges_recv[neigh]++;
	    }
	    if(j==0)
	      nc=(*localgrid)->grad[2*n+1];
	    else
	      nc=(*localgrid)->grad[2*n];
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4))
		num_edges_send[neigh]++;
	    }
	  }
	}
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    cell_send[neigh] = (int *)malloc(num_cells_send[neigh]*sizeof(int));
    cell_recv[neigh] = (int *)malloc(num_cells_recv[neigh]*sizeof(int));
    edge_send[neigh] = (int *)malloc(num_edges_send[neigh]*sizeof(int));
    edge_recv[neigh] = (int *)malloc(num_edges_recv[neigh]*sizeof(int));
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    kcellsend=0;
    kcellrecv=0;
    neighproc=(*localgrid)->myneighs[neigh];
    for(n=0;n<(*localgrid)->Nc;n++) {
      if(IsCellNeighborProc(n,maingrid,(*localgrid),myproc,neighproc)) {
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==2)
	  cell_send[neigh][kcellsend++]=n;
	if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==3)
	  cell_recv[neigh][kcellrecv++]=n;
      }
    }
    kedgesend=0;
    kedgerecv=0;
    for(n=0;n<(*localgrid)->Ne;n++)
      if((*localgrid)->mark[n]==5) 
	for(j=0;j<2;j++) {
	  nc = (*localgrid)->grad[2*n+j];
	  if(maingrid->part[(*localgrid)->mnptr[nc]]==neighproc) {
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4) &&
		 IsMember(ne,edge_recv[neigh],kedgerecv)==-1)
		edge_recv[neigh][kedgerecv++]=ne;
	    }
	    if(j==0)
	      nc=(*localgrid)->grad[2*n+1];
	    else
	      nc=(*localgrid)->grad[2*n];
	    for(nf=0;nf<NFACES;nf++) {
	      ne = (*localgrid)->face[nc*NFACES+nf];
	      if(n != ne && ((*localgrid)->mark[ne]<1 || (*localgrid)->mark[ne]>4) &&
		 IsMember(ne,edge_send[neigh],kedgesend)==-1)
		edge_send[neigh][kedgesend++]=ne;
	    }
	  }
	}
    num_edges_send[neigh]=kedgesend;
    num_edges_recv[neigh]=kedgerecv;
    edge_send[neigh]=ReSize(edge_send[neigh],num_edges_send[neigh]);
    edge_recv[neigh]=ReSize(edge_recv[neigh],num_edges_recv[neigh]);
  }

  /*
  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    printf("%d <-> %d, Cells to Send: %d, Receive: %d\n",
	   myproc,(*localgrid)->myneighs[neigh],
	   num_cells_send[neigh],num_cells_recv[neigh]);
    printf("%d <-> %d, Edges to Send: %d, Receive: %d\n",
	   myproc,(*localgrid)->myneighs[neigh],
	   num_edges_send[neigh],num_edges_recv[neigh]);
  }

  for(neigh=0;neigh<(*localgrid)->Nneighs;neigh++) {
    printf("%d -> %d (Cells): ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_cells_send[neigh];n++)
      printf("%d ",cell_send[neigh][n]);
    printf("\n");
    printf("%d <- %d (Cells) ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_cells_recv[neigh];n++)
      printf("%d ",cell_recv[neigh][n]);
    printf("\n");
    printf("%d -> %d (Edges): ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_edges_send[neigh];n++)
      printf("%d ",edge_send[neigh][n]);
    printf("\n");
    printf("%d <- %d (Edges): ",myproc,(*localgrid)->myneighs[neigh]);
    for(n=0;n<num_edges_recv[neigh];n++)
      printf("%d ",edge_recv[neigh][n]);
    printf("\n");
  }
  */

  (*localgrid)->cell_send=cell_send;
  (*localgrid)->cell_recv=cell_recv;
  (*localgrid)->edge_send=edge_send;
  (*localgrid)->edge_recv=edge_recv;
  (*localgrid)->num_cells_send=num_cells_send;
  (*localgrid)->num_cells_recv=num_cells_recv;
  (*localgrid)->num_edges_send=num_edges_send;
  (*localgrid)->num_edges_recv=num_edges_recv;
}

static void ReOrder0(gridT *maingrid, gridT **localgrid, int myproc)
{
  int *eorder, *eorderp, *corder, *corderp, *eflag, *cflag, *bctypeind, *marktypeind;
  int neigh, neighproc, j, k, n, nf, bctype, marktype, Nc=(*localgrid)->Nc, Ne=(*localgrid)->Ne;
  int nstart, nend, Nc_comp, Ne_comp;
  REAL *tmp, temp;

  corder=(int *)malloc(Nc*sizeof(int));
  eorder=(int *)malloc(Ne*sizeof(int));
  corderp=(int *)malloc(Nc*sizeof(int));
  eorderp=(int *)malloc(Ne*sizeof(int));
  tmp = (REAL *)malloc(2*(NFACES-1)*Ne*sizeof(REAL));

  // Place all of the computational cells before the boundary cells
  k=0;
  for(n=0;n<Nc;n++)
    corder[n]=0;
  for(n=0;n<Nc;n++) 
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)<MAXBCTYPES-1) 
      corder[k++]=n;
  (*localgrid)->Nc_comp = k;
  for(n=0;n<Nc;n++) 
    if(IsBoundaryCell((*localgrid)->mnptr[n],maingrid,myproc)==MAXBCTYPES-1) 
      corder[k++]=n;

  for(n=0;n<Nc;n++) 
    corderp[corder[n]]=n;

  // Place all of the computational edges before the boundary edges
  k=0;
  for(n=0;n<Ne;n++)
    eorder[n]=0;
  for(n=0;n<Ne;n++)
    if((*localgrid)->mark[n]==0 || (*localgrid)->mark[n]==5)
      eorder[k++]=n;
  (*localgrid)->Ne_comp = k;
  for(n=0;n<Ne;n++)
    if((*localgrid)->mark[n]!=0 && (*localgrid)->mark[n]!=5)
      eorder[k++]=n;

  /*  Sort(&(eorder[(*localgrid)->Ne_comp]),
       &((*localgrid)->eptr[(*localgrid)->Ne_comp]),
       (*localgrid)->Ne-(*localgrid)->Ne_comp);*/
  //  Sort(&eorder,&((*localgrid)->eptr[0]),Ne);

  for(n=0;n<Ne;n++)
    eorderp[eorder[n]]=n;

  // Reorder the data corresponding to cells with the corder array
  ReOrderRealArray((*localgrid)->Ac,corder,tmp,(*localgrid)->Nc,1);
  ReOrderRealArray((*localgrid)->xv,corder,tmp,(*localgrid)->Nc,1);
  ReOrderRealArray((*localgrid)->yv,corder,tmp,(*localgrid)->Nc,1);
  ReOrderRealArray((*localgrid)->dv,corder,tmp,(*localgrid)->Nc,1);
  ReOrderIntArray((*localgrid)->cells,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->face,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->normal,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->neigh,corder,(int *)tmp,(*localgrid)->Nc,NFACES);
  ReOrderIntArray((*localgrid)->vwgt,corder,(int *)tmp,(*localgrid)->Nc,1);
  ReOrderIntArray((*localgrid)->mnptr,corder,(int *)tmp,(*localgrid)->Nc,1);
  ReOrderIntArray((*localgrid)->Nk,corder,(int *)tmp,(*localgrid)->Nc,1);

  // Reorder the data corresponding to edges with the eorder array
  ReOrderRealArray((*localgrid)->df,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->dg,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->n1,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->n2,eorder,tmp,(*localgrid)->Ne,1);
  ReOrderRealArray((*localgrid)->xi,eorder,tmp,(*localgrid)->Ne,2*(NFACES-1));
  ReOrderIntArray((*localgrid)->grad,eorder,(int *)tmp,(*localgrid)->Ne,2);
  ReOrderIntArray((*localgrid)->eneigh,eorder,(int *)tmp,(*localgrid)->Ne,2*(NFACES-1));
  ReOrderIntArray((*localgrid)->edges,eorder,(int *)tmp,(*localgrid)->Ne,NUMEDGECOLUMNS);
  ReOrderIntArray((*localgrid)->mark,eorder,(int *)tmp,(*localgrid)->Ne,1);
  ReOrderIntArray((*localgrid)->eptr,eorder,(int *)tmp,(*localgrid)->Ne,1);
  ReOrderIntArray((*localgrid)->Nke,eorder,(int *)tmp,(*localgrid)->Ne,1);

  // Now adjust the pointers to point to the new locations
  // Face and eneigh arrays point to edges
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<NFACES;nf++) 
      tmp[n*NFACES+nf]=(*localgrid)->face[n*NFACES+nf];
  for(n=0;n<Nc;n++) 
    for(nf=0;nf<NFACES;nf++) 
      (*localgrid)->face[n*NFACES+nf]=eorderp[(int)(tmp[n*NFACES+nf])];

  for(n=0;n<Ne;n++)
    for(nf=0;nf<2*(NFACES-1);nf++)
      tmp[n*2*(NFACES-1)+nf]=(*localgrid)->eneigh[n*2*(NFACES-1)+nf];
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2*(NFACES-1);nf++)
      if((int)(tmp[n*2*(NFACES-1)+nf])!=-1)
	(*localgrid)->eneigh[n*2*(NFACES-1)+nf]=eorderp[(int)(tmp[n*2*(NFACES-1)+nf])];
      else
	(*localgrid)->eneigh[n*2*(NFACES-1)+nf]=-1;

  // Grad and neigh arrays point to cells
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2;nf++)
      tmp[2*n+nf]=(*localgrid)->grad[2*n+nf];
  for(n=0;n<Ne;n++)
    for(nf=0;nf<2;nf++)
      if((int)(tmp[2*n+nf])!=-1)
	(*localgrid)->grad[2*n+nf]=corderp[(int)(tmp[2*n+nf])];
      else
	(*localgrid)->grad[2*n+nf]=-1;

  for(n=0;n<Nc;n++)
    for(nf=0;nf<NFACES;nf++)
      tmp[n*NFACES+nf]=(*localgrid)->neigh[n*NFACES+nf];
  for(n=0;n<Nc;n++)
    for(nf=0;nf<NFACES;nf++)
      if((int)(tmp[n*NFACES+nf])!=-1)
	(*localgrid)->neigh[n*NFACES+nf]=corderp[(int)(tmp[n*NFACES+nf])];
      else
  	(*localgrid)->neigh[n*NFACES+nf]=-1;

  free(corder);
  free(eorder);
  free(corderp);
  free(eorderp);
  free(tmp);
}

void SendRecvCellData0(REAL *celldata, gridT *maingrid, gridT *localgrid, 
		       int myproc, MPI_Comm comm)
{
  int n, neigh, neighproc, receivesize, sendsize, noncontig, flag;
  int senstart, senend, recstart, recend;
  REAL **recv, **send;
  MPI_Status status;

  recv = (REAL **)malloc(maingrid->numneighs[myproc]*sizeof(REAL *));
  send = (REAL **)malloc(maingrid->numneighs[myproc]*sizeof(REAL *));

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    neighproc = maingrid->neighs[myproc][neigh];

    send[neigh] = (REAL *)malloc(localgrid->num_cells_send[neigh]*sizeof(REAL));
    recv[neigh] = (REAL *)malloc(localgrid->num_cells_recv[neigh]*sizeof(REAL));

    for(n=0;n<localgrid->num_cells_send[neigh];n++)
      send[neigh][n]=celldata[localgrid->cell_send[neigh][n]];

    MPI_Send((void *)(send[neigh]),localgrid->num_cells_send[neigh],MPI_DOUBLE,neighproc,1,comm); 
    MPI_Recv((void *)(recv[neigh]),localgrid->num_cells_recv[neigh],MPI_DOUBLE,neighproc,1,comm,&status);
  }
  MPI_Barrier(comm);

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    for(n=0;n<localgrid->num_cells_recv[neigh];n++)
      celldata[localgrid->cell_recv[neigh][n]]=recv[neigh][n];
  }

  for(neigh=0;neigh<maingrid->numneighs[myproc];neigh++) {
    free(send[neigh]);
    free(recv[neigh]);
  }
  free(send);
  free(recv);
}

/* From phys.c */
static REAL InterpolateEdge0(REAL xd, REAL yd, REAL zd,  
			    int j0, int k0, REAL z0,  int jnear, int knear, 
			    gridT *grid, physT *phys, propT *prop, char type) {
  int n, j, k, je, nk, status, nc1, nc2, sign;
  REAL U, V, un1, un2, ui, dz,  
    *uf = (REAL *)malloc(grid->Nnearestedges*sizeof(REAL)), 
    *x = (REAL *)malloc(grid->Nnearestedges*sizeof(REAL)), 
    *y = (REAL *)malloc(grid->Nnearestedges*sizeof(REAL));

  // U is the x component of velocity at the face
  // V is the y component of velocity at the face.
  // U^2 + V^2 = un^2 + ut^2

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

  for(n=0;n<grid->Nnearestedges;n++) {
    je = grid->nearestedges[j0][n];

    nc1=grid->grad[2*je];
    nc2=grid->grad[2*je+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    if(type=='u') {
      U=phys->u[je][knear]*grid->n1[je]-phys->ut[je][knear]*grid->n2[je];
      V=phys->u[je][knear]*grid->n2[je]+phys->ut[je][knear]*grid->n1[je];
      un1=U*grid->n1[j0]+V*grid->n2[j0];

      U=phys->u[je][knear+nk]*grid->n1[je]-phys->ut[je][knear+nk]*grid->n2[je];
      V=phys->u[je][knear+nk]*grid->n2[je]+phys->ut[je][knear+nk]*grid->n1[je];
      un2=U*grid->n1[j0]+V*grid->n2[j0];
    } else {
      un1=phys->wf[je][knear+nk];
      un2=phys->wf[je][knear+nk];
    }

    dz = 0.25*(grid->dzz[nc1][knear]+grid->dzz[nc2][knear]+
	       grid->dzz[nc1][knear+nk]+grid->dzz[nc2][knear+nk]);

    // If dz=0 this means that the interpolation is being performed in
    // a dry zone, so set uf[n]=un1 if this is so.
    if(dz==0)
      uf[n]=un1;
    else 
      uf[n] = un1 + sign*(zd-z0)*(un2-un1)/dz;

    /*
    if((un2>un1 && (uf[n]>un2 || uf[n]<un1)) ||
       (un2<un1 && (uf[n]<un2 || uf[n]>un1)))
      printf("Error! Vertical interpolation is not monotonic!: \n\
              knear=%d 1)%.4e 2)%.4e i)%.4e\n",
	     knear,un1,un2,uf[n]);
    */
    x[n] = grid->xe[je];
    y[n] = grid->ye[je];
  }

  kriging(xd,yd,&ui,x,y,uf,prop->kriging_cov,grid->Nnearestedges,&status);

  free(x);
  free(y);
  free(uf);

  return ui;
}
 
// Radiation condition
static void Radiate(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {
  REAL c, phin[NFACES], phin0, se, phinu;
  int i, iptr, j, jptr, k, n, nf, ne, ns[NFACES], count, nc1, nc2;
  
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    count=0;
    ns[0] = grid->grad[2*j];
    for(nf=0;nf<NFACES;nf++) {
      if(grid->neigh[ns[0]*NFACES+nf]!=-1)
	ns[count++]=grid->neigh[ns[0]*NFACES+nf];
    }

    for(k=grid->etop[j];k<grid->Nke[j];k++) {

      phin0=0;
      for(n=0;n<count;n++) {
	phin[n]=0;
	for(nf=0;nf<NFACES;nf++) {
	  ne=grid->face[ns[n]*NFACES+nf];
	  nc1=grid->grad[2*ne];
	  nc2=grid->grad[2*ne+1];
	  if(nc1==-1) nc1=nc2;
	  if(nc2==-1) nc2=nc1;
	  se = 0.5*(phys->s[nc1][k]+phys->s[nc2][k]);
	  phin[n]+=(grid->n1[ne]*grid->n1[j]+grid->n2[ne]*grid->n2[j])*se;
	  if(n==0)
	    phinu+=(grid->n1[ne]*grid->n1[j]+grid->n2[ne]*grid->n2[j])*
	      (phys->u[ne][k]*grid->n1[ne]*grid->n1[j]+phys->ut[ne][k]*grid->n2[ne]*grid->n2[j])*grid->df[ne];
	}
	phin0+=phin[n];
      }
      if(phin0!=0)
	c = -(phys->s[ns[0]][k]-phys->stmp[ns[0]][k])*count/prop->dt/phin0;
      else
	c = 0;
      //      printf("c = %e, phin0=%e\n",c,phin0);
      if(c<0)
      	c=0;
      if(c>.1)
	c=.1;
      phys->u[j][k]-=c*prop->dt/grid->Ac[ns[0]]*phinu;
      //      phys->u[j][k] = -(phys->h[ns[0]]+.5)*sqrt(GRAV/(grid->dv[ns[0]]+phys->h[ns[0]]));
    }
  }
}

  // Now correct the velocity change as a result of the change in volume in the upper layer
  // This does not work with wetting and drying!
  ComputeVelocityVector(phys->u,phys->uc,phys->vc,grid);
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    umom=0;
    vmom=0;
    dz=0;
    if(grid->ctopold[i]<grid->ctop[i]) {
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++) {
	umom+=phys->uc[i][k]*grid->dzzold[i][k];
	vmom+=phys->vc[i][k]*grid->dzzold[i][k];
      }
      phys->uc[i][grid->ctop[i]]=umom/grid->dzz[i][grid->ctop[i]];
      phys->vc[i][grid->ctop[i]]=vmom/grid->dzz[i][grid->ctop[i]];
    } else {
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++) 
	dz+=grid->dzzold[i][k];
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++) {
	phys->uc[i][k]=grid->dzzold[i][grid->ctopold[i]]*phys->uc[i][grid->ctopold[i]]/dz;
	phys->vc[i][k]=grid->dzzold[i][grid->ctopold[i]]*phys->vc[i][grid->ctopold[i]]/dz;
      }
    }
  }
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(grid->ctop[nc1]<grid->ctop[nc2]) 
      for(k=grid->ctop[nc1];k<=grid->ctop[nc2];k++)
	phys->u[j][k]=0.5*(grid->n1[j]*phys->uc[nc1][k]+grid->n2[j]*phys->vc[nc1][k]+
			   grid->n1[j]*phys->uc[nc2][grid->ctop[nc2]]+grid->n2[j]*phys->vc[nc2][grid->ctop[nc2]]);
    else
      for(k=grid->ctop[nc2];k<=grid->ctop[nc1];k++)
	phys->u[j][k]=0.5*(grid->n1[j]*phys->uc[nc1][grid->ctop[nc1]]+grid->n2[j]*phys->vc[nc1][grid->ctop[nc1]]+
			   grid->n1[j]*phys->uc[nc2][k]+grid->n2[j]*phys->vc[nc2][k]);
  }

static void UpdateScalars_new(gridT *grid, physT *phys, propT *prop, REAL **scal, REAL **Cn, 
			  REAL kappa, REAL kappaH, REAL **kappa_tv, REAL theta)
{
  int i, iptr, k, kmin, nf, ktop;
  int Nc=grid->Nc, normal, nc1, nc2, ne;
  REAL df, dg, Ac, dt=prop->dt, fab, *a, *b, *c, *d, *r, *ap, *am, dznew, mass, alpha;

  ap = phys->ap;
  am = phys->am;
  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;
  r = phys->bp;

  // 0 for central differencing, 1 for first-order upwinding
  alpha = 1;

  if(prop->n==1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
	Cn[i][k]=0;
  } else
    fab=1.5;
  fab=1;

  for(i=0;i<Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=scal[i][k];

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

    // Advective fluxes
    for(k=ktop;k<grid->Nk[i]+1;k++) {
      ap[k] = 0.5*(phys->w[i][k]+alpha*fabs(phys->w[i][k]));
      am[k] = 0.5*(phys->w[i][k]-alpha*fabs(phys->w[i][k]));
    }
    // Diffusive fluxes
    d[ktop]=0;
    for(k=ktop+1;k<grid->Nk[i];k++) 
      d[k] = 0*(2*kappa+kappa_tv[i][k-1]+kappa_tv[i][k])/(grid->dzz[i][k-1]+grid->dzz[i][k]);

    for(k=ktop+1;k<grid->Nk[i]-1;k++) {
      a[k]=-theta*dt*(d[k]-am[k]);
      b[k]=grid->dzz[i][k]+theta*dt*(d[k]+d[k+1]-am[k+1]+ap[k]);
      c[k]=-theta*dt*(d[k+1]+ap[k+1]);
    }

    // Top cell
    a[ktop]=0;
    b[ktop]=dznew+theta*dt*(d[ktop+1]-am[ktop+1]);
    c[ktop]=-theta*dt*(d[ktop+1]+ap[ktop+1]);

    // Bottom cell
    a[grid->Nk[i]-1]=-theta*dt*(d[grid->Nk[i]-1]-am[grid->Nk[i]-1]);
    b[grid->Nk[i]-1]=grid->dzz[i][grid->Nk[i]-1]+theta*dt*(d[grid->Nk[i]-1]+ap[grid->Nk[i]-1]);
    c[grid->Nk[i]-1]=0;

    // Advective fluxes
    for(k=0;k<grid->Nk[i]+1;k++) {
      ap[k] = 0.5*(phys->wtmp2[i][k]+alpha*fabs(phys->wtmp2[i][k]));
      am[k] = 0.5*(phys->wtmp2[i][k]-alpha*fabs(phys->wtmp2[i][k]));
    }

    // Create the right hand side
    for(k=ktop+1;k<grid->Nk[i]-1;k++) 
      r[k]=grid->dzzold[i][k]*phys->stmp[i][k]+
	(1-theta)*dt*((d[k]-am[k])*phys->stmp[i][k-1]-
		      (d[k]+d[k+1]-am[k+1]+ap[k])*phys->stmp[i][k]+
		      (d[k+1]+ap[k+1])*phys->stmp[i][k+1]);

    r[ktop]=0;
    if(grid->ctopold[i]<=grid->ctop[i])
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
	r[ktop]+=grid->dzzold[i][k]*phys->stmp[i][k];
    else
      r[ktop]+=grid->dzzold[i][ktop]*phys->stmp[i][ktop];

    if(ktop<grid->Nk[i]-1) {
      //Flux through bottom of top cell
      k=ktop;
      r[k]+=(1-theta)*dt*(-(d[k+1]-am[k+1])*phys->stmp[i][k]+
			 (d[k+1]+ap[k+1])*phys->stmp[i][k+1]);

      // Through top of bottom cell
      k=grid->Nk[i]-1;
      r[k]=grid->dzzold[i][k]*phys->stmp[i][k]+(1-theta)*dt*((d[k]-am[k])*phys->stmp[i][k-1]-
			 (d[k]+ap[k])*phys->stmp[i][k]);
    }

    // First add on the source term from the previous time step.
    for(k=ktop;k<grid->Nk[i];k++)
      r[k]+=(1-fab)*Cn[i][k];

    for(k=0;k<grid->Nk[i];k++)
      Cn[i][k]=0;

    // Now create the source term for the current time step
    for(nf=0;nf<NFACES;nf++) {
    
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      // Horizontal advection
      for(k=0;k<grid->etopold[ne];k++)
	ap[k]=0;
      for(k=grid->etopold[ne];k<grid->Nkc[ne];k++) 
	ap[k] = dt*df*normal/Ac*(0.5*(phys->utmp2[ne][k]+alpha*fabs(phys->utmp2[ne][k]))*
				 phys->stmp[nc2][k]*grid->dzzold[nc2][k]
				 +0.5*(phys->utmp2[ne][k]-alpha*fabs(phys->utmp2[ne][k]))*
				 phys->stmp[nc1][k]*grid->dzzold[nc1][k]);
      
      for(k=ktop+1;k<grid->Nk[i];k++) 
      	Cn[i][k]-=ap[k];

      for(k=0;k<=ktop;k++)
	Cn[i][ktop]-=ap[k];
    }

    // Add on the source from the current time step to the rhs.
    for(k=ktop;k<grid->Nk[i];k++)
      r[k]+=fab*Cn[i][k];

    if(grid->Nk[i]-ktop>1) 
      TriSolve(&(a[ktop]),&(b[ktop]),&(c[ktop]),&(r[ktop]),&(scal[i][ktop]),grid->Nk[i]-ktop);
    else if(b[ktop]!=0)
      scal[i][ktop]=r[ktop]/b[ktop];

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
