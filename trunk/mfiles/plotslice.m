%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: plotslice.m
% Description: Create a surface plot of 2d suntans data.
%
% Oliver Fringer
% Stanford University
% 15 Jun 04
%
% $Id: plotslice.m,v 1.1 2004-06-16 02:29:36 fringer Exp $
% $Log: not supported by cvs2svn $
%
% Surface plot determined by:
% PLOT: 0: q, 1: s, 2: u, 3: w
%
% Example: plotslice(1,'/home/data/suntans_data',10,0);
% will plot the salinity data without the free surface at time step 10
% for processor 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotslice(PLOT,datadir,n,proc)

  EMPTY=999999;         % Empty cells are defined by this
  dbl=8;                % Size of double precision in bytes
  precision='float64';  % Precision for reading in data
  
  % cellcentered data contains the voronoi points and the depths
  % at those points.
  cellcentereddata=load([datadir,'/celldata.dat.',num2str(proc)],'-ascii');
  xv = cellcentereddata(:,1);
  yv = cellcentereddata(:,2);
  dv = cellcentereddata(:,4);
  dz = load([datadir,'/vertspace.dat']);
  
  % Total number of cells in the horizontal Nc and vertical Nk
  Nc = length(xv);
  Nk = length(dz);
  
  % Length and depth of domain
  L = max(xv);
  dmax = max(dv);
  
  % Set up the Cartesian grid
  z = -dmax+cumsum(dz);
  [xs,is]=sort(xv);
  [X,Z]=meshgrid(xs,z);
  
  % Empty the cells below the bottom
  D = ones(Nk,1)*dv(is)';
  Z(find(Z<-D))=nan;
  dv = dv(is);
  
  % Open up file descriptors for binary files
  qfile = [datadir,'/q.dat.',num2str(proc)];
  qfid = fopen(qfile,'rb');
  ufile = [datadir,'/u.dat.',num2str(proc)];
  ufid = fopen(ufile,'rb');
  sfile = [datadir,'/s.dat.',num2str(proc)];
  sfid = fopen(sfile,'rb');
  hfile = [datadir,'/fs.dat.',num2str(proc)];
  hfid = fopen(hfile,'rb');
  
  % Total number of time steps to plot is obtained by checking
  % size of file
  status=fseek(sfid,0,'eof');
  nsteps = ftell(sfid)/Nk/Nc/dbl;
  status=fseek(sfid,0,'bof');
  
  % Determine hmax for plotting
  hd = getcdata(hfid,Nc*nsteps,0,1,precision);
  hmax = max(hd);
  hmin = min(hd);
  dtop = max(max(Z));  % Top of grid
  HFACT = 0.1;
  rH = 4;
  
  ax = [0 L -dmax dtop];
  
  U = zeros(Nk,Nc);
  W = zeros(Nk,Nc);
  
  component = 1;
  numcomponents = 3;

  ud = getcdata(ufid,numcomponents*Nc*Nk,0,n,precision);
  sd = getcdata(sfid,Nc*Nk,0,n,precision);
  qd = getcdata(qfid,Nc*Nk,0,n,precision);
  hd = getcdata(hfid,Nc,0,n,precision);
  
  S = reshape(sd,Nc,Nk);
  S = S(is,end:-1:1)';
  
  Q = reshape(qd,Nc,Nk);
  Q = Q(is,end:-1:1)';
  Q(find(S==EMPTY))=nan;
  H = hd(is);
  
  S(find(S==EMPTY))=nan;
  
  for k=1:Nk
    component = 1;
    nstart = 1+numcomponents*Nc*(k-1)+(component-1)*Nc;
    nend = nstart+Nc-1;
    us = ud(nstart:nend)';
    U(Nk-k+1,:) = us(is);
    
    component = 3;
    nstart = 1+numcomponents*Nc*(k-1)+(component-1)*Nc;
    nend = nstart+Nc-1;
    us = ud(nstart:nend)';
    W(Nk-k+1,:) = us(is);
  end
  
  switch(PLOT)
   case 0
    pcolor(X,Z,Q);
   case 1
    pcolor(X,Z,S);
   case 2
    pcolor(X,Z,U);
   case 3
    pcolor(X,Z,W);
   otherwise
  end
  shading flat;
  
