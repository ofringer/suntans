%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: sunmovie.m
% Description: Create a movie of the 2-d data created in suntans
% in the directory described by datadir.
%
% Oliver Fringer
% Stanford University
% 1 Jun 04
%
% $Id: sunmovie.m,v 1.3 2005-07-27 17:55:18 fringer Exp $
% $Log: not supported by cvs2svn $
% Revision 1.2  2004/06/01 04:24:56  fringer
% Changed daspect([1 3*dmax/L 1]) to daspect([1 5*dmax/L 1]) so that
% a 5:1 aspect ratio will be 1 to 1.
%
% Revision 1.1  2004/06/01 01:56:50  fringer
% This m-file takes a directory containining 2-d suntans (x-z) data
% and plots the output as a surface plot, then creates a movie.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datadir = '../main/examples/iwaves/data';

% Surface plot determined by:
% PLOT: 0: q, 1: s, 2: u, 3: w
PLOT=3;

% Whether or not to create a movie
% Requires the linux convert command and printmovie.m.
MOVIE=1;

% Show the free-surface profile
FS=1;

EMPTY=999999;         % Empty cells are defined by this
dbl=8;                % Size of double precision in bytes
precision='float64';  % Precision for reading in data

% cellcentered data contains the voronoi points and the depths
% at those points.
cellcentereddata=load([datadir,'/celldata.dat.0'],'-ascii');
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
z = -dmax+cumsum(dz(end:-1:1));
[xs,is]=sort(xv);
[X,Z]=meshgrid(xs,z);

% Empty the cells below the bottom
D = ones(Nk,1)*dv(is)';
Z(find(Z<-D))=nan;
dv = dv(is);

% Open up file descriptors for binary files
qfile = [datadir,'/q.dat.0'];
qfid = fopen(qfile,'rb');
ufile = [datadir,'/u.dat.0'];
ufid = fopen(ufile,'rb');
sfile = [datadir,'/s.dat.0'];
sfid = fopen(sfile,'rb');
hfile = [datadir,'/fs.dat.0'];
hfid = fopen(hfile,'rb');

% Total number of time steps to plot is obtained by checking
% size of file
status=fseek(sfid,0,'eof');
nsteps = ftell(sfid)/Nk/Nc/dbl;
status=fseek(sfid,0,'bof');

% Determine hmax for plotting
hd = getcdata(hfid,Nc*nsteps,1,precision);
hmax = max(hd);
hmin = min(hd);
dtop = max(max(Z));  % Top of grid
HFACT = 0.1;
rH = 4;

% For axes scaling
if(FS)
  ax = [0 L -dmax dtop+(1+rH)*HFACT*(dmax-dtop)];
else
  ax = [0 L -dmax dtop];
end

figure(1);
clf;
axis(ax);
axis image;
daspect([1 5*dmax/L 1]);
axis off;
hold on;

U = zeros(Nk,Nc);
W = zeros(Nk,Nc);

component = 1;
numcomponents = 3;
for n=1:nsteps,
  fprintf('On %d of %d\n',n,nsteps);
  ud = getcdata(ufid,numcomponents*Nc*Nk,n,precision);
  sd = getcdata(sfid,Nc*Nk,n,precision);
  qd = getcdata(qfid,Nc*Nk,n,precision);
  hd = getcdata(hfid,Nc,n,precision);

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
  
  cla;
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
  if(FS)
    plot(xv(is),dtop+HFACT*(dmax-dtop)+rH*HFACT*(dmax-dtop)*(H-hmin)/(hmax-hmin),'k-');
  end
  shading flat;
  
  pause(0);
  if(MOVIE)
    if(~exist('movies','dir'))
      qq=input('Movies directory does not exist.  Create? (y or n): ','s');
      switch(qq)
       case 'y'
        mkdir movies
       otherwise
        fprintf('Need to create the directory movies to create a movie.\n');
        return;
      end
    end
    printmovie(n,'movies/movie',1,50)
  end
end
