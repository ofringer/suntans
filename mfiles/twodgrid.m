%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: twogrid.m
% Description: Create a 2-d grid composed of equilateral triangles.
%
% Oliver Fringer
% Stanford University
% 27 July 04
%
% $Id: twodgrid.m,v 1.1 2004-07-27 21:51:12 fringer Exp $
% $Log: not supported by cvs2svn $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File name for the planar straight line graph
fname = '/home/data/suntans_data/zhonghua/twod.dat';

% Length and width of domain
L = 50000;
W0 = 200000;

% Number of cells along length
Nx = 20;

% Boundary condition types at specific boundaries
WestBoundary = 2;
NorthBoundary = 2;
EastBoundary = 2;
SouthBoundary = 2;

% Offset of lower leftmost cell
xoffset = 0;
yoffset = 0;

% Base length
dx = L/(Nx-1);

% Triangle height
W = dx/2*tan(pi/3);

% Number of rows determined by triangle height
numrows = fix(W0/W);

% Create the triangle vertices along lines
x1 = [-dx/2:dx:Nx*dx];
x2 = [Nx*dx-dx/2,Nx*dx-dx:-dx:0];
y0 = ones(size(x1));
x = [x1,x2];
y = [0*y0,W*y0];
for n=1:numrows-1,
  if(mod(n,2))
    x = [x,x1];
  else
    x = [x,x2];
  end
  y = [y,(n+1)*W*y0];
end

x = x + xoffset;
y = y + yoffset;

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

N = length(x);
ofile = fopen(fname,'w');

% Output in planar straigh line graph format for SUNTANS
fprintf(ofile,'%d\n',N+4);
for n=1:N, 
  fprintf(ofile,'%f %f 0\n',x(n),y(n));
end
fprintf(ofile,'%f %f 0\n',xmin,ymin);
fprintf(ofile,'%f %f 0\n',xmin,ymax);
fprintf(ofile,'%f %f 0\n',xmax,ymax);
fprintf(ofile,'%f %f 0\n',xmax,ymin);
fprintf(ofile,'4\n');
fprintf(ofile,'%d %d %d\n',N,N+1,WestBoundary);
fprintf(ofile,'%d %d %d\n',N+1,N+2,SouthBoundary);
fprintf(ofile,'%d %d %d\n',N+2,N+3,EastBoundary);
fprintf(ofile,'%d %d %d\n',N+3,N,NorthBoundary);
fprintf(ofile,'0\n');
fprintf(ofile,'%f\n',0);

fclose(ofile);