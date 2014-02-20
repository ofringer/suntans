%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: onedrid.m
% Description: Create a 1-d planar straight line graph composed of 
% equilateral triangles.
%
% Oliver Fringer
% Stanford University
% 27 July 04
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File name for the planar straight line graph
fname = 'oned.dat';

% Length of domain
L = 100000;

% Number of cells
Nx = 200;

% Boundary condition types
WestBoundary = 2;
EastBoundary = 1;

% Triangle height
W = L/Nx*tan(pi/3);

% Offset of lower leftmost cell
yoffset=0;
xoffset=0;

dx = 2*W/tan(pi/3);
x1 = [-dx/2:dx:L];
x2 = [0:dx:L+dx/2];

% Uncomment these to make the end triangles right triangles.
%x1(end)=L;
%x1(1)=0;

x = [x1,x2(end:-1:1)];
y = [ones(size(x1))*W,zeros(size(x2))];
x=x+xoffset;
y=y+yoffset;
amin = dx*W;

N = length(x);
ofile = fopen(fname,'w');

N = length(x);

fprintf(ofile,'%d\n',N);
for n=1:N,
  fprintf(ofile,'%f %f 0\n',x(n),y(n));
end
fprintf(ofile,'%d\n',N);
for n=1:N-1,
  if(n==N/2)
    fprintf(ofile,'%d %d %d\n',n-1,n,EastBoundary);
  else
    fprintf(ofile,'%d %d 0\n',n-1,n);
  end
end
fprintf(ofile,'%d %d %d\n',N-1,0,WestBoundary);
fprintf(ofile,'0\n');
fprintf(ofile,'%f\n',amin);

fclose(ofile);


