%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: twod_equilateral.m
% Description: Create a 2-d planar straight line graph composed of 
% equilateral triangles on all edges, unlike twodgrid.m,
% which places right triangles at the east and west boundaries.
%
% Oliver Fringer
% Stanford University
% 18 Nov 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File name for the planar straight line graph
fname = 'twod.dat';

% Length and width of domain
L = 64000;
W0 = 64000;

% Number of cells along length
Nx = 50;

% Boundary condition types at specific boundaries
WestBoundary = 1;
NorthBoundary = 1;
EastBoundary = 1;
SouthBoundary = 1;

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
x1 = [0:dx:L];
x2 = [dx/2:dx:L-dx/2];
y1 = zeros(size(x1));
y2 = zeros(size(x2));
x = [x1,x2];
y = [y1,y2+W];
b_ind = [1:length(x1)];
marks = SouthBoundary*ones(size(b_ind));
left = find(x==min(x2) & y==max(y));
marks_left = WestBoundary*ones(size(left));
right = find(x==max(x2) & y==max(y));
marks_right = EastBoundary*ones(size(left));

for n=2:numrows,
  if(~rem(n,2))
    x_last = [length(x)+1:length(x)+length(x1)];
    x = [x,x1];
    y = [y,y1+n*W];
    left = [left,find(x==min(x1) & y==max(y))];
    marks_left = [marks_left,WestBoundary];
    right = [right,find(x==max(x1) & y==max(y))];
    marks_right = [marks_right,EastBoundary];
  else
    x_last = [length(x)+1:length(x)+length(x2)];
    x = [x,x2];
    y = [y,y2+n*W];
    left = [left,find(x==min(x2) & y==max(y))];
    marks_left = [marks_left,WestBoundary];
    right = [right,find(x==max(x2) & y==max(y))];
    marks_right = [marks_right,EastBoundary];
  end
end
b_ind = [b_ind,right,x_last(end-1:-1:2),left(end:-1:1)];
marks = [marks,marks_right,NorthBoundary*ones(1,length(x_last)-2),marks_left];

x = x + xoffset;
y = y + yoffset;

N = length(x);
ofile = fopen(fname,'w');

% Output in planar straigh line graph format for SUNTANS
fprintf(ofile,'%d\n',N);
for n=1:N, 
  fprintf(ofile,'%f %f 0\n',x(n),y(n));
end
fprintf(ofile,'%d\n',length(b_ind));
for n=1:length(b_ind)-1
  fprintf(ofile,'%d %d %d\n',b_ind(n)-1,b_ind(n+1)-1,marks(n));
end
fprintf(ofile,'%d %d %d\n',b_ind(end)-1,b_ind(1)-1,marks(end));
fprintf(ofile,'0\n');
fprintf(ofile,'%f\n',0);

fclose(ofile);
