%
% HALF_CIRCLE
% HALF_CIRCLE creates the pslg for a rectangular domain with
%   a semi-circular bump at the southern boundary.  The pslg
%   is written into the file bump.dat.  Both east and west
%   boundaries are of type 2.
%
% Oliver Fringer
% Stanford University
% 30 Aug 05
%
fname = 'bump.dat';  % Output pslg
L = 50000;           % Domain length
W = 10000;           % Doman width
N = 20;              % Number of points in semi-circle
R = 5000;            % Radius of the semi-circle
amin = 5e5;          % Minimum area for triangle

WestBoundary = 2;    % Western boundary type
EastBoundary = 2;    % Eastern boundary type

thetas = [pi:-pi/N:0];
xp = L/2+R*cos(thetas);
yp = sqrt(R^2-(xp-L/2).^2);
y = [0,yp,0,W,W];
x = [0,xp,L,L,0];
marks = ones(size(x));
marks(end) = WestBoundary;
marks(end-2) = EastBoundary;

ofile = fopen(fname,'w');
fprintf(ofile,'%d\n',length(x));
for n=1:length(x)
  fprintf(ofile,'%d %d 0\n',x(n),y(n));
end
fprintf(ofile,'%d\n',length(x));
for n=1:length(x)-1
  fprintf(ofile,'%f %f %d\n',n-1,n,marks(n));
end
fprintf(ofile,'%f %f %d\n',length(x)-1,0,marks(end));
fprintf(ofile,'0\n');
fprintf(ofile,'%f\n',amin);
fclose(ofile);
