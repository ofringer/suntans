%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: quadrid.m
% Description: Create grid files for a quad grid.
%
% Oliver Fringer
% Stanford University
% 18 October 2012
%
% may have error when dx & dy is very small
% add ampfactor 
% revised by Yun Zhang
% 2/12/2013 @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length and width of domain
L = 200000;
W = 1732.1;

% Number of cells
Nx = 100;
Ny = 1;
ampfacx=Nx/L;
ampfacy=Ny/W;
L=Nx;
W=Ny;
% Boundary condition types
WestBoundary = 2;
EastBoundary = 1;
NorthBoundary = 1;
SouthBoundary = 1;

dx = 1;
dy = 1;


N = (Nx+1)*(Ny+1);
[xp,yp] = ndgrid([0:dx:L],[0:dy:W]);
xp = xp(:);
yp = yp(:);
Np = length(xp);

[xv,yv] = ndgrid([dx/2:dx:L-dx/2],[dy/2:dy:W-dy/2]);
xv = xv(:);
yv = yv(:);
Nc = length(xv);
 
cells = zeros(Nc,4);
for n=1:Nc
  xcell = xv(n) + [dx/2,dx/2,-dx/2,-dx/2];
  ycell = yv(n) + [-dy/2,dy/2,dy/2,-dy/2];
  for m=1:4
    cells(n,m) = find(xp==xcell(m) & yp==ycell(m));
  end
end  

[xeu,yeu] = ndgrid([0:dx:L],[dy/2:dy:W-dy/2]);
[xev,yev] = ndgrid([dx/2:dx:L-dx/2],[0:dy:W]);
xeu = xeu(:);
yeu = yeu(:);
xev = xev(:);
yev = yev(:);

Ne = Nx*(Ny+1) + Ny*(Nx+1);
mark = zeros(Ne,1);
edges = zeros(Ne,2);
grad = -ones(Ne,2);

k=1;
for n=1:length(xeu)
  edges(k,1) = find(xp==xeu(n) & yp==yeu(n)-0.5*dy);
  edges(k,2) = find(xp==xeu(n) & yp==yeu(n)+0.5*dy);
  if(xeu(n)==0)
    mark(k)=WestBoundary;
  elseif(xeu(n)==L)
    mark(k)=EastBoundary;
  else
    mark(k)=0;
  end

  xv1 = xeu(n)-0.5*dx;
  xv2 = xeu(n)+0.5*dx;
  yv1 = yeu(n);
  yv2 = yeu(n);
  ind1 = find(xv1==xv & yv1==yv);
  ind2 = find(xv2==xv & yv2==yv);
  if(~isempty(ind1))
    grad(k,1) = ind1;
  end
  if(~isempty(ind2))
    grad(k,2) = ind2;
  end

  k=k+1;
end
for n=1:length(xev)
  edges(k,1) = find(xp==xev(n)-0.5*dx & yp==yev(n));
  edges(k,2) = find(xp==xev(n)+0.5*dx & yp==yev(n));
  if(yev(n)==0)
    mark(k)=SouthBoundary;
  elseif(yev(n)==W)
    mark(k)=NorthBoundary;
  else
    mark(k)=0;
  end

  xv1 = xev(n);
  xv2 = xev(n);
  yv1 = yev(n)+0.5*dy;
  yv2 = yev(n)-0.5*dy;
  ind1 = find(xv1==xv & yv1==yv);
  ind2 = find(xv2==xv & yv2==yv);
  if(~isempty(ind1))
    grad(k,1) = ind1;
  end
  if(~isempty(ind2))
    grad(k,2) = ind2;
  end

  k=k+1;
end

neigh = -ones(Nc,4);
for n=1:Nc
  % clockwise get neigh to make we can get share node for each boundary
  ycell = yv(n);
  xcell = xv(n) - dx;
  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,1)=ind;
  end
  xcell = xv(n);
  ycell = yv(n) + dy;
  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,2)=ind;
  end
  ycell = yv(n);
  xcell = xv(n) + dx;
  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,3)=ind;
  end
  xcell = xv(n);
  ycell = yv(n) - dy;
  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,4)=ind;
  end
end
% 
% figure(1);
% clf;
% hold on;
% axis([-dx/2 L+dx/2 -dy/2 W+dy/2]);
% axis off;
% 
% for n=1:Ne
%   if(mark(n)==0)
%     color='k.-';
%   else
%     color='r.-';
%   end
%   
%   plot(xp(edges(n,:))',yp(edges(n,:))',color);
%   if(grad(n,1)~=-1 & grad(n,2)~=-1)
%     plot(xv(grad(n,:)),yv(grad(n,:)),'m-');
%   end
% end
%   
% for n=1:Nc
%   for m=1:4
%     if(neigh(n,m)~=-1)
%       plot([xv(n),xv(neigh(n,m))],[yv(n),yv(neigh(n,m))],'b.-');
%     end
%   end
% end

% output
celloutput=zeros(Nc,11);
celloutput(:,1)=4;
celloutput(:,2)=xv/ampfacx;
celloutput(:,3)=yv/ampfacy;
celloutput(:,4:7)=cells;
celloutput(:,8:11)=neigh;

edgeoutput=zeros(Ne,5);
%find open BC to make sure grad[2*j]~=-1
loc=find(mark~=0 & mark~=1);
edgeoutput(:,1:2)=edges(:,1:2);
edgeoutput(:,3)=mark;
edgeoutput(:,4:5)=grad(:,1:2);
for i=1:length(loc)
    if grad(i,1)==-1
        edgeoutput(loc,5)=-1;
        edgeoutput(loc,4)=grad(i,2);
    end
end

pointoutput=zeros(Np,3);
pointoutput(:,1)=xp/ampfacx;
pointoutput(:,2)=yp/ampfacy;

%check results (should be edge length)
dis=3*ones(Ne,1);
for i=1:Ne
    nc1=grad(i,1);
    nc2=grad(i,2);
    if(nc1~=-1 & nc2~=-1)
        dis(i)=((xv(nc1)-xv(nc2))^2+(yv(nc1)-yv(nc2))^2)^0.5;
    end
end



celloutput(:,4:7)=celloutput(:,4:7)-1;
for i=8:11
loc=find(celloutput(:,i)~=-1);
celloutput(loc,i)=celloutput(loc,i)-1;
end

edgeoutput(:,1:2)=edgeoutput(:,1:2)-1;

for i=4:5
loc=find(edgeoutput(:,i)~=-1);
edgeoutput(loc,i)=edgeoutput(loc,i)-1;
end


datadir='.';
cells_file = [datadir,'/cells.dat'];
cellsf = fopen(cells_file,'w');
fprintf(cellsf, '%1.0f %12.10e %12.10e %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f\n', celloutput');
status = fclose(cellsf);   


datadir='.';
points_file = [datadir,'/points.dat'];
pointsf = fopen(points_file,'w');
fprintf(pointsf, '%12.10e %12.10e %1.0f\n', pointoutput');
status = fclose(pointsf);   


datadir='.';
edges_file = [datadir,'/edges.dat'];
edgesf = fopen(edges_file,'w');
fprintf(edgesf, '%1.0f %1.0f %1.0f %1.0f %1.0f\n', edgeoutput');
status = fclose(edgesf);   




