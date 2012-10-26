%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: quadrid.m
% Description: Create grid files for a quad grid.
%
% Oliver Fringer
% Stanford University
% 18 October 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length and width of domain
L = 100;
W = 100;

% Number of cells
Nx = 10;
Ny = 10;

% Boundary condition types
WestBoundary = 1;
EastBoundary = 1;
NorthBoundary = 1;
SouthBoundary = 1;

dx = L/Nx;
dy = W/Ny;

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
  m=1;
  ycell = yv(n);
  for dxs=[-dx,dx]
    xcell = xv(n) + dxs;
    ind = find(xcell==xv & ycell==yv);
    if(~isempty(ind))
      neigh(n,m)=ind;
    end
    m=m+1;
  end
  xcell = xv(n);
  for dys=[-dy,dy]
    ycell = yv(n) + dys;
    ind = find(xcell==xv & ycell==yv);
    if(~isempty(ind))
      neigh(n,m)=ind;
    end
    m=m+1;
  end
end

figure(1);
clf;
hold on;
axis([-dx/2 L+dx/2 -dy/2 W+dy/2]);
axis off;

for n=1:Ne
  if(mark(n)==0)
    color='k.-';
  else
    color='r.-';
  end
  
  plot(xp(edges(n,:))',yp(edges(n,:))',color);
  if(grad(n,1)~=-1 & grad(n,2)~=-1)
    plot(xv(grad(n,:)),yv(grad(n,:)),'m-');
  end
end
  
for n=1:Nc
  for m=1:4
    if(neigh(n,m)~=-1)
      plot([xv(n),xv(neigh(n,m))],[yv(n),yv(neigh(n,m))],'b.-');
    end
  end
end
