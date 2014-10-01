addpath '../../../mfiles';

datadir='../data';
c = load([datadir,'/cells.dat']);
dz = load([datadir,'/vertspace.dat']);
z = getz(dz);

[Nc,cols] = size(c);
if(cols==8)
  xv = c(:,1);
  yv = c(:,2);
end

nsteps=20;
Nkmax=20;

ufid=fopen([datadir,'/u.dat'],'rb');
udata = reshape(fread(ufid,'float64'),Nc,Nkmax,3,nsteps+1);
u = squeeze(udata(:,:,1,:));
v = squeeze(udata(:,:,2,:));
w = squeeze(udata(:,:,3,:));

nplot = 21;
up = squeeze(u(:,:,nplot));
vp = squeeze(v(:,:,nplot));
wp = squeeze(w(:,:,nplot));

clear x0 y0 z0

Np = 10;
x0 = 0.2*ones(Np,1);
y0 = 0.25*ones(Np,1);
z0 = linspace(-0.24,-0.1,Np)';

X = xv*ones(1,Nkmax);
Y = yv*ones(1,Nkmax);
Z = ones(Nc,1)*z';

m=1;
dt=0.01;
ni = 20;
while(x0(p,m)<2)
  for p=1:Np
    Rh = (xv-x0(p,m)).^2+(yv-y0(p,m)).^2;
    [rh,ih] = sort(Rh);
    Rv = (z0(p,m)-z).^2;
    [rv,is] = sort(Rv);
    
    k1 = is(1);
    k2 = is(2);
    z1 = z(k1);
    z2 = z(k2);
    
    u1 = griddata(xv(ih(1:ni)),yv(ih(1:ni)),up(ih(1:ni),k1),x0(p,m),y0(p,m));
    u2 = griddata(xv(ih(1:ni)),yv(ih(1:ni)),up(ih(1:ni),k2),x0(p,m),y0(p,m));
    ui = u1 + (z0(p,m)-z1)*(u2-u1)/(z2-z1);
    
    v1 = griddata(xv(ih(1:ni)),yv(ih(1:ni)),vp(ih(1:ni),k1),x0(p,m),y0(p,m));
    v2 = griddata(xv(ih(1:ni)),yv(ih(1:ni)),vp(ih(1:ni),k2),x0(p,m),y0(p,m));
    vi = v1 + (z0(p,m)-z1)*(v2-v1)/(z2-z1);
    
    w1 = griddata(xv(ih(1:ni)),yv(ih(1:ni)),wp(ih(1:ni),k1),x0(p,m),y0(p,m));
    w2 = griddata(xv(ih(1:ni)),yv(ih(1:ni)),wp(ih(1:ni),k2),x0(p,m),y0(p,m));
    wi = w1 + (z0(p,m)-z1)*(w2-w1)/(z2-z1);
    
    x0(p,m+1) = x0(p,m) + dt*ui;
    y0(p,m+1) = y0(p,m) + dt*vi;
    z0(p,m+1) = z0(p,m) + dt*wi;
  end
  m=m+1;
end

figure(1);
clf;
plot3(x0',y0',z0');
axis image;
grid;
axis([0 2 0 0.4 -0.25 0]);
