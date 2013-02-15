function [xvs,yvs,xps,yps]=stretchgrid(xpg,ypg,Lr,Nr,rmax)

[Nx,Ny]=size(xpg);
Nx = Nx-1;
Ny = Ny-1;

if(Nr>=Nx)
  error(sprintf('Nr=%d (refined region) must satisfy Nr<=Nx (=%d)\n',Nr,Nx));
end

L = max(max(xpg));
Ls = (L-Lr)/2;
alpha = Ls*(Nr/Lr);

Ns = Nx - Nr;
if(rem(Ns,2)==0)
  Ns = Ns/2;
else
  Nr = Nr+1;
  Ns = (Nx-Nr)/2;
end
dxr = Lr/Nr;

fsolveopts = optimset('Diagnostics','off');
r = fsolve(@(r,Ns,alpha) r^Ns - alpha*r - 1 + alpha,1.1,...
	   fsolveopts,Ns,alpha);
if(r<1 | r>rmax)
  error(sprintf('Stretching factor is r=%.3f not in range [1,%.2f]',r,rmax));
else
  fprintf('Stretching factor is %.2f\n',r);
end

xps(Ns+1:Ns+Nr+1,:) = L/2+[-Lr/2:dxr:Lr/2]'*ones(1,Ny+1);
dx=dxr;
for j=Ns+1:-1:2
  xps(j-1,:) = (xps(j)-dx)*ones(1,Ny+1);
  dx=r*dx;
end

dx=dxr;
for j=Ns+Nr+1:Nx
  xps(j+1,:) = (xps(j)+dx)*ones(1,Ny+1);
  dx=r*dx;
end

yps = ypg;
xvs = 0.25*(xps(1:end-1,1:end-1)+xps(2:end,1:end-1)+...
	    xps(1:end-1,2:end)+xps(2:end,2:end));
yvs = 0.25*(yps(1:end-1,1:end-1)+yps(2:end,1:end-1)+...
	    yps(1:end-1,2:end)+yps(2:end,2:end));

xps = xps(:);
yps = yps(:);
xvs = xvs(:);
yvs = yvs(:);

