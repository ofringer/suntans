fontsize = 12;

L = 100;
k = pi/L;
D = 100;
delta = 20;
alpha = 0.99;
a = 1;
alpha0 = alpha;

deltarho = 1;
Nx = 100;
Nz = 100;
[x,z]=meshgrid(linspace(0,L,Nx),linspace(-D,0,Nz));

rho = -1/2*tanh(2*atanh(alpha)/delta*(z+D/2-a*cos(k*x)));

figure(1);
clf;

ax1 = axes('position',[.2 .3 .4 .4],'fontsize',fontsize);
contour(x/L,z/D,rho,linspace(-alpha*deltarho/2,alpha*deltarho/2,7),'k-');
text(.5,-.75,'+ \Delta\rho/2\rho_0','fontsize',fontsize,'horizontalalignment','center');
text(.5,-.25,'- \Delta\rho/2\rho_0','fontsize',fontsize,'horizontalalignment','center');
colormap gray;
shading interp;
axis image;
axis([0 1 -1 0]);
xlabel('x/L');
ylabel('z/D');
title('(a)');

ax2 = axes('position',[.6 .3 .2 .4],'fontsize',fontsize);
plot(rho(:,Nx/2),z(:,1)/D,'k-',[0 0],[-1 0],'k--');
set(gca,'xtick',[-.5 0 .5]);
set(gca,'ytick',[]);
axis([-.55 .55 -1 0]);
xlabel('\rho/\Delta\rho(x=L/2)');
title('(b)');
