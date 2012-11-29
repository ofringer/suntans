fontsize = 12;

L = 100;
k = pi/L;
D = 20;
delta = 5;
alpha = 0.99;
a = 1;
alpha0 = alpha;

deltarho = 1;
Nx = 200;
Nz = 20;
[x,z]=meshgrid(linspace(0,L,Nx),linspace(-D,0,Nz));

rho = -1/2*tanh(2*atanh(alpha)/delta*(x-L/2));

figure(1);
clf;
contour(x/L,z/D,rho,linspace(-alpha*deltarho/2,alpha*deltarho/2,7),'k-');
text(.25,-.5,'+ \Delta\rho/2\rho_0','fontsize',fontsize,'horizontalalignment','center');
text(.75,-.5,'- \Delta\rho/2\rho_0','fontsize',fontsize,'horizontalalignment','center');
colormap gray;
shading interp;
axis image;
axis([0 1 -1 0]);
daspect([.2 1 1]);
xlabel('x/L');
ylabel('z/D');
