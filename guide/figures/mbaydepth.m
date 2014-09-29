addpath ../../../../mfiles

p = load('mbaygrid/points.dat');
c = load('mbaygrid/cells.dat');
d = load('mbaygrid/mbay_bathy.dat-voro');

t = c(:,[3:5]);
x = p(:,1)/1000;
y = p(:,2)/1000;

figure(1)
unsurf(t,x,y,d(:,3));
colorbar
axis image;
set(gca,'box','on');
xlabel('x (km)');
ylabel('y (km)');
