p = load('plumegrid/points.dat');
c = load('plumegrid/cells.dat');
e = load('plumegrid/edges.dat');

x = p(:,1)/1000;
y = p(:,2)/1000;
ip = 1+c(:,[3:5,3]);
ie = 1+e(:,[1,2]);
ig = find(e(:,3)==2);


figure(1);
clf;
hold on;
plot(x(ip)',y(ip)','k-');
plot(x(ie(ig,:))',y(ie(ig,:))','k-','linewidth',5);
set(gca,'fontsize',12);
axis image;
axis([-.005 3.05 -.001 1.001]);
xlabel('x (km)')
ylabel('y (km)');
