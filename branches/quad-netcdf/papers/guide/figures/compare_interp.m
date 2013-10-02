dt = 90;
h1 = fread(fopen('h1.dat','rb'),'float64');
h2 = fread(fopen('h2.dat','rb'),'float64');
h3 = fread(fopen('h3.dat','rb'),'float64');

n = length(h3);

Tday = 86400;
tstart = 7*Tday;
t = tstart+[dt:dt:n*dt];

figure(1)
set(gca,'fontsize',14);
plot(t/Tday,h1,'r-',t/Tday,h2,'b-',t/Tday,h3,'k-');
axis tight;
xlabel('Days in 2006');
ylabel('h (m)');
daspect([4 1 1]);

