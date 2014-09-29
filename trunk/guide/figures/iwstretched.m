Nk = 100;
Nc = 100;
L = 100000;
D0 = 3000;
Ds = 500;
Ls = 20000;
xmid = 65000;

dx = L/Nc;
dz = load('iwstretched.dat');
z = cumsum(dz);
z = -z(end:-1:1);
x = [0:dx:L];
[x,z]=meshgrid(x,z);

d = D0 - (D0-Ds)*((x-xmid)/Ls+1/2);
d(find(x<=xmid-Ls/2))=D0;
d(find(x>xmid+Ls/2))=Ds;

iskip = 1;
kskip = 1;
is = [1:iskip:Nc];
ks = [1:kskip:Nk];

z(find(z<-d-1))=nan;
plot(x(is,ks)/1000,z(is,ks),'k-',x(is,ks)'/1000,z(is,ks)','k-');
axis image;
daspect([1 50 1]);
xlabel('x (km)');
ylabel('z (m)');
