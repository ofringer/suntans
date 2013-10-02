% $Id: nkenkc.m,v 1.1 2004-06-14 07:52:51 fringer Exp $
% $Log: not supported by cvs2svn $

open('nkenkc.fig');

Nkc = 5;

dz = .2;
dx = 2*dz;
z = [0:dz:(Nkc+1)*dz];
x = zeros(size(z));
z1 = z(3:end);
x1 = zeros(size(z1))+dx;
plot([x' x'+dx]',[z',z']','k',...
     [x' x'+dx]',[z',z']','k',...
     [x1' x1'+dx]',[z1',z1']','k',...
     [x1' x1'+dx]',[z1',z1']','k',...
     x,z,'k-',x+dx,z,'k-',x1+dx,z1,'k-');
axis image;
axis([-.1 .9 -.1 1.2*max(z)]);
axis off;

print -deps2 nkenkc