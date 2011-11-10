%
% ELLIPSES
%  ELLIPSES(X,Y,UAMP,UPHASE,VAMP,VPHASE,NP,SCALE) will plot the tidal ellipses
%  at the locations specified by the vectors X and Y with the specified
%  amplitudes and phases.  NP is the number of points used to plot the
%  ellipse and SCALE is a plotting scaling factor.  See CHECK_SUNTIDES
%  for an example.

% 
% Oliver Fringer
% Stanford University
% 9 Sep 07
%
function ellipses(x,y,uamp,uphase,vamp,vphase,np,scale)

  dy = max(y)-min(y);
  us = uamp.*cos(uphase);
  vs = vamp.*cos(vphase);
  umax = max(sqrt(us.^2+vs.^2));

  if(scale~=0)
    uamp = scale*uamp/umax;
    vamp = scale*vamp/umax;
  end

  thetas = linspace(0,2*pi,np);
  N = length(x);
  for n=1:N
    xe0 = x(n) + uamp(n)*cos(uphase(n));
    ye0 = y(n) + vamp(n)*cos(vphase(n));
    xe1 = x(n) + uamp(n)*cos(uphase(n)+pi/2);
    ye1 = y(n) + vamp(n)*cos(vphase(n)+pi/2);
    xe = x(n) + uamp(n)*cos(uphase(n)+thetas);    
    ye = y(n) + vamp(n)*cos(vphase(n)+thetas);        
    plot([xe,xe(1)],[ye,ye(1)],'k-',[x(n),xe0],[y(n),ye0],'k-',[x(n),xe1],[y(n),ye1],'k-');
  end

  if(scale~=0)
    quiver(mean(x),mean(y),scale,0,0,'k-');
    text(mean(x),mean(y)-0.1*dy,sprintf('%.2f cm/s',umax));
  end