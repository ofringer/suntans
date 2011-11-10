function ind = nearestx(x,y,x0,y0,N)
  
  Np = length(x);
  dist = sqrt((x-x0).^2+(y-y0).^2);
  [dist,is]=sort(dist);
  ind = is(1:N)';

