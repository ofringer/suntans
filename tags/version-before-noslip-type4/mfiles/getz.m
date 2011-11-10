function z = getz(dz)
  
  N = length(dz);
  dzh = 0.5*(dz(1:end-1)+dz(2:end));
  
  z = zeros(N,1);
  z(1) = -dz(1)/2;
  z(2:N) = z(1)-cumsum(dzh);