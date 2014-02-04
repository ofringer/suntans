function slicedata = getslice(xall,yall,data,xplot,yplot,Nkmax)

  Nslice = length(xplot);
  Nnear = 3;
  Ntotal = max(size(xall));
  slicedata = zeros(Nslice,Nkmax);

  ind = [1:Ntotal];

  inds = zeros(length(xplot),Nnear);
  coefs = zeros(length(xplot),Nnear);
  for n=1:length(xplot)
    inds(n,:) = nearestx(xall,yall,xplot(n),yplot(n),Nnear);
  end

  for n=1:length(xplot)
    A = [xall(inds(n,:)) , yall(inds(n,:)) , ones(3,1) ];
    f0 = data(inds(n,:),:);
    coefs = inv(A)*f0;
    slicedata(n,:) = coefs(1,:)*xplot(n)+coefs(2,:)*yplot(n)+coefs(3,:);
  end

