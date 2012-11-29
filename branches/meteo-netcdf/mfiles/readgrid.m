function [i,repeats,xall,yall,dall,z,triall,Nc,Nkall] = load_grid(dirname_g,numprocs)

  z = getz(load([dirname_g,'/vertspace.dat']));

  nn=1;
  for n=0:numprocs-1
    cp = load([dirname_g,'/cells.dat.',num2str(n)]);
    cdp = load([dirname_g,'/celldata.dat.',num2str(n)]);
    
    tri{nn} = cp(:,3:5);
    d{nn} = cdp(:,4);
    xv{nn} = cdp(:,1);
    yv{nn} = cdp(:,2);
    Nc(nn) = length(d{nn});
    Nk{nn} = cdp(:,5);
    
    nn=nn+1;
  end

  for n=0:numprocs-1
    if(n==0)
      xall = xv{1};
      yall = yv{1};
      dall = d{1};
      triall = tri{1};
      Nkall = Nk{1};
    else
      xall(end+1:end+Nc(n+1)) = xv{n+1};
      yall(end+1:end+Nc(n+1)) = yv{n+1};
      dall(end+1:end+Nc(n+1)) = d{n+1};
      triall(end+1:end+Nc(n+1),:) = tri{n+1};
      Nkall(end+1:end+Nc(n+1)) = Nk{n+1};
    end
  end

  R = sqrt(xall.^2+yall.^2);
  [dR,i] = sort(R);
  dR = dR(2:end)-dR(1:end-1);
  repeats = find(dR==0);
  xall(i(repeats)) = [];
  yall(i(repeats)) = [];
  dall(i(repeats)) = [];
  triall(i(repeats),:) = [];
  Nkall(i(repeats)) = [];
