%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: sunsurf.m
% Description: Create a surface (x-y) plot of 2d unstructured suntans data.
%
% Oliver Fringer
% Stanford University
% 30 Sep 04
%
% $Id: sunsurf.m,v 1.1 2004-09-30 21:37:17 fringer Exp $
% $Log: not supported by cvs2svn $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flag = sunsurf(PLOT,datadir,n,klevel,procs,varargin)

  flag = 1;
  EMPTY=999999;         % Empty cells are defined by this
  dbl=8;                % Size of double precision in bytes
  precision='float64';  % Precision for reading in data
  
  % cellcentered data contains the voronoi points and the depths
  % at those points.
  N=length(procs);
  
  held=ishold;
  if(N>1)
    if(~held)
      hold on;
    end
  end
  for proc=procs
    
    switch(PLOT)
     case 'q'
      filename = [datadir,'/q.dat.',num2str(proc)];
     case 's'
      filename = [datadir,'/s.dat.',num2str(proc)];
     case {'u','v','w'}
      filename = [datadir,'/u.dat.',num2str(proc)];
     case 'h'
      filename = [datadir,'/fs.dat.',num2str(proc)];
     case 'T'
      filename = [datadir,'/T.dat.',num2str(proc)];
     otherwise
    end

    fid = fopen(filename,'rb');
    if(fid<=0)
      fprintf('Error opening %s\n',filename);
      flag = 0;
      return;
    end
    
    pd = load([datadir,'/points.dat']);
    cp = load([datadir,'/cells.dat.',num2str(proc)]);
    cdp = load([datadir,'/celldata.dat.',num2str(proc)]);
    
    xp = pd(:,1);
    yp = pd(:,2);
    tri = cp(:,3:5);
    d = cdp(:,4);
    xv = cdp(:,1);
    yv = cdp(:,2);
    dz = load([datadir,'/vertspace.dat']);
    
    % Total number of cells in the horizontal Nc and vertical Nk
    Nc = length(xv);
    Nk = length(dz);
    z = -(sum(dz(1:klevel))-dz(klevel)/2)*ones(size(xv));
    
    switch(PLOT)
     case 'q'
      arraysize = Nc*Nk;
     case 's'
      arraysize = Nc*Nk;
     case {'u','v','w'}
      arraysize = Nc*Nk;
     case 'h'
      arraysize = Nc;
     case 'T'
      arraysize = Nc*Nk;
     otherwise
    end
    
    switch(PLOT)
     case 'u'
      fseek(fid,8*(3*arraysize*(n-1)+3*(klevel-1)*Nc),0);
     case 'v'
      fseek(fid,8*(3*arraysize*(n-1)+3*(klevel-1)*Nc+Nc),0);
     case 'w'
      fseek(fid,8*(3*arraysize*(n-1)+3*(klevel-1)*Nc+2*Nc),0);
     otherwise
      fseek(fid,8*(arraysize*(n-1)+(klevel-1)*Nc),0);
    end
    phi = fread(fid,Nc,'float64');
    
    phi(find(z<-d))=nan;
    
    unsurf(tri,xp,yp,phi,'edgecolor','none');
    
  end
  
  if(~held)
    hold off;
  end