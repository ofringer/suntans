%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: plotslice.m
% Description: Extract data for a surface plot of 2d suntans data.
%
% [x,z,data]=plotslice(plottype,dirname,timestep,processor);
% 
% plottype: 
%    'q': nonhydrostatic pressure
%    's': salinity
%    'u': u-velocity
%    'w': w-velocity
%    's0': background salinity
%    'h': free surface
%    'nut': eddy-viscosity
%    'kappat': scalar-diffusivity
%
% Example: [x,z,data]=plotslice(1,'/home/data/suntans_data',10,0);
% will extract the salinity data at output step 10
% for processor 0 and will return the x,z,data matrices for that plot,
% which can then be used with:  pcolor(x,z,data);
% If the free surface is selected, it can be viewed with
% plot(x(1,:),data);  
%
% Oliver Fringer
% Stanford University
% 15 Jun 04
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,z,data] = plotslice(PLOT,datadir,n,proc)

  EMPTY=999999;         % Empty cells are defined by this
  precision='float64';  % Precision for reading in data
  
  % cellcentered data contains the voronoi points and the depths
  % at those points.
  cellcentereddata=load([datadir,'/celldata.dat.',num2str(proc)],'-ascii');
  xv = cellcentereddata(:,1);
  yv = cellcentereddata(:,2);
  dv = cellcentereddata(:,4);
  dz = load([datadir,'/vertspace.dat']);
  
  % Total number of cells in the horizontal Nc and vertical Nk
  Nc = length(xv);
  Nk = length(dz);
  
  % Length and depth of domain
  L = max(xv);
  dmax = max(dv);
  
  % Set up the Cartesian grid
  z = -cumsum(dz);
  [xs,is]=sort(xv);
  [x,z]=meshgrid(xs,z);

  % Empty the cells below the bottom
  D = ones(Nk,1)*dv(is)';
  z(find(z<-D))=nan;
  dv = dv(is);

  % Open up file descriptors for binary files
  switch(PLOT)
   case 'q'
    file = [datadir,'/q.dat.',num2str(proc)];
   case 's'
    file = [datadir,'/s.dat.',num2str(proc)];    
   case {'u','w'}
    file = [datadir,'/u.dat.',num2str(proc)];    
   case 's0'
    file = [datadir,'/s0.dat.',num2str(proc)];    
   case 'h'
    file = [datadir,'/fs.dat.',num2str(proc)];    
   case 'nut'
    file = [datadir,'/nut.dat.',num2str(proc)];    
   case 'kappat'
    file = [datadir,'/kappat.dat.',num2str(proc)];    
   otherwise
    fprintf('Unrecognized plot variable.\n');
    fprintf('Use one of ''q'',''s'',''u'',''w'',''s0'',''h'',''nut'',''kappat''.\n');
    data=zeros(Nk,Nc);
    return;
  end

  fid = fopen(file,'rb');

  switch(PLOT)
   case {'u','w'}
    data = reshape(getcdata(fid,Nk*3*Nc,n,precision),Nk,3,Nc);
    if(PLOT=='u')
      data = squeeze(data(:,1,:));
    else
      data = squeeze(data(:,3,:));      
    end
    data = data(is,:)';
    data(find(data==EMPTY))=nan;
   case 'h'
    data = getcdata(fid,Nc,n,precision);
    data = data(is);
   otherwise
    data = reshape(getcdata(fid,Nc*Nk,n,precision),Nc,Nk);
    data = data(is,:)';
    data(find(data==EMPTY))=nan;
  end
  
