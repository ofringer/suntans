%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: plottransect.m
% Description: Plot a transect of suntans data.
%
% Oliver Fringer
% Stanford University
% 08/10/07
%
% This m-file plots a 2d x-z transect along a line defined by
% the end points (x1,y1) (x2,y2). 
%
% Instructions:
% 1) Set the desired end points (or specify an arbitrary set of
%    points that make up the transect in xplot,yplot).  The default
%    is to interpolate the data onto Nslice equispaced points on
%    the line with the given end points.  For multiscale grids it
%    is probably a good idea to adjust the spacing between points
%    appropriately.
% 2) Set the time step desired to plot in the data file
% 3) Set the number of processors.
% 4) Set the location of the data file in the variable data_directory
%    (i.e. the location of s.dat.*, u.data.*, etc...
% 5) Set the location of the grid files (cells.dat.*, celldata.dat.*, 
%    vertspace.dat) in the variable grid_directory.
% 6) Set the number of horizontal points that make up the transect.
% 7) To hide verbose output, set VERBOSE=0.
% 8) Change the 'readalldata' entry to change the variable being plotted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transect to plot from the point (x1,y1) to (x2,y2)
x1 = 0;
y1 = 0;
x2 = 1;
y2 = 1;

% Plot this time step in the given data file
nstep=1;            

% Number of processors
numprocs=1;           

% Location of data files
data_directory = './data';

% Directory where grid information is stored:
% vertspace.dat, cells.dat.*, celldata.dat.*
grid_directory = './grid';

% Number of points along transect
Nslice = 200;

% Display verbose information
VERBOSE = 1;           

% Read in the grid on all processors
if(~exist('xall','var'))
  if(VERBOSE)
    fprintf('Reading grid data...\n');
  end
  [ind,repeats,xall,yall,dall,z,triall,Nc,Nkall] = readgrid(grid_directory,numprocs);
  Nkmax = max(Nkall);
end

% Read in data from all processors
if(~exist('data','var'))

  % s-data:
  data = readalldata(data_directory,'s.dat',0,'s',numprocs,nstep,Nkmax,Nc,ind,repeats,VERBOSE);
  % u-data:
  %data = readalldata(data_directory,'u.dat',1,'u',numprocs,nstep,Nkmax,Nc,ind,repeats,VERBOSE);
  % v-data:
  % data = readalldata(data_directory,'u.dat',2,'u',numprocs,nstep,Nkmax,Nc,ind,repeats,VERBOSE);
  % w-data:
  % data = readalldata(data_directory,'u.dat',3,'u',numprocs,nstep,Nkmax,Nc,ind,repeats,VERBOSE);

  if(data==-1)
    return;
  end
else
  if(VERBOSE)
    fprintf(['Warning...data variable already exists. Use ''clear data'' '...
             'to reload the data to plot.\n']);
  end
end

% Compute the transects
xplot = linspace(x1,x2,Nslice);
yplot = linspace(y1,y2,Nslice);
if(VERBOSE)
  fprintf('Calculating slices...\n');
end
data_slice = getslice(xall,yall,data,xplot,yplot,Nkmax);
depth_slice = getslice(xall,yall,dall,xplot,yplot,Nkmax);
  
rplot = sqrt((xplot-xplot(1)).^2+(yplot-yplot(1)).^2);

[X0,Z0]=meshgrid(rplot,z);
X = X0';
Z = Z0';
Z(find(Z<-depth_slice))=nan;

% Create surface plots of the results.
figure(1);
clf;
pcolor(X,Z,data_slice);
shading flat;
hold on;
plot(X(:,end),-depth_slice(:,end),'k-');

