%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: compare_results.m
%
% Description: Create comparison plot of SUNTANS cavity flow simulation
% with data from Zang et al. (1994) "A Non-staggered Grid, Fractional Step 
% Method for Time-Dependent Incompressible Navier-Stokes Equations in 
% Curvilinear Coordinates", J. Comp. Physics, 114, pg 18-33, Figure 6
%
% To compare to the XY simulations from "make testXY" or "make
% testXY-quad", set plot_type='XY', otherwise, to compare to results from
% "make testXZ" or "make testXZ-quad", set plot_type='XZ';
%
% Phillip Wolfram
% Stanford University
% 14 Nov 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath '../../../mfiles';

datadir='../data';
lit_datadir='../LitData';

% XY plot, otherwise XZ plot
plot_type='XY';

% load in cell center locations 
cdp = load([datadir,'/celldata.dat.0']);

% Account for the first column in the hybrid-grid format
[R,C]=size(cdp);
if(C==20)
  vx = cdp(:,1); 
  vy = cdp(:,2);
else
  vx = cdp(:,2); 
  vy = cdp(:,3);
end

% Vertical grid
dz = load([datadir,'/vertspace.dat']);
vz = getz(dz);

% Profiles for plotting use a boundary-layer refined Chebychev grid
N = 30; 
i_cheby = [1:N];
x_cheby = 0.5*(1-cos(i_cheby*pi/N));

% Set up profiles for interpolation.  Note that there are two
% profiles needed, i.e. u(y) and v(x) for the XY case and u(z) and
% w(x) for the XZ case.  
if(strcmp(plot_type,'XY'))
  xprof_1 = x_cheby;
  yprof_1 = 0.5*ones(1,N); 
  zprof_1 = -0.5*ones(1,N);

  xprof_2 = 0.5*ones(1,N); 
  yprof_2 = x_cheby;
  zprof_2 = -0.5*ones(1,N);
else
  xprof_1 = x_cheby;
  yprof_1 = 0.5*ones(1,N); 
  zprof_1 = -0.5*ones(1,N);

  xprof_2 = 0.5*ones(1,N); 
  yprof_2 = 0.5*ones(1,N);
  zprof_2 = -x_cheby(end:-1:1);
end

% Load velocity data, but check to make sure how many steps are
% available first.  This allows the mfile to be executed before the
% job comes to completion.
fid = fopen([datadir,'/u.dat.0'],'rb');
Nkmax = str2num(getvalue([datadir,'/suntans.dat'],'Nkmax'));
Nc = length(vx); 
arraysize = Nc*3*Nkmax;
nsteps = length(fread(fid,'float64'))/arraysize;
fclose(fid);

% Number of steps is known so can load latest data
fid = fopen([datadir,'/u.dat.0'],'rb');
fseek(fid,(nsteps-1)*8*arraysize,0);
data = fread(fid, arraysize,'float64'); 
fclose(fid);

data = reshape(data,Nc,3,Nkmax);
U = squeeze(data(:,1,:)); 
V = squeeze(data(:,2,:)); 
W = squeeze(data(:,3,:));

% Interpolate onto the profiles.  Note that the profiles need to be
% shifted to match the results in the Zang et al. paper, which
% forced the lid at the y=1 surface of an x-y box with u(y=1)=1.  In the present
% cases, the y=0 surface is forced with u(y=0)=+1 for the x-y box,
% while the x=0 surface is forced with w(x=0)=-1.
if(strcmp(plot_type,'XY'))
  vprof = griddata(vx,vy,squeeze(V),xprof_1,yprof_1);
  uprof = griddata(vx,vy,squeeze(U),xprof_2,yprof_2(end:-1:1));

  % v(x)
  xplot = 2*xprof_1-1;
  vplot = -vprof;

  % u(y)
  yplot = 2*yprof_2-1;
  uplot = uprof;
else
  [Vx,Vz]=ndgrid(vx,vz);
  uplot = griddata(Vx',Vz',U',xprof_2,zprof_2);
  wplot = griddata(Vx',Vz',W',xprof_1,zprof_1);

  % v(x)
  xplot = 2*(zprof_2+0.5);
  vplot = -uplot(end:-1:1);

  % u(y)
  yplot = 2*xprof_1-1;
  uplot = -wplot(end:-1:1);
end

% This is the data digitized from the Zang et al. paper.
load([lit_datadir,'/6diamonds']);
xdata_1 = Xdata; 
vdata_1 = Ydata;

load([lit_datadir,'/6circles']);
udata_2 = Ydata; 
ydata_2 = Xdata;

figure(1)
clf;
hold on;

ax1 = gca;
set(ax1, 'XColor','k','YColor','k');
axis([0 1 -1 1])
ax2 = axes('Position',get(ax1,'Position'), ...
    'XAxisLocation', 'top',...
    'YAxisLocation','right',...
    'Color','none');

axis(ax1);
line(xplot,vplot,'color','k');
line(2*xdata_1-1,vdata_1,'color','k','linestyle','none','marker','x');
axis([-1 1 -1 1]);

axis(ax2);
line(uplot,yplot,'color','k');
line(udata_2,2*ydata_2-1,'color','k','linestyle','none','marker','o');
axis([-1 1 -1 1]);

xlabel(ax1,'u/U')
xlabel(ax2,'x/B')
ylabel(ax1,'v/U')
ylabel(ax2,'y/D')
