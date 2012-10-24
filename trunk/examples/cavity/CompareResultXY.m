%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: MakePlots.m
% Description: Create comparison plot of SUNTANS cavity flow simulation
% with data from Zang et al. (1994) "A Non-staggered Grid, Fractional Step 
% Method for Time-Dependent Incompressible Navier-Stokes Equations in 
% Curvilinear Coordinates", J. Comp. Physics, 114, pg 18-33, Figure 6
%
% Phillip Wolfram
% Stanford University
% 14 Nov 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CompareResultXY
%% this function assumes that there is only one processor and is not general
clear all; close all; clc;

% XY plot, otherwise XZ plot
XY = 1;

% load in cell center locations 
cdp = load('./dataXY/celldata.dat.0'); vx = cdp(:,1); vy = cdp(:,2);
z = getz(load('./dataXY/vertspace.dat')); dz = max(z(1:end-1)-z(2:end));
zmax = max(z); zmin = min(z);

% setup profiles
N = 70; dx=1/N;
xprof = linspace(dx/2,1-dx/2,N);
yprof = 0.5*ones(1,N); 
zprof = 0.5*(zmin+zmax)*ones(1,N);
Iset1 = [1:N];

xprof = [xprof dx/2];
yprof = [yprof yprof(end)];
zprof = [zprof 0.5*(zmin+zmax)];

xprof = [xprof 0.5*ones(1,N)];
yprof = [yprof linspace(dx/2,1-dx/2,N)];
zprof = [zprof 0.5*(zmin+zmax)*ones(1,N)];
Iset2 = [N+2:length(xprof)];

% load in velocity information corresponding to (vx, vy)
fid = fopen('./dataXY/u.dat.0','rb');
Nkmax = 1; Nc = length(vx); arraysize = Nc*3*Nkmax;
fseek(fid,8*arraysize,0);
data = fread(fid, arraysize,'float64'); 
fclose(fid);

data = reshape(data,Nc,3,Nkmax);
U = squeeze(data(:,1,:)); 
V = squeeze(data(:,2,:)); 
W = squeeze(data(:,3,:));

% interpolate values onto desired profiles
if(XY)
    F = TriScatteredInterp(vx, vy, U); G =  TriScatteredInterp(vx, vy, V);
    uprof = F(xprof,yprof); vprof = G(xprof,yprof);
else
    % form xy interpolant over the data for the entire level
    for level = 1: Nkmax
        F = TriScatteredInterp(vx, vy, U(:,level));
        G = TriScatteredInterp(vx, vy, V(:,level));
        H = TriScatteredInterp(vx, vy, W(:,level));
        Utempdata(:,level) = F(xprof,yprof);
        Vtempdata(:,level) = G(xprof,yprof);
        Wtempdata(:,level) = H(xprof,yprof);
    end

    % interpolate in z from xy interpolated z data
    if(Nkmax > 1)
        for prof = 1:length(xprof)
            uprof(prof) = interp1(z, Utempdata(prof,:), zprof(prof));
            vprof(prof) = interp1(z, Vtempdata(prof,:), zprof(prof));
            wprof(prof) = interp1(z, Wtempdata(prof,:), zprof(prof));
        end
    else
        % XY case (1 level)
        uprof = Utempdata; vprof = Vtempdata; wprof = Wtempdata;
    end 
    
end

% plot values and compare with data from literature
figure()
hold on
h11 = line(xprof(Iset1), -vprof(Iset1),'Color','k');
load ./LitData/6diamonds; Vdata = Ydata; Xdata = Xdata;
h12 = line(Xdata, Vdata,'Color','k','LineStyle','none','Marker','x');
ax1 = gca;
set(ax1, 'XColor','k','YColor','k');
axis([0 1 -1 1])
ax2 = axes('Position',get(ax1,'Position'), ...
    'XAxisLocation', 'top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','b', 'Ycolor','b');

h21 =  line(uprof(Iset2), yprof(Iset2(end:-1:1)),'Color','b','Parent',ax2);
load ./LitData/6circles; Udata = Ydata; Ydata = Xdata;
h22 = line(Udata,Ydata,'Color','b','LineStyle','none','Marker','x');
axis([-1 1 0 1])
xlabel(ax1,'u/U')
xlabel(ax2,'x/B')
ylabel(ax1,'v/U')
ylabel(ax2,'y/D')
saveas(gcf,'XY.eps','epsc2')
end

function z = getz(dz)
  
  N = length(dz);
  dzh = 0.5*(dz(1:end-1)+dz(2:end));
  
  z = zeros(N,1);
  z(1) = -dz(1)/2;
  z(2:N) = z(1)-cumsum(dzh);
end
