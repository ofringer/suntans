%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lineplot.m 
% usage: 1. line plot for specific location for different fields
%           for the results in tri or quad grids
%        2. use linear interpolation to get the data 
%        3. can use in comparison between tri and quad results
% Parameters:
% 1) type = one of:
% i) 'q' nonhydrostatic pressure
% ii) 's' salinity
% iii) 'u','v','w' velocity fields
% iv) 'h free surface
% 2) datadir=location of data, i.e. '/home/data'
% 3) nstep= time step to plot
% 4) kevel=vertical level to plot
% 5) procs: total number of procs 
% 6) x,y: the line you want interpolate to 
% 7) METHOD: method for interpolation      
%    'nearest'   - Nearest neighbor interpolation
%    'linear'    - Linear interpolation (default)
%    'natural'   - Natural neighbor interpolation
%    'cubic'     - Cubic interpolation (2D only)
%    'v4'        - MATLAB 4 griddata method (2D only)
% 8) LINECOLOR: the color of the line
% 9) DIMENSION: 2 FOR 2D VECTOR(for simple line or y(x) is constant)
%               3 FOR 3D VECTOR
% 10) X: only for DIMENSION=2, 0 for x, 1 for y
% input: suntans results 
% output: figure and data along the line you assign
% Made by Yun Zhang
% 01/25/2013 @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi1 = lineplot(PLOT,griddir,datadir,n,klevel,procs,x,y,METHOD,LINECOLOR,DIMENSION,X)
% read data
flag = 1;
EMPTY=999999;         % Empty cells are defined by this
dbl=8;                % Size of double precision in bytes
precision='float64';  % Precision for reading in data
%  procs=0;              % how many processors 
%  klevel=1;             % which layer you want to plot
%  PLOT='u';
%  n=1;                  % time step
%  x=3:6:120;
%  y=30*ones(size(x));
% cellcentered data contains the voronoi points and the depths
% at those points.
N=length(procs);  
% held=ishold;
% if(~held)
%     hold on;
% end

% % directory of data
%  datadir='tridata/trinodragamp0.05';
%  griddir='tridata/trinodragamp0.05';


% get filename
proc=0;

for proc=procs
switch(PLOT)
    case 'q'
        filename = [datadir,'/q.dat.',num2str(proc)];
    case 's'
        filename = [datadir,'/s.dat.',num2str(proc)];
    case {'u','v','w','Q'}
        filename = [datadir,'/u.dat.',num2str(proc)];
    case 'T'
        filename = [datadir,'/T.dat.',num2str(proc)];
    case {'h','D','d'}
        filename = [datadir,'/fs.dat.',num2str(proc)];
    case{'n'}
        filename = [datadir,'/nut.dat.',num2str(proc)];
    case{'tk'}  
        filename = [datadir,'/tk.dat.',num2str(proc)];
    otherwise
end
hfilename = [datadir,'/fs.dat.',num2str(proc)];
fid = fopen(filename,'rb');
hfid = fopen(hfilename,'rb');
if(fid<=0)
    fprintf('Error opening %s\n',filename);
    flag = 0;
    return;
end
if(hfid<=0 & proc==procs(1))
    fprintf('Error opening fs data\n continue with no fs info');
end

% get grid data
pd = load([griddir,'/points.dat']);
cp = load([griddir,'/cells.dat.',num2str(proc)]);
cdp = load([griddir,'/celldata.dat.',num2str(proc)]);
edp= load([griddir,'/edges.dat.',num2str(proc)]);
xp = (pd(:,1));%+558442)/1000;
yp = (pd(:,2));%+5318790)/1000;
nfaces=cp(:,1);
tri2=edp(:,1:2);
mark=edp(:,3);
d = cdp(:,5);
xv = cdp(:,2);
yv = cdp(:,3);
dz = load([griddir,'/vertspace.dat']);
% Total number of cells in the horizontal Nc and vertical Nk
Nc = length(xv);
Nk = length(dz);
tri=zeros(Nc,max(nfaces));
for i=1:Nc
   tri(i,1:nfaces(i))=cp(i,4:(4+nfaces(i)-1));
end
z = -(sum(dz(1:klevel))-dz(klevel)/2)*ones(size(xv));

% get array size
switch(PLOT)
    case 'q'
        arraysize = Nc*Nk;
    case {'s','T','n','tk'}
        arraysize = Nc*Nk;
    case {'u','v','w','Q'}
        arraysize = Nc*Nk;
    case {'h','D','d'}
        arraysize = Nc;
    otherwise
end

 phi = zeros(Nc,1);
    tmp =phi;
% read data we want    
switch(PLOT)
    case 'u'
        fseek(fid,8*(3*arraysize*(n-1)+3*(klevel-1)*Nc),0);
        phi = fread(fid,Nc,'float64');
   case 'v'
        fseek(fid,8*(3*arraysize*(n-1)+3*(klevel-1)*Nc+Nc),0);
        phi = fread(fid,Nc,'float64');
   case 'w'
        fseek(fid,8*(3*arraysize*(n-1)+3*(klevel-1)*Nc+2*Nc),0);
        phi = fread(fid,Nc,'float64');
   case 'Q'
        fseek(fid,8*(3*arraysize*(n-1)+3*(klevel-1)*Nc),0);
        tmp = fread(fid,Nc,'float64');
        phi=fread(fid,Nc,'float64');
        phi=sqrt(phi.^2+tmp.^2); % same sign as u
        %phi=tmp*cos(0.865)+phi*sin(0.865);
   otherwise % s T nut
        fseek(fid,8*(arraysize*(n-1)+(klevel-1)*Nc),0);
        phi = fread(fid,Nc,'float64');
end

%%%%% set the data to nan for cells below depth %%%%
phi(find(phi==EMPTY))=nan   ; 2; -1;nan;
%%%%% set the data to nan for cells above free surface (value=zero) %%%%
phi(find(phi==0))=nan;

%%%%% when there is free surface data, check drying condition
if(hfid>0) 
    hphi=zeros(Nc,1);
    fseek(hfid,8*Nc*(n-1),0);
    hphi = fread(hfid,Nc,'float64');
       
    switch(PLOT)
        case 'D'
            phi=hphi;
            phi=phi+d;
            phi(find(phi<1e-3))=0;
        case 'h'
            phi=hphi;
            phi(find(hphi+d<1e-3))=0;
            phi=phi;
        case {'s','T','n','tk'}
            phi(find(hphi+d<1e-3))=0;
        case 'd'
            phi=-d+4; % elevation
            phi(find(hphi+d<0.1))=10;
        otherwise
            phi(find(hphi+d<5e-2))=10;
    end
end

% interpolation 
result=griddata(xv,yv,phi,x,y,METHOD);
if DIMENSION==3
   plot3(x,y,result,LINECOLOR)
   aaa=num2str(n);
   aaa2=strcat('Interpolation results along the line you assign in 3D vector at step ',aaa);
   title(aaa2);
   xlabel('x');
   ylabel('y');
   zlabel(PLOT);
end
if DIMENSION==2
    if X==0;
        plot(x,result,LINECOLOR);
        aaa=num2str(n);
        aaa2=strcat('Interpolation results along x axis at step ',aaa);
        title(aaa2)
        xlabel('x');
        ylabel(PLOT);
    elseif X==1
        plot(y,result,LINECOLOR);
        aaa=num2str(n);
        aaa2=strcat('Interpolation results along y axis at step ',aaa);
        title(aaa2)
        xlabel('y');
        ylabel(PLOT);
    end
end
end

phi1=result;




