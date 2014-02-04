%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: profplot.m
% Description: Plot the data produced by the profile output in suntans.
%
% Oliver Fringer
% Stanford University
% 07/15/06
%
% This mfile will plot the data output along the transect specified 
% by the parameters in suntans.dat.  An example of the parameters
% in suntans.dat is given below:
%
% ProfileVariables	default          # Output u, s, s0, and h
% DataLocations	        dataxy.dat       # dataxy.dat contains column x-y data
% ProfileDataFile	profdata.dat     # Information about profiles is in profdata.dat
% ntoutProfs		10               # Output profile data every 10 time steps
% NkmaxProfs		100              # Only output the top 10 z-levels
% numInterpPoints	1                # Output data at the three nearest neighbors to each input point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirname = '../main/examples/iwaves/data';
EMPTY = 999999;

fname = [dirname,'/profdata.dat'];

fid = fopen(fname,'rb');
numTotalDataPoints = fread(fid,1,'int32');
numInterpPoints = fread(fid,1,'int32');
Nkmax = fread(fid,1,'int32');
nsteps = fread(fid,1,'int32');
ntoutProfs = fread(fid,1,'int32');  
dt = fread(fid,1,'float64');
dz = fread(fid,Nkmax,'float64');
dataIndices = fread(fid,numTotalDataPoints,'int32');
dataXY = fread(fid,2*numTotalDataPoints,'float64');
xv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
yv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
fclose(fid);

dataX = dataXY(1:2:end);
dataY = dataXY(2:2:end);

dzh = 0.5*(dz(1:end-1)+dz(2:end));
z = zeros(Nkmax,1);
z(1) = -dz(1)/2;
z(2:Nkmax) = z(1) - cumsum(dzh);
dZ = dz*ones(1,length(dataX));

x = sqrt((dataX-dataX(1)).^2+(dataY-dataY(1)).^2)';
X = ones(Nkmax,1)*x;
Z = z*ones(1,numTotalDataPoints);

nout=nsteps/ntoutProfs;

% These locations will be stored at each time step and then
% plotted as a time series.
iplot=fix(numTotalDataPoints/2);
kplot=fix(Nkmax/2);

% How often to plot
nplot=20;

fname = [dirname,'/fs.dat.prof'];
hfid = fopen(fname,'rb');
fname = [dirname,'/u.dat.prof'];
ufid = fopen(fname,'rb');
fname = [dirname,'/s.dat.prof'];
sfid = fopen(fname,'rb');
fname = [dirname,'/s0.dat.prof'];
s0fid = fopen(fname,'rb');

s0data = reshape(fread(s0fid,Nkmax*numInterpPoints*numTotalDataPoints,'float64'),...
                 Nkmax,numInterpPoints,numTotalDataPoints);
  
for n=1:nout
  rtime = n*ntoutProfs*dt;
    
  fprintf('On %d of %d\n',n,nsteps/ntoutProfs);

  hdata = reshape(fread(hfid,numInterpPoints*numTotalDataPoints,'float64'),...
                  numInterpPoints,numTotalDataPoints);
  
  udata = reshape(fread(ufid,3*Nkmax*numInterpPoints*numTotalDataPoints,'float64'),...
                  Nkmax,numInterpPoints,numTotalDataPoints,3);
  udata(find(udata==EMPTY))=nan;
    
  sdata = reshape(fread(sfid,Nkmax*numInterpPoints*numTotalDataPoints,'float64'),...
                  Nkmax,numInterpPoints,numTotalDataPoints);
  sdata(find(sdata==EMPTY))=nan;
  
  h = squeeze(hdata(1,:));
  s = squeeze(sdata(:,1,:));
  s0 = squeeze(s0data(:,1,:));      
  dZ(find(isnan(s)))=0;
  depth = sum(dZ);      
  
  u = squeeze(udata(:,1,:,1));
  v = squeeze(udata(:,1,:,2));
  w = squeeze(udata(:,1,:,3));

  if(~rem(n,nplot))
    figure(1)
    subplot(3,1,1)
    pcolor(X,Z,u);
    shading flat;
    colorbar
    set(gca,'xticklabel','');
    title('u');
    
    subplot(3,1,2)
    pcolor(X,Z,w)
    shading flat;
    colorbar
    set(gca,'xticklabel','');
    title('w');
    
    subplot(3,1,3)  
    pcolor(X,Z,s);
    shading flat;
    colorbar  
    title('\rho''');
  end
  
  hplot(n) = h(iplot);
  uplot(n) = u(kplot,iplot);
  wplot(n) = w(kplot,iplot);
  splot(n) = s(kplot,iplot);
  
  pause(0);
end

dT = ntoutProfs*dt;
t = dT*[1:nout];

figure(2)

subplot(4,1,1)
plot(t,hplot);
set(gca,'xticklabel','');
ylabel('h(t)');

subplot(4,1,2)
plot(t,uplot);
set(gca,'xticklabel','');
ylabel('u(t)');

subplot(4,1,3)
plot(t,wplot);
set(gca,'xticklabel','');
ylabel('w(t)');

subplot(4,1,4)
plot(t,splot);
xlabel('t');
ylabel('\rho(t)');