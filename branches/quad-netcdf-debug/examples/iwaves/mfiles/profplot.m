%
% Plot profiles of u(z,t) at the different stations specified in
% ../rundata/dataxy.dat
%
dirname = '../data';
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

z = getz(dz);

fname = [dirname,'/u.dat.prof'];
fid = fopen(fname,'rb');

data = fread(fid,'float64');
data(find(data==EMPTY))=nan;

nout = length(data)/(3*Nkmax*numInterpPoints*numTotalDataPoints);
udata = reshape(data,Nkmax,numInterpPoints,numTotalDataPoints,3,nout);
udata = squeeze(udata);

Time = ones(Nkmax,1)*[1:nout]*dt*ntoutProfs;
Z = z*ones(1,nout);
Tday = 86400;

u = squeeze(udata(:,:,1,:));
v = squeeze(udata(:,:,2,:));
w = squeeze(udata(:,:,3,:));

figure(1);
for loc=1:3
  
  uplot = squeeze(u(:,loc,:));
  kmax = max(find(~isnan(uplot(:,1))));
  d0 = sum(dz(1:kmax));
  u0 = sum((dz(1:kmax)*ones(1,nout)).*uplot(1:kmax,:))/d0;
  
  subplot(3,1,loc)
  pcolor(Time/Tday,Z,uplot-ones(Nkmax,1)*u0);
  title(sprintf('Baroclinic u(z,t) at Location %d (x = %.0f km)',loc,dataX(loc)/1000));
  shading flat;
  axis([0 max(Time(1,:)/Tday) -d0 0]);
  colorbar;
  
  if(loc==3) 
    xlabel('Time (days)'); 
  else
    set(gca,'xticklabel','');
  end
  ylabel('Depth (m)');
end
