clear
%
lats=0.0;
late=65.5;
lat=lats:(1/2):late;
lat=lat';
lons=100;
lone=310.5;
lon=lons:(1/2):lone; 
%
latu=lat;
lonu=lon-0.5*(1/2);
latv=lat-0.5*(1/2);
lonv=lon;
lon=lon-360;
lonu=lonu-360;
lonv=lonv-360;
%
lat=lat*ones(1,422);
lon=ones(132,1)*lon;
latu=latu*ones(1,422);
lonu=ones(132,1)*lonu;
latv=latv*ones(1,422);
lonv=ones(132,1)*lonv;

save latlons.mat
