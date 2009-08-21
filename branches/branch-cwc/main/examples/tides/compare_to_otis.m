%
% Compare suntans result to OTIS result.
% 

%
%   Oliver Fringer
%   Stanford University
%   9 Sep 07
%

% datadir where suntans data is stored
datadir = './data';
suntansdatafile = [datadir,'/suntans.dat'];

% Load in required data and default values
setup_tides

% Load points desired from dataxy.dat
dataxy = load([datadir,'/dataxy.dat']);
x = dataxy(1);
y = dataxy(2);
[lon,lat]=m_xy2ll(x,y);
ddata = load([datadir,'/mbay_bathy.dat-voro'],'ascii');

% Obtain values for simulation from suntans.dat
dt = getvalue(suntansdatafile,'dt');
ntoutProfs = getvalue(suntansdatafile,'ntoutProfs');
nsteps = getvalue(suntansdatafile,'nsteps');
toffSet = getvalue(suntansdatafile,'toffSet');

% Get tidal signals from OTIS
Tday = 86400;
daystart = toffSet;
numdays = dt*nsteps/Tday;
totis = [daystart:dt*ntoutProfs/Tday:daystart+numdays]';

% Get data from OTIS; note that OTIS data is in cm/s, so it
% needs to be converted
[omegas,uamp,uphase,vamp,vphase,hamp,hphase]=get_tides(year,lat,lon);
omegas = omegas*Tday;
q = ones(size(totis));
hotis = sum(((q*hamp).*cos(totis*omegas+q*hphase))')'/100;
uotis = sum(((q*uamp).*cos(totis*omegas+q*uphase))')'/100;
votis = sum(((q*vamp).*cos(totis*omegas+q*vphase))')'/100;

% Get data from SUNTANS run
tsuntans = [daystart:dt*ntoutProfs/Tday:daystart+numdays];
hsuntans = fread(fopen([datadir,'/fs.dat.prof'],'rb'),'float64');
Usuntans = fread(fopen([datadir,'/u.dat.prof'],'rb'),'float64');
usuntans = Usuntans(1:3:end);
vsuntans = Usuntans(2:3:end);
wsuntans = Usuntans(3:3:end);

figure(1);
subplot(3,1,1)
plot(totis,uotis,'k-',tsuntans(1:length(hsuntans)),usuntans,'r');
axis tight;
ylabel('u (m s^{-1})');

subplot(3,1,2)
plot(totis,votis,'k-',tsuntans(1:length(hsuntans)),vsuntans,'r');
axis tight;
ylabel('v (m s^{-1})');

subplot(3,1,3)
plot(totis,hotis,'k-',tsuntans(1:length(hsuntans)),hsuntans,'r');
axis tight;
ylabel('h (m)');
xlabel(sprintf('Days in %d (SUNTANS: red, OTIS: black)\n',year));




