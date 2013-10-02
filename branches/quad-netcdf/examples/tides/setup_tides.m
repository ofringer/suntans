addpath ../../mfiles
addpath ../../mfiles/suntides
addpath packages/m_map

% Load the m_map files to convert from x,y to lon,lat using UTM
m_proj('utm','lon',[-128 -119],'lat',[34 40],'zon',[10],'hem',[0],'ell',['wgs84']);

% Location of xy points of suntans tidal boundaries
tideprefix = './data/tidexy.dat';

% Output location of tidal component files for input into
% suntans.
outputprefix = './data/tidecomponents.dat';

% Location of OTIS data; This example is for the Southern Pacific Ocean
tidespath = './tide_data/pac';

% Year of simulation (so that phases are relative to start of year)
year = 2006;

% Number of processors
numprocs = 1;

% Data is in lon/lat (not xy projection)
% XY = [];
% Otherwise specify conversion function to convert from xy to lon/lat
XY = 'm_xy2ll';

% Verbose output
VERBOSE = 1;

% Initialize tidal files
initialize_tides(tidespath);