%
% SUNTIDES_EXAMPLE Example use of SUNTIDES (also called by SUNTIDES
% without arguments).

%
%   Oliver Fringer
%   Stanford University
%   9 Sep 07
%
function suntides_example

  % Location of xy points of suntans tidal boundaries
  tideprefix = './data/tidal_files/tidexy.dat';
  
  % Output location of tidal component files for input into
  % suntans.
  outputprefix = './data/tidal_files/tidecomponents.dat';
  
  % Location of OTIS data; This example is for the Southern Pacific Ocean
  % Only needs allconsts.mat and latlons.mat in this directory.
  tidespath = 'data/so';
  
  % Year of simulation (so that phases are relative to start of year)
  year = 2004;
  
  % Number of processors
  numprocs = 32;
  
  % Data is in lon/lat (not xy projection)
  XY = [];
  % Otherwise specify conversion function to convert from xy to lon/lat
  % XY = @xy2lonlat;
  
  % Verbose output
  VERBOSE = 1;
  
  % Create files for input into suntans
  suntides(tideprefix,outputprefix,tidespath,year,numprocs,XY,VERBOSE);
  
  % Check these files
  check_suntides(tideprefix,outputprefix,tidespath,year,numprocs,XY,VERBOSE);
