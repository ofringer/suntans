%
% SUNTIDES_EXAMPLE Example use of SUNTIDES (also called by SUNTIDES
% without arguments).

%
%   Oliver Fringer
%   Stanford University
%   9 Sep 07
%

% Load in required data and default values
setup_tides
  
% Create files for input into suntans
suntides(tideprefix,outputprefix,tidespath,year,numprocs,XY,VERBOSE);
