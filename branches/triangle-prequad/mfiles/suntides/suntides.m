% 
% SUNTIDES Write tidal component input files for SUNTANS using OTIS.
%   SUNTIDES(TIDEPREFIX,OUTPUTPREFIX,TIDESPATH,YEAR,NUMPROCS,'XY2LL',VERBOSE)
%   takes as its input a list of lat-lon coordinates in the
%   file defined by TIDEPREFIX.n, where n is the processor number (0 to
%   NUMPROCS-1).  It then writes the tidal component files in a format
%   suitable for input into SUNTANS into the file OUTPUTPREFIX.  Note that tidal component files for
%   processors with no boundary edges are still written, but limited
%   data is written to these files.
%
%   The path to the OTIS data must be specified by the TIDESPATH
%   variable so that the initialize_tides script (in this directory)
%   can load the following files:
% 
%     load([TIDESPATH,'/so/allconsts.mat']);
%     load([TIDESPATH,'/so/latlons.mat']);
%     load([TIDESPATH,'/lun_node.dat']);
%     load([TIDESPATH,'/v0u.dat']);
%
%   The YEAR must also be specified.  Also note that if the lat/lon
%   coordinates are outside the range of the data in the specified path
%   then this script will return with an error.
%
%   It is assumed that the coordinates in the TIDESPREFIX file are in
%   lon,lat; however, if they are in x-y coordinates, SUNTIDES must be
%   called with a function that returns the lon-lat coordinates of the
%   given x-y coordinates in the variable XY2LL, where XY2LL is a
%   string containing the function that returns lon/lat given x,y, i.e.
%
%   [LON,LAT]=XY2LL(X,Y);
%
%   Otherwise leave this variable empty with [].
%
%   The last argument is interpreted as the VERBOSE option, which is either 1, for verbose
%   output, or 0, for silence.
%
%   SUNTIDES without any arguments runs the example in suntides_example.

%
%   Oliver Fringer
%   Stanford University
%   15 Dec 06
%
function suntides(varargin)

  if(nargin==0)
    suntides_example;
    return;
  end
  
  if(nargin~=7)
    fprintf('Error: Wrong number of arguments; Argument list must be:\n');
    fprintf('\tsuntides(tideprefix,outputprefix,tidespath,year,numprocs,''xy2ll'',verbose)\n');
    return;
  end
  tideprefix = varargin{1};
  outputprefix = varargin{2};
  tidespath = varargin{3};
  year = varargin{4};
  numprocs = varargin{5};
  if(isempty(varargin{6}))
    XY=0;
  else
    XY=1;
    XY2LL = varargin{6};
  end
  VERBOSE = varargin{7};
  
  INT = 'int';
  REAL = 'float64';
  
  % Initialize the tidal data for use in the interpolation of OTIS tidal
  % constituents.
  initialize_tides(tidespath);
  
  % Get omegas
  [omegas,uamp,uphase,vamp,vphase,hamp,hphase]=get_tides(0,0,0);
  numtides = length(omegas);
  
  for proc=0:numprocs-1
    
    fname = sprintf('%s.%d',tideprefix,proc);
    ofile = sprintf('%s.%d',outputprefix,proc);
    
    fid = fopen(ofile,'wb');
    numboundaryedges = 0;
    if(exist(fname,'file'))
      data = load(fname,'-ascii');
      if(~isempty(data))
        if(XY)
          eval(['[xlon,ylat]=',XY2LL,'(data(:,1),data(:,2));']);
          data = [xlon,ylat];
        end
        numboundaryedges = length(data(:,1));
      end
    end
    
    if(numboundaryedges>0)
      if(VERBOSE)
        fprintf('File %s, found %d boundary points.\n',fname,numboundaryedges);
      end      
      
      for n=1:numboundaryedges
        lon = data(n,1);
        lat = data(n,2);
        [omegas,uamp,uphase,vamp,vphase,hamp,hphase]=get_tides(year,lat,lon);
        if(isempty(uamp))
          return;
        end
        if(n==1)
          fwrite(fid,numtides,INT);
          fwrite(fid,numboundaryedges,INT);
          fwrite(fid,omegas,REAL);
        end
        fwrite(fid,uamp,REAL);
        fwrite(fid,uphase,REAL);
        fwrite(fid,vamp,REAL);
        fwrite(fid,vphase,REAL);
        fwrite(fid,hamp,REAL);
        fwrite(fid,hphase,REAL);
      end
    else
      fwrite(fid,numtides,INT);
      fwrite(fid,numboundaryedges,INT);
      fwrite(fid,omegas,REAL);
      
      if(VERBOSE)
        fprintf('No data for processor %d.\n',proc);
      end
    end
  end
  