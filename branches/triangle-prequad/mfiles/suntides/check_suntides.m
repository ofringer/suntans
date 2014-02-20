% 
% CHECK_SUNTIDES Check output from SUNTIDES.
%   CHECK_SUNTIDES(TIDEPREFIX,OUTPUTPREFIX,TIDESPATH,YEAR,NUMPROCS,'XY2LL',VERBOSE)
%   takes as its input a list of lat-lon coordinates in the
%   file defined by TIDEPREFIX.n, where n is the processor number (0 to
%   NUMPROCS-1), and plots the tidal ellipses at each one of these points in
%   the file OUTPUTPREFIX.n.
%
%   Although the path to the OTIS data must be specified by the TIDESPATH
%   variable, it is not used but the variable list is chose to match that
%   of SUNTIDES.
%
%   It is assumed that the coordinates in the TIDESPREFIX file are in
%   lon,lat; however, if they are in x-y coordinates, this function must be
%   called using a function that returns the lon-lat coordinates of the
%   given x-y coordinates.  This is done with
%
%   SUNTIDES(TIDEPREFIX,OUTPUTPREFIX,TIDESPATH,YEAR,NUMPROCS,'XY2LL')
%
%   Where XY2LL is the name of a function (a string) that returns lon/lat given x,y, i.e.
%   [LON,LAT]=XY2LL(X,Y);
%
%   If the last argument is an integer (either argument 6 or 7), it
%   is interpreted as the VERBOSE option, which is either 1, for verbose
%   output, or 0, for silence.

%
%   Oliver Fringer
%   Stanford University
%   9 Sep 07
%
function check_suntides(varargin)
  
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
  
  uamps = [];
  uphases = [];
  vamps = [];
  vphases = [];
  hamps = [];
  hphases = [];
  boundarypoints = [];
  
  tidecomponents = {'M2','S2','N2','K2','O1','K1','P1','Q1'};
  
  for proc=0:numprocs-1
    fname = sprintf('%s.%d',outputprefix,proc);
    
    if(exist(fname,'file'))
      fid = fopen(fname,'rb');
      numtides = fread(fid,1,INT);
      numboundaryedges = fread(fid,1,INT);
      omegas = fread(fid,numtides,REAL);
      
      if(numboundaryedges>0)
        data = load(sprintf('%s.%d',tideprefix,proc));
        boundarypoints = [boundarypoints ; data];
        
        if(VERBOSE)
          fprintf('File %s, found %d boundary points.\n',fname,numboundaryedges);
        end      
        
        for n=1:numboundaryedges
          if(n==1)          
            uamp = zeros(numboundaryedges,numtides);
            uphase = zeros(numboundaryedges,numtides);
            vamp = zeros(numboundaryedges,numtides);
            vphase = zeros(numboundaryedges,numtides);
            hamp = zeros(numboundaryedges,numtides);
            hphase = zeros(numboundaryedges,numtides);
          end
          uamp(n,:) = fread(fid,[1,numtides],REAL);
          uphase(n,:) = fread(fid,[1,numtides],REAL);
          vamp(n,:) = fread(fid,[1,numtides],REAL);
          vphase(n,:) = fread(fid,[1,numtides],REAL);
          hamp(n,:) = fread(fid,[1,numtides],REAL);
          hphase(n,:) = fread(fid,[1,numtides],REAL);
        end
        
        uamps = [uamps ; uamp];
        uphases = [uphases ; uphase];
        vamps = [vamps ; vamp];
        vphases = [vphases ; vphase];
        hamps = [hamps ; hamp];
        hphases = [hphases ; hphase];
      else
        if(VERBOSE)
          fprintf('File %s has no boundary points.\n',fname);
        end
      end
      
      fclose(fid);
    else
      if(VERBOSE)
        fprintf('No data for processor %d.\n',proc);
      end
    end
  end
  
  numboundaryedges = length(boundarypoints);
  lon = boundarypoints(:,1);
  lat = boundarypoints(:,2);
  if(XY)
    eval(['[lon,lat]=',XY2LL,'(lon,lat);']);
  end
  
  lonmin = min(lon);
  lonmax = max(lon);
  latmin = min(lat);
  latmax = max(lat);
  dlon = lonmax-lonmin;
  dlat = latmax-latmin;
  
  % For plotting - extra point for scale vector
  lonp = [lon;mean(lon)];
  latp = [lat;mean(lat)];
  
  AX = [lonmin-0.1*dlon lonmax+0.1*dlon latmin-0.1*dlat latmax+0.1*dlat];
  skip = 2;
  
  figure(1);
  clf;

  NUM_PLOT=40;
  SCALE_FACTOR=0.1;
  for o=1:numtides
    figure(o);
    clf;
    hold on;
    axis image;
    axis(AX);
    title(sprintf('%s (%.2f hr)',tidecomponents{o},2*pi/omegas(o)/3600));
    set(gca,'box','on');
    xlabel('^o E');
    ylabel('^o N');
    
    u = uamps(:,o).*cos(uphases(:,o));
    v = vamps(:,o).*cos(vphases(:,o));
    uscale = max(sqrt(u.^2+v.^2));
    
    up = [u;uscale];
    vp = [v;0];
    ellipses(lon([1:skip:end-1]),lat([1:skip:end-1]),...
             uamps([1:skip:end-1],o),uphases([1:skip:end-1],o),...
             vamps([1:skip:end-1],o),vphases([1:skip:end-1],o),NUM_PLOT,SCALE_FACTOR);
    if(o<numtides)
      fprintf('Showing %s; hit enter to continue...\n',tidecomponents{o});
      pause;
    else
      fprintf('Showing %s.\n',tidecomponents{o});
    end
  end
  