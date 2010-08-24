%  
%  SHOWTIDE Plot tidal time series of u, v, and h at a given location.
%    SHOWTIDE(tidespath,lon,lat,daystart,numdays) Plots the tidal time
%    series at the location (lon,lat) over the period of length numdays
%    starting at daystart, which must be in a format that datestr
%    recognizes (i.e. '01/01/2000');  The location of the OTIS data must
%    be specified by tidespath.  Time series are plotted using eight tidal
%    constituents, namely m2,s2,n2,k2,o1,k1,p1,q1.
%
%    SHOWTIDE(tidespath,lon,lat,daystart,numdays,constituents) prints the
%    desired tidal amplitudes to the screen and only plots those that are
%    specified.  For example, if constituents is {'m2','s2'}, the tidal
%    amplitudes for m2 and s2 are printed to the screen, and the tidal time
%    series only contains the effects of those two constituents.
%    Specifying 'all' displays all constituents.
%
%    SHOWTIDE without any arguments displays an example.

%
%  Oliver Fringer
%  Stanford University
%  08/18/07
%
function showtide(varargin)

    if(nargin==0)
        exampledirectory = '../sopred';
        if(~exist(exampledirectory,'dir'))
            fprintf('Error: OTIS tidal data for southern ocean must be in %s to run the example.\n',exampledirectory);
            return;
        end
        lat_test = -19.57;
        lon_test = -243.8;
        daystart_test = '3/1/2004';
        numdays_test = 60;
        showtide(exampledirectory,lon_test,lat_test,daystart_test,numdays_test,'all');
        return;
    end
    
    allcomponents = {'m2','s2','n2','k2','o1','k1','p1','q1'};
    tconst = zeros(size(allcomponents));

    if(nargin<5 | nargin>6)        
        fprintf('Error: Wrong number of arguments; Argument list must be one of:\n');
        fprintf('\tshowtide\n');
        fprintf('\tshowtide(tidespath,lon,lat,daystart,numdays)\n');
        fprintf('\tshowtide(tidespath,lon,lat,daystart,numdays,constituents)\n');
        return;
    end
    tidespath=varargin{1};
    lon=varargin{2};
    lat=varargin{3};
    daystart=varargin{4};
    numdays=varargin{5};
    if(nargin>5)
        constituents=varargin{6};
    else
        constituents='';
    end
    
    year = str2num(datestr(daystart,10));
    initialize_tides(tidespath);
    [omegas,uamp,uphase,vamp,vphase,hamp,hphase]=get_tides(year,lat,lon);
    if(isempty(uamp))
        return;
    end
    
    if(strcmp(constituents,'all'))
        constituents = allcomponents;
    end
    disp_consts=1;
    if(strcmp(constituents,''))
        tconst = ones(size(omegas));
        constituents = allcomponents;
    else
        fprintf('Tidal components (phase relative to year %d, u,v in cm/s and h in cm):\n',year);
        fprintf('\t\t\t\tuamp\tuphase\tvamp\tvphase\thamp\thphase\n');
        for nc=1:length(allcomponents)
            if(~isempty(find(strcmp(allcomponents{nc},constituents))))
                tconst(nc)=1;
                fprintf('%s(%.2f hr):\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',...
                    allcomponents{nc},2*pi/omegas(nc)/3600,...
                    uamp(nc),uphase(nc),vamp(nc),vphase(nc),hamp(nc),hphase(nc));
            end
        end
    end
    omegas = omegas.*tconst;
    
    Tday = 3600*24;
    dthour = 3600;
    nstart = datenum(daystart) - datenum(['01/01/',num2str(year)]);
    t = [nstart:1/24:nstart+numdays]'*Tday;

    u = zeros(size(t));
    v = zeros(size(t));
    h = zeros(size(t));

    q = ones(size(t));
    U = (q*uamp).*cos(t*omegas+q*uphase);
    V = (q*vamp).*cos(t*omegas+q*vphase);
    H = (q*hamp).*cos(t*omegas+q*hphase);
 
    u = sum(U')';
    v = sum(V')';
    h = sum(H')';

    % Add constituents to title
    cstr = 'Constituents: ';
    if(isa(constituents,'cell'))
        for nc=1:length(constituents)
            cstr = [cstr,' ',constituents{nc}];
        end
    else
        cstr = [cstr,' ',constituents];
    end 
    
    subplot(3,1,1)
    plot(t/Tday,u,'k-');
    set(gca,'xticklabel','');
    ylabel('u (cm/s)');
    title(cstr);
    
    subplot(3,1,2)
    plot(t/Tday,v,'k-');
    set(gca,'xticklabel','');
    ylabel('v (cm/s)');

    subplot(3,1,3)
    plot(t/Tday,h,'k-');
    xlabel(sprintf('Days in %d',year));
    ylabel('h (cm)');