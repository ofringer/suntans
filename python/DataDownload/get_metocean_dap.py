"""
Tools for downloading specfific ocean/atmosphere/climate
datasets from an opendap server
"""

from mythredds import GetDAP, Dataset
from datetime import datetime,timedelta

###
# Dictionary containing opendap server specific
# information for each model
###
metoceandict = {\
    'HYCOM':{\
        'ncurl':'http://tds.hycom.org/thredds/dodsC/glb_analysis',\
        'type':'ocean',\
        'u':'u',\
        'v':'v',\
        'temp':'temperature',\
        'salt':'salinity',\
        'ssh':'ssh',\
    },\
    'GFS':{\
        'ncurl':[],\
        'type':'atmosphere',\
        'multifile':True,\
        'uwind':'ugrd10m',\
        'vwind':'vgrd10m',\
        'tair':'tmp2m',\
        'pair':'prmslmsl',\
        'rh':'rh2m',\
        'cloud':'tcdcclm',\
        'rain':'pratesfc',\
    },\
} # End of dictionary

#############################
# Dataset specific classes 
#############################

# These are for doing thigs like retrieving file names,
# setting coordinate projections etc
class GFSFiles:
    """
    Class that returns the filenames of global forecast system output
    for a given time range
    """

    resolution = 4 # 3 = one degree, 4 = 0.5 degree
    #dt = 3 # time interval between files in hours
    #baseurl =\
        #'http://nomads.ncdc.noaa.gov/thredds/dodsC/gfs-%03d/%s/%s/gfs_%d_%s_0000_%03d.grb2'
    dt = 24 
    baseurl =\
        'http://nomads.ncep.noaa.gov:9090/dods/gfs_0p50/gfs%s/gfs_0p50_00z'

    def __init__self(**kwargs):
        self.__dict__.update(**kwargs)

    def __call__(self,trange):

        self.basetime = datetime(1900,1,1)

        trange = [datetime.strptime(trange[0],'%Y%m%d.%H%M%S'),\
            datetime.strptime(trange[1],'%Y%m%d.%H%M%S')]

        # Generate a list of time variables
        dt = timedelta(hours=self.dt)
        t1 = trange[0]
        time = [t1]
        while t1 <= trange[1]:
            t1+=dt
            time.append(t1)
            print t1

        return time,[self.generate_url(tt) for tt in time]

    def generate_url(self,time):
        """
        Generate a url for a given time
        """
        def _gen_url(yymmdd,yyyymm,hours):
            #return self.baseurl%(self.resolution,\
            #    yyyymm,yymmdd,self.resolution,\
            #    yymmdd,hours)
            return self.baseurl%(yymmdd)


        yymmdd = datetime.strftime(time,'%Y%m%d')
        basetime = datetime.strptime(yymmdd,'%Y%m%d')

        # Generate the string
        yyyymm = datetime.strftime(time,'%Y%m')
        hours = (time-basetime).total_seconds()/3600

        url = _gen_url(yymmdd,yyyymm,hours)

        # Check if the url exists
        if not basetime == self.basetime:
            print 'Checking if url exists...\n\t%s'%url
            try:
                # Update to a new data
                #f = urllib2.urlopen('%s.html'%url)
                nc = Dataset(url)
                self.basetime = basetime
                print 'yes'
                nc.close()
                return url
            except:
                print 'File does not exist - we are in the forecast\
                    stage...(%s)'%(yymmdd)
                # Generate a string from the old basetime 
                yymmdd = datetime.strftime(self.basetime,'%Y%m%d')
                yyyymm = datetime.strftime(self.basetime,'%Y%m')
                hours = (time-self.basetime).total_seconds()/3600
                url = _gen_url(yymmdd,yyyymm,hours)
                return url


##############################
# Main functions for calling
##############################
def get_metocean_dap(xrange,yrange,zrange,trange,outfile,\
        gridfile=None,name='HYCOM'):
    """
    Extracts any met-ocean model that doesn't have server specific needs

    Inputs:
        name is dataset name in metoceandict (HYCOM by default)
        xrange/yrange/zrange: lists with i.e. [min(x),max(x)]

        trange: lists with start and end time in string format i.e.
            ['20101201.000000','20101202.000000']
    """
    oceandict=metoceandict[name]

    # Construct the dap class
    TDS = GetDAP(gridfile=gridfile,**oceandict)
    # Call the object
    TDS(xrange,yrange,trange,zrange=zrange,outfile=outfile)

    print 'Done.'


def get_metocean_local(ncfile,varname,name='HYCOM',TDS=None,\
        xrange=None,yrange=None,zrange=None,trange=None):
    """
    Retrieves variable data from a local file

    Use common name for variable i.e. u,v,temp,tair,uwind, etc

    Returns: data, ncobject (containing coordinates, etc)

    This is sort of a long way around but is general
    """
    oceandict=metoceandict[name]

    # 
    ncvar = oceandict[varname]
    oceandict['ncurl']=ncfile
    if oceandict.has_key('multifile'):
        oceandict['multifile']=False

    # Construct the dap class
    if TDS==None:
        TDS = GetDAP(vars=[varname],**oceandict)

    ncobj = getattr(TDS,'_%s'%varname)

    # Call the object
    data = ncobj.get_data(ncvar,xrange,yrange,zrange,trange)

    return data, ncobj


def get_gfs_tds(xrange,yrange,zrange,trange,outfile):
    """
    Extract Global forecasting system data
    """
    gfsdict = metoceandict['GFS']

    # Get the file names for the given time range from the class
    gfs = GFSFiles()
    time,files = gfs(trange)

    # Update the dictionary
    gfsdict['ncurl']=files

    # Create the thredds object
    TDS = GetDAP(**gfsdict)
    
    # Call the object
    TDS(xrange,yrange,trange,zrange=zrange,outfile=outfile)

