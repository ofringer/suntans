"""
Various tools for interacting with opendap thredds datasets

Matt Rayson
Stanford University
October 2014
"""

import os
import numpy as np
from netCDF4 import Dataset, num2date,date2num,date2index
from datetime import datetime
import othertime

import pdb

####
# Utility functions
####
def hasdim(nc,dim):
    return nc.dimensions.has_key(dim)
def hasvar(nc,var):
    return nc.variables.has_key(var)
def gettime(nc,timename):
    t = nc.variables[timename]
    return num2date(t[:],t.units)



class GridDAP(object):
    """
    General class for extracting gridded opendap data
    """

    timecoord = 'time'
    xcoord = 'lon'
    ycoord = 'lat'
    zcoord = None

    gridfile = None
    multifile = False

    def __init__(self,url,**kwargs):

        self.__dict__.update(kwargs)

        # Open the opendap dataset
        if isinstance(url,str):
            self._nc = Dataset(url)
        elif isinstance(url,Dataset):
            # we're passing Dataset object
            self._nc = url
        elif isinstance(url,list) and self.multifile:
            # We are passing it a list for multifile variables
            self._nc = Dataset(url[0])
            self._ncfiles = url # List of file names

        else:   
            raise Exception, 'Unknown input url type'

        # Open the grid netcdf (if it exists)
        if not self.gridfile==None:
            self._gridnc = Dataset(self.gridfile)
        else:
            self._gridnc = self._nc

        self.get_coords()

    def get_coords(self):
        """ Download the coordinates"""

        self.X = self._gridnc.variables[self.xcoord][:]
        self.Y = self._gridnc.variables[self.ycoord][:]

        if not self.zcoord == None:
            self.Z = self._gridnc.variables[self.zcoord][:]

    def get_time_indices(self,trange):
        """
        Find the time indices
        """

        if not self.multifile:
            self.time = gettime(self._nc,self.timecoord)

        else:
            self._MF = MFncdap(self._ncfiles,timevar=self.timecoord)
            self.time = self._MF.time

        # Time
        if trange==None:
            self.t1=0
            self.t2 = self.time.shape[0]
        else:
            self.t1 = othertime.findNearest(trange[0],self.time)
            self.t2 = othertime.findNearest(trange[1],self.time)

        if self.t1==self.t2:
            self.t2+=1

        self.nt = self.t2 - self.t1


    def get_indices(self,xrange,yrange,zrange):
        """ Find the domain indices"""

        # Check if X/Y are 1D or gridded

        if self.X.ndim==1:
            if xrange==None:
                self.x1=0
                self.x2=self.X.shape[0]
            else:
                self.x1 = self.find_nearest_1d(self.X,xrange[0])
                self.x2 = self.find_nearest_1d(self.X,xrange[1])

            if yrange==None:
                self.y1=0
                self.y2 = self.Y.shape[0]
            else:
                y1 = self.find_nearest_1d(self.Y,yrange[0])
                y2 = self.find_nearest_1d(self.Y,yrange[1])
                self.y1 = min([y1,y2])
                self.y2 = max([y1,y2])
        
        elif self.X.ndim==2:
            if xrange or yrange == None:
                self.x1 = 0
                self.x2 = self.X.shape[1]
                self.y1 = 0
                self.y2 = self.Y.shape[0]
            else:
                ind = self.find_nearest_2d(\
                    self.X,self.Y,[xrange[0],yrange[0]])
                self.y1,self.x1 = ind[0][0],ind[0][1]
                ind = self.find_nearest_2d(\
                    self.X,self.Y,[xrange[1],yrange[1]])
                self.y2,self.x2 = ind[0][0],ind[0][1]

        self.nx = self.x2-self.x1
        self.ny = self.y2-self.y1

        if not self.zcoord == None:
            if zrange==None:
                self.z1=0
                self.z2=self.Z.shape[0]
            else:
                self.z1 = self.find_nearest_1d(zrange[0],self.Z)
                self.z2 = self.find_nearest_1d(zrange[1],self.Z)
                if self.z1==self.z2:
                    self.z2+=1
            self.nz = self.z2 - self.z1


    def get_data(self,varname,xrange,yrange,zrange,trange):
        """
        Download the actual data

        Set xrange/yrange/zrange/trange=None to get all
        
        """
        
        self.get_indices(xrange,yrange,zrange)
        self.get_time_indices(trange)

        # Store the local coordinates
        self.localtime = self.time[self.t1:self.t2]

        if not self.multifile:
            return self.get_data_singlefile(varname,self._nc,self.t1,self.t2)
        else:
            return self.get_data_multifile(varname) 

    def get_data_multifile(self,varname):
        """
        Download the data from multiple files and stick in one array

        """
        if self.zcoord==None:
            sz = (self.nt,self.ny,self.nx)
        else:
            sz = (self.nt,nz,self.ny,self.nx)

        data = np.zeros(sz)
        tindex,ncurl = self._MF(self.localtime.tolist())

        ncold = ncurl[0]
        nc = Dataset(ncold)
        ii=-1
        for t,name in zip(tindex,ncurl):
            ii+=1
            if not ncold == name:
                nc.close()
                nc = Dataset(name)
            else:
                ncold = name

            data[ii,...] = self.get_data_singlefile(varname,nc,t,t+1)

        nc.close()
        return data
            

    def get_data_singlefile(self,varname,nc,t1,t2):  
        """
        Download the data from a single file
        """
        ndim = nc.variables[varname].ndim

        if ndim==3:
            data = nc.variables[varname]\
                [t1:t2,self.y1:self.y2,self.x1:self.x2]

        elif ndim==4:
            data = nc.variables[varname]\
                [t1:t2,self.z1:self.z2,\
                self.y1:self.y2,self.x1:self.x2]

            self.localZ = self.Z[self.z1:self.z2]

        if self.X.ndim==1:
            self.localX = self.X[self.x1:self.x2]
            self.localY = self.Y[self.y1:self.y2]
        else:
            self.localX = self.X[self.y1:self.y2,self.x1:self.x2]
            self.localY = self.Y[self.y1:self.y2,self.x1:self.x2]

        return data

    def find_nearest_1d(self,data,value):
        dist = np.abs(data-value)
        return np.argwhere(dist==dist.min())[0][0]

    def find_nearest_2d(self,xdata,ydata,xy):
        dist = np.sqrt( (xy[0]-xdata)**2. + (xy[1]-ydata)**2.)
        return np.argwhere(dist==dist.min())


class GetDAP(object):
    """
    High level class for downloading data
    """
    ncurl = None # file location (can be url or local)
    type = 'ocean' # Ocean or atmoshpere
    multifile = False

    # Ocean variable names (name of variables in file)
    u = 'u'
    v = 'v'
    temp = 'temp'
    salt = 'salt'
    ssh = 'ssh'
    # Default variables to extract
    oceanvars = ['ssh','u','v','temp','salt']

    # Atmosphere variable names
    uwind = 'uwind'
    vwind = 'vwind'
    tair = 'tair'
    pair = 'pair'
    rh = 'rh'
    cloud = 'cloud'
    rain = 'rain'

    # Default variables to extract
    airvars = ['uwind','vwind','tair','pair','rh','cloud','rain']

    tformat = '%Y%m%d.%H%M%S'

    # Location of cached grid file
    gridfile = None

    def __init__(self,vars=None,**kwargs):
        """
        Initialise the variables to outfile
        """
        self.__dict__.update(kwargs)

        # Open the file
        if not self.multifile:
            self._nc = Dataset(self.ncurl)
        else:
            self._nc = Dataset(self.ncurl[0])
            self._ncfiles = self.ncurl

        # Set the variables that need downloading
        if vars==None and self.type=='ocean':
            self.vars = self.oceanvars
        elif vars==None and self.type=='atmosphere':
            self.vars = self.airvars
        else:
            self.vars = vars

        # For each variable:
        # Get the coordinate variables
        
        # Initiate the griddap class for each variable
        #   (this does most of the work)
        self.init_var_grids()

    def __call__(self,xrange,yrange,trange,zrange=None,outfile=None):    
        """
        Download the data, save if outfile is specified
        """

        trange=self.check_trange(trange)

        # Load the data for all variable and write
        for vv in self.vars:
            data = self.load_data(vv,xrange,yrange,zrange,trange)

            # Write the data
            if not outfile == None:
                self.write_var(outfile,vv,data)


    def load_data(self,varname,xr,yr,zr,tr):
        ncvar = getattr(self,varname)
        dap = getattr(self,'_%s'%varname)
        print 'Loading variable: %s...'%ncvar

        return dap.get_data(ncvar,xr,yr,zr,tr)

    def init_var_grids(self):
        """
        Initialise each variable into an object with name "_name"
        """
        for vv in self.vars:
            print 'loading grid data for variable: %s...'%getattr(self,vv)
            attrname = '_%s'%vv
            timecoord,xcoord,ycoord,zcoord = \
                self.get_coord_names(getattr(self,vv))
            
            if self.multifile:
                nc = self._ncfiles
            else:
                nc = self._nc
                
            dap = GridDAP(nc,xcoord=xcoord,ycoord=ycoord,\
                timecoord=timecoord,zcoord=zcoord,\
                gridfile=self.gridfile,multifile=self.multifile)

            setattr(self,attrname,dap)
 

    def check_trange(self,trange):
        if isinstance(trange[0],str):
            t1 = datetime.strptime(trange[0],self.tformat)
            t2 = datetime.strptime(trange[1],self.tformat)
            return [t1,t2]
        elif isinstance(trange[0],datetime):
            return trange
        else:
            raise Exception, 'unknown time format'

    def get_coord_names(self,varname):
        """
        Try to automatically workout the coordinate names

        This is probably not very robust
        """
        dims = self._nc.variables[varname].dimensions

        ndims = len(dims)
        try:
            coordinates = self._nc.variables[varname].coordinates

            # Dimensions are usually (time,depth,y,x)
            # Coordinates are generally (lon,lat,date)
            coordlist = coordinates.split()
            coordlist.reverse()

            timecoord,zcoord,ycoord,xcoord = \
                dims[0],dims[1],coordlist[-2],coordlist[-1]
        except:
            if ndims==4:
                timecoord,zcoord,ycoord,xcoord = dims
            else:
                 timecoord,ycoord,xcoord = dims
        if ndims==3:
            zcoord=None

        return timecoord,xcoord,ycoord,zcoord


    def write_var(self,outfile,var,data):

        # Check if the file is open and/or exists
        if os.path.isfile(outfile):
            print 'File exists - appending...'
            self._outnc = Dataset(outfile,'a')
        else:
            if not self.__dict__.has_key('_outnc'):
                print '\tOpening %s'%outfile
                self._outnc = Dataset(outfile,'w')

        localvar = var
        remotevar = getattr(self,var)
        dapobj = getattr(self,'_%s'%var)

        # Get the coordinate variables
        timecoord,xcoord,ycoord,zcoord = \
            self.get_coord_names(remotevar)
     

        # Write the coordinate variables (will not be overwritten)
        self.create_ncvar_fromvarname(xcoord,dimsizes=dapobj.localX.shape)
        self._outnc.variables[xcoord][:]=dapobj.localX
        self.create_ncvar_fromvarname(ycoord,dimsizes=dapobj.localY.shape)
        self._outnc.variables[ycoord][:]=dapobj.localY

        is3d = False
        if not zcoord==None:
            self.create_ncvar_fromvarname(zcoord,dimsizes=dapobj.localZ.shape)
            self._outnc.variables[zcoord][:]=dapobj.localZ
            is3d = True

        # Write the time variable
        self.write_time_var(dapobj.localtime,timecoord)

        # Create the actual variable
        if is3d:
            dimsize = (None, dapobj.nz, dapobj.ny, dapobj.nx)
        else:
            dimsize = (None, dapobj.ny, dapobj.nx)
        self.create_ncvar_fromvarname(remotevar,dimsizes=dimsize)


        # Write to the variable
        self._outnc.variables[remotevar]\
            [self.tindex[0]:self.tindex[1],...] = data


    def write_time_var(self,tout,timecoord):
        """
        Write the time variale into the local file

        if it exists append the data into the right spot
        """
        # Check if local file has time variable
        if hasvar(self._outnc,timecoord):
            # Load the local time variable
            localtime = gettime(self._outnc,timecoord)

            # Get the index for insertion
            t = self._outnc.variables[timecoord]

            try:
                # This doesn't work if there is only one time value
                t1 = date2index(tout[0],t)
            except:
                if tout[0]>localtime[0]:
                    t1 = 1
                else:
                    t1 = 0
            self.tindex = [t1,t1+tout.shape[0]]
            print self.tindex

            # Check to see if we need to write
            if localtime.shape[0]<self.tindex[1]:
                self._outnc.variables[timecoord][self.tindex[0]:self.tindex[1]] = \
                    date2num(tout,t.units)

        else:
            # Create the variable (note unlimited dimension size)
            self.create_ncvar_fromvarname(timecoord,dimsizes=(None,))

            # Write the data (converts to netcdf units)
            t = self._outnc.variables[timecoord]
            self._outnc.variables[timecoord][:] = \
                date2num(tout,t.units)

            self.tindex = [0,tout.shape[0]]


    def create_gridfile_fromvar(self,outfile,varname):
        """
        Create a local file with the grid information from the variable
        """
        if os.path.isfile(outfile):
            print 'File exists - appending...'
            self._outnc = Dataset(outfile,'a')
        else:
            if not self.__dict__.has_key('_outnc'):
                print '\tOpening %s'%outfile
                self._outnc = Dataset(outfile,'w')

        print 'Creating local netcdf file with grid data...'

        # Get the variables to store in the grid file
        timecoord,xcoord,ycoord,zcoord = self.get_coord_names(varname)

        # Create the grid netcdf file
        self.create_ncfile(outfile)

        # Create the variables
        # Download the coordinate variable data
        if not zcoord==None:
            self.create_ncvar_fromvarname(zcoord)
            print 'Downloading Z...'
            self._outnc.variables[zcoord][:] = self._nc.variables[zcoord][:]

        #X
        self.create_ncvar_fromvarname(xcoord)
        print 'Downloading X...'
        self._outnc.variables[xcoord][:] = self._nc.variables[xcoord][:]

        #Y
        self.create_ncvar_fromvarname(ycoord)
        print 'Downloading Y...'
        self._outnc.variables[ycoord][:] = self._nc.variables[ycoord][:]
        

    def create_ncvar_fromvarname(self,varname,dimsizes=None):
        """
        Create a netcdf variable using a remote variable name

        dimsizes can be passed as a list which will override the native grid

        """

        if hasvar(self._outnc,varname):
            print 'Variable exists - exiting.'
            return

        # List the dimensions
        dims = self._nc.variables[varname].dimensions
        for ii,dim in enumerate(dims):
            # Create the dimension if it doesn't exist
            if not hasdim(self._outnc,dim):
                if dimsizes==None:
                    dimsize = self._nc.dimensions[dim].__len__()
                else:
                    dimsize = dimsizes[ii]

                self._outnc.createDimension(dim,dimsize)

        # Create the variable
        print varname
        V = self._nc.variables[varname]
        tmp= self._outnc.createVariable(varname, V.dtype, V.dimensions)

        # Copy the attributes
        for aa in V.ncattrs():
            tmp.setncattr(aa,getattr(V,aa))

    def create_ncfile(self,ncfile):
        nc = Dataset(ncfile,'w')
        nc.Title = '%s model data'%(self.type)
        nc.url = '%s'%(self.ncurl)
        
        self._outnc=nc


class MFncdap(object):
    """
    Multi-file class for opendap netcdf files
    
    MFDataset module is not compatible with opendap data 
    """
    
    timevar = 'time'
    tformat = '%Y%m%d.%H%M%S'
    
    def __init__(self,ncfilelist,**kwargs):
        
        self.__dict__.update(kwargs)
        print 'Retrieving the time information from files...'
        
        self._timelookup = {}
        timeall = []
        #self.time = np.zeros((0,))
        
        for f in ncfilelist:
            print f
            nc = Dataset(f)
            t = nc.variables[self.timevar]
            time = num2date(t[:],t.units).tolist()
            nc.close()
            
            #self.timelookup.update({f:time})
            timeall.append(time)

            # Create a dictionary of time and file
            for ii,tt in enumerate(time):
                tstr = datetime.strftime(tt,self.tformat)
                if self._timelookup.has_key(tstr):
                    # Overwrite the old key
                    self._timelookup[tstr]=(f,ii)
                else:
                    self._timelookup.update({tstr:(f,ii)})

            
        #self.time = np.asarray(self.time)
        timeall = np.array(timeall)
        self.time = np.unique(timeall)
            
    def __call__(self,time):
        """
        Return the filenames and time index of the closest time
        """
        
        fname = []
        tind =[]
        for t in time:
            tstr = datetime.strftime(t,self.tformat)
            f,ii = self._timelookup[tstr]
            fname.append(f)
            tind.append(ii)

#        for t in time:
#            flag=1
#            for f in self.timelookup.keys():
#
#                if t >= self.timelookup[f][0] and t<=self.timelookup[f][-1]:
#                    #print 'Found tstep %s'%datetime.strftime(t,'%Y-%m-%d %H:%M:%S')
#                    tind.append(othertime.findNearest(t,self.timelookup[f][:]))
#                    fname.append(f)
#                    flag=0
#                    continue
#
##            if flag:
##                print 'Warning - could not find matching file for time:%s'%datetime.strptime(t,'%Y-%m-%d %H:%M:%S')
##                tind.append(-1)
##                fname.append(-1)
        
        return tind, fname
 
