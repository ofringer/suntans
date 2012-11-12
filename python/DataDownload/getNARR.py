# -*- coding: utf-8 -*-
"""
Tool for downloading data from the North American Regional Reanalysis (NARR) product

Note that this requires the netCDF4 module to be compiled with opendap support (libcurl)

Created on Thu Nov 08 17:22:40 2012

@author: mrayson

"""

import numpy as np
from netCDF4 import Dataset
from osgeo import osr
from datetime import datetime, timedelta

from pydap.client import open_url

class getNARR(object):
    """
    General class for initialising a NARR data object
    
    Example Usage:
        tstart = '20100101'
        tend = '20100103'
        bbox = [-95.40,-94.49,28.8,29.9]
        varname = 'Relative_humidity'
        
        narr = getNARR(tstart,tend,bbox)
        RH = narr(varname)
        
    """
    basedir = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/'
    basefile = 'narr-a_221_'
    verbose = True
    
    def __init__(self,tstart,tend,bbox,**kwargs):
        
        self.__dict__.update(kwargs)
        self.tstart=tstart
        self.tend=tend
        self.bbox=bbox
        
        # Get the filenames and the coordinates
        self.getFileNames()
        
        self.getNARRlatlon()
        
        self.findXYindices()
        
        self.getArraySize()
        
    def __call__(self,varname):
        """
        Call function to extract variable, "varname"
        """
        
        self.getDimInfo(varname)
        
        return self.getDataPyDap(varname)
    
    def getData(self,vv):
        """
        Downloads the variable data using the NetCDF4 module
        
        This module has memory leaks if not compiled with the correct netcdf/hdf/curl libraries
        
        Safer to use pydap.
        """
        
        if self.verbose:
            print 'Retrieving variable: %s...'%vv
        
        data = np.zeros((self.nt,self.ny,self.nx))
        
        tt=-1
        for ff in self.grbfiles:
            tt+=1
            if self.verbose:
                print '     File: %s...'%ff
                
            nc = Dataset(ff,'r')
            
            # get the height coordinate for 4D arrays
            if self.ndim == 4:
                 # dimension order [time, z, y, x]
                data[tt,:,:] = nc.variables[vv][:,0,self.y1:self.y2,self.x1:self.x2]
            elif self.ndim == 3:
                data[tt,:,:] = nc.variables[vv][:,self.y1:self.y2,self.x1:self.x2]
                
            nc.close()    
        
        return data
        
    def getDataPyDap(self,vv):
        """
        Downloads the variable data using the pydap library
        """
        
        if self.verbose:
            print 'Retrieving variable: %s...'%vv
        
        data = np.zeros((self.nt,self.ny,self.nx))
        
        tt=-1
        for ff in self.grbfiles:
            tt+=1
            if self.verbose:
                print '     File: %s...'%ff
                
            nc = open_url(ff)
            
            # get the height coordinate for 4D arrays
            if self.ndim == 4:
                 # dimension order [time, z, y, x]
                data[tt,:,:] = nc[vv][:,0,self.y1:self.y2,self.x1:self.x2]
            elif self.ndim == 3:
                data[tt,:,:] = nc[vv][:,self.y1:self.y2,self.x1:self.x2]
                   
        
        return data
        
    def getDimInfo(self,vv):
        """
        Gets the variable dimension data
        """
        nc = Dataset(self.grbfiles[0],'r')
        self.dimdata = nc.variables[vv].dimensions
        self.ndim  = nc.variables[vv].ndim
        
        # get the height coordinate for 4D arrays
        if self.ndim == 4:
             # dimension order [time, z, y, x]
            self.z = nc.variables[self.dimdata[1]][0]
        else:
            self.z=0.0
        nc.close()
    
    def getArraySize(self):
        
        self.nt = len(self.time)
        self.nx = self.x2-self.x1
        self.ny = self.y2-self.y1
        
    def getFileNames(self):
        """
        Return the filenames for each time step between the two dates
        """

        # build a list of timesteps
        t1 = datetime.strptime(self.tstart,'%Y%m%d')
        t2 = datetime.strptime(self.tend,'%Y%m%d')
        
        self.time = []
        t0=t1
        while t0 < t2:
            self.time.append(t0)
            t0 += timedelta(hours=3)
            
        # Build the list of filenames
        self.grbfiles=[]
        for tt in self.time:
            filestr='%s%s/%s/%s%s_000.grb'%(self.basedir,datetime.strftime(tt,'%Y%m'),\
            datetime.strftime(tt,'%Y%m%d'),self.basefile,datetime.strftime(tt,'%Y%m%d_%H%M'))
            self.grbfiles.append(filestr)
            
    def findXYindices(self):
        """
        Find the indices of the lower left and upper right corners from the grid
        """
        
        dist = np.sqrt( (self.lon - self.bbox[0])**2 + (self.lat - self.bbox[2])**2)
        j = np.argwhere(dist==dist.min())
        self.y1=j[0,0]
        self.x1=j[0,1]
        
        dist = np.sqrt( (self.lon - self.bbox[1])**2 + (self.lat - self.bbox[3])**2)
        j = np.argwhere(dist==dist.min())
        self.y2=j[0,0]
        self.x2=j[0,1]
        
        # Resize the lat lon arrays
        self.lat=self.lat[self.y1:self.y2,self.x1:self.x2]
        self.lon=self.lon[self.y1:self.y2,self.x1:self.x2]

    def getNARRlatlon(self):
        """
        Returns the NARR grid in lat/lon coordinates
        
        Need to convert from their Lambert Conformal projection to WGS84
        
        *** THIS NEEDS TO BE CHECKED THOROUGHLY!!! ***
        """
        
        nc = Dataset(self.grbfiles[0],'r')
        
        x = nc.variables['x'][:]*1000.0
        y = nc.variables['y'][:]*1000.0
        
        nc.close()
        
        # NARR grid projection information:
        # -----------------------------------    
        #grid_mapping_name: lambert_conformal_conic
        #standard_parallel: 50.0
        #longitude_of_central_meridian: -107.0
        #latitude_of_projection_origin: 50.0
        #earth_shape: spherical
        #earth_radius: 6367470.21484375
        #GRIB_param_Dx: 32463.0
        #GRIB_param_Dy: 32463.0
        #GRIB_param_GDSkey: 55295
        #GRIB_param_La1: 1.0
        #GRIB_param_Latin1: 50.0
        #GRIB_param_Latin2: 50.0
        #GRIB_param_Lo1: -145.5
        #GRIB_param_LoV: -107.0
        #GRIB_param_NpProj: true
        #GRIB_param_Nx: 349
        #GRIB_param_Ny: 277
        #GRIB_param_ProjFlag: 0
        #GRIB_param_ResCompFlag: 8
        #GRIB_param_SpLat: 0.0
        #GRIB_param_SpLon: 0.0
        #GRIB_param_VectorComponentFlag: gridRelative
        #GRIB_param_Winds: Relative
        #GRIB_param_grid_name: Lambert_Conformal
        #GRIB_param_grid_radius_spherical_earth: 6367.47
        #GRIB_param_grid_shape: spherical
        #GRIB_param_grid_shape_code: 0
        #GRIB_param_grid_type: 3
        #GRIB_param_grid_units: m
        #GRIB_param_scanning_mode: 64
        #DODS:
        #  strlen: 0
        
        # Define the input coordinate system
        # See this website:
            # http://trac.osgeo.org/proj/wiki/GenParms
        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=lcc +lat_0=50.0 +lon_0=-107.0 +lon_1=-145.5 +lat_1=50.0 +a=6367470.21484375 +b=6367470.21484375 +units=m')    
        
        # set the output coordinate system
        srsout = osr.SpatialReference()
        srsout.SetWellKnownGeogCS( "WGS84" );
        
        # define the transformation object
        ct = osr.CoordinateTransformation(srs, srsout)
        
        nx = len(x)
        ny = len(y)
        self.lon = np.zeros((ny,nx))
        self.lat = np.zeros((ny,nx)) 
        for ii in range(nx):
            for jj in range(ny):
                X,Y,z =  ct.TransformPoint(x[ii],y[jj])
                self.lon[jj,ii]=X
                self.lat[jj,ii]=Y
        
    
         

#############
# Testing stuff

#tstart = '20100101'
#tend = '20100102'
#bbox = [-95.40,-94.49,28.8,29.9]
##
#
#narr=getNARR(tstart,tend,bbox)
#
##RH = narr('Relative_humidity')
#cloud = narr('Total_cloud_cover')