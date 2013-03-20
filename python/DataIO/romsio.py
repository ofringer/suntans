# -*- coding: utf-8 -*-
"""
Tools for dealing with ROMS model output

See Octant project as well

Created on Fri Mar 08 15:09:46 2013

@author: mrayson
"""

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from scipy import interpolate
from interpXYZ import interpXYZ
import othertime

try:
    from octant.slice import isoslice
except:
    print 'Warning - could not import octant package.'

import pdb

class roms_grid(object):
    """
    Class for ROMS grid
    """
    def __init__(self,ncfile):
        self.grdfile = ncfile
        
        self.readGrid()
        
    def readGrid(self):
        """
        Read in the main grid variables from the grid netcdf file
        """
        
        try: 
            nc = MFDataset(self.grdfile, 'r')
        except:
            nc = Dataset(self.grdfile, 'r')  
        
        self.angle = nc.variables['angle'][:]
        self.lon_rho =  nc.variables['lon_rho'][:]
        self.lat_rho =  nc.variables['lat_rho'][:]
        self.lon_psi =  nc.variables['lon_psi'][:]
        self.lat_psi =  nc.variables['lat_psi'][:]
        self.lon_u =  nc.variables['lon_u'][:]
        self.lon_v =  nc.variables['lon_v'][:]
        self.lat_u =  nc.variables['lat_u'][:]
        self.lat_v =  nc.variables['lat_v'][:]
        self.h = nc.variables['h'][:]
        self.f = nc.variables['f'][:]
        
        self.mask_rho =  nc.variables['mask_rho'][:]
        self.mask_psi =  nc.variables['mask_psi'][:]
        self.mask_u =  nc.variables['mask_u'][:]
        self.mask_v =  nc.variables['mask_v'][:]

        nc.close()
        
    def findNearset(self,x,y,grid='rho'):
        """
        Return the J,I indices of the nearst grid cell to x,y
        """
        
        if grid == 'rho':
            lon = self.lon_rho
            lat = self.lat_rho
        elif grid == 'u':
            lon = self.lon_u
            lat = self.lat_u
        elif grid =='v':
            lon = self.lon_v
            lat = self.lat_v
        elif grid =='psi':
            lon = self.lon_psi
            lat = self.lat_psi
            
        dist = np.sqrt( (lon - x)**2 + (lat - y)**2)
        
        return np.argwhere(dist==dist.min())
        
    def utmconversion(self,lon,lat,utmzone,isnorth):
        """
        Convert the ROMS grid to utm coordinates
        """
        from maptools import ll2utm
        
        M,N = lon.shape
        
        xy = ll2utm(np.hstack((np.reshape(lon,(M*N,1)),np.reshape(lat,(M*N,1)))),utmzone,north=isnorth)
        
        return np.reshape(xy[:,0],(M,N)), np.reshape(xy[:,1],(M,N)) 
        
        
class ROMS(roms_grid):
    """
    General class for reading and plotting ROMS model output
    """
    
    varname = 'zeta'
    JRANGE = None
    IRANGE = None
    
    zlayer = False # True load z layer, False load sigma layer
    K = [0] # Layer to extract, 0 bed, -1 surface, -99 all
    tstep = [0] # - 1 last step, -99 all time steps
    
    clim = None # Plot limits
    
        
    def __init__(self,romsfile,**kwargs):
        
        self.__dict__.update(kwargs)
        
        self.romsfile = romsfile
        
        # Load the grid        
        roms_grid.__init__(self,self.romsfile)
        
        # Open the netcdf object
        self._openNC()
        
        # Load the time information
        self._loadTime()
        
        # Check the spatial indices of the variable
        self._loadVarCoords()
        
        self._checkCoords(self.varname)
        
        # Check the vertical coordinates                
        self._readVertCoords()
        
        self._checkVertCoords(self.varname)
        
    
    def listCoordVars(self):
        """
        List all of the variables that have the 'coordinate' attribute
        """
        
        for vv in self.nc.variables.keys():
            if hasattr(self.nc.variables[vv],'coordinates'):
                print '%s - %s'%(vv,self.nc.variables[vv].long_name)
            
    def loadData(self,varname=None,tstep=None):
        """
        Loads model data from the netcdf file
        """
        
        if varname == None:
            varname=self.varname
            self._checkCoords(varname)
        else:
            self._checkCoords(varname)
            
            if self.ndim == 4:
                self._checkVertCoords(varname)
            
        if tstep == None:
            tstep = self.tstep
            
            
        if self.ndim == 2:
            data = self.nc.variables[varname][self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
        elif self.ndim == 3:
            data = self.nc.variables[varname][tstep,self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
        elif self.ndim == 4:
            data = self.nc.variables[varname][tstep,self.K,self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
        
        if self.ndim == 4 and self.zlayer==True:
            # Slice along z layers
            print 'Extracting data along z-coordinates...'
            dataz = np.zeros((len(tstep),)+self.X.shape)
            
            for ii,tt in enumerate(tstep):
                Z = self.calcDepth(zeta=self.loadData(varname='zeta',tstep=[tt]))
                dataz[ii,:,:] = isoslice(data[ii,:,:,:].squeeze(),Z,self.Z)
                
            data = dataz
        
        self._checkCoords(self.varname)            
        # Reduce rank
        return data.squeeze()
        
    def calcDepth(self,zeta=None):
        """
        Calculates the depth array for the current variable
        """
        
        return get_depth(self.S,self.C,self.hc,self.h,zeta=zeta, Vtransform=self.Vtransform).squeeze()
     
    def pcolor(self,data=None,titlestr=None,**kwargs):
        """
        Pcolor plot of the data in variable
        """
        
        if data==None:
            data=self.loadData()
            
        if self.clim==None:
            clim=[data.min(),data.max()]
        else:
            clim=self.clim
            
        fig = plt.gcf()
        ax = fig.gca()
        
        p1 = plt.pcolor(self.X,self.Y,data,vmin=clim[0],vmax=clim[1],**kwargs)
        
        ax.set_aspect('equal')
        plt.colorbar(p1)
        
        if titlestr==None:
            plt.title(self._genTitle(self.tstep[0]))
        else:
            plt.title(titlestr)
        
        return p1
    
    def contourf(self, data=None, VV=20, titlestr=None,**kwargs):
        """
        contour plot of the data in variable
        """
        
        if data==None:
            data=self.loadData()
            
        if self.clim==None:
            clim=[data.min(),data.max()]
        else:
            clim=self.clim
            
        fig = plt.gcf()
        ax = fig.gca()
        
        p1 = plt.contourf(self.X,self.Y,data,VV,vmin=clim[0],vmax=clim[1],**kwargs)
        
        ax.set_aspect('equal')
        plt.colorbar(p1)
        
        if titlestr==None:
            plt.title(self._genTitle(self.tstep[0]))
        else:
            plt.title(titlestr)
        
        return p1
        
    def contourbathy(self,clevs=np.arange(0,3000,100),**kwargs):
                
        p1 = plt.contour(self.lon_rho,self.lat_rho,self.h,clevs,**kwargs)
        
        return p1
    
    def _genTitle(self,tstep):
        """
        Generates a title for plots
        """
        if self.zlayer:
            titlestr = '%s [%s]\nz: %6.1f m, %s'%(self.long_name,self.units,self.Z,datetime.strftime(self.time[tstep],'%d-%b-%Y %H:%M:%S'))            
        else:
            titlestr = '%s [%s]\nsigma[%d], %s'%(self.long_name,self.units,self.K[0],datetime.strftime(self.time[tstep],'%d-%b-%Y %H:%M:%S'))            
            
        return titlestr

    def _checkCoords(self,varname):
        """
        Load the x and y coordinates of the present variable, self.varname
        """
        #print 'updating coordinate info...'
        C = self.varcoords[varname].split()        
        self.ndim = len(C)
        
        self.xcoord = C[0]
        self.ycoord = C[1]
          
        if self.JRANGE==None:
            self.JRANGE = [0,self[self.xcoord].shape[0]]
        if self.IRANGE==None:
            self.IRANGE = [0,self[self.xcoord].shape[1]]
            
        # Check the dimension size
        if self.JRANGE[1] > self[self.xcoord].shape[0]:
            print 'Warning JRANGE outside of size range. Setting equal size.'
            self.JRANGE[1] = self[self.xcoord].shape[0]
            
        if self.IRANGE[1] > self[self.xcoord].shape[1]:
            print 'Warning JRANGE outside of size range. Setting equal size.'
            self.IRANGE[1] = self[self.xcoord].shape[1]
            
        if not self.__dict__.has_key('X'):
            self.X = self[self.xcoord][self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
            self.Y = self[self.ycoord][self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
            
        # Load the long_name and units from the variable
        self.long_name = self.nc.variables[varname].long_name
        self.units = self.nc.variables[varname].units
    
    def _checkVertCoords(self,varname):
        """
        Load the vertical coordinate info
        """
        
        # First put K into a list
        #if not type(self.K)=='list':
        #    self.K = [self.K]
        try:
            K = self.K[0] #  a list
            self.K = self.K
        except:
            # not a list
            self.K = [self.K]
         
        C = self.varcoords[varname].split() 
        
        ndim = len(C)
        
        if ndim == 4:
            self.zcoord = C[2] 
            self.Nz = len(self[self.zcoord])
            
            if self.K[0] == -99:
                self.K = range(0,self.Nz)
                
            if self.zlayer==True: # Load all layers when zlayer is true
                self.Z = self.K[0]
                self.K = range(0,self.Nz)
                
        if self.zcoord == 's_rho':
            self.S = self.s_rho[self.K]
            self.C = self.Cs_r[self.K]
        elif self.zcoord == 's_w':
            self.S = self.s_w[self.K]
            self.C = self.Cs_w[self.K]
            
        
    def _readVertCoords(self):
        """
        Read the vertical coordinate information
        """
        nc = self.nc
        
        self.Cs_r = nc.variables['Cs_r'][:]
        self.Cs_w = nc.variables['Cs_w'][:]
        self.s_rho = nc.variables['s_rho'][:]
        self.s_w = nc.variables['s_w'][:]
        self.hc = nc.variables['hc'][:]
        self.Vstretching = nc.variables['Vstretching'][:]
        self.Vtransform = nc.variables['Vtransform'][:]

        
    def _loadVarCoords(self):
        """
        Load the variable coordinates into a dictionary
        """
        self.varcoords={}
        for vv in self.nc.variables.keys():
            if hasattr(self.nc.variables[vv],'coordinates'):
                self.varcoords.update({vv:self.nc.variables[vv].coordinates})
                
    def _openNC(self):
        """
        Load the netcdf object
        """
        try: 
            self.nc = MFDataset(self.romsfile, 'r')
        except:
            self.nc = Dataset(self.romsfile, 'r')
            
    def _loadTime(self):
        """
        Load the netcdf time as a vector datetime objects
        """
        #nc = Dataset(self.ncfile, 'r', format='NETCDF4') 
        nc = self.nc
        t = nc.variables['ocean_time']
        self.time = num2date(t[:],t.units)  
        
    def __getitem__(self,y):
        x = self.__dict__.__getitem__(y)
        return x
        
    def __setitem__(self,key,value):
        
        if key == 'varname':
            self.varname=value
            self._checkCoords(value)
        else:
            self.__dict__[key]=value
            
            

class roms_subset(roms_grid):
    """
    Class for subsetting ROMS output
    """
    gridfile = None
    
    def __init__(self,ncfiles,bbox,timelims,**kwargs):
        self.__dict__.update(kwargs)
        
        if self.gridfile==None:
            self.gridfile=ncfiles[0]
        
        self.ncfiles = ncfiles
        self.x0 = bbox[0]
        self.x1 = bbox[1]
        self.y0 = bbox[2]
        self.y1 = bbox[3]
        
        # Step 1) Find the time steps
        self.t0 = datetime.strptime(timelims[0],'%Y%m%d%H%M%S')
        self.t1 = datetime.strptime(timelims[1],'%Y%m%d%H%M%S')
        
        # Multifile object        
        ftime = MFncdap(ncfiles,timevar='ocean_time')
        
        ind0 = othertime.findNearest(self.t0,ftime.time)
        ind1 = othertime.findNearest(self.t1,ftime.time)
        
        self.time = ftime.time[ind0:ind1]
        self.tind,self.fname = ftime(self.time) # list of time indices and corresponding files
        
        self.Nt = len(self.tind)
        
        # Step 2) Subset the grid variables
        roms_grid.__init__(self,self.gridfile)
        
        self.SubsetGrid()
        
        # Step 3) Read the vertical coordinate variables
        self.ReadVertCoords()
        

    def SubsetGrid(self):
        """
        Subset the grid variables
        """
        #Find the grid indices
        ind = self.findNearset(self.x0,self.y0)
        
        self.J0=ind[0][0]
        self.I0=ind[0][1]
        
        ind = self.findNearset(self.x1,self.y1)
        self.J1=ind[0][0]
        self.I1=ind[0][1]
        
        # Define the dimensions
        M = self.J1-self.J0
        N = self.I1-self.I0
        
        self.eta_rho = M
        self.xi_rho = N
        self.eta_psi = M-1
        self.xi_psi = N-1
        self.eta_u = M-1
        self.xi_u = N
        self.eta_v = M
        self.xi_v = N-1
        
        # Subset the horizontal coordinates
        self.lon_rho = self.lon_rho[self.J0:self.J1,self.I0:self.I1]
        self.lat_rho = self.lat_rho[self.J0:self.J1,self.I0:self.I1]
        self.mask_rho = self.mask_rho[self.J0:self.J1,self.I0:self.I1]

        
        self.lon_psi = self.lon_psi[self.J0:self.J1-1,self.I0:self.I1-1]
        self.lat_psi = self.lat_psi[self.J0:self.J1-1,self.I0:self.I1-1]
        self.mask_psi = self.mask_psi[self.J0:self.J1-1,self.I0:self.I1-1]
        
        self.lon_u = self.lon_u[self.J0:self.J1-1,self.I0:self.I1]
        self.lat_u = self.lat_u[self.J0:self.J1-1,self.I0:self.I1]
        self.mask_u = self.mask_u[self.J0:self.J1-1,self.I0:self.I1]
        
        self.lon_v = self.lon_v[self.J0:self.J1,self.I0:self.I1-1]
        self.lat_v = self.lat_v[self.J0:self.J1,self.I0:self.I1-1]
        self.mask_v = self.mask_v[self.J0:self.J1,self.I0:self.I1-1]
        
        self.h = self.h[self.J0:self.J1,self.I0:self.I1]
        self.angle = self.angle[self.J0:self.J1,self.I0:self.I1]

    def ReadVertCoords(self):
        """
        
        """
        nc = Dataset(self.fname[0])
        
        self.Cs_r = nc.variables['Cs_r'][:]
        #self.Cs_w = nc.variables['Cs_w'][:]
        self.s_rho = nc.variables['s_rho'][:]
        #self.s_w = nc.variables['s_w'][:]
        self.hc = nc.variables['hc'][:]
        self.Vstretching = nc.variables['Vstretching'][:]
        self.Vtransform = nc.variables['Vtransform'][:]

        nc.close()
        
    def ReadData(self,tstep):
        """
        Reads the data from the file for the present time step
        """
           
        fname = self.fname[tstep]
        t0 = self.tind[tstep]
        
        print 'Reading data at time: %s...'%datetime.strftime(self.time[tstep],'%Y-%m-%d %H:%M:%S')        
        
        nc = Dataset(fname)
        
        self.ocean_time = nc.variables['ocean_time'][t0]
        
        self.zeta = nc.variables['zeta'][t0,self.J0:self.J1,self.I0:self.I1]
        self.temp = nc.variables['temp'][t0,:,self.J0:self.J1,self.I0:self.I1]
        self.salt = nc.variables['salt'][t0,:,self.J0:self.J1,self.I0:self.I1]
        self.u = nc.variables['u'][t0,:,self.J0:self.J1-1,self.I0:self.I1]
        self.v = nc.variables['v'][t0,:,self.J0:self.J1,self.I0:self.I1-1]
        
        nc.close()
        
    def Writefile(self,outfile,verbose=True):
        """
        Writes subsetted grid and coordinate variables to a netcdf file
        
        Code modified from roms.py in the Octant package
        """
        self.outfile = outfile
        
        Mp, Lp = self.lon_rho.shape
        M, L = self.lon_psi.shape
        
        N = self.s_rho.shape[0] # vertical layers
        
        xl = self.lon_rho[self.mask_rho==1.0].ptp()
        el = self.lat_rho[self.mask_rho==1.0].ptp()
        
        # Write ROMS grid to file
        nc = Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        nc.Description = 'ROMS subsetted history file'
        nc.Author = ''
        nc.Created = datetime.now().isoformat()
        nc.type = 'ROMS HIS file'
        
        nc.createDimension('xi_rho', Lp)
        nc.createDimension('xi_u', Lp)
        nc.createDimension('xi_v', L)
        nc.createDimension('xi_psi', L)
        
        nc.createDimension('eta_rho', Mp)
        nc.createDimension('eta_u', M)
        nc.createDimension('eta_v', Mp)
        nc.createDimension('eta_psi', M)
        
        nc.createDimension('s_rho', N)
        nc.createDimension('ocean_time', None)      
        
        nc.createVariable('xl', 'f8', ())
        nc.variables['xl'].units = 'meters'
        nc.variables['xl'] = xl
        
        nc.createVariable('el', 'f8', ())
        nc.variables['el'].units = 'meters'
        nc.variables['el'] = el
        
        nc.createVariable('spherical', 'S1', ())
        nc.variables['spherical'] = 'F'
        
        def write_nc_var(var, name, dimensions, units=None):
            nc.createVariable(name, 'f8', dimensions)
            if units is not None:
                nc.variables[name].units = units
            nc.variables[name][:] = var
            if verbose:
                print ' ... wrote ', name
                
        def create_nc_var(name, dimensions, units=None):
            nc.createVariable(name, 'f8', dimensions)
            if units is not None:
                nc.variables[name].units = units
            if verbose:
                print ' ... wrote ', name
        
        # Grid variables
        write_nc_var(self.angle, 'angle', ('eta_rho', 'xi_rho'))
        write_nc_var(self.h, 'h', ('eta_rho', 'xi_rho'), 'meters')
        
        write_nc_var(self.mask_rho, 'mask_rho', ('eta_rho', 'xi_rho'))
        write_nc_var(self.mask_u, 'mask_u', ('eta_u', 'xi_u'))
        write_nc_var(self.mask_v, 'mask_v', ('eta_v', 'xi_v'))
        write_nc_var(self.mask_psi, 'mask_psi', ('eta_psi', 'xi_psi'))
        
        write_nc_var(self.lon_rho, 'lon_rho', ('eta_rho', 'xi_rho'), 'meters')
        write_nc_var(self.lat_rho, 'lat_rho', ('eta_rho', 'xi_rho'), 'meters')
        write_nc_var(self.lon_u, 'lon_u', ('eta_u', 'xi_u'), 'meters')
        write_nc_var(self.lat_u, 'lat_u', ('eta_u', 'xi_u'), 'meters')
        write_nc_var(self.lon_v, 'lon_v', ('eta_v', 'xi_v'), 'meters')
        write_nc_var(self.lat_v, 'lat_v', ('eta_v', 'xi_v'), 'meters')
        write_nc_var(self.lon_psi, 'lon_psi', ('eta_psi', 'xi_psi'), 'meters')
        write_nc_var(self.lat_psi, 'lat_psi', ('eta_psi', 'xi_psi'), 'meters')
        
        # Vertical coordinate variables
        write_nc_var(self.s_rho, 's_rho', ('s_rho'))
        write_nc_var(self.Cs_r, 'Cs_r', ('s_rho'))

        write_nc_var(self.hc, 'hc', ())
        write_nc_var(self.Vstretching, 'Vstretching', ())
        write_nc_var(self.Vtransform, 'Vtransform', ())
        
        # Create the data variables
        create_nc_var('ocean_time',('ocean_time'),'seconds since 1970-01-01 00:00:00')
        create_nc_var('zeta',('ocean_time','eta_rho','xi_rho'),'meter')
        create_nc_var('salt',('ocean_time','s_rho','eta_rho','xi_rho'),'psu')
        create_nc_var('temp',('ocean_time','s_rho','eta_rho','xi_rho'),'degrees C')
        create_nc_var('u',('ocean_time','s_rho','eta_u','xi_u'),'meter second-1')
        create_nc_var('v',('ocean_time','s_rho','eta_v','xi_v'),'meter second-1')
        
        nc.close()
        
    def Writedata(self, tstep):
        
        nc = Dataset(self.outfile, 'a')
        
        nc.variables['ocean_time'][tstep]=self.ocean_time
        nc.variables['zeta'][tstep,:,:]=self.zeta
        nc.variables['salt'][tstep,:,:,:]=self.salt
        nc.variables['temp'][tstep,:,:,:]=self.temp
        nc.variables['u'][tstep,:,:,:]=self.u
        nc.variables['v'][tstep,:,:,:]=self.v

        nc.close()
        
    def Go(self):
        """
        Downloads and append each time step to a file
        """
        for ii in range(0,self.Nt):
            self.ReadData(ii)
            self.Writedata(ii)
            
        print '##################\nDone!\n##################'
        
class roms_interp(roms_grid):
    """
    Class for intperpolating ROMS output in space and time
    """
    
    utmzone = 15
    isnorth = True
    
    # Interpolation options
    interpmethod='idw' # 'nn', 'idw', 'kriging', 'griddata'
    NNear=3
    p = 1.0 #  power for inverse distance weighting
    # kriging options
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0

    def __init__(self,romsfile, xi, yi, zi, timei, **kwargs):
        
        self.__dict__.update(kwargs)
        
        self.romsfile = romsfile
        self.xi = xi
        self.yi = yi
        self.zi = zi
        self.timei = timei
        
        # Step 1) Find the time steps
        self.t0 = timei[0]
        self.t1 = timei[-1]
        
        # Multifile object        
        ftime = MFncdap(self.romsfile,timevar='ocean_time')
        
        ind0 = othertime.findNearest(self.t0,ftime.time)
        ind1 = othertime.findNearest(self.t1,ftime.time)
        
        self.time = ftime.time[ind0:ind1+1]
        self.tind,self.fname = ftime(self.time) # list of time indices and corresponding files
        
        # Step 2) Prepare the grid variables for the interpolation class
        roms_grid.__init__(self,self.romsfile[0])
        
        # rho points
        x,y = self.utmconversion(self.lon_rho,self.lat_rho,self.utmzone,self.isnorth)
        self.xy_rho = np.vstack((x[self.mask_rho==1],y[self.mask_rho==1])).T
        
        # uv point (averaged onto interior rho points)
        self.mask_uv = self.mask_rho[0:-1,0:-1]
        x = x[0:-1,0:-1]
        y = y[0:-1,0:-1]
        self.xy_uv = np.vstack((x[self.mask_uv==1],y[self.mask_uv==1])).T
        
        # Step 3) Build the interpolants for rho and uv points
        self.xy_out = np.hstack((xi,yi))  

        
        self.Frho = interpXYZ(self.xy_rho,self.xy_out,method=self.interpmethod,NNear=self.NNear,\
            p=self.p,varmodel=self.varmodel,nugget=self.nugget,sill=self.sill,vrange=self.vrange)
        
        self.Fuv = interpXYZ(self.xy_uv,self.xy_out,method=self.interpmethod,NNear=self.NNear,\
            p=self.p,varmodel=self.varmodel,nugget=self.nugget,sill=self.sill,vrange=self.vrange)
        
        # Read the vertical coordinate
        self.ReadVertCoords()
        # Dimesions sizes
        self.Nx = self.xy_out.shape[0]
        self.Nz = self.zi.shape[0]
        self.Nt = len(self.timei)
        
        self.Nz_roms = self.s_rho.shape[0]
        self.Nt_roms = self.time.shape[0]
        
    def interp(self,zinterp='linear',tinterp=3):
        """
        Performs the interpolation in this order:
            1) Interpolate onto the horizontal coordinates
            2) Interpolate onto the vertical coordinates
            3) Interpolate onto the time coordinates
        """
        
        # Initialise the output arrays @ roms time step
        zetaroms, temproms, saltroms, uroms, vroms = self.initArrays(self.Nt_roms,self.Nx,self.Nz)
        
        tempold = np.zeros((self.Nz_roms,self.Nx))
        saltold = np.zeros((self.Nz_roms,self.Nx))
        uold = np.zeros((self.Nz_roms,self.Nx))
        vold = np.zeros((self.Nz_roms,self.Nx))

        # Interpolate h
        h = self.Frho(self.h[self.mask_rho==1])
        
        # Loop through each time step            
        for tstep in range(0,self.Nt_roms):
        
            # Read all variables
            self.ReadData(tstep)
                    
            # Interpolate zeta
            zetaroms[tstep,:] = self.Frho(self.zeta[self.mask_rho==1])
            
            # Interpolate other 3D variables
            for k in range(0,self.Nz_roms):
                tmp = self.temp[k,:,:]
                tempold[k,:] = self.Frho(tmp[self.mask_rho==1])
                
                tmp = self.salt[k,:,:]
                saltold[k,:] = self.Frho(tmp[self.mask_rho==1])
                
                tmp = self.u[k,:,:]
                uold[k,:] = self.Fuv(tmp[self.mask_uv==1])
                
                tmp = self.v[k,:,:]
                vold[k,:] = self.Fuv(tmp[self.mask_uv==1])
    
            # Calculate depths (zeta dependent)
	    #zroms = get_depth(self.s_rho,self.Cs_r,self.hc, h, zetaroms[tstep,:], Vtransform=self.Vtransform)
	    zroms = get_depth(self.s_rho,self.Cs_r,self.hc, h, zeta=zetaroms[tstep,:], Vtransform=self.Vtransform)
    
            # Interpolate vertically
            for ii in range(0,self.Nx):
                y = tempold[:,ii]
                Fz = interpolate.interp1d(zroms[:,ii],y,kind=zinterp,bounds_error=False,fill_value=y[0])
                temproms[tstep,:,ii] = Fz(self.zi)
                
                y = saltold[:,ii]
                Fz = interpolate.interp1d(zroms[:,ii],y,kind=zinterp,bounds_error=False,fill_value=y[0])
                saltroms[tstep,:,ii] = Fz(self.zi)
                
                y = uold[:,ii]
                Fz = interpolate.interp1d(zroms[:,ii],y,kind=zinterp,bounds_error=False,fill_value=y[0])
                uroms[tstep,:,ii] = Fz(self.zi)
                
                y = vold[:,ii]
                Fz = interpolate.interp1d(zroms[:,ii],y,kind=zinterp,bounds_error=False,fill_value=y[0])
                vroms[tstep,:,ii] = Fz(self.zi)
            
        # End time loop
        
        # Initialise the output arrays @ output time step
        
        # Interpolate temporally
        if self.Nt_roms > 1:
            troms = othertime.SecondsSince(self.time)
            tout = othertime.SecondsSince(self.timei)
            
            Ft = interpolate.interp1d(troms,zetaroms,axis=0,kind=tinterp,bounds_error=False)
            zetaout = Ft(tout)
            
            Ft = interpolate.interp1d(troms,temproms,axis=0,kind=tinterp,bounds_error=False)
            tempout = Ft(tout)
            
            Ft = interpolate.interp1d(troms,saltroms,axis=0,kind=tinterp,bounds_error=False)
            saltout = Ft(tout)
            
            Ft = interpolate.interp1d(troms,uroms,axis=0,kind=tinterp,bounds_error=False)
            uout = Ft(tout)
            
            Ft = interpolate.interp1d(troms,vroms,axis=0,kind=tinterp,bounds_error=False)
            vout = Ft(tout)
        else:
            zetaout = zetaroms
            tempout = temproms
            saltout = saltroms
            uout = uroms
            vout = vroms
        
        return zetaout, tempout, saltout, uout, vout
        
    def initArrays(self,Nt,Nx,Nz):
        
        zetaout = np.zeros((Nt,Nx))
        tempout = np.zeros((Nt,Nz,Nx))
        saltout = np.zeros((Nt,Nz,Nx))
        uout = np.zeros((Nt,Nz,Nx))
        vout = np.zeros((Nt,Nz,Nx))
        
        return zetaout, tempout, saltout, uout, vout
            
            
    def ReadData(self,tstep):
        """
        Reads the data from the file for the present time step
        """
           
        fname = self.fname[tstep]
        t0 = self.tind[tstep]
        
        print 'Interpolating data at time: %s of %s...'%(datetime.strftime(self.time[tstep],'%Y-%m-%d %H:%M:%S'),\
        datetime.strftime(self.time[-1],'%Y-%m-%d %H:%M:%S'))
        
        nc = Dataset(fname)
        
        self.ocean_time = nc.variables['ocean_time'][t0]
        
        self.zeta = nc.variables['zeta'][t0,:,:]
        self.temp = nc.variables['temp'][t0,:,:,:]
        self.salt = nc.variables['salt'][t0,:,:,:]
        u = nc.variables['u'][t0,:,:,:]
        v = nc.variables['v'][t0,:,:,:]
    
        nc.close()
        
        # Rotate the vectors
        self.u,self.v = rotateUV( (u[...,:,0:-1]+u[...,:,1::])*0.5,(v[...,0:-1,:]+v[...,1::,:])*0.5,self.angle[0:-1,0:-1])
    
    def ReadVertCoords(self):
        """
        
        """
        nc = Dataset(self.romsfile[0])
        
        self.Cs_r = nc.variables['Cs_r'][:]
        #self.Cs_w = nc.variables['Cs_w'][:]
        self.s_rho = nc.variables['s_rho'][:]
        #self.s_w = nc.variables['s_w'][:]
        self.hc = nc.variables['hc'][:]
        self.Vstretching = nc.variables['Vstretching'][:]
        self.Vtransform = nc.variables['Vtransform'][:]

        nc.close()
    
    
class MFncdap(object):
    """
    Multi-file class for opendap netcdf files
    
    MFDataset module is not compatible with opendap data 
    """
    
    timevar = 'time'
    
    def __init__(self,ncfilelist,**kwargs):
        
        self.__dict__.update(kwargs)
        
        self.timelookup = {}
        self.time = np.zeros((0,))
        for f in ncfilelist:
            print f
            nc = Dataset(f)
            t = nc.variables[self.timevar]
            time = num2date(t[:],t.units)
            nc.close()
            
            self.timelookup.update({f:time})
            self.time = np.hstack((self.time,np.asarray(time)))
            
        self.time = np.asarray(self.time)
    
            
    def __call__(self,time):
        """
        Return the filenames and time index of the closest time
        """
        
        fname = []
        tind =[]
        for t in time:
            flag=1
            for f in self.timelookup.keys():

                if t >= self.timelookup[f][0] and t<=self.timelookup[f][-1]:
#                    print 'Found tstep %s'%datetime.strptime(t,'%Y-%m-%d %H:%M:%S')
                    tind.append(othertime.findNearest(t,self.timelookup[f][:]))
                    fname.append(f)
                    flag=0

#            if flag:
#                print 'Warning - could not find matching file for time:%s'%datetime.strptime(t,'%Y-%m-%d %H:%M:%S')
#                tind.append(-1)
#                fname.append(-1)
        
        return tind, fname
            
    
def get_depth(S,C,hc,h,zeta=None, Vtransform=1):
    """
    Calculates the sigma coordinate depth
    """
    if zeta == None:
        zeta = 0.0*h
        
    N = len(S)
    
    #Nj,Ni = np.size(h)
    shp = (N,)+h.shape    
    
    z = np.zeros(shp)
    
    if Vtransform == 1:
        for k in range(0,N):
            z0 = (S[k]-C[k])*hc + C[k]*h
            z[k,...] = z0 + (zeta *(1.0 + z0/h))
    elif Vtransform == 2:
        for k in range(0,N):
            z0 = (hc*S[k]+C[k]*h)/(hc+h)
            z[k,...] = zeta + (zeta+h)*z0
    
    return z
        
def rotateUV(uroms,vroms,ang):
    """
    Rotates ROMS output vectors to cartesian u,v
    """
    
    u = uroms*np.cos(ang) - vroms*np.sin(ang)
    v = uroms*np.sin(ang) + vroms*np.cos(ang)
    
    return u,v
    
###############        
## Testing
##grdfile = 'http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc'
#grdfile = 'C:\\Projects\\GOMGalveston\\MODELLING\\ROMS\\txla_grd_v4_new.nc'
##grd = roms_grid(grdfile)
#
##ncfiles = ['http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_%04d.nc'%i for i in range(1,3)]
##MF = MFncdap(ncfiles,timevar='ocean_time')
##
##tsteps = [datetime(2003,2,16)+timedelta(hours=i*4) for i in range(0,24)]
##tind,fname = MF(tsteps)
#
#ncfiles = ['http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_%04d.nc'%i for i in range(100,196)]
#timelims = ('20090501000000','20090701000000')
##timelims = ('20090501000000','20090502000000')
#bbox = [-95.53,-94.25,28.3,30.0]
#
#roms = roms_subset(ncfiles,bbox,timelims,gridfile=grdfile)
#outfile = 'C:\\Projects\\GOMGalveston\\MODELLING\\ROMS\\txla_subset_HIS_MayJun2009.nc'
#roms.Writefile(outfile)
#roms.Go()
#
##roms2 = roms_subset([outfile],bbox,timelims)
