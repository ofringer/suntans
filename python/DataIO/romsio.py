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
# Private modules
from interpXYZ import interpXYZ
import othertime
from timeseries import timeseries

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
        
        varnames = ['angle','lon_rho','lat_rho','lon_psi','lat_psi','lon_u','lat_u',\
        'lon_v','lat_v','h','f','mask_rho','mask_psi','mask_u','mask_v','pm','pn']
        
        for vv in varnames:
            try:
                setattr(self,vv,nc.variables[vv][:])
            except:
                print 'Cannot find variable: %s'%vv

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
        nc.createDimension('xi_u', L)
        nc.createDimension('xi_v', Lp)
        nc.createDimension('xi_psi', L)
        
        nc.createDimension('eta_rho', Mp)
        nc.createDimension('eta_u', Mp)
        nc.createDimension('eta_v', M)
        nc.createDimension('eta_psi', M)
        
        nc.createDimension('s_rho', N)
        nc.createDimension('s_w', N+1)
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
            
        
        # Grid variables
        write_nc_var(self.angle, 'angle', ('eta_rho', 'xi_rho'))
        write_nc_var(self.h, 'h', ('eta_rho', 'xi_rho'), 'meters')
        write_nc_var(self.f, 'f', ('eta_rho', 'xi_rho'), 'seconds-1')
        
        write_nc_var(self.mask_rho, 'mask_rho', ('eta_rho', 'xi_rho'))
        write_nc_var(self.mask_u, 'mask_u', ('eta_u', 'xi_u'))
        write_nc_var(self.mask_v, 'mask_v', ('eta_v', 'xi_v'))
        write_nc_var(self.mask_psi, 'mask_psi', ('eta_psi', 'xi_psi'))
        
        write_nc_var(self.lon_rho, 'lon_rho', ('eta_rho', 'xi_rho'), 'degrees')
        write_nc_var(self.lat_rho, 'lat_rho', ('eta_rho', 'xi_rho'), 'degrees')
        write_nc_var(self.lon_u, 'lon_u', ('eta_u', 'xi_u'), 'degrees')
        write_nc_var(self.lat_u, 'lat_u', ('eta_u', 'xi_u'), 'degrees')
        write_nc_var(self.lon_v, 'lon_v', ('eta_v', 'xi_v'), 'degrees')
        write_nc_var(self.lat_v, 'lat_v', ('eta_v', 'xi_v'), 'degrees')
        write_nc_var(self.lon_psi, 'lon_psi', ('eta_psi', 'xi_psi'), 'degrees')
        write_nc_var(self.lat_psi, 'lat_psi', ('eta_psi', 'xi_psi'), 'degrees')
        
        write_nc_var(self.pm, 'pm', ('eta_rho', 'xi_rho'), 'degrees')
        write_nc_var(self.pn, 'pn', ('eta_rho', 'xi_rho'), 'degrees')
        
        # Vertical coordinate variables
        write_nc_var(self.s_rho, 's_rho', ('s_rho'))
        write_nc_var(self.s_w, 's_w', ('s_w'))
        write_nc_var(self.Cs_r, 'Cs_r', ('s_rho'))
        write_nc_var(self.Cs_w, 'Cs_w', ('s_w'))

        write_nc_var(self.hc, 'hc', ())
        write_nc_var(self.Vstretching, 'Vstretching', ())
        write_nc_var(self.Vtransform, 'Vtransform', ())

        
    def nc_add_dimension(self,outfile,name,length):
        """
        Add a dimension to an existing netcdf file
        """
        nc = Dataset(outfile, 'a')
        nc.createDimension(name, length)
        nc.close()
        
    def nc_add_var(self,outfile,data,name,dimensions,units=None,long_name=None,coordinates=None):
        """
        Add a new variable and write the data
        """        
        nc = Dataset(outfile, 'a')
        nc.createVariable(name, 'f8', dimensions)
        if units is not None:
            nc.variables[name].units = units
        if coordinates is not None:
            nc.variables[name].coordinates = coordinates
        if long_name is not None:
            nc.variables[name].long_name = long_name
        
        nc.variables[name][:] = data
            
        nc.close()
        
    def nc_add_varnodata(self,outfile,name,dimensions,units=None,long_name=None,coordinates=None):
        """
        Add a new variable and doesn't write the data
        """        
        nc = Dataset(outfile, 'a')
        nc.createVariable(name, 'f8', dimensions)
        if units is not None:
            nc.variables[name].units = units
        if coordinates is not None:
            nc.variables[name].coordinates = coordinates
        if long_name is not None:
            nc.variables[name].long_name = long_name
        
            
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
        try:
            self._loadTime()
        except:
            print 'No time variable.'
                
        
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
            
        if self.ndim==1:
            data = self.nc.variables[varname][tstep]
        elif self.ndim == 2:
            data = self.nc.variables[varname][self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
        elif self.ndim == 3:
            data = self.nc.variables[varname][tstep,self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
        elif self.ndim == 4:
            data = self.nc.variables[varname][tstep,self.K,self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
        
        if self.ndim == 4 and self.zlayer==True:
            # Slice along z layers
            print 'Extracting data along z-coordinates...'
            dataz = np.zeros((len(tstep),)+self.Z.shape+self.X.shape)
            
            for ii,tt in enumerate(tstep):
                #Z = self.calcDepth(zeta=self.loadData(varname='zeta',tstep=[tt]))
                Z = self.calcDepth()
                if len(Z.shape) > 1:
                    dataz[ii,:,:] = isoslice(data[ii,:,:,:].squeeze(),Z,self.Z)
                else:
                    # Isoslice won't work on 1-D arrays
                    F = interpolate.interp1d(Z,data[ii,:,:,:].squeeze(),bounds_error=False)
                    dataz[ii,:,:] = F(self.Z)[:,np.newaxis,np.newaxis]
                
            data = dataz
        
        self._checkCoords(self.varname)            
        # Reduce rank
        return data.squeeze()
        
    def loadTimeSeries(self,x,y,z=None,varname=None,trange=None):
        """
        Load a time series at point x,y
        
        Set z=None to load all layers, else load depth
        """
        
        if varname == None:
            self.varname = self.varname
        else:
            self.varname = varname
            self._checkCoords(self.varname)

                     
        if self.ndim == 4:
            self._checkVertCoords(self.varname)
            if z == None:
                self.zlayer=False
                self.K = [-99]
            else:
                self.zlayer=True
                self.K = [z]
                
        if trange==None:
            tstep=np.arange(0,self.Nt)
        
        # Set the index range to grab            
        JI = self.findNearset(x,y,grid=self.gridtype)
        
        self.JRANGE = [JI[0][0], JI[0][0]+1]
        self.IRANGE = [JI[0][1], JI[0][1]+1]
                
        if self.zlayer:
            Zout = z
        else:
            # Return the depths at each time step
            h = self.h[self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]].squeeze()
            zeta=self.loadData(varname='zeta',tstep=tstep)
            h = h*np.ones(zeta.shape)
            Zout = get_depth(self.S,self.C,self.hc,h,zeta=zeta, Vtransform=self.Vtransform).squeeze()
     
        return self.loadData(varname=varname,tstep=tstep), Zout
        
    def calcDepth(self,zeta=None):
        """
        Calculates the depth array for the current variable
        """
        #h = self.h[self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]].squeeze()
        if self.gridtype == 'rho':
            h = self.h
        elif self.gridtype == 'psi':
            h = 0.5 * (self.h[1:,1:] + self.h[0:-1,0:-1])
        elif self.gridtype == 'u':
            h = 0.5 * (self.h[:,1:] + self.h[:,0:-1])
        elif self.gridtype == 'v':
            h = 0.5 * (self.h[1:,:] + self.h[0:-1,:])
         
        return get_depth(self.S,self.C,self.hc,h,zeta=zeta, Vtransform=self.Vtransform).squeeze()
        
    def depthInt(self,var,grid='rho',cumulative=False):
        """
        Depth-integrate data in variable, var (array [Nz, Ny, Nx])
        
        Set cumulative = True for cumulative integration i.e. for pressure calc.
        """
        
        sz = var.shape
        if not sz[0] == self.Nz:
            raise Exception, 'length of dimension 0 must equal %d (currently %d)'%(self.Nz,sz[0])
        
        if not len(sz)==3:
            raise Exception, 'only 3-D arrays are supported.'
          
        if grid == 'rho':
            h = self.h
        elif grid == 'psi':
            h = 0.5 * (self.h[1:,1:] + self.h[0:-1,0:-1])
        elif grid == 'u':
            h = 0.5 * (self.h[:,1:] + self.h[:,0:-1])
        elif grid == 'v':
            h = 0.5 * (self.h[1:,:] + self.h[0:-1,:])
            
        z_w = get_depth(self.s_w,self.Cs_w,self.hc,h,Vtransform=self.Vtransform).squeeze()
        
        dz = np.diff(z_w,axis=0)
        
        if cumulative:
            return np.cumsum(dz*var,axis=0)
        else:
            return np.sum(dz*var,axis=0)
        
    def depthAvg(self,var,grid='rho'):
        """
        Depth-average data in variable, var (array [Nz, Ny, Nx])
        """
        
        sz = var.shape
        if not sz[0] == self.Nz:
            raise Exception, 'length of dimension 0 must equal %d (currently %d)'%(self.Nz,sz[0])
        
        if not len(sz)==3:
            raise Exception, 'only 3-D arrays are supported.'
          
        if grid == 'rho':
            h = self.h
        elif grid == 'psi':
            h = 0.5 (self.h[1:,1:] + self.h[0:-1,0:-1])
        elif grid == 'u':
            h = 0.5 (self.h[:,1:] + self.h[:,0:-1])
        elif grid == 'v':
            h = 0.5 (self.h[1:,:] + self.h[0:-1,:])
            
        z_w = get_depth(self.s_w,self.Cs_w,self.hc,h,Vtransform=self.Vtransform).squeeze()
        
        dz = np.diff(z_w,axis=0)
        
        return np.sum(dz*var,axis=0) / h
        
    def areaInt(self,var,grid='rho'):
        """
        Calculate the area integral of var
        """
        if grid == 'rho':
            dx = 1.0/self.pm
            dy = 1.0/self.pn
        elif grid == 'psi':
            dx = 1.0/(0.5*(self.pm[1:,1:] + self.pm[0:-1,0:-1]))
            dy = 1.0/(0.5*(self.pn[1:,1:] + self.pn[0:-1,0:-1]))
        elif grid == 'u':
            dx = 1.0/(0.5 * (self.pm[:,1:] + self.pm[:,0:-1]))
            dy = 1.0/(0.5 * (self.pn[:,1:] + self.pn[:,0:-1]))
        elif grid == 'v':
            dx = 0.5 * (self.pm[1:,:] + self.pm[0:-1,:])
            dy = 0.5 * (self.pn[1:,:] + self.pn[0:-1,:])
            
        A = dx*dy
        return np.sum(var*A)
        
    def gradZ(self,var,grid='rho',cumulative=False):
        """
        Depth-gradient of data in variable, var (array [Nz, Ny, Nx])
        
        """
        
        sz = var.shape
        #print sz
        if not sz[0] == self.Nz:
            raise Exception, 'length of dimension 0 must equal %d (currently %d)'%(self.Nz,sz[0])
        
#        if not len(sz)==3:
#            raise Exception, 'only 3-D arrays are supported.'
          
#        if grid == 'rho':
#            h = self.h
#        elif grid == 'psi':
#            h = 0.5 (self.h[1:,1:] + self.h[0:-1,0:-1])
#        elif grid == 'u':
#            h = 0.5 (self.h[:,1:] + self.h[:,0:-1])
#        elif grid == 'v':
#            h = 0.5 (self.h[1:,:] + self.h[0:-1,:])

        
        h = self.h[self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]].squeeze()            
        z_r = get_depth(self.s_rho,self.Cs_r,self.hc,h,Vtransform=self.Vtransform).squeeze()
                
        dz = np.diff(z_r,axis=0)
        dz_mid = 0.5 * (dz[1:,...] + dz[0:-1,...]) # N-2
        
        var_mid =  0.5 * (var[1:,...] + var[0:-1,...])
        
        dv_dz = np.zeros(sz)
        # 2-nd order mid-points
        dv_dz[1:-1,...] = (var_mid[1:,...] - var_mid[0:-1,...]) / dz_mid
        # 1st order end points
        dv_dz[0,...] = (var[1,...] - var[0,...]) / dz[0,...]
        dv_dz[-1,...] = (var[-1,...] - var[-2,...]) / dz[-1,...]
        
        return dv_dz

    def MLD(self,tstep,thresh=-0.006,z_max=-20.0):
        """
        Mixed layer depth calculation
        
        thresh is the density gradient threshold
        z_max is the min mixed layer depth 
        """
        
        # Load the density data
        self.K=[-99]
        
        drho_dz=self.gradZ(self.loadData(varname='rho',tstep=tstep))
        
        # Mask drho_dz where z >= z_max
        z = self.calcDepth()
        mask = z >= z_max
        
        drho_dz[mask] = 0.0
        #
        
        mld_ind = np.where(drho_dz <= thresh)
        
        zout = -99999.0*np.ones(z.shape)
        zout[mld_ind[0],mld_ind[1],mld_ind[2]] = z[mld_ind[0],mld_ind[1],mld_ind[2]]
        
        mld = np.max(zout,axis=0)
        
        # Isoslice averages when there is more than one value
        #mld = isoslice(z,drho_dz,thresh)
        
        mld = np.max([mld,-self.h],axis=0)
        
        return mld

    def MLDmask(self,mld,grid='rho'):
        """
        Compute a 3D mask for variables beneath the mixed layer
        """
        if grid == 'rho':
            h = self.h
        elif grid == 'psi':
            h = 0.5 * (self.h[1:,1:] + self.h[0:-1,0:-1])
        elif grid == 'u':
            h = 0.5 * (self.h[:,1:] + self.h[:,0:-1])
        elif grid == 'v':
            h = 0.5 * (self.h[1:,:] + self.h[0:-1,:])
            
        z = get_depth(self.s_rho,self.Cs_r,self.hc,h,Vtransform=self.Vtransform).squeeze()
        
        mask = np.zeros(z.shape)
        for jj in range(mld.shape[0]):
            for ii in range(mld.shape[1]):
                ind = z[:,jj,ii] >= mld[jj,ii]
                if np.size(ind)>0:
                    mask[ind,jj,ii]=1.0
        
        return mask
        
        

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
    
    def contourf(self, data=None, clevs=20, titlestr=None,**kwargs):
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
        
        p1 = plt.contourf(self.X,self.Y,data,clevs,vmin=clim[0],vmax=clim[1],**kwargs)
        
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
    
    def getTstep(self,tstart,tend,timeformat='%Y%m%d.%H%M'):
        """
        Returns a vector of the time indices between tstart and tend
        
        tstart and tend can be string with format=timeformat ['%Y%m%d.%H%M' - default]
        
        Else tstart and tend can be datetime objects
        """
        
        try:
            t0 = datetime.strptime(tstart,timeformat)
            t1 = datetime.strptime(tend,timeformat)
        except:
            # Assume the time is already in datetime format
            t0 = tstart
            t1 = tend
                        
        n1 = othertime.findNearest(t0,self.time)
        n2 = othertime.findNearest(t1,self.time)
        
        return np.arange(n1,n2)
    
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
        
        if self.ndim==1:
            return
            
        self.xcoord = C[0]
        self.ycoord = C[1]
          
        if self.JRANGE==None:
            self.JRANGE = [0,self[self.xcoord].shape[0]+1]
        if self.IRANGE==None:
            self.IRANGE = [0,self[self.xcoord].shape[1]+1]
            
        # Check the dimension size
        if self.JRANGE[1] > self[self.xcoord].shape[0]+1:
            print 'Warning JRANGE outside of size range. Setting equal size.'
            self.JRANGE[1] = self[self.xcoord].shape[0]+1
            
        if self.IRANGE[1] > self[self.xcoord].shape[1]+1:
            print 'Warning JRANGE outside of size range. Setting equal size.'
            self.IRANGE[1] = self[self.xcoord].shape[1]+1
            
        self.X = self[self.xcoord][self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
        self.Y = self[self.ycoord][self.JRANGE[0]:self.JRANGE[1],self.IRANGE[0]:self.IRANGE[1]]
            
        # Load the long_name and units from the variable
        try:
            self.long_name = self.nc.variables[varname].long_name
        except:
            self.long_name = varname
            
        self.units = self.nc.variables[varname].units
        
        # Set the grid type
        if self.xcoord[-3:]=='rho':
            self.gridtype='rho'
            self.mask=self.mask_rho
        elif self.xcoord[-3:]=='n_u':
            self.gridtype='u'
            self.mask=self.mask_u
        elif self.xcoord[-3:]=='n_v':
            self.gridtype='v'
            self.mask=self.mask_v
    
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
                self.Z = np.array(self.K)
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
            self.nc = MFDataset(self.romsfile)
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
        self.Nt = np.size(self.time)
        
    def __getitem__(self,y):
        x = self.__dict__.__getitem__(y)
        return x
        
    def __setitem__(self,key,value):
        
        if key == 'varname':
            self.varname=value
            self._checkCoords(value)
        else:
            self.__dict__[key]=value

            
class roms_timeseries(ROMS, timeseries):
    """
    Class for loading a timeseries object from ROMS model output
    """
    
    IJ = False
    varname = 'u'
    zlayer=False
    
    def __init__(self,ncfile,XY,z=None,**kwargs):
        """
        Loads a time series from point X,Y. Set z = None (default) to load all layers
        
        if self.IJ = True, loads index X=I, Y=J directly
        """
        self.__dict__.update(kwargs)
        self.XY = XY
        self.z = z
        
        # Initialise the class
        ROMS.__init__(self,ncfile,varname=self.varname,K=[-99])
        
        self.tstep = range(0,self.Nt) # Load all time steps
        
        self.update()
        
        
    def update(self):
        """
        Updates the class
        """
        # 
        self._checkCoords(self.varname)
        # Load  I and J indices from the coordinates
        self.setIJ(self.XY)
        
        # Load the vertical coordinates
        if not self.z == None:
            self.zlayer = True
            
        if self.zlayer == False:
            self.Z = self.calcDepth()
        else:
            self.Z = self.z
            
        # Load the data into a time series object
        timeseries.__init__(self,self.time[self.tstep],self.loadData())
        
    def contourf(self,clevs=20,**kwargs):
        """
        z-t contour plot of the time series
        """
        h1 = plt.contourf(self.time[self.tstep],self.Z,self.y.T,clevs,**kwargs)
        
        #plt.colorbar()
        
        plt.xticks(rotation=17)
        return h1
        
        
    def setIJ(self,xy):
        if self.IJ:
            I0 = xy[0]
            J0 = xy[1]
        else:
            ind = self.findNearset(xy[0],xy[1],grid=self.gridtype)
            J0=ind[0][0]
            I0=ind[0][1]
            
        self.JRANGE = [J0,J0+1]
        self.IRANGE = [I0,I0+1]
        
    def __setitem__(self,key,value):
        
        if key == 'varname':
            self.varname=value
            self.update()
            
        elif key == 'XY':
            self.XY = value
            self.update()            
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
        
    def interp(self,zinterp='linear',tinterp='linear',setUV=True,seth=True):
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
            if seth:
                zetaroms[tstep,:] = self.Frho(self.zeta[self.mask_rho==1])
            
            # Interpolate other 3D variables
            for k in range(0,self.Nz_roms):
                tmp = self.temp[k,:,:]
                tempold[k,:] = self.Frho(tmp[self.mask_rho==1])
                
                tmp = self.salt[k,:,:]
                saltold[k,:] = self.Frho(tmp[self.mask_rho==1])
                
                if setUV:
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
                
                if setUV:
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
	    print 'Temporally interpolating ROMS variables...'
            troms = othertime.SecondsSince(self.time)
            tout = othertime.SecondsSince(self.timei)
            if seth:
                print '\tzeta...'
                Ft = interpolate.interp1d(troms,zetaroms,axis=0,kind=tinterp,bounds_error=False)
                zetaout = Ft(tout)
            else:
                zetaout=-1

            print '\ttemp...'
            Ft = interpolate.interp1d(troms,temproms,axis=0,kind=tinterp,bounds_error=False)
            tempout = Ft(tout)
            print '\tsalt...'
            Ft = interpolate.interp1d(troms,saltroms,axis=0,kind=tinterp,bounds_error=False)
            saltout = Ft(tout)
            if setUV:
                print '\tu...'
                Ft = interpolate.interp1d(troms,uroms,axis=0,kind=tinterp,bounds_error=False)
                uout = Ft(tout)
                print '\tv...'
                Ft = interpolate.interp1d(troms,vroms,axis=0,kind=tinterp,bounds_error=False)
                vout = Ft(tout)
            else:
                uout = vout = -1
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
