# -*- coding: utf-8 -*-
"""
    Tools for interpolating irregularly spaced data onto irregularly spaced points
    
    Largely a wrapper for other interpolation methods such as scipy.griddata and
    scipy.Rbf (radial basis functions).
    
    
    
"""

import gzip
from scipy import spatial
import numpy as np
from maptools import ll2utm, readShpBathy, readraster
from kriging import kriging
from netCDF4 import Dataset

from scipy.interpolate import griddata
import time
import matplotlib.pyplot as plt

# testing stuff
import pdb
#import sunpy

class interpXYZ(object):
    "Class for interpolating xyz data"""
    
    ## Properties ##
   
    
    # Interpolation options
    method = 'nn' # One of 'nn', 'idw', 'kriging', 'griddata'
    maxdist=np.inf
    NNear = 3 # Number of points to include in interpolation (only applicable to idw and kriging)
    p = 1.0 # power for inverse distance weighting
    
    # kriging options
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0
    
    def __init__(self,XY,XYout,**kwargs):
        
        self.__dict__.update(kwargs)
        self.XY = XY
        self.XYout = XYout
        
        if self.method=='nn':
            #print 'Building DEM with Nearest Neighbour interpolation...'
            self._nearestNeighbour()
        
        elif self.method=='idw':
            #print 'Building DEM with Inverse Distance Weighted Interpolation...'
            self._invdistweight()
            
        elif self.method=='kriging':
            #print 'Building DEM with Kriging Interpolation...'
            self._krig()
            
        #elif self.method=='griddata':
            #print 'Building DEM using griddata...'
            # Do nothing
            
        else:
            print 'Error - Unknown interpolation type: %s.'%self.method
        
    def __call__(self,Zin):
        """
        
        """
        
        self.Zin = Zin
        
        if self.method in ['nn','idw','kriging']:
            #print 'Building DEM with Nearest Neighbour interpolation...'
            self.Z = self.Finterp(Zin)
            
        elif self.method=='griddata':
            #print 'Building DEM using griddata...'
            self._griddata()
            
        else:
            print 'Error - Unknown interpolation type: %s.'%self.method
        
        
        return self.Z
        
    def _nearestNeighbour(self):
        """ Nearest neighbour interpolation algorithm
            Sets any points outside of maxdist to NaN
        """ 

        self.Finterp = nn(self.XY,self.XYout,maxdist=self.maxdist)
      

    def _griddata(self):
        """Wrapper for griddata"""
        print 'Interpolating %d data points'%self.npt
        self.Z = griddata((self.XY[:,0],self.XY[:,1]), self.Zin, (self.grd.X, self.grd.Y), method='linear')
        
        
    
    def _invdistweight(self):
        """ Inverse distance weighted interpolation """
        
        self.Finterp=idw(self.XY,self.XYout,maxdist=self.maxdist,NNear=self.NNear,p=self.p)
        
    
    def _krig(self):    
        """ Kriging interpolation"""
       
        self.Finterp = kriging(self.XY,self.XYout,maxdist=self.maxdist,NNear=self.NNear)
                
    
    def clipPoints(self,LL):
        """ Clips points outside of the bounding box"""
        X = LL[:,0]
        Y = LL[:,1]
        ind = np.all([X>=self.bbox[0],X<=self.bbox[1],Y>=self.bbox[2],Y<=self.bbox[3]],axis=0)                    
        
        print 'Clipped %d points.'%(self.npt-sum(ind))
        self.Zin = self.Zin[ind]
        self.npt = len(self.Zin)
        return np.concatenate((np.reshape(X[ind],(self.npt,1)),np.reshape(Y[ind],(self.npt,1))),axis=1)
        
    def save(self,outfile='DEM.nc'):
        """ Saves the DEM to a netcdf file"""
        
        # Create the global attributes
        if self.isnorth:
            proj = "UTM %d (%s) in northern hemisphere."%(self.utmzone,self.CS)
        else:
            proj = "UTM %d (%s) in southern hemisphere."%(self.utmzone,self.CS)
        
        intparamstr = 'Interpolation Type: %s, Number of neighbours: %d, Maximum search distance: %3.1f m'%(self.method,self.NNear,self.maxdist)
        if self.method=='idw':
            intparamstr += ', IDW power: %2.1f'%self.p
        elif self.method=='kriging':
            intparamstr += ', Variogram model: %s, sill: %3.1f, nugget: %3.1f, range: %3.1f'%(self.varmodel,self.sill,self.nugget,self.vrange)
            
        globalatts = {'title':'DEM model',\
        'history':'Created on '+time.ctime(),\
        'Input dataset':self.infile,\
        'Projection':proj,\
        'Interpolation Parameters':intparamstr}
        
        
        nc = Dataset(outfile, 'w', format='NETCDF4')
        # Write the global attributes
        for gg in globalatts.keys():
            nc.setncattr(gg,globalatts[gg])
            
        # Create the dimensions
        dimnamex = 'nx'
        dimlength = self.grd.nx
        nc.createDimension(dimnamex,dimlength)
        dimnamey = 'ny'
        dimlength = self.grd.ny
        nc.createDimension(dimnamey,dimlength)
        
        # Create the lat lon variables
        tmpvarx=nc.createVariable('X','f8',(dimnamex,))
        tmpvary=nc.createVariable('Y','f8',(dimnamey,))
        tmpvarx[:] = self.grd.X[0,:]
        tmpvary[:] = self.grd.Y[:,0]
        # Create the attributes
        tmpvarx.setncattr('long_name','Easting')
        tmpvarx.setncattr('units','metres')
        tmpvary.setncattr('long_name','Northing')
        tmpvary.setncattr('units','metres')
        
        # Write the topo data
        tmpvarz=nc.createVariable('topo','f8',(dimnamey,dimnamex),zlib=True,least_significant_digit=1)
        tmpvarz[:] = self.Z
        tmpvarz.setncattr('long_name','Topographic elevation')
        tmpvarz.setncattr('units','metres')
        tmpvarz.setncattr('coordinates','X, Y')
        tmpvarz.setncattr('positive','up')
        tmpvarz.setncattr('datum',self.vdatum)
        
        nc.close()
        
        print 'DEM save to %s.'%outfile
        
        
        
    def scatter(self,**kwargs):
        fig= plt.figure(figsize=(9,8))
        #h.imshow(np.flipud(self.Z),extent=[bbox[0],bbox[1],bbox[3],bbox[2]])
        plt.scatter(np.ravel(self.grd.X),np.ravel(self.grd.Y),c=np.ravel(self.Z),s=10,**kwargs)
        plt.colorbar()
        return fig


class Inputs(object):
    """
        Class for handling input data from different file formats
        
    """
    
    # Projection information
    convert2utm=True
    CS='NAD83'
    utmzone=15
    isnorth=True
    vdatum = 'MSL'
     
    def __init__(self,infile,**kwargs):
    
        self.infile = infile        
        self.__dict__.update(kwargs)
        
        
        # Read in the array
        print 'Reading data from: %s...'%self.infile
        if self.infile[-3:]=='.gz':
            LL,self.Zin = read_xyz_gz(self.infile)
        elif self.infile[-3:] in ['txt','dat']:
            LL,self.Zin = read_xyz(self.infile)
            self.Zin = np.ravel(self.Zin)
        elif self.infile[-3:]=='shp':
            LL,self.Zin = readShpBathy(self.infile,FIELDNAME='contour')
        elif self.infile[-3:]=='.nc':
            self.loadnc()
            LL = self._returnXY(self.xgrd,self.ygrd)
            self.Zin = np.ravel(self.Zin)
        elif self.infile[-3:] in ['dem','asc']:
            xgrd,ygrd,self.Zin = readraster(self.infile)
            LL = self._returnXY(xgrd,ygrd)
            self.Zin = np.ravel(self.Zin)
        
        self.npt = len(LL)
        
        if self.convert2utm:                     
            # Convert the coordinates
            print 'Transforming the coordinates to UTM...'
            self.XY=ll2utm(LL,self.utmzone,self.CS,self.isnorth)
        else:
            self.XY=LL
   
        self._returnNonNan()
        
    def _returnXY(self, x, y):
        """
        Returns gridded points as a vector
        """
        X,Y = np.meshgrid(x,y)
        nx =  np.prod(np.shape(X))
        return np.hstack((np.reshape(np.ravel(X),(nx,1)),np.reshape(np.ravel(Y),(nx,1))))
        
    def _returnNonNan(self):
        
        ind = np.isnan(self.Zin)
        ind = ind==False
        
        self.Zin=self.Zin[ind]
        self.XY = self.XY[ind,:]
    
    def loadnc(self):
        """ Load the DEM data from a netcdf file"""        
        nc = Dataset(self.infile, 'r')
        try:
            self.xgrd = nc.variables['X'][:]
            self.ygrd = nc.variables['Y'][:]
            self.Zin = nc.variables['topo'][:]
        except:
            self.xgrd = nc.variables['x'][:]
            self.ygrd = nc.variables['x'][:]
            self.Zin = nc.variables['z'][:]
                
        nc.close()

class idw(object):
    """Inverse distance weighted interpolation function"""
    
    maxdist=300
    NNear=3
    p=1
    
    def __init__(self,XYin,XYout,**kwargs):
        self.__dict__.update(kwargs)

        
        # Compute the spatial tree
        kd = spatial.cKDTree(XYin)
        
        # Perform query on all of the points in the grid
        dist,self.ind=kd.query(XYout,distance_upper_bound=self.maxdist,k=self.NNear)
        
        # Calculate the weights
        self.W = 1/dist**self.p
        Wsum = np.sum(self.W,axis=1)
        
        for ii in range(self.NNear):
            self.W[:,ii] = self.W[:,ii]/Wsum
            
        # create the mask
        mask = (dist==np.inf)
        self.ind[mask]=1
    
    def __call__(self,Zin):
        # Fill the array and resize it
        Zin = np.squeeze(Zin[self.ind])
        
        # Compute the weighted sums and mask the blank points
        return np.sum(Zin*self.W,axis=1)
        #Z[mask]=np.nan
        
class nn(object):
    """ 
    Nearest neighbour interpolation algorithm
    Sets any points outside of maxdist to NaN
    """
    maxdist = 1000.0
    
    def __init__(self,XYin,XYout,**kwargs):
        self.__dict__.update(kwargs)
        
        # Compute the spatial tree
        kd = spatial.cKDTree(XYin)
        # Perform query on all of the points in the grid
        dist,self.ind=kd.query(XYout,distance_upper_bound=self.maxdist)
        # create the mask
        self.mask = (dist==np.inf)
        self.ind[self.mask]=1
    
    def __call__(self,Zin):
        # Fill the array and resize it
        Z = Zin[self.ind]
        Z[self.mask]=np.nan
        
        return Z            
## Other functions that don't need to be in a class ##

    
def read_xyz_gz(fname):
    # Read the raw data into an array
    f = gzip.open(fname,'r')
    
    npts = line_count(f)-1
    XY = np.zeros((npts,2))
    Z = np.zeros((npts,1))
    ii=-1
    for line in f:
        ii+=1
        if ii > 0:
            xyz = line.split(', ')
            XY[ii-1,0] = float(xyz[0]) 
            XY[ii-1,1] = float(xyz[1])
            Z[ii-1,0] = float(xyz[2])
            
    f.close()
      
    return XY,Z
    
def read_xyz(fname):
    # Read the raw data into an array
    f =open(fname,'r')
    
    npts = line_count(f)-1
    XY = np.zeros((npts,2))
    Z = np.zeros((npts,1))
    ii=-1
    for line in f:
        ii+=1
        if ii > 0:
            xyz = line.split()
            XY[ii-1,0] = float(xyz[0]) 
            XY[ii-1,1] = float(xyz[1])
            Z[ii-1,0] = float(xyz[2])
#            try: # comma delimeted
#                xyz = line.split(', ')
#                XY[ii-1,0] = float(xyz[0]) 
#                XY[ii-1,1] = float(xyz[1])
#                Z[ii-1,0] = float(xyz[2])
#            except: # space delimitede
#                xyz = line.split(' ')
#                XY[ii-1,0] = float(xyz[0]) 
#                XY[ii-1,1] = float(xyz[1])
#                print xyz[2]
#                Z[ii-1,0] = float(xyz[2])
                
            
    f.close()
      
    return XY,Z
    

def line_count(f):
    for i, l in enumerate(f):
        pass
#    try:
#        f.rewind()
#    except ValueError:
#        f.seek(0,0) 
    f.seek(0,0)
    return i + 1

def tile_vector(count,chunks):
    rem = np.remainder(count,chunks)
    
    cnt2 = count-rem
    dx = cnt2/chunks
    
    if count != cnt2:
        pt1 = range(0,cnt2,dx)
        pt2 = range(dx,cnt2,dx) + [count]
    else:
        pt1 = range(0,count-dx,dx)
        pt2 = range(dx,count,dx)  
    return pt1,pt2
    

    
################
# Testing sections

## Initialise the Input class
#infile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAA_25m_UTM_DEM.nc'
#indata = Inputs(infile,convert2utm=False)
#
## Initialise the Interpolation class
#print 'Building interpolant class...'
#F = interpXYZ(indata.XY,indata.Zin,method='idw',NNear=3)
#
## Initialise the interpolation points
#
#print 'Loading suntans grid points...'
#filename = 'C:/Projects/GOMGalveston/MODELLING/GRIDS/GalvestonCoarse'
#grd = sunpy.Grid(filename)
#xy = np.column_stack((grd.xv,grd.yv))
#
## Interpolate the data
#print 'Interpolating data...'
#grd.dv = F(xy)
#
#grd.clim=[-20.0,0.0]
#plt.figure()
#grd.plot(cmap=plt.cm.gist_earth)
#plt.show()
#
#print 'Smoothing the data...'
#Fsmooth = interpXYZ(xy,grd.dv,method='kriging',NNear=5)
#grd.dv = Fsmooth(xy)
#
#plt.figure()
#grd.plot(cmap=plt.cm.gist_earth)
#plt.show()
