# -*- coding: utf-8 -*-
"""
    Interpolate sounding data onto a regular grid
"""

import gzip
from scipy import spatial
import numpy as np
from maptools import ll2utm, readShpBathy, readDEM
from kriging import kriging
from netCDF4 import Dataset

from scipy.sparse import coo_matrix
from scipy.interpolate import griddata
import time
import matplotlib.pyplot as plt

# testing stuff
import pdb

class demBuilder(object):
    "Class for building a DEM on a regular cartesian grid from different inputs"""
    
    ## Properties ##
    infile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/NOSHyrdrographicSurveys/xyz/NOS_Galveston_Survey_Depths.csv.gz'
    
    # Projection information
    convert2utm=True
    CS='NAD83'
    utmzone=15
    isnorth=True
    vdatum = 'MSL'
    
    # DEM Grid params
    dx=25
    bbox = [-95.45,-94.44,28.8,29.8]
    
    # Interpolation options
    interptype = 'nn' # One of 'nn', 'blockavg', 'idw'
    maxdist=200
    NNear = 3 # Number of points to include in interpolation (only applicable to idw and kriging)
    p = 1.0 # power for inverse distance weighting
    
    # kriging options
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0
    
    def __init__(self,**kwargs):
        
        self.__dict__.update(kwargs)
        
        # Check if the input file is not a list
        T = type(self.infile)
        
        if T!=list:
            self.multifile=False
            # Read in the array
            print 'Reading data from: %s...'%self.infile
            if self.infile[-3:]=='.gz':
                LL,self.Zin = read_xyz_gz(self.infile)
            if self.infile[-3:]=='txt':
                LL,self.Zin = read_xyz(self.infile)
            elif self.infile[-3:]=='shp':
                LL,self.Zin = readShpBathy(self.infile)
            elif self.infile[-3:]=='dem':
                LL,self.Zin = readDEM(self.infile,True)
            if self.infile[-3:]=='.nc':
                self.loadnc(fv=2)
                LL = self._returnXY()
                self.Zin = np.ravel(self.Zin)
            
            self.npt = len(self.Zin)
            if self.convert2utm:
                if self.bbox==None:
                    # Work out the domain limits from the input file
                    pdb.set_trace()
                    
                    self.bbox = [LL[:,0].min(),LL[:,0].max(),LL[:,1].min(),LL[:,1].max()]
                else:
                    # Clip the points outside of the domain
                    print 'Clipping points outside of the bounding box...'
                    LL=self.clipPoints(LL)   
                # Convert the coordinates
                print 'Transforming the coordinates to UTM...'
                self.XY=ll2utm(LL,self.utmzone,self.CS,self.isnorth)
            else:
                self.XY=LL
        else: # Multiple files
            self.multifile=True
        
        # Create the grid object
        self.grd = Grid(self.bbox,self.dx,self.dx,utmzone=self.utmzone,CS=self.CS,isnorth=self.isnorth)
        
    def build(self):
        
        tic=time.clock()
        if self.multifile==False:
            if self.interptype=='nn':
                print 'Building DEM with Nearest Neighbour interpolation...'
                self.nearestNeighbour()
                
            elif self.interptype=='blockavg':
                print 'Building DEM with Block Averaging...'
                self.blockAvg()
            
            elif self.interptype=='idw':
                print 'Building DEM with Inverse Distance Weighted Interpolation...'
                self.invdistweight()
                
            elif self.interptype=='kriging':
                print 'Building DEM with Kriging Interpolation...'
                self.krig()
            elif self.interptype=='griddata':
                print 'Building DEM using griddata...'
                self.griddata()
            else:
                print 'Error - Unknown interpolation type: %s.'%self.interptype
        else: # Multiple file interpolation
            print 'Multiple input files detected - setting "interptype" to "blockavg".'
            self.interptype = 'blockavg'
            self.Z = np.zeros((self.grd.ny,self.grd.nx))
            self.N = np.zeros((self.grd.ny,self.grd.nx))
            ctr=0
            for f in self.infile:
                ctr+=1
                # Read in the array
                print 'Reading data file (%d of %d): %s...'%(ctr,len(self.infile),f)
                if f[-3:]=='.gz':
                    LL,self.Zin = read_xyz_gz(f)
                if f[-3:]=='txt':
                    LL,self.Zin = read_xyz(f)
                elif f[-3:]=='shp':
                    LL,self.Zin = readShpBathy(f)
                elif f[-3:]=='dem':
                    LL,self.Zin = readDEM(f,True)
                
                self.npt = len(self.Zin)
                
                if self.convert2utm:
            
                    # Clip the points outside of the domain
                    #print 'Clipping points outside of the bounding box...'
                    #LL=self.clipPoints(LL)   
                    
                    # Convert the coordinates
                    print 'Transforming the coordinates to UTM...'
                    self.XY=ll2utm(LL,self.utmzone,self.CS,self.isnorth)
                else:
                    self.XY=LL

                del LL
                # Interpolate
                print 'Building DEM with Block Averaging...'
                self.blockAvgMulti()
                
                # Memory cleanup
                del self.XY
                del self.Zin
                
            
            # Compute the block average for all of the files
            self.Z = np.divide(self.Z,self.N)


        toc=time.clock()
        print 'Elapsed time %10.3f seconds.'%(toc-tic)
    
    def _returnXY(self):
        """
        Returns gridded points as a vector
        """
        X,Y = np.meshgrid(self.xgrd,self.ygrd)

        return np.column_stack((np.ravel(X),np.ravel(Y)))
        
    def clipPoints(self,LL):
        """ Clips points outside of the bounding box"""
        X = LL[:,0]
        Y = LL[:,1]
        ind = np.all([X>=self.bbox[0],X<=self.bbox[1],Y>=self.bbox[2],Y<=self.bbox[3]],axis=0)                    
        
        print 'Clipped %d points.'%(self.npt-sum(ind))
        self.Zin = self.Zin[ind]
        self.npt = len(self.Zin)
        return np.concatenate((np.reshape(X[ind],(self.npt,1)),np.reshape(Y[ind],(self.npt,1))),axis=1)
        
    def nearestNeighbour(self):
        """ Nearest neighbour interpolation algorithm
            Sets any points outside of maxdist to NaN
        """ 
        
        MAXSIZE = 10e6
        nchunks = np.ceil(self.grd.npts*self.NNear/MAXSIZE)
        
        if nchunks == 1:
            Z = nn(self.XY,self.Zin,self.grd.ravel(),maxdist=self.maxdist)
        else:
           pt1,pt2=tile_vector(int(self.grd.npts),int(nchunks))
           Z = np.zeros((self.grd.npts,))
           XYout = self.grd.ravel()
           for p1,p2 in zip(pt1,pt2):
               print 'Interpolating tile %d to %d of %d...'%(p1,p2,self.grd.npts)
               Z[p1:p2] = nn(self.XY,self.Zin,XYout[p1:p2,:],maxdist=self.maxdist)
        
        self.Z = np.reshape(Z,(self.grd.ny,self.grd.nx))
    
    def griddata(self):
        """Wrapper for griddata"""
        print 'Interpolating %d data points'%self.npt
        self.Z = griddata((self.XY[:,0],self.XY[:,1]), self.Zin, (self.grd.X, self.grd.Y), method='linear')
        
        
    def blockAvg(self):
        
        """Block averaging interpolation"""
        
        # Get the grid indices
        J,I = self.grd.returnij(self.XY[:,0],self.XY[:,1])
        
        # Average onto the grid
        Z = np.zeros((self.grd.ny,self.grd.nx))
        N = np.zeros((self.grd.ny,self.grd.nx))
        ctr=-1
        for jj,ii in zip(J,I):
            ctr+=1
            if jj != -1 and ii != -1:
                Z[jj,ii] += self.Zin[ctr]
                N[jj,ii] += 1.0
            
        self.Z = np.divide(Z,N)
        self.N = N
    
    def blockAvgMulti(self):
        
        """Block averaging interpolation"""
        print 'Interpolating %d data points'%self.npt
        # Get the grid indices
        J,I = self.grd.returnij(self.XY[:,0],self.XY[:,1])
        
        # Zero out of bound points
        sumpts = np.ones((self.npt,))
        ind = I==-1
        self.Zin[ind]=0.0
        I[ind]=0
        sumpts[ind]=0
        ind = J==-1
        self.Zin[ind]=0.0
        J[ind]=0
        sumpts[ind]=0
        
        #  Use the sparse matrix library for accumulation
        self.Z += coo_matrix((np.ravel(self.Zin),(J,I)),\
            shape=(self.grd.ny,self.grd.nx)).todense()   
        self.N += coo_matrix((sumpts,(J,I)),\
            shape=(self.grd.ny,self.grd.nx)).todense()   

##        # Average onto the grid
##        ctr=-1
##        for jj,ii in zip(J,I):
##            ctr+=1
##            if jj != -1 and ii != -1:
##                self.Z[jj,ii] += self.Zin[ctr]
##                self.N[jj,ii] += 1.0
            
    
    def invdistweight(self):
        """ Inverse distance weighted interpolation """
        
        # Break it down into smaller chunks
        MAXSIZE = 10e6
        nchunks = np.ceil(self.grd.npts*self.NNear/MAXSIZE)
        
        if nchunks == 1:
            Z=idw(self.XY,self.Zin,self.grd.ravel(),maxdist=self.maxdist,NNear=self.NNear,p=self.p)
        else:
           pt1,pt2=tile_vector(int(self.grd.npts),int(nchunks))
           Z = np.zeros((self.grd.npts,))
           XYout = self.grd.ravel()
           for p1,p2 in zip(pt1,pt2):
               print 'Interpolating tile %d to %d of %d...'%(p1,p2,self.grd.npts)
               Z[p1:p2]=idw(self.XY,self.Zin,XYout[p1:p2,:],maxdist=self.maxdist,NNear=self.NNear,p=self.p)
        
        self.Z = np.reshape(Z,(self.grd.ny,self.grd.nx))
    
    def krig(self):    
        """ Kriging interpolation"""
         # Break it down into smaller chunks
        MAXSIZE = 15e6
        nchunks = np.ceil(self.grd.npts*self.NNear/MAXSIZE)
        
        if nchunks == 1:
            self.Finterp = kriging(self.XY,self.grd.ravel(),maxdist=self.maxdist,NNear=self.NNear)
            Z = self.Finterp(self.Zin)
        else:
            pt1,pt2=tile_vector(int(self.grd.npts),int(nchunks))
            Z = np.zeros((self.grd.npts,))
            XYout = self.grd.ravel()
            for p1,p2 in zip(pt1,pt2):
                print 'Interpolating tile %d to %d of %d...'%(p1,p2,self.grd.npts)
                self.Finterp = kriging(self.XY,XYout[p1:p2,:],maxdist=self.maxdist,NNear=self.NNear)
                Z[p1:p2] = self.Finterp(self.Zin)
                
        self.Z = np.reshape(Z,(self.grd.ny,self.grd.nx))
    
    def loadnc(self,fv=1):
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
        
        self.xgrd=self.xgrd[::fv]
        self.ygrd=self.ygrd[::fv]
        self.Zin=self.Zin[::fv,::fv]
           
    def save(self,outfile='DEM.nc'):
        """ Saves the DEM to a netcdf file"""
        
        # Create the global attributes
        if self.isnorth:
            proj = "UTM %d (%s) in northern hemisphere."%(self.utmzone,self.CS)
        else:
            proj = "UTM %d (%s) in southern hemisphere."%(self.utmzone,self.CS)
        
        intparamstr = 'Interpolation Type: %s, Number of neighbours: %d, Maximum search distance: %3.1f m'%(self.interptype,self.NNear,self.maxdist)
        if self.interptype=='idw':
            intparamstr += ', IDW power: %2.1f'%self.p
        elif self.interptype=='kriging':
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
        
        
    def plot(self,**kwargs):
        h= plt.figure(figsize=(9,8))
        #h.imshow(np.flipud(self.Z),extent=[bbox[0],bbox[1],bbox[3],bbox[2]])
        plt.imshow(np.flipud(self.Z),extent=[self.grd.x0,self.grd.x1,self.grd.y0,self.grd.y1],**kwargs)
        plt.colorbar()
        return h
        
    def scatter(self,**kwargs):
        fig= plt.figure(figsize=(9,8))
        #h.imshow(np.flipud(self.Z),extent=[bbox[0],bbox[1],bbox[3],bbox[2]])
        plt.scatter(np.ravel(self.grd.X),np.ravel(self.grd.Y),c=np.ravel(self.Z),s=10,**kwargs)
        plt.colorbar()
        return fig

    def contourf(self,vv=range(-10,0),**kwargs):
        fig= plt.figure(figsize=(9,8))
        #h.imshow(np.flipud(self.Z),extent=[bbox[0],bbox[1],bbox[3],bbox[2]])
        plt.contourf(self.grd.X,self.grd.Y,self.Z,vv,**kwargs)
        plt.colorbar()
        plt.axis('equal')
        return fig
        
class Grid(object):
    """ Cartesian grid object"""
    
    CS='NAD83'
    utmzone=15
    isnorth=True
    
    def __init__(self,bbox,dx,dy,**kwargs):
        self.__dict__.update(kwargs)

        # Generate the grid
        xy0 = ll2utm([bbox[0],bbox[2]],self.utmzone,self.CS,self.isnorth)
        xy1 = ll2utm([bbox[1],bbox[3]],self.utmzone,self.CS,self.isnorth)    
        self.x0 = xy0[0,0]
        self.y0 = xy0[0,1]
        self.x1 = xy1[0,0]
        self.y1 = xy1[0,1]
        self.dx=dx
        self.dy=dy
        
        xgrd = np.arange(self.x0,self.x1,dx)
        ygrd = np.arange(self.y0,self.y1,dy)
        self.nx = len(xgrd)
        self.ny = len(ygrd)
        self.npts = self.nx*self.ny
        
        self.X,self.Y = np.meshgrid(xgrd,ygrd)
        
    def ravel(self):
        """ Returns the grid coordinates as a vector"""
        return np.concatenate( (np.reshape(np.ravel(self.X),(self.npts,1)),\
            np.reshape(np.ravel(self.Y),(self.npts,1))),axis=1)
            
    def returnij(self,x,y):
        """
        Returns the grid cell indices that points x,y reside inside of.
        
        """
        I = np.ceil( (x-self.x0)/self.dx)
        J =np.ceil( (y-self.y0)/self.dy)
        
        J = np.array(J,dtype=int)
        I = np.array(I,dtype=int)
        
        # blank out bad cells
        J[J<0]=-1
        J[J>self.ny-1]=-1
        I[I<0]=-1
        I[I>self.nx-1]=-1
        
        return J,I
                
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
            xyz = line.split(', ')
            XY[ii-1,0] = float(xyz[0]) 
            XY[ii-1,1] = float(xyz[1])
            Z[ii-1,0] = float(xyz[2])
            
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
    
def idw(XYin,Zin,XYout,maxdist=300,NNear=3,p=1):
    """Inverse distance weighted interpolation function"""
    # Compute the spatial tree
    kd = spatial.cKDTree(XYin)
    
    # Perform query on all of the points in the grid
    dist,ind=kd.query(XYout,distance_upper_bound=maxdist,k=NNear)
    
    # Calculate the weights
    W = 1/dist**p
    Wsum = np.sum(W,axis=1)
    
    for ii in range(NNear):
        W[:,ii] = W[:,ii]/Wsum
        
    # create the mask
    mask = (dist==np.inf)
    ind[mask]=1
    
    # Fill the array and resize it
    Zraw = np.squeeze(Zin[ind])
    
    # Compute the weighted sums and mask the blank points
    return np.sum(Zraw*W,axis=1)
    #Z[mask]=np.nan
def nn(XYin,Zin,XYout,maxdist=300):
    """ Nearest neighbour interpolation algorithm
        Sets any points outside of maxdist to NaN
    """
    # Compute the spatial tree
    kd = spatial.cKDTree(XYin)
    # Perform query on all of the points in the grid
    dist,ind=kd.query(XYout,distance_upper_bound=maxdist)
    # create the mask
    mask = (dist==np.inf)
    ind[mask]=1
    
    # Fill the array and resize it
    Z = Zin[ind]
    Z[mask]=np.nan
    
    return(Z)
################
# Testing sections
#dem = demBuilder(dx=100,interptype='kriging',maxdist=500,NNear=6,vrange=200)
#dem = demBuilder(dx=50,interptype='idw',maxdist=500,NNear=3)
#dem.build()
#
#ncfile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/testDEM.nc'
#dem.save(ncfile)
#
#f=dem.contourf(vv=range(-15,1),vmin=-15,vmax=0,cmap=plt.cm.gist_earth)
#f.savefig(ncfile[:-2]+'pdf')
#dem.plot(vmin=-15,vmax=0,cmap=plt.cm.gist_earth)
#dem.scatter(vmin=-15,vmax=0,cmap=plt.cm.gist_earth)

#infile = 'E:/Projects/GOMGalveston/DATA/Bathymetry/LIDAR/USACE_2009/29094_48_43_raw.txt'
#LL,Z = read_xyz(infile)
