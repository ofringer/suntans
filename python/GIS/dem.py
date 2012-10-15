# -*- coding: utf-8 -*-
"""

General class for parsing different DEMs
Created on Mon Sep 10 14:43:26 2012

@author: mrayson
"""

import numpy as np
from netCDF4 import Dataset
from scipy import spatial
import matplotlib.pyplot as plt
from interp_xyz import tile_vector
import time
import shutil
import gdal
from gdalconst import * 

import pdb

class DEM(object):
    """
        General DEM class
    """    
    
    infile = None
    W = 1.0 # Weight
    maxdist = 250.0
    
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)
        if self.infile[-3:]=='.nc':
            xgrd,ygrd,self.Z = self.loadnc()
        elif self.infile[-3:] in ['dem','asc']:
            xgrd,ygrd,self.Z = self.readraster()
        
        # Generate the grid
        self.x0 = xgrd.min()
        self.y0 = ygrd.min()
        self.x1 = xgrd.max()
        self.y1 = ygrd.max()
        self.dx= xgrd[1]-xgrd[0]
        self.dy= ygrd[1]-ygrd[0]
#        

        self.nx = len(xgrd)
        self.ny = len(ygrd)
        self.npts = self.nx*self.ny
#        
        self.X,self.Y = np.meshgrid(xgrd,ygrd)
        
    def loadnc(self):
        """ Load the DEM data from a netcdf file"""        
        nc = Dataset(self.infile, 'r')
        try:
            X = nc.variables['X'][:]
            Y = nc.variables['Y'][:]
            Z = nc.variables['topo'][:]
        except:
            X = nc.variables['x'][:]
            Y = nc.variables['x'][:]
            Z = nc.variables['z'][:]
                
        nc.close()
        return X,Y,Z
        
    def readraster(self):
        """ Loads the data from a DEM raster file"""
        # register all of the drivers
        gdal.AllRegister()
        # open the image
        ds = gdal.Open(self.infile, GA_ReadOnly)
        
        # Read the x and y coordinates
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        bands = ds.RasterCount
        
        geotransform = ds.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        
        x = originX + np.linspace(0,cols-1,cols)*pixelWidth
        y = originY + np.linspace(0,rows-1,rows)*pixelHeight
        
        # Read the actual data
        data = ds.ReadAsArray(0,0,cols,rows)
        
        # Remove missing points
        data[data==-32767]=np.nan
        
        return x, y, data
    
        
    def ravel(self):
        """ Returns the grid coordinates as a vector"""
        return np.concatenate( (np.reshape(np.ravel(self.X),(self.npts,1)),\
            np.reshape(np.ravel(self.Y),(self.npts,1))),axis=1)
            
    def nanxy(self):
        """
            Returns the x,y locations of the nan points
        """
        ind = np.isnan(self.Z)
        nc = np.sum(ind)
        xy = np.zeros((nc,2)) 
        n = -1
        for jj in range(0,self.ny):  
            for ii in range(0,self.nx):  
                if ind[jj,ii]:
                    n+=1
                    xy[n,0]=self.X[jj,ii]
                    xy[n,1]=self.Y[jj,ii]
        
        return xy
        
#        ind = np.isnan(np.ravel(self.Z))
#        nc = np.sum(ind)
#        
#        x=np.ravel(self.X)
#        y=np.ravel(self.Y)
#        
#        return np.concatenate((np.reshape(x[ind],(nc,1)),np.reshape(y[ind],(nc,1))),axis=1)

        
    def nonnanxy(self):
        """
            Returns the x,y locations of the non-nan points
        """
        ind = np.isnan(self.Z)
        ind = ind==False
        nc = np.sum(ind)
        xy = np.zeros((nc,2)) 
        n = -1
        for jj in range(0,self.ny):  
            for ii in range(0,self.nx):  
                if ind[jj,ii]:
                    n+=1
                    xy[n,0]=self.X[jj,ii]
                    xy[n,1]=self.Y[jj,ii]
        
        return xy
        
#        ind = np.isnan(np.ravel(self.Z))
#        ind = ind==False
#        nc = np.sum(ind)
#        print nc
#        x=np.ravel(self.X)
#        y=np.ravel(self.Y)
#        
#        return np.concatenate((np.reshape(x[ind],(nc,1)),np.reshape(y[ind],(nc,1))),axis=1)

    def returnij(self,x,y):
        """def contourf(self,Z,vv=range(-10,0),**kwargs):
        fig= plt.figure(figsize=(9,8))
        plt.contourf(self.X,self.Y,Z,vv,**kwargs)
        plt.colorbar()
        plt.axis('equal')
        return fig
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
        
    def calcWeight(self):
        
        """ Calculate the weight at each point """
        MAXPOINTS=20e6
        weight = np.zeros((self.ny,self.nx))
        
        # Calculate the distance from each point to a nan point
        xy = self.nonnanxy()
        xynan = self.nanxy()
        
        # Compute the spatial tree
        kd = spatial.cKDTree(xynan)
        
        nxy = len(xy)
        
        if nxy <= MAXPOINTS:
            # Perform query on all of the points in the grid
            dist,ind=kd.query(xy)
            
            # Compute the actual weight
            w = dist/self.maxdist
            w[dist>self.maxdist]=1.0
            w=self.W*w
            
            # Map onto the grid
            J,I=self.returnij(xy[:,0],xy[:,1])
            weight[J,I]=w
        else:
            print 'Dataset too large - calculating weights for chunks...'
            nchunks = np.ceil(len(xy)/MAXPOINTS)
            pt1,pt2=tile_vector(len(xy),int(nchunks))
            for p1,p2 in zip(pt1,pt2):
                print 'Calculating points %d to %d of %d...'%(p1,p2,nxy)
                dist,ind=kd.query(xy[p1:p2,:])
                # Compute the actual weight
                w = dist/self.maxdist
                w[dist>self.maxdist]=1.0
                w=self.W*w
                
                # Map onto the grid
                J,I=self.returnij(xy[p1:p2,0],xy[p1:p2,1])
                weight[J,I]=w
        
        return weight   
        
    def contourf(self,Z,vv=range(-10,0),**kwargs):
        fig= plt.figure(figsize=(9,8))
        plt.contourf(self.X,self.Y,Z,vv,**kwargs)
        plt.colorbar()
        plt.hold(True)
        plt.contour(self.X,self.Y,Z,[0.0,0.0],colors='k',linewidths=0.02)
        plt.axis('equal')
        return fig
        
    def plot(self,Z,**kwargs):
        h= plt.figure(figsize=(9,8))
        #h.imshow(np.flipud(self.Z),extent=[bbox[0],bbox[1],bbox[3],bbox[2]])
        plt.imshow(Z,extent=[self.x0,self.x1,self.y0,self.y1],**kwargs)
        plt.colorbar()
        return h
        
    def savenc(self,outfile='DEM.nc'):
        """ Saves the DEM to a netcdf file"""
        
        # Create the global attributes
        
        globalatts = {'title':'DEM model',\
        'history':'Created on '+time.ctime(),\
        'Input dataset':self.infile}
        
        
        nc = Dataset(outfile, 'w', format='NETCDF4')
        # Write the global attributes
        for gg in globalatts.keys():
            nc.setncattr(gg,globalatts[gg])
            
        # Create the dimensions
        dimnamex = 'nx'
        dimlength = self.nx
        nc.createDimension(dimnamex,dimlength)
        dimnamey = 'ny'
        dimlength = self.ny
        nc.createDimension(dimnamey,dimlength)
        
        # Create the lat lon variables
        tmpvarx=nc.createVariable('X','f8',(dimnamex,))
        tmpvary=nc.createVariable('Y','f8',(dimnamey,))
        tmpvarx[:] = self.X[0,:]
        tmpvary[:] = self.Y[:,0]
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
        
        nc.close()
        
        print 'DEM save to %s.'%outfile

    
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
    
def blendDEMs(ncfile,outfile,W,maxdist):
    ### Combine multiple files ###   
    
    #Calculate the weights for each file
    nfiles = len(ncfile)
    ii=-1
    for nc in ncfile:
        ii+=1
        d = DEM(infile=nc,W=W[ii],maxdist=maxdist[ii])
        print 'Calculating weights for %s...'%nc
        print 'Weight = %6.3f, maxdist = %f'%(W[ii],maxdist[ii])
        w=d.calcWeight()
        ny = d.ny
        nx = d.nx
        
#        if ii == 1:
#            f=d.contourf(w,vv=np.linspace(0,W[ii],10))
#            f.savefig('%s_Weights.pdf'%outfile[:-2])
#            del f
        del d
        if ii == 0:
            Wall = np.zeros((ny,nx,nfiles))
        Wall[:,:,ii]=w
        del w
        
    # Normalise the weights
    print 'Normalising the weights...'
    Wsum = np.sum(Wall,axis=2)
    for ii in range(0,nfiles):
        Wall[:,:,ii] = np.squeeze(Wall[:,:,ii]) / Wsum
        
        
        
        
    # Re-load in the depths from each file and sum
    print 'Writing to an output file...'
    Zout = np.zeros((ny,nx))
    filestr = ''
    ii=-1
    for infile in ncfile:
        ii+=1
        nc = Dataset(infile, 'r')
        Zin = nc.variables['topo'][:]
        nc.close()
        
        Zin[np.isnan(Zin)]=0.0
        Zout +=  np.squeeze(Wall[:,:,ii]) * Zin 
        filestr +='%s, '%infile
    
    # Copy the data to a new netcdf file
    shutil.copyfile(ncfile[-1],outfile)
    nc = Dataset(outfile, 'r+')
    nc.variables['topo'][:]=Zout
    
    globalatts = {'title':'DEM model',\
        'history':'Created on '+time.ctime(),\
        'Input datasets':filestr}
    # Write the global attributes
    for gg in globalatts.keys():
        nc.setncattr(gg,globalatts[gg])
        
    nc.close()
    
    print 'Completed write to %s.'%outfile


#ncfile = [\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/USACELIDAR_dx25_blockavg.nc',\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAADEM_dx25_IDW_dist100_NNEar3.nc',\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAASoundingsDEM_dx25_KRIG_dist500_Nnear3_range200.nc',\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/TNRIS_dx25_GridData.nc'\
#]
#W = [50.0,1.0,10.0,0.1]
#maxdist = [100.0,100.,1500.,1000.0]
#outfile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/Blended/NOAA_Blended_All.nc'
#
#blendDEMs(ncfile,outfile,W,maxdist)

#d = DEM(infile=outfile)
#f=d.contourf(d.Z,vv=range(-15,3),vmin=-15,vmax=4,cmap=plt.cm.gist_earth)
#f.savefig(outfile[:-2]+'pdf')
#plt.show()

## Load in other formats
#infile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/NGDCCoastalRelief/galveston_tx.asc'
#print 'Loading %s...'%infile
#d = DEM(infile=infile)
#
#print 'Saving to an image...'
##f=d.contourf(d.Z,vv=range(-15,3),vmin=-15,vmax=4,cmap=plt.cm.gist_earth)
#f=d.plot(d.Z,vmin=-20,vmax=5,cmap=plt.cm.gist_earth)
#f.savefig(infile[:-3]+'png',dpi=1200)
#d.savenc(outfile='C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAA_10m_DEM.nc')

infile='C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAA_10m_DEM.nc'
print 'Loading %s...'%infile
d = DEM(infile=infile)

print 'Saving to an image...'
f=d.contourf(d.Z+0.14,vv=range(-20,3),vmin=-20,vmax=4,cmap=plt.cm.gist_earth)
f.savefig(infile[:-3]+'_MSL.'+'pdf')