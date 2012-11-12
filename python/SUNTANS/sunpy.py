# -*- coding: utf-8 -*-
"""

Tools for handling SUNTANS output data

Created on Mon Sep 24 16:55:45 2012

@author: mrayson
"""

from netCDF4 import num2date
from netCDF4 import MFDataset, Dataset
import numpy as np
from datetime import datetime
import os, time, getopt, sys

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.animation as animation


import pdb

class Spatial(object):
    
    """ Class for reading SUNTANS spatial netcdf output files """
    
    def __init__(self,ncfile, **kwargs):
        
        self.ncfile = ncfile
        
        # Open the netcdf file
        self.__openNc()
        
        # Load the grid
        self.grid = Grid(self.ncfile)  
        self.xy = self.grid.cellxy()
        
        # Load the time variable
        self.loadTime()
        
        # Set some default parameters
        self.tstep=0
        self.klayer=0 # -1 get seabed
        self.j=np.arange(0,self.grid.Nc) # Grid cell number for time series
        self.variable='eta'
        
        # Plotting parmaters
        self.clim=None
        
        self.__dict__.update(kwargs)
        
        # Update tstep 
        self.__updateTstep()
        

     
    def loadData(self):
        """ 
            Load the specified suntans variable data as a vector
            
        """
        nc = self.nc
        
        self.long_name = nc.variables[self.variable].long_name
        self.units= nc.variables[self.variable].units
        #        ndims = len(nc.variables[self.variable].dimensions)
        ndim = nc.variables[self.variable].ndim
        if ndim==1:
            self.data=nc.variables[self.variable][self.j]
        elif ndim==2:
            self.data=nc.variables[self.variable][self.tstep,self.j]
        else:
            if self.klayer==-1: # grab the seabed values
                klayer = np.arange(0,self.grid.Nkmax)

                if type(self.tstep)==int:
                    data=nc.variables[self.variable][self.tstep,klayer,self.j]
                    self.data = data[self.grid.Nk[self.j],self.j]
                else: # need to extract timestep by timestep for animations to save memory
                    self.data=np.zeros((len(self.tstep),len(self.j)))
                    i=-1
                    for t in self.tstep:
                        i+=1
                        data=nc.variables[self.variable][t,klayer,self.j]
                        self.data[i,:] = data[self.grid.Nk[self.j],self.j]
            else:
                self.data=nc.variables[self.variable][self.tstep,self.klayer,self.j]
        
        
        
        
    def loadTime(self):
         """
            Load the netcdf time as a vector datetime objects
         """
         #nc = Dataset(self.ncfile, 'r', format='NETCDF4') 
         nc = self.nc
         t = nc.variables['time']
         self.time = num2date(t[:],t.units)


    
    def plot(self,**kwargs):
        """
          Plot the unstructured grid data
        """
        # Load the data if it's needed
        if not self.__dict__.has_key('data'):
            self.loadData()
            
        if self.__dict__.has_key('patches'):
             self.patches.set_array(self.data)
             self.ax.add_collection(self.patches)
        else:
            #tic = time.clock()
            self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,self.data,xlim=self.grid.xlims,ylim=self.grid.ylims,\
                clim=self.clim,**kwargs)
            plt.title(self.__genTitle())
            #print 'Elapsed time: %f seconds'%(time.clock()-tic)
            
    def plotvtk(self,**kwargs):
        """
          Plot the unstructured grid data using vtk libraries
        """
        # Mayavi libraries

        
        # Load the data if it's needed
        if not self.__dict__.has_key('data'):
            self.loadData()
        
        points = np.column_stack((self.grid.xp,self.grid.yp,0.0*self.grid.xp))
                
        self.fig,h,ug,title = unsurfm(points,self.grid.cells,self.data,clim=self.clim,title=self.__genTitle())
        
            
    def plotTS(self,j=0,**kwargs):
        """
         Plots a time series of the data at a given grid cell given by:
             self.j, self.klayer
        """

        # Load the time-series
        self.tstep = np.arange(0,len(self.time))
        self.j=j
        self.loadData()
        
        plt.ioff()
        self.fig = plt.gcf()
        ax = self.fig.gca()
        h = plt.plot(self.time,self.data,'b') 
        c
        ax.set_title(titlestr)
        
       
        
    def savefig(self,outfile,dpi=150):
        mod=self.fig.__class__.__module__
        name=self.fig.__class__.__name__
        
        if mod+'.'+name == 'matplotlib.figure.Figure':
            self.fig.savefig(outfile,dpi=dpi)
        else:
            self.fig.scene.save(outfile)
            
        print 'SUNTANS image saved to file:%s'%outfile
    
    def animate(self,**kwargs):
        """
        Animates a spatial plot over all time steps
        
        Animation object is stored in the 'anim' property
        """
        #anim = unanimate(self.xy,self.data,self.tstep,xlim=self.grid.xlims,ylim=self.grid.ylims,clim=self.clim,**kwargs)
        #return anim
        
        
    
        # Create the figure and axes handles
        #plt.ion()
        fig = plt.gcf()
        ax = fig.gca()
        #ax.set_animated('True')
        
        # Find the colorbar limits if unspecified
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(self.data))
            self.clim.append(np.max(self.data))
           
        collection = PolyCollection(self.xy)
        collection.set_array(np.array(self.data[0,:]))
        collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        ax.add_collection(collection)    
        ax.set_xlim(self.grid.xlims)
        ax.set_ylim(self.grid.ylims)
        ax.axis('equal') 
        title=ax.set_title("")
        fig.colorbar(collection)
        
        def init():
            collection.set_array([])
            title.set_text("")
            return (collection,title)
               
        def updateScalar(i):
            collection.set_array(np.array(self.data[i,:]))
            collection.set_edgecolors(collection.to_rgba(np.array((self.data[i,:])))) 
            title.set_text(self.__genTitle(i))
            return (title,collection)
  
        self.anim = animation.FuncAnimation(fig, updateScalar, frames=len(self.tstep), interval=50, blit=True)

    def saveanim(self,outfile):
        """
        Save the animation object to an mp4 movie
        """
        
#        try:

        print 'Building animation sequence...'
        #self.anim.save(outfile, fps=15,bitrate=1800,extra_args=['-vcodec', 'libx264']) # plays in a web browser
        self.anim.save(outfile, fps=15,bitrate=1800)
        print 'Complete - animation saved to: %s'%outfile
#        except:
#            print 'Error with animation generation - check if either ffmpeg or mencoder are installed.'
            
    def animateVTK(self):
        """
        Animate a scene in the vtk window
        
        ## Not Working ##
        """
        
        # Load all the time steps into memory
        self.tstep=np.arange(0,len(self.time))
        nt = len(self.tstep)
        self.loadData()
        
        # Generate the initial plot
        points = np.column_stack((self.grid.xp,self.grid.yp,0.0*self.grid.xp))
        
        titlestr='%s [%s]\n Time: %s'%(self.long_name,self.units,\
                datetime.strftime(self.time[self.tstep[0]],'%d-%b-%Y %H:%M:%S'))
                
        self.fig,self.h,ug,title = unsurfm(points,self.grid.cells,np.array(self.data[0,:]),clim=self.clim,title=titlestr)        
        # Animate the plot by updating the scalar data in the unstructured grid object
        @mlab.animate
        def anim():
            ii=-1
            while 1:
                if ii<nt:
                    ii+=1
                else:
                    ii=0
                ug.cell_data.scalars = self.data[ii,:]
                titlestr='%s [%s]\n Time: %s'%(self.long_name,self.units,\
                datetime.strftime(self.time[self.tstep[0]],'%d-%b-%Y %H:%M:%S'))
                title.text=titlestr
                
                self.fig.scene.render()
                yield
        
        a = anim() # Starts the animation.

    def __del__(self):
        self.nc.close()
        
    def __openNc(self):
        #nc = Dataset(self.ncfile, 'r', format='NETCDF4') 
        try: 
            self.nc = MFDataset(self.ncfile, 'r')
        except:
            self.nc = Dataset(self.ncfile, 'r')
        
    def __genTitle(self,tt=0):
        
        if self.klayer>=0:
            zlayer = '%3.1f [m]'%self.grid.z_r[self.klayer]
        elif self.klayer==-1:
            zlayer = 'seabed'
        titlestr='%s [%s]\n z: %s, Time: %s'%(self.long_name,self.units,zlayer,\
                datetime.strftime(self.time[tt],'%d-%b-%Y %H:%M:%S'))
                
        return titlestr
        
    def __updateTstep(self):
        """
        Updates the tstep variable: -99 all steps, -1 last step
        """
        try:
            if self.tstep.any()==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep.any()==-1:
                self.tstep=len(self.time)-1
        except:
            if self.tstep==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep==-1:
                self.tstep=len(self.time)-1

            
###############################################################        
class Grid(object):
    
    """ Class for handling SUNTANS grid data"""
    
    def __init__(self,infile ,**kwargs):
               
        self.__dict__.update(kwargs)
        
        if os.path.isdir(infile):
            # Load ascii grid file
            self.infile = infile
            self.__loadascii()            
            
        else:
            # Load the grid fromm a netcdf file
            self.infile = infile
            self.__loadnc()

        # Find the grid limits
        self.xlims = [self.xp.min(),self.xp.max()]
        self.ylims = [self.yp.min(),self.yp.max()]
        
        # Cell polygon attribute for plotting     
        self.xy = self.cellxy()
        
    
    def __loadascii(self):
        """
        Load the grid variables from the ascii files: points.dat, edges.dat, cells.dat
        """
        pointdata = readTXT(self.infile+'/points.dat')
        celldata = readTXT(self.infile+'/cells.dat')
        edgedata = readTXT(self.infile+'/edges.dat')
        
        
        self.xp = pointdata[:,0]
        self.yp = pointdata[:,1]
        self.dv = pointdata[:,2] # zero to start
        
        self.xv = celldata[:,0]
        self.yv = celldata[:,1]
        self.cells = np.asarray(celldata[:,2:5],int)
        self.neigh = np.asarray(celldata[:,5:8])
        self.Nc = len(self.xv)
        
        self.edges = np.asarray(edgedata[:,0:2],int)
        self.mark = np.asarray(edgedata[:,2],int)
        self.grad = np.asarray(edgedata[:,3:5],int)
        
    def __loadnc(self):
        
        """Load the grid variables into the object from a netcdf file"""
        
        #print self.infile
        
        try: 
            nc = MFDataset(self.ncfile, 'r')
        except:
            nc = Dataset(self.ncfile, 'r')     
        
        self.xp = nc.variables['xp'][:]
        self.yp = nc.variables['yp'][:]
        self.xv = nc.variables['xv'][:]
        self.yv = nc.variables['yv'][:]
        self.dv = nc.variables['dv'][:]
        self.cells = nc.variables['cells'][:]
        self.Nc = len(self.xv)
        self.dz = nc.variables['dz'][:]
        self.Nk = nc.variables['Nk'][:]
        self.Nkmax = self.Nk.max()
        self.Nk -= 1 # needs to be zero based
        try:
            self.Ac = nc.variables['Ac'][:]
        except:
            print 'Warning no area variable, Ac, present...' 
        
        try:
            self.z_r = nc.variables['z_r'][:]
        except:
            print 'Warning no depth variable, z_r, present...' 
        
        
        
        #print nc.variables.keys()
        nc.close()
     
    def plot(self,**kwargs):
        """
          Plot the unstructured grid data
        """
        if self.__dict__.has_key('clim'):
            clim = self.clim
        else:
            clim = [self.dv.min(), self.dv.max()]
       
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,self.dv,xlim=self.xlims,ylim=self.ylims,\
            clim=clim,**kwargs)
        
        plt.title('SUNTANS Grid Bathymetry [m]')
        
    def plotvtk(self,**kwargs):
        """
          Plot the unstructured grid data using vtk libraries
        """
        if self.__dict__.has_key('clim'):
            clim = self.clim
        else:
            clim = [self.dv.min(), self.dv.max()]
        points = np.column_stack((self.xp,self.yp,0.0*self.xp))
        f,h = unsurfm(points,self.cells,self.dv,clim=clim,title='SUNTANS Grid Bathymetry [m]',\
            colormap='gist_earth')
        return f,h
        
            
    def cellxy(self):
        """ 
        Returns a list of Nx2 vectors containing the grid cell node coordinates
            
        Used by spatial ploting routines 
        """
        xynodes = []
        for ii in range(0,self.Nc):
            x=self.xp[self.cells[ii,:]]
            y=self.yp[self.cells[ii,:]]
            xynodes.append(closePoly(x,y))
        
        return xynodes
        
    def saveBathy(self,filename):
        """
            Saves the grid bathymetry to an xyz ascii file
        """
        f = open(filename,'w')
        
        for x,y,z in zip(self.xv,self.yv,self.dv):
            f.write('%10.6f %10.6f %10.6f\n'%(x,y,z))
            
        f.close()
      
class Profile(object):
    """
        Class for handling SUNTANS profile netcdf files
    """        
    
    def __init__(self,ncfile ,**kwargs):
        
       
       self.ncfile = ncfile
       
       self.__loadMeta()
       
       # Set defaults
       self.indices = np.arange(0,self.Np)
       self.tstep = 0 # -99 all steps, -1 last step
       self.klayer = np.arange(0,self.Nz)
       self.variable = 'u'
       self.xcoord = 'xp'
       self.ycoord = 'z'
       self.clim = None
       self.clevels = 12 # Number of contour levels
       
       # Update user-defined properties
       self.__dict__.update(kwargs)
        
       # Update tstep 
       self.__updateTstep()

    
    def __loadMeta(self):
        """
        Loads the metadata from the profile netcdf file
        """
        #nc = Dataset(self.ncfile, 'r', format='NETCDF4') 
        try: 
            nc = MFDataset(self.ncfile, 'r')
        except:
            nc = Dataset(self.ncfile, 'r')
        
        # Note that some of these variable names may change
        try:
            self.xp = nc.variables['x'][:]
        except:
            self.xp = nc.variables['xv'][:]
        try:    
            self.yp = nc.variables['y'][:]
        except:
            self.yp = nc.variables['yv'][:]
        try:
            self.dv = nc.variables['h'][:]
        except:
            self.dv = nc.variables['dv'][:]
        try:
            self.dz = nc.variables['vertspace'][:]
        except:
            self.dz = nc.variables['dz'][:]
        try:
            self.z = - nc.variables['vertdepth'][:]
        except:
            self.z =  - nc.variables['z_r'][:]
        try:
            self.Nk = nc.variables['klayers'][:]
        except:
            self.Nk = nc.variables['Nk'][:]
            
        self.Np = len(self.xp)
        self.Nz = len(self.dz)
        
        self.xlims = [self.xp.min(),self.xp.max()]
        self.ylims = [self.yp.min(),self.yp.max()]
        
        try:
            t = nc.variables['suntime']
        except:
            t = nc.variables['time']
        self.time = num2date(t[:],t.units)
        
        #print nc.variables.keys()
        nc.close()
        
    def loadData(self):
        """
        Loads the actual data for the given variable, indices, tsteps and zlayers
        """ 
        
        # "Higher order" variable stuff
        tmpvar = self.variable
        if tmpvar == 'ubar':
            self.variable = 'u'
        if tmpvar == 'vbar':
            self.variable = 'v'
        
        # Load the data
        
        #nc = Dataset(self.ncfile, 'r', format='NETCDF4')        
        try: 
            nc = MFDataset(self.ncfile, 'r')
        except:
            nc = Dataset(self.ncfile, 'r')
        
        self.long_name = nc.variables[self.variable].long_name
        self.units= nc.variables[self.variable].units
        #        ndims = len(nc.variables[self.variable].dimensions)
        #print nc.variables[self.variable].dimensions
        ndim = nc.variables[self.variable].ndim
        if ndim==1:
            self.data=nc.variables[self.variable][self.indices]
        elif ndim==2:
            self.data=nc.variables[self.variable][self.tstep,self.indices]
        else:
            self.data=nc.variables[self.variable][self.tstep,self.klayer,self.indices]
    
        nc.close()
        
        # Calculate the higher order variables from the raw data
        # To Add:
        #   uprime, vprime
        #   time mean variables
        #   seabed values of variables
        #   hydrostatic pressure perturbation?
        if tmpvar in ['ubar','vbar']:
            self.data=depthave(self.data,self.dz,np.abs(self.dv[self.indices]))
            
        self.variable = tmpvar
        
        
        # Update the colorbar limits if not set
        if self.clim==None:
            self.clim=[np.min(self.data),np.max(self.data)]
        
    
    def __loadXY(self):
        """
        Loads the X-Y coordinates used by 2-D plotting commands
        """
        if self.xcoord=='xp':
            self.xplot = self.xp[self.indices]
        elif self.xcoord=='yp':
            self.xplot = self.yp[self.indices]
        elif self.xcoord=='time':
            self.xplot = self.time[self.tstep]
            
        if self.ycoord=='z':
            self.yplot = self.z[self.klayer]
    
    def __updateTstep(self):
        """
        Updates the tstep variable: -99 all steps, -1 last step
        """
        try:
            if self.tstep.any()==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep.any()==-1:
                self.tstep=len(self.time)-1
        except:
            if self.tstep==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep==-1:
                self.tstep=len(self.time)-1
    def __checkDims(self):
        """
        Check that the dimensions sizes match for plotting
        
        If not transpose the data
        """        
        rc = np.shape(self.data)
        nx = len(self.xplot)
        ny = len(self.yplot)
        
        if ny!=rc[0] or nx !=rc[1]:
            self.data=np.transpose(self.data)
        
    def plotIndices(self):
        """
        Plots the locations of the points with the index numbers overlaid
        """
        offset=20
        plt.figure()
        plt.plot(self.xp,self.yp,'b.')
        plt.plot(self.xp[self.indices],self.yp[self.indices],'o',markeredgecolor='r',markerfacecolor=None)
        for s in range(self.Np):
            plt.text(self.xp[s]+offset,self.yp[s]+offset,'%d'%s)
            
        plt.axis('equal')
        plt.show()
        
    def pcolor(self,data=None,**kwargs):
        """
        Pcolor plot of the given variable (doesn't like time variable)
        """     
        if not self.__dict__.has_key('xplot'):
            self.__loadXY()
        if not self.__dict__.has_key('data'):
            self.loadData()
            self.__checkDims()  
        if data == None:
            data=self.data
          
        #plt.ion()
        self.fig=plt.gcf().s
        self.ax =plt.gca()
        self.h = plt.pcolor(self.xplot,self.yplot,data,**kwargs)
        self.cb = plt.colorbar()
        
        
    def contourf(self,data=None,V=None,**kwargs):
        """
        Filled contour plot of the given variable
        """
        if not self.__dict__.has_key('xplot'):
            self.__loadXY()
        if not self.__dict__.has_key('data'):
            self.loadData()
            self.__checkDims()  
        if data == None:
            data=self.data
        if V == None:
            V = np.linspace(self.clim[0],self.clim[1],num=self.clevels)    

        #plt.ion()
        self.fig=plt.gcf()
        self.ax =plt.gca()
        self.h = plt.contourf(self.xplot,self.yplot,data,V,**kwargs)
        self.cb = plt.colorbar()  
        
        
    def contour(self,data=None,V=None,**kwargs):
        """
        Contour plot of the given variable
        """
        if not self.__dict__.has_key('xplot'):
            self.__loadXY()
        if not self.__dict__.has_key('data'):
            self.loadData()
            self.__checkDims()  
        if data == None:
            data=self.data
        if V == None:
            V = np.linspace(self.clim[0],self.clim[1],num=self.clevels)
            
            
        #plt.ion()
        self.fig=plt.gcf()
        self.ax =plt.gca()
        self.h = plt.contour(self.xplot,self.yplot,data,V,**kwargs)
        self.cb = plt.colorbar()  


    def savefig(self,outfile,dpi=150):
        """ Saves a figure to file (matplotlib only)"""
        
        self.fig.savefig(outfile,dpi=dpi)
         
        print 'SUNTANS image saved to file:%s'%outfile
        
    def animate(self,fig=None,ax=None,h=None,cb=None,tsteps=None):
        """
        Method to animate the current method for multiple time steps
        """
         
        if tsteps == None:
            tsteps = np.arange(self.tstep,len(self.time))
            
        #pdb.set_trace()
        
        
        # Set the tsteps and grab the data
        tstepold = self.tstep
        self.tstep = tsteps
        self.loadData()
        
        if (ax==None) or (h ==None) or (fig==None):
            fig, ax, h, cb = self.pcolor(data=np.squeeze(self.data[0,:,:]))
         
        
        # Update the image
        ax.set_animated(True)
        h.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        cb.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        cb.draw_all()
        for ii in tsteps:
            zdata = np.squeeze(self.data[ii,:,:])
            h.set_array(np.ravel(zdata))
            ax.set_title('%d'%ii)
            fig.canvas.draw()
#        def updateanim(ii):
#            zdata = np.squeeze(self.data[ii,:,:])
#            h.set_array(np.ravel(zdata))
#            ax.set_title('%d'%ii)
#            
#        ani = animation.FuncAnimation(fig, updateanim, tsteps, interval=10)                
        
        # Set the data array to pre-animation array
        self.tstep = tstepold
        self.loadData()
        
    
    

####################################################################
#
# General functions to be used by all classes
#
####################################################################        
def closePoly(x,y):

    """ 
    Returns an Nx2 closed polygon for vectors x and y
        
    This output is required for plotting by unsurf. 
    """
    
    nx = len(x)
    ny = len(y)
    if nx != ny:
        print "Error: The lengths of vector x and y must be equal"
        return

    x = np.reshape(x,(nx,1))
    y = np.reshape(y,(ny,1))

    x = np.vstack((x,x[0]))
    y = np.vstack((y,y[0]))
    
    return np.hstack((x,y))

def depthave(data,dz,h):
    """ Calculate the depth average of a variable
    Variable should have dimension: [nz*nx] or [nt*nz*nx]
    
    """
    ndim = np.ndim(data)
    if ndim == 2:
        return depthint(data,dz) / h
    elif ndim == 3:
        nt = np.size(data,0)
        h = np.tile(h,(nt,1))
        return depthint(data,dz) / h
            
def depthint(data,dz):
    """
    Calculates the depth integral of a variable: data
    Variable should have dimension: [nz*nx] or [nt*nz*nx]
    
    """
    ndim = np.ndim(data)
    
    nz = np.size(dz)
    
    if ndim == 2:
        nx = np.size(data,1)
        dz2=np.reshape(dz,(nz,1))
        dz2 = np.tile(dz2,(1,nx))
        return np.sum(data*dz2,axis=0)
    if ndim == 3:
        nt = np.size(data,0)
        nx = np.size(data,2)
        dz3=np.reshape(dz,(1,nz,1))
        dz3 = np.tile(dz3,(nt,1,nx))
        return np.sum(data*dz3,axis=0)
        
def readTXT(fname,sep=None):
    """
    Reads a txt file into an array of floats
    """
    
    fp = file(fname,'rt')
    data = np.array([map(float,line.split(sep)) for line in fp])
    fp.close()
    
    return data

def unsurfm(points, cells, z,clim=None,title=None,**kwargs):
    """
    Plot cell-centred data using the mayavi/tvtk libraries
    
    """        
    if clim==None:
        clim=[]
        clim.append(np.min(z))
        clim.append(np.max(z))
    
    try:    
        tri_type = tvtk.Triangle().cell_type
    except:
        # Load tvtk libraries here as they slow the script down
        from tvtk.api import tvtk
        from mayavi import mlab
        tri_type = tvtk.Triangle().cell_type
        
    ug = tvtk.UnstructuredGrid(points=points)
    ug.set_cells(tri_type, cells)
    
    ug.cell_data.scalars = z
    ug.cell_data.scalars.name = 'suntans_scalar'
    
    d = mlab.pipeline.add_dataset(ug)
    h=mlab.pipeline.surface(d,vmin=clim[0],vmax=clim[1],**kwargs)
    f=mlab.gcf()
    f.scene.background = (0.,0.,0.)
    mlab.colorbar(object=h,orientation='vertical')
    mlab.view(0,0)
    title=mlab.title(title,height=0.95,size=0.15)
    
    return f, h, ug,title
    
def unsurf(xy,z,xlim=[0,1],ylim=[0,1],clim=None,**kwargs):
    """
    Plot cell-centred data on an unstructured grid as a series of patches
        
    Similar functionality to the suntans matlab unsurf function.
        
    Inputs:
        xy -list[Nc] of N*2 arrays of closed polygons
        z - scalar vector[Nc]
            
    """     
    
    # Create the figure and axes handles
    plt.ioff()
    fig = plt.gcf()
    ax = fig.gca()
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Find the colorbar limits if unspecified
    if clim==None:
        clim=[]
        clim.append(np.min(z))
        clim.append(np.max(z))
    
    collection = PolyCollection(xy)
    collection.set_array(np.array(z))
    collection.set_array(np.array(z))
    collection.set_clim(vmin=clim[0],vmax=clim[1])
    #collection.set_linewidth(0)
    collection.set_edgecolors(collection.to_rgba(np.array(z)))    
    
    ax.add_collection(collection)
    ax.axis('equal')
    axcb = fig.colorbar(collection)
    
    return fig, ax, collection, axcb
    
def unanimate(xy,z,tsteps,xlim=[0,1],ylim=[0,1],clim=None,**kwargs):
    """
    Plot cell-centred data on an unstructured grid as a series of patches
        
    Similar functionality to the suntans matlab unsurf function.
        
    Inputs:
        xy -list[Nc] of N*2 arrays of closed polygons
        z - scalar array [nt x Nc]
        t - vector of time steps
            
    """     
    from matplotlib import animation
    
    # Create the figure and axes handles
    #plt.ion()
    fig = plt.gcf()
    ax = fig.gca()
    #ax.set_animated('True')
    
    # Find the colorbar limits if unspecified
    if clim==None:
        clim=[]
        clim.append(np.min(z))
        clim.append(np.max(z))
    

        
    collection = PolyCollection(xy)
    collection.set_array(np.array(z[0,:]))
    collection.set_clim(vmin=clim[0],vmax=clim[1])
    ax.add_collection(collection)    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.axis('equal') 
    title=ax.set_title("")
    fig.colorbar(collection)
    
    def init():
        collection.set_array([])
        title.set_text("")
        return (collection,title)

        
    def updateScalar(i):
        ts = i
        collection.set_array(np.array(z[i,:]))
        collection.set_edgecolors(collection.to_rgba(np.array((z[i,:]))))    
        title.set_text('%d'%ts)
        return (title,collection)

    
    anim = animation.FuncAnimation(fig, updateScalar, frames=200, interval=50, blit=True)
    return anim
#    print 'Building animation sequence...'
#    anim.save('C:/Projects/GOMGalveston/CODE/PYTHON/SUNTANS/test_animation.mp4', fps=15)


def usage():
    """
    Command line usage output
    """
    print "--------------------------------------------------------------"
    print "sunpy.py   -h                 # show this help message      "
    print "         -f suntans.nc        # SUNTANS output netcdf file  "
    print "         -v varname           # variable name to plot [default: u]      "
    print "         -k  N                # vertical layer to plot [default: 0]"
    print "         -t  N                # time step to plot. -1 = last step. [default: 0]"
    print "         -j  N                # grid cell index to plot (timeseries only) [default: 0]"
    print '         -c "N N"             # Color bar limits !! IN DOUBLE QUOTES !! [default: None]'
    print "         -s figure.png        # Save to a figure"
    print "         --animate            # Animate plot through all time steps"
    print "         --timeseries         # Plot time series of individual point"
    print "         --profile            # time-depth contour plot at cell: j"
    print "         --vtk                # Use the vtk plotting libraries"
    print "\n Example Usage:"
    print "--------"
    print " python sunpy.py -f suntans.nc -v temp -k 5 -t 10 -c '10 29' "
    print ""

    
if __name__ == '__main__':
    """
    Command line call to sunpy
    """        
    
    # Defaults
    k = 0
    t = 0
    j = 0
    varname = 'u'
    plottype = 0
    usevtk = False
    clim = None
    save = False
    
    try:
            opts,rest = getopt.getopt(sys.argv[1:],'hf:v:t:k:j:c:s:',
                                      ['animate',
                                       'timeseries',
                                       'profile',
                                       'vtk'])
    except getopt.GetoptError,e:
        print e
        print "-"*80
        usage()
        exit(1)

    for opt,val in opts:
        if opt == '-h':
            usage()
            exit(1)
        elif opt == '-f':
            ncfile = val
        elif opt == '-v':
            varname = val
        elif opt == '-t':
            t = int(val)
        elif opt == '-k':
            k = int(val)
        elif opt == '-j':
            j = int(val)
        elif opt == '-c':
            clim = np.asarray(val.split(),dtype=float)
        elif opt == '-s':
            outfile = val
            save = True
        elif opt == '--animate':
            plottype = 1
        elif opt == '--timeseries':
            plottype = 2
        elif opt == '--profile':
            plottype = 3
        elif opt == '--vtk':
            usevtk = True
            
    # Load the class and plot
    if plottype in [0,1,2]:
        sun = Spatial(ncfile,klayer=k,tstep=t,variable=varname,clim=clim)
    
    if plottype == 0:
        # Spatial Plot
        if usevtk:
            #mlab.figure(size=(640,480))
            sun.plotvtk()
            if save:
                sun.savefig(outfile)
            #mlab.show()

        else:
            plt.figure(figsize=(10,7))
            sun.plot() 
            if save:
                sun.savefig(outfile)
            #plt.show()
        
    elif plottype == 1:
        # Animation
        sun.tstep=np.arange(0,len(sun.time))
        sun.loadData()
        plt.figure(figsize=(10,7))
        sun.animate()
        if save:
            sun.saveanim(outfile)
    
    elif plottype == 2:
        # Time series plot
        plt.figure(figsize=(10,7))
        sun.plotTS(j=j)
        if save:
            sun.savefig(outfile)
        #plt.show()
    
    elif plottype == 3:
        sun = Profile(ncfile,indices=j,tstep=-99,xcoord='time',variable=varname,clim=clim)

        plt.figure(figsize=(10,5))
        sun.contourf()
        titlestr='%s [%s]'%(sun.long_name,sun.units)
        plt.title(titlestr)
        sun.ax.set_xlabel('Time')
        sun.ax.set_ylabel('z [m]')

        if save:
            sun.savefig(outfile)


############################
# Testing stuff -spatial plots and animation
#ncfile = 'C:/Projects/GOMGalveston/MODELLING/GalvestonSquare/rundata/suntan_output.nc.0'
#sun = Spatial(ncfile,klayer=0,tstep=20)
#plottype = 1 # 0 - spatial, 1 - animation, 2 - time series
#sun.variable='Tair'
#
#if plottype == 0:
#    # Spatial Plot
##    plt.figure(figsize=(10,7))
##    sun.plot()
##    plt.show()
#    
#    mlab.figure(size=(640,480))
#    sun.plotvtk()
#    
#    sun.savefig('C:/Projects/GOMGalveston/MODELLING/GalvestonSquare/rundata/test.png')
#    
#elif plottype == 1:
#    # Animation
##    sun.tstep=np.arange(0,len(sun.time))
##    sun.loadData()
##    ani=unanimate(sun.xy,sun.data,sun.tstep,xlim=sun.grid.xlims,ylim=sun.grid.ylims)
#    sun.animateVTK()
#
#elif plottype == 2:
#    # Time series plot
#    sun.plotTS(j=80)
    

############################
# Profile testing
#ncfile =  'E:/Projects/ScottReef/MODELLING/SUNTANS/netcdf/ScottRfTVD3_prof.nc';
#
#sun = Profile(ncfile,indices=np.arange(26,528),tstep=100,variable='ubar')
#sun.loadData()
#u=sun.data
#fig, ax, h, cb 2= sun.contourf()
#sun.plotIndices()