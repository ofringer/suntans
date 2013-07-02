# -*- coding: utf-8 -*-
"""
SUNTANS particle tracking 

Created on Wed Apr 17 09:54:48 2013

@author: mrayson
"""

from sunpy import Spatial, Grid
import othertime
from trisearch import TriSearch

from datetime import datetime,timedelta
from scipy import spatial
import numpy as np
import matplotlib.pyplot as plt


from time import clock
import pdb


class SunTrack(Spatial):
    """
    Particle tracking class
    """
    
    verbose = True
    
    # Interpolation method
    interp_method = 'idw' # 'idw' or 'nearest' or 'mesh'
    
    advect_method = 'rk2' # 'euler' or 'rk2'
    
    # 
    
    def __init__(self,ncfile,**kwargs):
        """
        Initialize the 3-D grid, etc
        """
        
        self.__dict__.update(kwargs)

        Spatial.__init__(self,ncfile,klayer=[-99],**kwargs)
            
        #self.klayer=np.arange(0,self.Nkmax)
        
        # Initialise the 3-D grid
        self.init3Dgrid()
            
        # Step 2) Initialise the interpolation function
        if self.interp_method in ('idw','nearest'):
            self.UVWinterp = interp3D(self.xv3d,self.yv3d,self.zv3d,method=self.interp_method)
            # Interpolation function to find free surface and seabed
            self.Hinterp = interp3D(self.xv,self.yv,0*self.xv,method='nearest')
            
        elif self.interp_method == 'mesh':
            self.UVWinterp = interp3Dmesh(self.xp,self.yp,-self.z_w,self.cells,self.mask3D)
            



    def __call__(self,x,y,z, timeinfo, runmodel=True, **kwargs):
        """
        Run the particle model 
        
        Inputs:
            xyz - release locations of particles [n_parts x 3]
            timeinfo - tuple (starttime,endtime,dt) 
        """

        self.__dict__.update(kwargs)
        
        self.getTime(timeinfo)
        
        # Initialise the particle dictionary
        self.particles={'X':x,'Y':y,'Z':z}
        
        if self.verbose:
            print '#######################################################'
            print 'Running SUNTANS particle tracking\n'
            print 'Releasing %d particles.\n'%self.particles['X'].shape[0]
            print '######################################################'
            
        
        # Initialise the currents
        self.initCurrents()
        
        # Start time stepping
        if runmodel:
            for ii,time in enumerate(self.time_track):
                # Step 1) Advect particles
                self.advectParticles(time,self.time_track_sec[ii])
                # Step 2) Check bounds
                
    def animate(self,x,y,z,timeinfo,xlims=None,ylims=None,outfile=None):
        """
        Animate the particles on the fly using matplotlib plotting routines
        """
        import matplotlib.animation as animation
        
        # Set the xy limits
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        self.__call__(x,y,z,timeinfo,runmodel=False)
        # Plot a map of the bathymetry
        self.fig = plt.figure(figsize=(10,8))
        ax=plt.gca()
        self.contourf(z=-self.dv,clevs=30,titlestr='',xlims=xlims,ylims=ylims,cmap='bone')
        
        # Plot the particles at the first time step
        h1 = plt.plot([],[],'y.',markersize=1.0)
        h1 = h1[0]
        
        title=ax.set_title("")
        
        def init():
            h1.set_xdata(self.particles['X'])
            h1.set_ydata(self.particles['Y'])
            title=ax.set_title(self.genTitle(0))
            return (h1,title)

        def updateLocation(ii):
            if ii==0:
                # Re-initialise the particle location
                self.particles={'X':x,'Y':y,'Z':z}
                
            self.advectParticles(self.time_track[ii],self.time_track_sec[ii])
            h1.set_xdata(self.particles['X'])
            h1.set_ydata(self.particles['Y'])
            title.set_text(self.genTitle(ii))
            return (h1,title)
        
        self.anim = animation.FuncAnimation(self.fig, updateLocation, frames=len(self.time_track), interval=50, blit=True)
        
        if outfile==None:
            plt.show()
        else:
            self.saveanim(outfile)
            
    def animate3D(self,x,y,z,timeinfo,outfile=None):
        """
        3D animation of the particles using mayavi
        """
        
        from mayavi import mlab
        from suntvtk import SunTvtk
        
        # Initiate the particle model        
        self.__call__(x,y,z,timeinfo,runmodel=False)
        
        # Initiate the suntans tvtk class
        self.vtk = SunTvtk(self.ncfile)
        
        # Plot the bathymetry
        self.vtk.plotbathy3d(colormap='bone')
        
        # Plot the particles
        self.vtkobj = mlab.points3d(self.particles['X'],self.particles['Y'],self.particles['Z']*self.vtk.zscale,\
            color=(1.0,1.0,0.0),scale_mode='none',scale_factor=100.0,opacity=0.8)
                      
        nt = len(self.time_track)         
    
        @mlab.animate
        def anim():
            ii=-1
            while 1:
                if ii<nt-1:
                    ii+=1
                else:
                    ii=0
                
                # Advect the particles
                self.advectParticles(self.time_track[ii],self.time_track_sec[ii])
                
                # Update the plot object
                self.vtkobj.mlab_source.x = self.particles['X']
                self.vtkobj.mlab_source.y = self.particles['Y']
                self.vtkobj.mlab_source.z = self.particles['Z']*self.vtk.zscale
               
                #titlestr=self._Spatial__genTitle(tt=ii)
                #self.title.text=titlestr
                
                self.vtk.fig.scene.render()
                yield
        
        if outfile==None:
            anim() # Starts the animation.
    
    def saveanim3D(self,outfile,frate=15):
        """
        Saves an animation of the current scene
        """
        from mayavi import mlab
        import os
        # This works for mp4 and mov...
        cmdstring ='ffmpeg -r %d -i ./.tmpanim%%04d.png -y -loglevel quiet -c:v libx264 -crf 23 -pix_fmt yuv420p %s'%(frate,outfile)
        
        self.vtk.fig.scene.anti_aliasing_frames = 0
        self.vtk.fig.scene.disable_render = True
        
        nt = len(self.time_track) 
        png_list=[]
        for ii in range(nt):

            self.advectParticles(self.time_track[ii],self.time_track_sec[ii])
                
            # Update the plot object
            self.vtkobj.mlab_source.x = self.particles['X']
            self.vtkobj.mlab_source.y = self.particles['Y']
            self.vtkobj.mlab_source.z = self.particles['Z']*self.vtk.zscale
            
            self.vtk.fig.scene.render()
            
            # This bit saves each img
            outimg='./.tmpanim%04d.png'%ii
            png_list.append(outimg)
#            self.fig.scene.save_png(outimg)
            mlab.savefig(outimg,figure=self.vtk.fig)

            print 'Saving image %s of %d...'%(outimg,nt)
            
        # Call ffmpeg within python
        os.system(cmdstring)
        
        print '####\n Animation saved to: \n %s\n####' % outfile
        # Delete the images
        print 'Cleaning up temporary images...'
        for ff in png_list:
            os.remove(ff)
        print 'Complete.'
        
    def genTitle(self,ii):
            return 'Particle time step: %s'%datetime.strftime(self.time_track[ii],'%d-%b-%Y %H:%M:%S')
    
    def advectParticles(self,timenow,tsec):
        """
        Advect the particles
        """          
        t0 = clock()
        if self.verbose:
            print '\tTime step: ',timenow
         
        
        if self.advect_method=='euler':
            self.euler(timenow,tsec)
            
        elif self.advect_method=='rk2':
            self.rk2(timenow,tsec)
            
        else:
            raise Exception, 'unknown advection scheme: %s. Must be "euler" or "rk2"'%self.advect_method
            
        
        t1 = clock()
        if self.verbose:
            print '\t\tElapsed time: %s seconds.'%(t1-t0)
            
    def euler(self,timenow,tsec):
        """
        euler time integration
        """
        self.updateCurrents(timenow,tsec)
        
        # Interpolate the currents
        u = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.u)
        v = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.v)
        w = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.w)
        
        self.particles['X'] += u*self.dt
        self.particles['Y'] += v*self.dt
        self.particles['Z'] += w*self.dt
        
        # Check the vertical bounds of a particle 
        self.particles['Z'] = self.checkVerticalBounds(self.particles['X'],self.particles['Y'],self.particles['Z'])
        
    def rk2(self,timenow,tsec):
        """
        2nd order Runge-Kutta advection scheme
        """
        
        self.updateCurrents(timenow,tsec)
        
        # Interpolate the currents
        u = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.u)
        v = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.v)
        w = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.w)

        x1 = self.particles['X'] + 0.5*self.dt*u
        y1 = self.particles['Y'] + 0.5*self.dt*v
        z1 = self.particles['Z'] + 0.5*self.dt*w
        
        # Check the vertical bounds of a particle 
        z1 = self.checkVerticalBounds(x1,y1,z1)
        
        # Update the currents again
        self.updateCurrents(timenow+timedelta(seconds=self.dt*0.5),tsec+self.dt*0.5)
        
        u = self.UVWinterp(x1,y1,z1,self.u)
        v = self.UVWinterp(x1,y1,z1,self.v)
        w = self.UVWinterp(x1,y1,z1,self.w)
        
        self.particles['X'] += u*self.dt
        self.particles['Y'] += v*self.dt
        self.particles['Z'] += w*self.dt
        
        # Check the vertical bounds of a particle again 
        self.particles['Z'] = self.checkVerticalBounds(self.particles['X'],self.particles['Y'],self.particles['Z'])
        
    def checkVerticalBounds(self,x,y,z):
        """
        Checks that particles are not above the surface or below the seabed
        
        """
        SMALL = 0.001
        
        # find the free surface and seabed at the particle locations
        if not self.interp_method == 'mesh':
            eta_P = self.Hinterp(x,y,z,self.eta)
            h_P = self.Hinterp(x,y,z,-self.dv)
        else:
            ind = self.UVWinterp.cellind
            mask=ind==-1
            ind[mask]=0.0
            eta_P = self.eta[ind]
            h_P = -self.dv[ind]
            eta_P[mask]=0.0
            h_P[mask]=0.0
        
        np.where(z > eta_P, eta_P-SMALL, z)
        np.where(z < h_P, h_P+SMALL, z)
        
        return z
        
    
    def initCurrents(self):
        """
        Initialise the forward and backward currents and free-surface height
        """
        
        tindex = othertime.findGreater(self.time_track[0],self.time)
        
        self.time_index = tindex
        
        # Check the time index here
        if tindex == None:
            raise Exception, 'start time less than model time: ',self.time[0]
        elif self.time_index==0:
            self.time_index=1
            
        self.uT = np.zeros((self.nActive,2))
        self.vT = np.zeros((self.nActive,2))
        self.wT = np.zeros((self.nActive,2))
        self.etaT = np.zeros((self.Nc,2))
        
        self.uT[:,0], self.vT[:,0], self.wT[:,0], self.etaT[:,0] = self.getUVWh(tindex-1)
        self.uT[:,1], self.vT[:,1], self.wT[:,1], self.etaT[:,1] = self.getUVWh(tindex)
        
        self.timeInterpUVW(self.time_track_sec[0],self.time_index)
        
    def updateCurrents(self,timenow,tsec):
        """
        Checks to see if the currents need updating
        """
        
        tindex = othertime.findGreater(timenow,self.time)
        
        if not tindex == self.time_index:
            print 'Reading SUNTANS currents at time: ',timenow
            self.uT[:,0]=self.uT[:,1]
            self.vT[:,0]=self.vT[:,1]
            self.wT[:,0]=self.wT[:,1]
            self.etaT[:,0]=self.etaT[:,1]
            self.uT[:,1], self.vT[:,1], self.wT[:,1], self.etaT[:,1] = self.getUVWh(tindex)
            
            self.time_index = tindex
        
        # Temporally interpolate onto the model step
        self.timeInterpUVW(tsec,self.time_index) 
        
    def timeInterpUVW(self,tsec,tindex):
        """
        Temporally interpolate the currents  onto the particle time step, tsec.
        
        Linear interpolation used.
        """
        
        t0 = self.time_sec[tindex-1]
        t1 = self.time_sec[tindex]
        dt = t1 - t0
        
        w1 = (tsec-t0)/dt
        w2 = 1.0-w1
        
        self.u = self.uT[:,0]*w1 + self.uT[:,1]*w2 
        self.v = self.vT[:,0]*w1 + self.vT[:,1]*w2 
        self.w = self.wT[:,0]*w1 + self.wT[:,1]*w2 
        self.eta = self.etaT[:,0]*w1 + self.etaT[:,1]*w2 
        
    def getUVWh(self,tstep):
        """
        Load the velocity arrays
        """
        self.tstep = tstep
        
        u = self.loadData(variable='uc')
        v = self.loadData(variable='vc')
        w = self.loadData(variable='w')
        eta = self.loadData(variable='eta')

        return u[self.mask3D], v[self.mask3D], w[self.mask3D], eta    
        
    def getTime(self,timeinfo):
        """
        Get the particle time step info
        """
        self.dt = timeinfo[2]
        
        self.time_track = othertime.TimeVector(timeinfo[0],timeinfo[1],timeinfo[2])
        
        self.time_track_sec = othertime.SecondsSince(self.time_track)
        self.time_sec = othertime.SecondsSince(self.time)
        
        self.time_index = -9999
        
        
    def init3Dgrid(self):
        """
        Constructs the 3d cell unstructured grid object
        
        Includes only active vertical grid layers
        
        """
        nc = self.Nc
        nz = self.Nkmax+1
        nv = len(self.xp)
        
        self.returnMask3D()
        self.nActive = np.sum(self.mask3D) # Total number of active cells
        
        self.cells3d = np.zeros((self.nActive,6))
        self.xv3d = np.zeros((self.nActive,))
        self.yv3d = np.zeros((self.nActive,))
        self.zv3d = np.zeros((self.nActive,))
        pt1=0
        for k in range(1,nz):
            masklayer = self.mask3D[k-1,:]
            nc = np.sum(masklayer)
            pt2 = pt1+nc            
            self.cells3d[pt1:pt2,0] = self.cells[masklayer,0]+(k-1)*nv
            self.cells3d[pt1:pt2,1] = self.cells[masklayer,1]+(k-1)*nv
            self.cells3d[pt1:pt2,2] = self.cells[masklayer,2]+(k-1)*nv
            self.cells3d[pt1:pt2,3] = self.cells[masklayer,0]+k*nv
            self.cells3d[pt1:pt2,4] = self.cells[masklayer,1]+k*nv
            self.cells3d[pt1:pt2,5] = self.cells[masklayer,2]+k*nv
            
            self.xv3d[pt1:pt2] = self.xv[masklayer]
            self.yv3d[pt1:pt2] = self.yv[masklayer]
            self.zv3d[pt1:pt2] = -self.z_r[k-1]
            
            pt1=pt2
            #print k, nc
            
        self.verts = np.zeros((nv*nz,3))
        pv1 = 0
        for k in range(0,nz):
            self.verts[pv1:pv1+nv,0] = self.xp
            self.verts[pv1:pv1+nv,1] = self.yp
            self.verts[pv1:pv1+nv,2] = -self.z_w[k]
            pv1 += nv

    
    def returnMask3D(self):
        """
        Returns the 3D mask [Nk x Nc] True = Active, False = Ghost
        
        """
        self.mask3D = np.ones((self.Nkmax,self.Nc),dtype=bool)
        
        for i in range(self.Nc):
            if self.Nk[i] == self.Nkmax:
                Nk = self.Nk[i]
            else:
                Nk = self.Nk[i]+1
            self.mask3D[Nk:self.Nkmax,i]=False
            

class interp3Dmesh(TriSearch):
    """
    3D interpolation class for an unstructured grid
    
    Uses knowledge of the mesh to efficiently update the velocities based on 
    the particle location.
    """
    
    def __init__(self,x,y,z,cells,mask):
        
        # Initialise the trisearch array
        TriSearch.__init__(self,x,y,cells)
        
        self.z = np.sort(z)
        self.Nk = z.size-1
        
        self.mask3d = mask
        
        self.maskindex = np.zeros(self.mask3d.shape,dtype=np.int32)
        rr=0
        for ii in range(self.mask3d.shape[0]):
            for jj in range(self.mask3d.shape[1]):
                if self.mask3d[ii,jj]:
                    self.maskindex[ii,jj]=rr
                    rr+=1
                  
    def __call__(self,X,Y,Z,data):
        
        # The Update the cell index using TriSearch class
        if not self.__dict__.has_key('cellind'):
            print ' Finding initial particle index...'
            TriSearch.__call__(self,X,Y)
        else:
            if np.sum(np.abs(X-self.xpt))>0+1e-8:
                #print ' updating location index...'
                self.updatexy(X,Y)
                        
        # Find the k-index
        kind=self.z.searchsorted(Z)
        kind = self.Nk - kind 
        
        kind[kind>=self.Nk-1] = self.Nk-1
                
        # Calculate the linear index
        #ind = kind*self.mask.shape[1]+self.cellind
        
        ind = self.maskindex[kind,self.cellind]
        maskpts = self.mask3d[kind,self.cellind]
        
        ind[maskpts == False] = 0
        
        # Return the nearest data point (... for now)
        dataout = data[ind]
                
        # Mask bogey points
        dataout[maskpts==False]=0.0
        dataout[self.cellind==-1]=0.0
                
        return dataout
    

class interp3D(object):
    """
    3D interpolation class
    
    This is a "points-in-space" interpolation method and 
    has no knowledge of the mesh.
    """
    
    method = 'idw' # 'idw' or 'nearest'
    
    def __init__(self,Xin,Yin,Zin,**kwargs):
        """
        Initialise the interpolation object
        """
        self.__dict__.update(kwargs)
        
        self.XYZin = np.vstack((Xin,Yin,Zin)).T
        
        if self.method == 'idw':
            self.F = idw(self.XYZin)
            
        elif self.method == 'nearest':
            self.F = nearest(self.XYZin)
            
        else:
            raise Exception, 'Unknown interpolation method: %s'%self.method
            
    def __call__(self,X,Y,Z,data):
        
        XYZ = np.vstack((X,Y,Z)).T

        return self.F(XYZ,data)
            
        
        
class idw(object):
    """
    Inverse distance weighted interpolation function
    """
    
    maxdist=2000.0
    NNear=6
    p=1.0
    zaspect = 1000.0 # scales the vertical coordinate by this to avoid selecting all points in the vertical
    
    def __init__(self,XYZin,**kwargs):
        
        self.__dict__.update(kwargs)
        
        # Compute the spatial tree 
        XYZin[:,2] = self.zaspect * XYZin[:,2]
        self.kd = spatial.cKDTree(XYZin)
    
    def __call__(self,XYZout,Zin):
        
        EPS = 1e-9
        
        XYZout[:,2] = self.zaspect * XYZout[:,2]
         # Perform query on all of the points in the grid
        dist,self.ind = self.kd.query(XYZout,distance_upper_bound=self.maxdist,k=self.NNear)
        
        # Calculate the weights
        dist += EPS # Add small value to avoid divide by zero
        self.W = 1.0/dist**self.p
        Wsum = np.sum(self.W,axis=1)
        
        for ii in range(self.NNear):
            self.W[:,ii] = self.W[:,ii]/Wsum
            
        # create the mask
        mask = (dist==np.inf)
        self.ind[mask]=1
        
        # Fill the array and resize it
        Zin = np.squeeze(Zin[self.ind])
        
        # Compute the weighted sums and mask the blank points
        return np.sum(Zin*self.W,axis=1)
        #Z[mask]=np.nan
        
class nearest(object):
    """ 
    Nearest neighbour interpolation algorithm
    Sets any points outside of maxdist to NaN
    """
    maxdist = 2000.0
    
    def __init__(self,XYZin,**kwargs):
        self.__dict__.update(kwargs)
        
        # Compute the spatial tree
        self.kd = spatial.cKDTree(XYZin)

    def __call__(self,XYZout,Zin):
        
        # Perform query on all of the points in the grid
        dist,self.ind=self.kd.query(XYZout,distance_upper_bound=self.maxdist)
        # create the mask
        self.mask = (dist==np.inf)
        self.ind[self.mask]=1
        
        # Fill the array and resize it
        Z = Zin[self.ind]
        Z[self.mask]=np.nan
        
        return Z            

def GridParticles(grdfile,dx,dy,nz):
    """
    Returns the locations of particles on a regular grid inside of suntans grid
    """
    
    # Load the suntans grid
    sun = Grid(grdfile)

    # Load a trisearch object    
    tri = TriSearch(sun.xp,sun.yp,sun.cells,verbose=False)
    
    # Construct a 2D mesh of particles
    x = np.arange(sun.xlims[0],sun.xlims[1],dx)
    y = np.arange(sun.ylims[0],sun.ylims[1],dy)
    
    X,Y = np.meshgrid(x,y)
    
    X=X.ravel()
    Y=Y.ravel()
    # Check which particles are inside the grid 
    cellind = tri(X,Y)
    
    mask = cellind!=-1
    xout = X[mask]
    yout = Y[mask]
    
    nx = xout.shape[0]
    
    # Construct the 3D mesh
    xout = np.repeat(xout.reshape((nx,1)),nz,axis=1)
    yout = np.repeat(yout.reshape((nx,1)),nz,axis=1)
    
    zout = np.linspace(0.05,0.95,nz)
    zout = np.repeat(zout.reshape((nz,1)),nx,axis=1)
        
    zout *= -sun.dv[cellind[mask]]
    zout = zout.T
    
    return xout.ravel(), yout.ravel(), zout.ravel()
    
    
    
    
    
    