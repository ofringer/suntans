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
from scipy import spatial, sparse
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date
import matplotlib.nxutils as nxutils #inpolygon equivalent lives here
import operator

from time import clock
import pdb


class SunTrack(Spatial):
    """
    Particle tracking class
    """
    
    verbose = True
    
    # Interpolation method
    interp_method = 'mesh' # 'idw' or 'nearest' or 'mesh'
    # IF interp_method == 'mesh'
    interp_meshmethod = 'nearest' # 'linear' or 'nearest'
    
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
	    if self.interp_meshmethod == 'nearest':
		self.UVWinterp = interp3Dmesh(self.xp,self.yp,-self.z_w,self.cells,self.mask3D,method='nearest')
	    elif self.interp_meshmethod == 'linear':
            	self.UVWinterp = interp3Dmesh(self.xp,self.yp,-self.z_w,self.cells,self.mask3D,method='linear',grdfile=self.ncfile)



    def __call__(self,x,y,z, timeinfo,outfile=None,dtout=3600.0,\
    	tstart=None,agepoly=None,age=None,agemax=None,runmodel=True, **kwargs):
        """
        Run the particle model 
        
        Inputs:
            x, y, z - release locations of particles [n_parts x 3]
            timeinfo - tuple (starttime,endtime,dt) 

	    (optional)
	    tstart - start time of each particle (format seconds since 1990-01-01)
	    	if None defaults to starting all particles from starttime.
	    agepoly - x/y coordinates of a polygon which ages particles
	    age, agemax - vectors of length n_parts containing the initial age and agemax. (leave
	    	as None to set to zero)
        """

        self.__dict__.update(kwargs)
        
        self.getTime(timeinfo)
        
        # Initialise the particle dictionary
        self.particles={'X':x,'Y':y,'Z':z,'X0':x,'Y0':y,'Z0':z}
	
	# Initialise the age calculation
	self._calcage = False
	if not agepoly == None:
	    self._calcage = True
	    self.agepoly = agepoly
	    if age==None:
	    	age=np.zeros_like(x)
	    if agemax==None:
	    	agemax=np.zeros_like(x)
	    
	self.particles.update({'age':age,'agemax':agemax})

	# Set the particle time start and activate if necessary
	self.initParticleTime(tstart)

	# Initialse the outut netcdf file
	if not outfile==None:
            self.initParticleNC(outfile,x.shape[0],age=self._calcage)

        if self.verbose:
            print '#######################################################'
            print 'Running SUNTANS particle tracking\n'
            print 'Releasing %d particles.\n'%self.particles['X'].shape[0]
            print '######################################################'
            
        
        # Initialise the currents
        self.initCurrents()
        
        # Start time stepping
	tctr=dtout
	ctr=0
        if runmodel:
            for ii,time in enumerate(self.time_track):
	    	tctr+=self.dt
                # Step 1) Advect particles
                self.advectParticles(time,self.time_track_sec[ii])

                # Write output if needed 
		if not outfile==None and tctr>=dtout:
		    self.writeParticleNC(outfile,self.particles['X'],self.particles['Y'],self.particles['Z'],\
		    	self.time_track_sec[ii],ctr,age=self.particles['age'],agemax=self.particles['agemax'])
		    tctr=0
		    ctr+=1
    
    def animate(self,x,y,z,timeinfo,agepoly=None,agemax=7.0,xlims=None,ylims=None,outfile=None):
        """
        Animate the particles on the fly using matplotlib plotting routines
        """
        import matplotlib.animation as animation
        
	agescale = 1.0/86400.0

        # Set the xy limits
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        self.__call__(x,y,z,timeinfo,agepoly=agepoly,runmodel=False)
        # Plot a map of the bathymetry
        self.fig = plt.figure(figsize=(10,8))
        ax=plt.gca()
        self.contourf(z=-self.dv,clevs=30,titlestr='',xlims=xlims,ylims=ylims,cmap='bone')
	
        
	plotage=self._calcage

        # Plot the particles at the first time step
	if not plotage:
	    h1 = plt.plot([],[],'y.',markersize=1.0)
	    h1 = h1[0]

	else:
	    h1 = plt.scatter(self.particles['X'],self.particles['Y'],s=1.0,c=self.particles['age'],vmin=0,vmax=agemax,edgecolors=None)
	    self.fig.colorbar(h1)
        
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
		if self._calcage:
		    self.particles.update({'age':np.zeros_like(x),'agemax':np.zeros_like(x)})
                
            self.advectParticles(self.time_track[ii],self.time_track_sec[ii])
	    if plotage:
		# Update the scatter object
		h1.set_offsets(np.vstack([self.particles['X'],self.particles['Y']]).T)
		h1.set_array(self.particles['age'])*agescale
		h1.set_edgecolors(h1.to_rgba(np.array(self.particles['age'])))    
	    else:
		h1.set_xdata(self.particles['X'])
		h1.set_ydata(self.particles['Y'])

            title.set_text(self.genTitle(ii))
            return (h1,title)
        
        self.anim = animation.FuncAnimation(self.fig, updateLocation, frames=len(self.time_track), interval=50, blit=True)
        
        if outfile==None:
            plt.show()
        else:
            self.saveanim(outfile)
    
    def animateNC(self,ncfile,plotage=False,agemax=7.0,xlims=None,ylims=None,outfile=None,**kwargs):
        """
        Animate the particles from a previously generated netcdf file
        """
        import matplotlib.animation as animation
        
	agescale = 1.0/86400.0
        # Set the xy limits
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        # Open the netcdf file
        nc = Dataset(ncfile,'r')
        
        # Read the time data into 
        t = nc.variables['tp']
        self.time_track = num2date(t[:],t.units)
        
        # Plot a map of the bathymetry
        self.fig = plt.figure(figsize=(10,8))
        ax=plt.gca()
        self.contourf(z=-self.dv,clevs=30,titlestr='',xlims=xlims,ylims=ylims,cmap='bone')
        
        # Plot the particles at the first time step
        if not plotage:
#	    h1 = plt.plot([],[],'y.',markersize=1.0)
	    h1 = plt.plot([],[],'.',**kwargs)
	    h1 = h1[0]

	else:
	    h1 = plt.scatter(nc.variables['xp'][0,:],nc.variables['yp'][0,:],s=1.0,c=nc.variables['age'][0,:],vmin=0,vmax=agemax,edgecolors=None)
	    self.fig.delaxes(self.fig.axes[1])
	    self.cb = self.fig.colorbar(h1)
	    #self.cb.on_mappable_changed(h1) # Updates the colobar

        
        title=ax.set_title("")
        
	def updateLocation(ii):
	    xp = nc.variables['xp'][:,ii]
	    yp = nc.variables['yp'][:,ii]
	    if plotage:
		# Update the scatter object
		h1.set_offsets(np.vstack([xp,yp]).T)
		age=nc.variables['age'][:,ii]*agescale
		h1.set_array(age)
		h1.set_edgecolors(h1.to_rgba(np.array(age)))    
	    else:
		h1.set_xdata(xp)
		h1.set_ydata(yp)

            title.set_text(self.genTitle(ii))
            return (h1,title)
        
        self.anim = animation.FuncAnimation(self.fig, updateLocation, frames=len(self.time_track), interval=50, blit=True)
        
        if outfile==None:
            plt.show()
        else:
            self.saveanim(outfile) 
        
        nc.close()

    def animate_xz(self,x,y,z,timeinfo,agepoly=None,agemax=86400.0,xlims=None,ylims=None,outfile=None):
        """
        Animate the particles on the fly using matplotlib plotting routines
        """
        import matplotlib.animation as animation
        
           
        self.__call__(x,y,z,timeinfo,agepoly=agepoly,runmodel=False)

	# Set the xy limits
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=[-self.dv.max(),1.0]
         
        self.fig = plt.figure(figsize=(10,8))
        ax=plt.gca()
	ax.set_xlim(xlims)
	ax.set_ylim(ylims)
	
	plotage=self._calcage

        # Plot the particles at the first time step
	if not plotage:
	    h1 = plt.plot([],[],'y.',markersize=1.0)
	    h1 = h1[0]

	else:
	    h1 = plt.scatter(self.particles['X'],self.particles['Z'],s=1.0,c=self.particles['age'],vmin=0,vmax=agemax,edgecolors=None)
	    self.fig.colorbar(h1)
        
        title=ax.set_title("")
        
        def init():
            h1.set_xdata(self.particles['X'])
            h1.set_ydata(self.particles['Z'])
            title=ax.set_title(self.genTitle(0))
            return (h1,title)

        def updateLocation(ii):
            if ii==0:
                # Re-initialise the particle location
                self.particles={'X':x,'Y':y,'Z':z}
		if self._calcage:
		    self.particles.update({'age':np.zeros_like(x),'agemax':np.zeros_like(x)})
                
            self.advectParticles(self.time_track[ii],self.time_track_sec[ii])
	    if plotage:
		# Update the scatter object
		h1.set_offsets(np.vstack([self.particles['X'],self.particles['Z']]).T)
		h1.set_array(self.particles['age'])
		h1.set_edgecolors(h1.to_rgba(np.array(self.particles['age'])))    
	    else:
		h1.set_xdata(self.particles['X'])
		h1.set_ydata(self.particles['Z'])

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
         
	# Activate particles for current time step
	self.activateParticles(tsec)

        # Advection step 
        if self.advect_method=='euler':
            self.euler(timenow,tsec)
            
        elif self.advect_method=='rk2':
            self.rk2(timenow,tsec)
            
        else:
            raise Exception, 'unknown advection scheme: %s. Must be "euler" or "rk2"'%self.advect_method

	# Reset inactive particles
	self.resetInactiveParticles()

	# Call the age calculation
	if self._calcage:
	   self.CalcAge()
            
        
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
        v = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.v,update=False)
        w = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.w,update=False)
        #u,v,w = self.timeInterpUVWxyz(tsec,self.particles['X'],self.particles['Y'],self.particles['Z'])
        
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
        v = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.v,update=False)
        w = self.UVWinterp(self.particles['X'],self.particles['Y'],self.particles['Z'],self.w,update=False)
        #u,v,w = self.timeInterpUVWxyz(tsec,self.particles['X'],self.particles['Y'],self.particles['Z'])

        x1 = self.particles['X'] + 0.5*self.dt*u
        y1 = self.particles['Y'] + 0.5*self.dt*v
        z1 = self.particles['Z'] + 0.5*self.dt*w
        
        # Check the vertical bounds of a particle 
        z1 = self.checkVerticalBounds(x1,y1,z1)
        
        # Update the currents again
        self.updateCurrents(timenow+timedelta(seconds=self.dt*0.5),tsec+self.dt*0.5)
        
        u = self.UVWinterp(x1,y1,z1,self.u)
        v = self.UVWinterp(x1,y1,z1,self.v,update=False)
        w = self.UVWinterp(x1,y1,z1,self.w,update=False)
        #u,v,w = self.timeInterpUVWxyz(tsec,x1,y1,x1)
        
        self.particles['X'] += u*self.dt
        self.particles['Y'] += v*self.dt
        self.particles['Z'] += w*self.dt
        
        # Check the vertical bounds of a particle again 
        self.particles['Z'] = self.checkVerticalBounds(self.particles['X'],self.particles['Y'],self.particles['Z'])

	# Check the horizontal coordinates
	#    This is done by default in the messh interpolation class
	#self.particles['X'],self.particles['Y'] = self.checkHorizBounds(self.particles['X'],self.particles['Y'])
        
    def checkHorizBounds(self,x,y):
    	"""
	NOT USED
	Moves particles outside of the horizontal bounds of the grid back
	to the closet cell centre.
	"""

	ind = self.UVWinterp.cellind==-1
	if not self.__dict__.has_key('kd'):
	    self.kd = spatial.cKDTree(np.vstack((self.xv,self.yv)).T)

	xy = np.vstack((x[ind],y[ind])).T
	dist,cell = self.kd.query(xy)
	x[ind]=self.xv[cell]
	y[ind]=self.yv[cell]

	return x, y

    def checkVerticalBounds(self,x,y,z):
        """
        Checks that particles are not above the surface or below the seabed
        
	(Artificially moves them to the surface or bed if they are).
        """
        SMALL = 0.001
	zbed = -self.dv
	#zbed = -self.z_w[self.Nk-1] # Set the bottom one layer above the seabed
        
        ## find the free surface and seabed at the particle locations
        #if not self.interp_method == 'mesh':
        #    eta_P = self.Hinterp(x,y,z,self.eta)
        #    h_P = self.Hinterp(x,y,z,zbed)
        #else:
        #    ind = self.UVWinterp.cellind
        #    mask=ind==-1
        #    ind[mask]=0.0
        #    eta_P = self.eta[ind]
        #    h_P = zbed[ind]
        #    eta_P[mask]=0.0
        #    h_P[mask]=0.0
        #
	#pdb.set_trace()
	#indtop = np.where(z>eta_P)
	#indbot = np.where(z<h_P)
	#z[indtop[0]] = eta_P[indtop[0]]-SMALL
	#z[indbot[0]] = h_P[indbot[0]]+SMALL

	if not self.interp_method == 'mesh':
            eta_P = self.Hinterp(x,y,z,self.eta)
            h_P = self.Hinterp(x,y,z,zbed)
        else:
            ind = self.UVWinterp.cellind
            eta_P = self.eta[ind]
            h_P = zbed[ind]
        
	indtop = np.where( operator.and_(z>eta_P, ind!=-1) )
	indbot = np.where( operator.and_(z<h_P, ind!=-1) )

	z[indtop[0]] = eta_P[indtop[0]]-SMALL
	z[indbot[0]] = h_P[indbot[0]]+SMALL


        #np.where(z > eta_P, eta_P-SMALL, z)
        #np.where(z < h_P, h_P+SMALL, z)
        
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
            if self.verbose:
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
        
    def timeInterpUVWxyz(self,tsec,x,y,z):
        """
	NOT USED AT THIS STAGE

        Temporally interpolate the currents  onto the particle time step, tsec.
        
        Linear interpolation used.
        """
        tindex=self.time_index
        
        t0 = self.time_sec[tindex-1]
        t1 = self.time_sec[tindex]
        dt = t1 - t0
        
        w1 = (tsec-t0)/dt
        w2 = 1.0-w1
        
        uT0 = self.UVWinterp(x,y,z,self.uT[:,0],update=True)
        vT0 = self.UVWinterp(x,y,z,self.vT[:,0],update=False) # Locations of the particles are not checked for updates
        wT0 = self.UVWinterp(x,y,z,self.wT[:,0],update=False)
        
        uT1 = self.UVWinterp(x,y,z,self.uT[:,1],update=False)
        vT1 = self.UVWinterp(x,y,z,self.vT[:,1],update=False)
        wT1 = self.UVWinterp(x,y,z,self.wT[:,1],update=False)
        
        u = uT0*w1 + uT1*w2 
        v = vT0*w1 + vT1*w2 
        w = wT0*w1 + wT1*w2 
        
        return u,v,w
        
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
        
        self.time_track = othertime.TimeVector(timeinfo[0],timeinfo[1],timeinfo[2],timeformat ='%Y%m%d.%H%M%S')
        
        self.time_track_sec = othertime.SecondsSince(self.time_track)
        self.time_sec = othertime.SecondsSince(self.time)
        
        self.time_index = -9999

    def initParticleTime(self,tstart):
    	"""
	Initialise the particle start time
	"""
	if tstart==None:
	    self.particles['tstart'] = self.time_track_sec[0]*np.ones_like(self.particles['X'])
	else:
	    self.particles['tstart'] = tstart

	# Set all particles as inactive to start
	self.particles['isActive'] = np.zeros_like(self.particles['X'])

    def activateParticles(self,tsec):
	"""
	Activate the particles that start past the present time

	Note that inactive particles are still advected - their positions
	are simply reset after each time step
	"""
	ind = self.particles['tstart'] <= tsec

	self.particles['isActive'][ind]=True
	print '\t%d Active particles'%(np.sum(self.particles['isActive']))

    def resetInactiveParticles(self):
	"""
	Reset the positions of inactive particles to x0. 
	"""
	ind = self.particles['isActive']==False
	self.particles['X'][ind]=self.particles['X0'][ind]
	self.particles['Y'][ind]=self.particles['Y0'][ind]
	self.particles['Z'][ind]=self.particles['Z0'][ind]

        
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

    def CalcAge(self):
        """     
	Calculate the age of a particle inside of the age polygon
	"""
	#print '\t\tCalculating the particle age...'
	inpoly = nxutils.points_inside_poly(np.vstack((self.particles['X'],self.particles['Y'])).T,self.agepoly)

	self.particles['age'][inpoly] = self.particles['age'][inpoly] + self.dt
	self.particles['age'][inpoly==False]=0.0

	# Update the agemax attribute
	self.particles['agemax'] = np.max([self.particles['age'],self.particles['agemax']],axis=0)



    def initParticleNC(self,outfile,Np,age=False):
        """
        Export the grid variables to a netcdf file
        """
	import os

        if self.verbose:
            print '\nInitialising particle netcdf file: %s...\n'%outfile
            
	# Global Attributes
        nc = Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
        nc.Description = 'Particle trajectory file'
        nc.Author = os.getenv('USER')
        nc.Created = datetime.now().isoformat()

	#tseas = self.time_sec[1] - self.time_sec[0]
	#nsteps = np.floor(tseas/self.dt)
	#nc.nsteps = '%f (number of linear interpolation steps in time between model outputs)'%nsteps
	#nc.tseas = '%d (Time step (seconds) between model outputs'%tseas
	#nc.dt = '%f (Particle model time steps [seconds])'%self.dt
	nc.dataset_location = '%s'%self.ncfile

	# Dimensions
        nc.createDimension('ntrac', Np)
        nc.createDimension('nt', 0) # Unlimited
        
        # Create variables
        def create_nc_var( name, dimensions, attdict, dtype='f8'):
	    tmp=nc.createVariable(name, dtype, dimensions)
            for aa in attdict.keys():
                tmp.setncattr(aa,attdict[aa])
          
        create_nc_var('tp',('nt'),{'units':'seconds since 1990-01-01 00:00:00','long_name':"time at drifter locations"},dtype='f8')
        create_nc_var('xp',('ntrac','nt'),{'units':'m','long_name':"Easting coordinate of drifter",'time':'tp'},dtype='f8')
        create_nc_var('yp',('ntrac','nt'),{'units':'m','long_name':"Northing coordinate of drifter",'time':'tp'},dtype='f8')
        create_nc_var('zp',('ntrac','nt'),{'units':'m','long_name':"vertical position of drifter (negative is downward from surface)",'time':'tp'},dtype='f8')
	if age:
	    create_nc_var('age',('ntrac','nt'),{'units':'seconds','long_name':"Particle age",'time':'tp'},dtype='f8')
	    create_nc_var('agemax',('ntrac','nt'),{'units':'seconds','long_name':"Maximum particle age",'time':'tp'},dtype='f8')

        nc.close()
    
    def writeParticleNC(self,outfile,x,y,z,t,tstep,age=None,agemax=None):
        """
        Writes the particle locations at the output time step, 'tstep'
        """
        if self.verbose:
            print 'Writing netcdf output at tstep: %d...\n'%tstep
            
        nc = Dataset(outfile, 'a')
        
        nc.variables['tp'][tstep]=t
        nc.variables['xp'][:,tstep]=x
        nc.variables['yp'][:,tstep]=y
        nc.variables['zp'][:,tstep]=z

	if not age==None:
	    nc.variables['age'][:,tstep]=age
	if not agemax==None:
	    nc.variables['agemax'][:,tstep]=agemax

        nc.close()
#################
# "Post-processing 
#################
    def plotTrackNC(self,ncfile,xlims=None,ylims=None,**kwargs):	
	"""
	Plots the tracks of all particles from a particle netcdf file
	"""
	# Set the xy limits
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
	
	# Load all of the data
	nc = Dataset(ncfile,'r')
	xp = nc.variables['xp'][:]
	yp = nc.variables['yp'][:]
	nc.close()

	# Plot a map of the bathymetry
        self.fig = plt.figure(figsize=(10,8))
        ax=plt.gca()
        self.contourf(z=-self.dv,clevs=30,titlestr='',xlims=xlims,ylims=ylims,cmap='bone')

	# plot the tracks
        plt.plot(xp[:,0],yp[:,0],'o',color='g',markersize=3,label='_nolegend_',alpha=0.4)
	# Plot tracks
	plt.plot(xp.T,yp.T,'-',color='grey',linewidth=.2)

	plt.title('Particle Trajectory\n(file name: %s)'%ncfile)

#	plt.show()

    def plotAgeMax(self,ncfile,dx,dy,vmax=7,xlims=None,ylims=None,**kwargs):
    	"""
	Filled contour plot of the maximum age of particle based on their initial locations
	"""
	scalefac = 1/86400.
	# Set the xy limits
	if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims


	xp,yp,agemax = self.readAgeMax(ncfile)
	
	grd = RegGrid(xlims,ylims,dx,dy)

	agegrid = grd.griddata(xp,yp,agemax)*scalefac

	fig=plt.gcf()
	self.contourf(z=-self.dv,clevs=20,titlestr='',xlims=xlims,ylims=ylims,\
	    filled=False,colors='k',linestyles='solid',linewidths=0.2)
	plt.pcolor(grd.X,grd.Y,agegrid,vmax=vmax,**kwargs)
	plt.colorbar()
	plt.title('Particle Age [days]\n(file name: %s)'%ncfile)

    def readAgeMax(self,ncfile):
	"""
	Reads the maximum age from the last time step in a netcdf file
	"""

	# Load all of the data
	nc = Dataset(ncfile,'r')
	xp = nc.variables['xp'][:,0]
	yp = nc.variables['yp'][:,0]
	try:
	    # Load the age from the last time step
	    agemax = nc.variables['agemax'][:,-1]
	except:
	    raise Exception, ' "agemax" variable not present in file: %s'%ncfile
	nc.close()
	
	return xp, yp, agemax


    def calcAgeMaxProb(self,ncfiles,dx,dy,exceedance_thresh,plot=True,xlims=None,ylims=None,**kwargs):
    	"""
	Calculate the probability of the age exceeding a threshold time from
	a series of particle tracking runs. 
	"""
	scalefac = 1/86400.
	# Set the xy limits
	if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims

	# Create the probability grid
	grd = RegGrid(xlims,ylims,dx,dy)
	nruns = len(ncfiles)
	prob = np.zeros((grd.ny,grd.nx))

	# Loop through the files and add one to the probability matrix
	# where the age exceeds the threshold
	for ncfile in ncfiles:
	    print 'Reading file: %s'%ncfile
	    xp,yp,agemax = self.readAgeMax(ncfile)
	
	    agegrid = grd.griddata(xp,yp,agemax)

	    ind = agegrid >= exceedance_thresh

	    prob[ind] = prob[ind]+1

	# Return the probability as a percentage
	prob = prob/nruns*100.0

	if plot:
	    fig=plt.gcf()
	    self.contourf(z=-self.dv,clevs=20,titlestr='',xlims=xlims,ylims=ylims,\
		filled=False,colors='k',linestyles='solid',linewidths=0.2)
	    plt.pcolor(grd.X,grd.Y,prob,vmax=100,cmap=plt.cm.Spectral_r)
	    plt.colorbar()
	    plt.title('Probability (%%) of particle age exceeding %3.1f days'%(exceedance_thresh*scalefac))


class RegGrid(object):
    """
    Class for handling regular grids
    """
    def __init__(self,xlims,ylims,dx,dy):
        self.ox = xlims[0]
        self.oy = ylims[0]
        self.dx = dx
        self.dy = dy
	self.xlims=xlims
	self.ylims=ylims

	#  Construct a 2D mesh of particles
        x = np.arange(xlims[0],xlims[1],dx)
	y = np.arange(ylims[0],ylims[1],dy)

	self.X, self.Y = np.meshgrid(x,y)

	shp = self.X.shape
	self.ny = shp[0]
	self.nx = shp[1]
	
    def returnij(self,x,y):
        """
	Returns the i,j (cols,rows) of the points in i,j

	Returns NaNs for points that are out of bounds
	"""

	i = np.ceil( (x-self.ox)/self.dx )
	j = np.ceil( (y-self.oy)/self.dy )

	# Check the bounds
	ind = x > self.xlims[1]
	ind = i >= self.nx
	i[ind]=np.nan
	j[ind]=np.nan

	ind = x < self.ox
	ind = i < 0
	i[ind]=np.nan
	j[ind]=np.nan

	ind = y > self.ylims[1]
	ind = j >= self.ny
	i[ind]=np.nan
	j[ind]=np.nan

	ind = y < self.oy
	ind = j < 0
    	i[ind]=np.nan
	j[ind]=np.nan

	return i,j

    def griddata(self,x,y,z):
        """
        Grids data in the vectors x, y and z.

        Uses sparse matrices and therefore duplicate entries are averaged
        """
        # Get the coordinates
        i,j = self.returnij(x,y)
        ind = np.isfinite(i)

        # Build two sparse matrices - one for the data and one to count 
        # the number of entries
        ctr = np.ones(z.shape)

        data = sparse.coo_matrix( (z[ind],(j[ind],i[ind])),shape=(self.ny,self.nx) )
        count = sparse.coo_matrix( (ctr[ind],(j[ind],i[ind])),shape=(self.ny,self.nx) )

        data = data/count

        return np.array(data.todense())

class interp3Dmesh(TriSearch,Grid):
    """
    3D interpolation class for an unstructured grid
    
    Uses knowledge of the mesh to efficiently update the velocities based on 
    the particle location.
    """


    def __init__(self,x,y,z,cells,mask,method='nearest',grdfile=None):
        
	self.method=method
        # Initialise the trisearch array
        TriSearch.__init__(self,x,y,cells,force_inside=True)

	if self.method == 'linear':
	    Grid.__init__(self,grdfile)
	    self.datatmp = np.zeros(mask.shape,dtype=np.double)
        
        self.z = np.sort(z)
	self.z[-1]=10.0 # Set the surface layer to large
        self.Nkmax = z.size-1
        
        self.mask3d = mask
        
        self.maskindex = np.zeros(self.mask3d.shape,dtype=np.int32)
        rr=0
        for ii in range(self.mask3d.shape[0]):
            for jj in range(self.mask3d.shape[1]):
                if self.mask3d[ii,jj]:
                    self.maskindex[ii,jj]=rr
                    rr+=1
                  
    def __call__(self,X,Y,Z,data,update=True):
        
        if update:
            # The Update the cell index using TriSearch class
            if not self.__dict__.has_key('cellind'):
                #print ' Finding initial particle index...'
                TriSearch.__call__(self,X,Y)
            else:
                if np.sum(np.abs(X-self.xpt))>0+1e-8:
                    #print ' updating location index...'
                    self.updatexy(X,Y)
                        
        # Find the k-index
        kind=self.z.searchsorted(Z)
        kind = self.Nkmax - kind 
        
        kind[kind>=self.Nkmax-1] = self.Nkmax-1

        ind = self.maskindex[kind,self.cellind]
        maskpts = self.mask3d[kind,self.cellind]
        
        ind[maskpts == False] = 0
        
        # Return the nearest data point (... for now)
	if self.method == 'nearest':
	    dataout = data[ind]
	if self.method == 'linear':
	    dataout = self.lininterp(X,Y,Z,data,kind)
                
        # Mask bogey points
        dataout[maskpts==False]=0.0
        dataout[self.cellind==-1]=0.0
                
        return dataout

    def lininterp(self,X,Y,Z,data,k):
    	"""
	Wrapper for linear interpolation

	interpLinear method in sunpy needs to be called layer by layer
	"""
	# Put the input data back into a 2D array (Nz,Nc)
	self.datatmp[self.mask3d]=data

	# Loop through layer by layer
	dataout = np.zeros_like(X)
	for kk in range(0,self.Nkmax):
	    ind = operator.and_(k == kk,self.cellind!=-1)
	    if any(ind):
		dataout[ind] = self.interpLinear(self.datatmp[kk,:],X[ind],Y[ind],self.cellind[ind],k=kk)

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

def GridParticles(grdfile,dx,dy,nz,xypoly=None,splitvec=1):
    """
    Returns the locations of particles on a regular grid inside of suntans grid

    Inputs:
    	grdfile - netcdf filename containing the suntans grid
	dx,dy - resolution in x and y component respectively.
	nz - number of particles in the vertical. Particles are arranged in sigma layers.
	xypoly - [optional] coordinates of an additional bounding polygon [Nx2 array]
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

    # Check which particles are also inside of the polygon
    if not xypoly == None:
	inpoly = nxutils.points_inside_poly(np.vstack((X,Y)).T,xypoly)
	mask = operator.and_(mask,inpoly)

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
    
    xout = xout.ravel()
    yout = yout.ravel()
    zout = zout.ravel()

    # Rearrange the vectors to avoid clustering (this helps even the MPI workload)
    #xout_split=[]
    #yout_split=[]
    #zout_split=[]
    #for start in range(splitvec):
    #    xout_split.append(xout[start::splitvec])
    #    yout_split.append(yout[start::splitvec])
    #    zout_split.append(zout[start::splitvec])

    #xout = np.hstack(xout_split)
    #yout = np.hstack(yout_split)
    #zout = np.hstack(zout_split)
    

    return xout, yout, zout
    
