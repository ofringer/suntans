# -*- coding: utf-8 -*-
"""
Tools for interfacing SUNTANS with the tvtk libraries

Created on Tue Dec 18 16:12:27 2012

@author: mrayson
"""

import numpy as np
from sunpy import Spatial
from tvtk.api import tvtk
from mayavi import mlab

import pdb 

class SunTvtk(Spatial):
    """
    Wrapper for suntans grid on top of tvtk unstructured grid object
    
    The main object is stored in the 'ug' attribute
    
    The scalar data is stored in the 'data' attribute
    
    # TODO
        - Interpolation routine
    """
    
    is3D = False
    vector_overlay=False
    zscale = 500.0
    clim = None
    kstart=0 # Starting klayer - set > 0 to ignore top 'kstart' cells
    offscreen=False
    
    def __init__(self,infile,**kwargs):
        
        self.__dict__.update(kwargs)

        Spatial.__init__(self,infile,**kwargs)
            
        #
        if self.is3D:
            self.klayer=np.arange(self.kstart,self.Nkmax)
            self.Nkmax-=self.kstart
            self.Nk -= self.kstart
            self.data = np.zeros((self.Nc,self.Nkmax))
            self.data = np.ravel(self.data)
            
            self.initTvtk3D()
        else:
            # Initialize the 2D object
            self.data = np.zeros((self.Nc,))
            self.returnPoints()
            self.initTvtk2D()

	if self.offscreen:
	    print 'Using offscreen rendering.'
	    mlab.options.offscreen=True

    def initTvtk2D(self):
        """
        Initialise the actual 2 dimensional tvtk object
        """
        
        tri_type = tvtk.Triangle().cell_type
        
        self.ug = tvtk.UnstructuredGrid(points=self.points)
        self.ug.set_cells(tri_type, self.cells)
    
        self.ug.cell_data.scalars = self.data
        self.ug.cell_data.scalars.name = 'suntans_scalar'
        
    def returnPoints(self):
        """
        XYZ location of the points
        """        
        
        self.points = np.column_stack((self.xp,self.yp,0.0*self.xp))
    
    def initTvtk3Dold(self):
        """
        Constructs the 3d cell unstructured grid object
        
        This method includes all vertical grid layers
        
        # See this example
        https://github.com/enthought/mayavi/blob/master/examples/mayavi/advanced_visualization/unstructured_grid.py
        
        """
        nc = self.Nc
        nz = self.Nkmax+1
        nv = len(self.xp)
                
        # Create the index to vertices (nodes) and the coordinates of the vertices (verts) arrays
        nodes = np.zeros((nc*(nz-1),6))
        pt1=0
        for k in range(1,nz):
            nodes[pt1:nc*k,0] = self.cells[:,0]+(k-1)*nv
            nodes[pt1:nc*k,1] = self.cells[:,1]+(k-1)*nv
            nodes[pt1:nc*k,2] = self.cells[:,2]+(k-1)*nv
            nodes[pt1:nc*k,3] = self.cells[:,0]+k*nv
            nodes[pt1:nc*k,4] = self.cells[:,1]+k*nv
            nodes[pt1:nc*k,5] = self.cells[:,2]+k*nv
            pt1+=nc
            
        verts = np.zeros((nv*nz,3))
        pv1 = 0
        for k in range(0,nz):
            verts[pv1:pv1+nv,0] = self.xp
            verts[pv1:pv1+nv,1] = self.yp
            verts[pv1:pv1+nv,2] = -self.z_w[k] * self.zscale
            pv1 += nv
            
        wedge_type = tvtk.Wedge().cell_type
        self.ug = tvtk.UnstructuredGrid(points=verts)
        self.ug.set_cells(wedge_type, nodes)
        
        self.ug.cell_data.scalars = self.data
        self.ug.cell_data.scalars.name = 'suntans_scalar' 
    
    def initTvtk3D(self):
        """
        Constructs the 3d cell unstructured grid object
        
        This method includes only active vertical grid layers
        
        """
        nc = self.Nc
        nz = self.Nkmax+1
        nv = len(self.xp)
        
        self.returnMask3D()
        self.nActive = np.sum(self.mask3D) # Total number of active cells
        
        nodes = np.zeros((self.nActive,6))
        pt1=0
        for k in range(1,nz):
            masklayer = self.mask3D[k-1,:]
            nc = np.sum(masklayer)
            pt2 = pt1+nc            
            nodes[pt1:pt2,0] = self.cells[masklayer,0]+(k-1)*nv
            nodes[pt1:pt2,1] = self.cells[masklayer,1]+(k-1)*nv
            nodes[pt1:pt2,2] = self.cells[masklayer,2]+(k-1)*nv
            nodes[pt1:pt2,3] = self.cells[masklayer,0]+k*nv
            nodes[pt1:pt2,4] = self.cells[masklayer,1]+k*nv
            nodes[pt1:pt2,5] = self.cells[masklayer,2]+k*nv
            pt1=pt2
            #print k, nc
            
        self.verts = np.zeros((nv*nz,3))
        pv1 = 0
        for k in range(0,nz):
            self.verts[pv1:pv1+nv,0] = self.xp
            self.verts[pv1:pv1+nv,1] = self.yp
            self.verts[pv1:pv1+nv,2] = -self.z_w[k] * self.zscale
            pv1 += nv
            
        wedge_type = tvtk.Wedge().cell_type
        self.ug = tvtk.UnstructuredGrid(points=self.verts)
        self.ug.set_cells(wedge_type, nodes)
        
        self.ug.cell_data.scalars = self.data
        self.ug.cell_data.scalars.name = 'suntans_scalar' 
        
    
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
            
    def newscene(self,size=(800,700)):
        """
        Creates a new scene
        """
        
        self.fig=mlab.figure(bgcolor=(0.,0.,0.),size=size)
        self.fig.scene.z_plus_view()
        #mlab.view(0,0)
        #self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
        
    def colorbar(self):
        """
        Adds a colorbar for the object in 'h'
        """
        self.cb = mlab.colorbar(object=self.h,orientation='vertical')
    
    def surface(self,clim=None,**kwargs):
        """
        Surface plot of the scalar in the 'data' attribute with mayavi
        
        Works on the 2D and 3D data
        """

        if clim==None:
            clim = [self.data.min(), self.data.max()]
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
        
        src = mlab.pipeline.add_dataset(self.ug)
        self.h=mlab.pipeline.surface(src,vmin=clim[0],vmax=clim[1],**kwargs)
        
        # Add a colorbar if the isn't one
        if not self.__dict__.has_key('cb'):
            self.colorbar() 
            
        # Add a title if there isn't one
        if not self.__dict__.has_key('title'):
            self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
        
    def contour(self,vv=[10],clim=None,**kwargs):
        """
        Filled contour plot of scalar data
        """
        
        if clim==None:
            clim = [self.data.min(), self.data.max()]
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
        
        # Need to set use point (vertex) data
        src = mlab.pipeline.cell_to_point_data(self.ug)
        
        # Add the contour_surface module to the scene
        self.h=mlab.pipeline.contour_surface(src,contours=vv,line_width=1.0,vmax=clim[1],vmin=clim[0],**kwargs)
        self.h.contour.filled_contours=True # This is the trick to fill the contours
        
        # Add a colorbar if the isn't one
        if not self.__dict__.has_key('cb'):
            self.colorbar() 
            
        # Add a title if there isn't one
        if not self.__dict__.has_key('title'):
            self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
    
    def isosurface(self,vv=[4.0],clim=None,**kwargs):
        """
        3D isosurfaces of scalar data
        """
        if self.clim==None:
            self.clim = [self.data.min(), self.data.max()]
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
        
        # Convert the cell centred data into a scene source
        # Need to set use point (vertex) data
        src = mlab.pipeline.cell_to_point_data(self.ug)
        
        # Add the iso_surface module to the scene
        self.h=mlab.pipeline.iso_surface(src,contours=vv,line_width=1.0,**kwargs)
        
        # Add a colorbar if the isn't one
        if not self.__dict__.has_key('cb'):
            self.colorbar() 
            
        # Add a title if there isn't one
        if not self.__dict__.has_key('title'):
            self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
    
    def volume(self,clim=None,**kwargs):
        """
        3D volumetric plot of scalar data
        
        **Warning** This is really slow and memory hungry!
        """
        if self.clim==None:
            self.clim = [self.data.min(), self.data.max()]
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
        
        # Convert the cell centred data into a scene source
        # Need to set use point (vertex) data
        src = mlab.pipeline.cell_to_point_data(self.ug)
        
        # Add the volume module to the scene
        self.h=mlab.pipeline.volume(src,**kwargs)
        
        # Add a colorbar if the isn't one
        if not self.__dict__.has_key('cb'):
            self.colorbar() 
            
        # Add a title if there isn't one
        if not self.__dict__.has_key('title'):
            self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
            
    def sliceplane(self,plane_orientation='y_axes',**kwargs):
        """
        Applies the image plane widget to the dataset
        """        
            
        if self.clim==None:
            self.clim = [self.data.min(), self.data.max()]
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()

        src = mlab.pipeline.cell_to_point_data(self.ug)
        
        self.h=mlab.pipeline.scalar_cut_plane(src,plane_orientation=plane_orientation,view_controls=True,line_width=0.5)
        
       # Add a colorbar if the isn't one
        if not self.__dict__.has_key('cb'):
            self.colorbar() 
            
        # Add a title if there isn't one
        if not self.__dict__.has_key('title'):
            self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15) 
            
    def vector(self,color=(0,0,0),subsample=1,scale=1e-3,line_width=1.0,**kwargs):
        """
        Overlay vectors onto the current scence
        """
        
        self.vector_overlay=True
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
        
        self.loadVectorVTK()
        src = mlab.pipeline.cell_to_point_data(self.ug)
        self.h=mlab.pipeline.vectors(src,color=color,mask_points=subsample,scale_factor=1./scale,scale_mode='vector',line_width=line_width,**kwargs)
        
    def streamline(self,**kwargs):
        
        """
        Plots streamlines
        
        
        See here:
            http://mindseye.no/2010/09/25/using-mayavi-to-visualize-electric-fields/
        """
        
        self.vector_overlay=True
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
        
        self.loadVectorVTK()
        src = mlab.pipeline.add_dataset(self.ug)
        #src = mlab.pipeline.cell_to_point_data(self.ug)
        self.h=mlab.pipeline.streamline(src,**kwargs)
        
        self.h.stream_tracer.initial_integration_step = 0.1
        self.h.stream_tracer.maximum_propagation=10000.0
        self.h.stream_tracer.integration_direction = 'both'
        
    def streammesh(self,nx=30,ny=30,**kwargs):
        """
        Plots streamlines originating from points on a mesh
        """
        
        
        self.vector_overlay=True
        
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
        
        self.loadVectorVTK()
        src = mlab.pipeline.add_dataset(self.ug)
#        src = mlab.pipeline.cell_to_point_data(self.ug) # This colours the streamlines by the scalar field
             
        X = np.linspace(self.xv.min(),self.xv.max(),nx)
        y0 =self.yv.min()
        y1 = self.yv.max()
        streams=[]
        for x in X:
            stream=mlab.pipeline.streamline(src,seedtype='line',**kwargs)
            stream.stream_tracer.initial_integration_step = 0.1
            stream.stream_tracer.maximum_propagation=10000.0
            stream.stream_tracer.integration_direction = 'both'
            stream.seed.widget.point1 = [x,y1,0.1]
            stream.seed.widget.point2 = [x,y0,0.1]
            stream.seed.widget.resolution = ny
            streams.append(stream)

        # Go back through and turn them off
        for stream in streams:
            stream.seed.widget.enabled = False
         
        return streams


    # Plane method  # Won't work on 2D data
#        stream=mlab.pipeline.streamline(src,seedtype='plane')
#        stream.seed.widget.point1 = [self.xv.min(),self.yv.min(),0.0]
#        stream.seed.widget.point2 = [self.xv.max(),self.yv.max(),0.0]
#        stream.seed.widget.normal_to_z_axis=True
#        pdb.set_trace()
        
        
        # Point method - very slow way of creating mesh
#        X,Y = np.meshgrid(np.linspace(self.xv.min(),self.xv.max(),nx),np.linspace(self.yv.min(),self.yv.max(),ny))
#        X=X.ravel()
#        Y=Y.ravel()
#        Z=Y*0
#        stream=[]        
#        for x,y,z in zip(X,Y,Z):
#            h=mlab.pipeline.streamline(src,seedtype='point',**kwargs)
#            #pdb.set_trace()
#            h.seed.widget.position = [x,y,z] # set the stream widget to the same position as the charges
#            h.seed.widget.enabled = False
#            h.stream_tracer.initial_integration_step = 0.1
#            h.stream_tracer.maximum_propagation=10000.0
#            h.stream_tracer.integration_direction = 'both'
#            stream.append(h)
            
                    
        
    def animate(self,tstep=None):
        """
        Animates the current scene through all time steps (Interactive)
        """
        # Load all the time steps into memory (this can get large)
        #self.tstep=np.arange(0,len(self.time))
        #nt = len(self.tstep)
        #Spatial.loadData(self)
        
        # Load one time step at a time into memory (slower run time)
        if tstep==None:
            tstep=np.arange(0,len(self.time))
        nt = len(tstep)
        
        @mlab.animate
        def anim():
            ii=-1
            while 1:
                if ii<nt-1:
                    ii+=1
                else:
                    ii=0
                
                # Load all time steps
                #data=self.data[ii,:,:]
                #data=np.ravel(data[self.mask3D])
                #self.ug.cell_data.scalars = data
                #self.ug.cell_data.scalars.name = 'suntans_scalar'
                  
                # Loading one time step at a time
                if self.vector_overlay==True:
                    self.loadVectorVTK()
                    
                self.tstep = tstep[ii]
                self.loadData()
                self.ug.cell_data.scalars = self.data
                self.ug.cell_data.scalars.name = 'suntans_scalar'
                
                
                    
                titlestr=self._Spatial__genTitle(tt=ii)
                self.title.text=titlestr
                
                self.fig.scene.render()
                yield
        
        anim() # Starts the animation.
    
    def saveanimation(self,outfile,tstep=None,frate=10):
        """
        Saves an animation of the current scene to a file (non-interative)
        
        Saves a sequence of images and uses ffmpeg to convert to a video format.
        
        (I can't get it to pipe directly to ffmpeg as is done by matplotlib)
        
        See this:
            http://stackoverflow.com/questions/4092927/generating-movie-from-python-without-saving-individual-frames-to-files
        and this...
            http://stackoverflow.com/questions/13163106/ffmpeg-converting-image-sequence-to-video-results-in-blank-video
        """
#        
        import os
        # This works for mp4 and mov on Windows...
        #cmdstring ='ffmpeg -r %d -i ./.tmpanim%%04d.png -y -loglevel quiet -c:v libx264 -crf 23 -pix_fmt yuv420p %s'%(frate,outfile)
	# Linux command
        cmdstring ='ffmpeg -f image2 -r %d -y -loglevel quiet -b:v 7200k -i ./.tmpanim%%04d.png %s'%(frate,outfile)

        # Load one time step at a time into memory (slower run time)
        if tstep==None:
            tstep=np.arange(0,len(self.time))
        nt = len(tstep)
        
        png_list=[]
        for ii in range(nt):

            # Loading one time step at a time
            self.tstep = tstep[ii]
            
            if self.vector_overlay==True:
                self.loadVectorVTK()
                    
            self.loadData()
            self.ug.cell_data.scalars = self.data
            self.ug.cell_data.scalars.name = 'suntans_scalar'
                
            titlestr=self._Spatial__genTitle(tt=ii)
            self.title.text=titlestr
            
            self.fig.scene.render()
            
            # This bit saves each img
            outimg='./.tmpanim%04d.png'%ii
            png_list.append(outimg)
#            self.fig.scene.save_png(outimg)
            mlab.savefig(outimg,figure=self.fig)

            print 'Saving image %s of %d...'%(outimg,nt)
            
        # Call ffmpeg within python
        os.system(cmdstring)
        
        print '####\n Animation saved to: \n %s\n####' % outfile
        # Delete the images
        print 'Cleaning up temporary images...'
        for ff in png_list:
            os.remove(ff)
        print 'Complete.'
            
#    def savefig(self,outfile):
#        """
#        Saves the current scene to an image
#        """
#        self.fig.scene.save(outfile)
#            
#        print 'SUNTANS image saved to file:%s'%outfile
        
    def plotbathy3d(self,clims=None,**kwargs):
        """
        Adds 3D plot of the bathymetry to the current scene
        """
        # Create a new scene if there isn't one
        if not self.__dict__.has_key('fig'):
            self.newscene()
            
        depth = -self.dv
        if self.clim==None:
            clim = [depth.min(), depth.max()]
            
        #  Create an unstructured grid object to interpolate cells onto points
        points = np.column_stack((self.xp,self.yp,0.0*self.xp))
        tri_type = tvtk.Triangle().cell_type
        
        ug = tvtk.UnstructuredGrid(points=points)
        ug.set_cells(tri_type, self.cells)
    
        ug.cell_data.scalars = depth
        ug.cell_data.scalars.name = 'suntans_depth'
        
        # Interpolate the cell data onto the points
        F = mlab.pipeline.cell_to_point_data(ug)
        dp = mlab.pipeline.probe_data(F,self.xp,self.yp,0.0*self.xp)
        
        # Now set up a new object with the 3D points
        points = np.column_stack((self.xp,self.yp,dp*self.zscale))
        ug = tvtk.UnstructuredGrid(points=points)
        ug.set_cells(tri_type, self.cells)
        ug.cell_data.scalars = depth
        ug.cell_data.scalars.name = 'suntans_depth'
        
        # Plot as a 3D surface
        src = mlab.pipeline.add_dataset(ug)
        h=mlab.pipeline.surface(src,vmin=clim[0],vmax=clim[1],**kwargs)
        return h
        
    def loadData(self):
        """
        Overloaded loadData function - updates the unstructured grid object
        """
        Spatial.loadData(self)
        if self.is3D:
            self.data=np.ravel(self.data[self.mask3D])
        else:
            self.data=np.ravel(self.data)
            
        self.ug.cell_data.scalars = self.data
        self.ug.cell_data.scalars.name = 'suntans_scalar' 
        self.ug.modified()      
        
    def loadVectorVTK(self):
        """
        Loads vector data into the unstructured grid object
        """
        u,v,w = self.getVector() # This seems to do the masking (sub-classing??)
        
       
        if self.is3D==False:
            # Point data
            up = self.cell2node(u)
            vp = self.cell2node(v)
            wp = self.cell2node(w)
            wp *= 0.0 # Zero w for 2-D plots
        
            velocity = np.array((up,vp,wp)).T
            self.ug.point_data.vectors =  velocity
            self.ug.point_data.vectors.name = 'suntans_vector' 

        else: # 3D
            velocity = np.array((u,v,w)).T
            self.ug.cell_data.vectors =  velocity
            self.ug.cell_data.vectors.name = 'suntans_vector' 
        
        self.ug.modified()  
        
    def __setitem__(self,key,value):
        """
        Updates the unstructured grid object, ug, when the 'data' key is updated
        
        This is required when e.g., a variable is modified and re-inserted into the object
        """
        if key == "data":
            self.data=value
            self.ug.cell_data.scalars = value
            self.ug.cell_data.scalars.name = 'suntans_scalar' 
            self.ug.modified()
        elif key == 'is3D':
            self.is3D=value
            if self.is3D==True:
                self.klayer=np.arange(0,self.Nkmax)
                
########                      
# Testing stuff here

#ncfile = 'C:/Projects/GOMGalveston/MODELLING/GalvestonCoarse/rundata/GalvCoarse_TidesRivers3D*nc'
#ncfile = 'E:/Projects/GOMGalveston/MODELLING/SCENARIOS/CoarseHarmonicTideRivers/GalvCoarse_TidesRivers3D*nc'
#ncfile = 'C:\\Projects\\GOMGalveston\\MODELLING\\GalvestonMedium\\rundata\\GalvMedium_TidesRivers_Aug_049.nc'

#sunvtk = SunTvtk(ncfile,variable='salt',tstep=1,is3D=True,zscale=500.,kstart=1)
#sunvtk = SunTvtk(ncfile,variable='w',tstep=100,is3D=False,klayer=[-1])
#sunvtk.loadData()

#sunvtk.plotbathy3d(colormap='bone')
#sunvtk.isosurface(vv=[4.0,12.0,20.0,28.0],transparent=False,opacity=0.7)
#sunvtk.surface(representation='wireframe',opacity=0.3)
#sunvtk.surface()
#sunvtk.contour(vv=[4.0,8.0,12.0,16.0,20.0,24.0,28.0,32.0, 33.0],opacity=1.0)
#sunvtk.sliceplane()
#sunvtk.volume()

#sunvtk.animate()
#
#sunvtk.saveanimation('C:\Projects\GOMGalveston\MOVIES\CoarsePlume2D.mov')

