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
    """
    
    is3D = False
    zscale = 500.0
    
    def __init__(self,infile,**kwargs):
        
        self.__dict__.update(kwargs)

        Spatial.__init__(self,infile,**kwargs)
            
        #
        if self.is3D:
            self.klayer=np.arange(0,self.Nkmax)
            self.data = np.zeros((self.Nc,self.Nkmax))
            self.data = np.ravel(self.data)
            
            self.initTvtk3D()
        else:
            # Initialize the 2D object
            self.data = np.zeros((self.Nc,))
            self.returnPoints()
            self.initTvtk2D()

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
            
    def plot(self,clim=None,**kwargs):
        """
        Plots the scalar in the 'data' attribute with mayavi
        """

        if clim==None:
            clim = [self.data.min(), self.data.max()]
            
        self.fig=mlab.gcf()
        self.fig.scene.background = (0.,0.,0.)
        
        src = mlab.pipeline.add_dataset(self.ug)
        h=mlab.pipeline.surface(src,vmin=clim[0],vmax=clim[1],**kwargs)
        mlab.colorbar(object=h,orientation='vertical')
        mlab.view(0,0)
        self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
    
    def contour(self,vv=[20],clim=None,**kwargs):
        """
        Filled contour plot of scalar data
        """
            
        if clim==None:
            clim = [self.data.min(), self.data.max()]
            
        self.fig=mlab.gcf()
        self.fig.scene.background = (0.,0.,0.)
        
        # Convert the cell centred data into a scene source
        # Need to set use point (vertex) data
        src = mlab.pipeline.cell_to_point_data(self.ug)
        
        # Add the contour_surface module to the scene
        h=mlab.pipeline.contour_surface(src,contours=vv,line_width=1.0)
        h.contour.filled_contours=True # This is the trick to fill the contours
        
        mlab.colorbar(object=h,orientation='vertical')
        mlab.view(0,0)
        self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
    
    def isosurface(self,vv=[4.0],clim=None,**kwargs):
        
        if not self.is3D:
            raise RuntimeError('isosurface can only be used when is3D=True')
        
        if clim==None:
            clim = [self.data.min(), self.data.max()]
            
        self.fig=mlab.gcf()
        self.fig.scene.background = (0.,0.,0.)
        
        # Convert the cell centred data into a scene source
        # Need to set use point (vertex) data
        src = mlab.pipeline.cell_to_point_data(self.ug)
        
        # Add the iso_surface module to the scene
        h=mlab.pipeline.iso_surface(src,contours=vv,line_width=1.0,**kwargs)
        
        mlab.colorbar(object=h,orientation='vertical')
        mlab.view(0,0)
        self.title=mlab.title(self._Spatial__genTitle(),height=0.95,size=0.15)
        
    def loadData(self):
        """
        Overloaded loadData function - updates the unstructured grid object
        """
        Spatial.loadData(self)
        self.data=np.ravel(self.data[self.mask3D])
        self.ug.cell_data.scalars = self.data
        self.ug.cell_data.scalars.name = 'suntans_scalar' 
        self.ug.modified()        
        
    def __setitem__(self,key,value):
        """
        Updates the unstructured grid object, ug, when the 'data' key is updated
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
                       
# Testing stuff here

ncfile = 'C:/Projects/GOMGalveston/MODELLING/GalvestonCoarse/rundata/GalvCoarse_TidesRivers3D*nc'

sunvtk = SunTvtk(ncfile,variable='salt',tstep=4,is3D=True)

sunvtk.loadData()

sunvtk.isosurface(vv=[4.0,12.0,20.0,28.0],transparent=False,opacity=0.5)

#sunvtk.plot(representation='wireframe')

#sunvtk.contour(vv=40)

#sunvtk.contour(vv=map(None,np.linspace(-0.3 ,0.3,50)))