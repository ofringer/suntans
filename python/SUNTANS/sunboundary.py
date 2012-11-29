# -*- coding: utf-8 -*-
"""
Tools for modifying suntans boundary condition files


Examples:
---------

Example 1) Modify the boundary condition markers with a shapefile:
------------------------------------------------------------------
    
    >>from sunboundary import modifyBCmarker
    >># Path to the grid
    >>suntanspath = 'C:/Projects/GOMGalveston/MODELLING/GRIDS/GalvestonCoarseBC'
    >># Name of the shapefile
    >>bcfile = '%s/Galveston_BndPoly.shp'%suntanspath 
    >>modifyBCmarker(suntanspath,bcfile)
    
Created on Fri Nov 02 15:24:12 2012

@author: mrayson
"""


import sunpy
import numpy as np
import matplotlib.pyplot as plt
from maptools import readShpPoly
import matplotlib.nxutils as nxutils #inpolygon equivalent lives here

class Boundary(object):
    """
    Generic SUNTANS boundary class
    """
    def __init__(self,suntanspath):
        self.grd = sunpy.Grid(suntanspath)
        
    def loadBoundary(self):
        """
        Load the coordinates and indices for type 2 and 3 BC's
        """
        ind2 = np.argwhere(self.grd.mark==2)
        ind3 = np.argwhere(self.grd.mark==3)
        
        # Edge index of type 2 boundaries
        self.edgep = ind2
        self.N2 = len(self.edgep)
        
        # Cell index of type 3 boundaries
        cellp1 = self.grd.grad[ind3,0]
        cellp2 = self.grd.grad[ind3,1]
        self.cellp=[]
        for c1,c2 in zip(cellp1,cellp2):
            if c1==-1:
                self.cellp.append(c2)
            elif c2==-1:
                self.cellp.append(c1)
        
        self.N3 = self.cellp
        
        # Store the coordinates of the type 2 and 3 boundaries
        self.xv = self.grd.xv[self.cellp]
        self.yv = self.grd.yv[self.cellp]
        
        # Find the edge points
        xe = np.mean(self.grd.xp[self.grd.edges],axis=1)
        ye = np.mean(self.grd.yp[self.grd.edges],axis=1)
        self.xe = xe[self.edgep]
        self.ye = ye[self.edgep]
        
        
        
def modifyBCmarker(suntanspath,bcfile):
    """
    Modifies SUNTANS boundary markers with a shapefile

    The shapefile must contain polygons with the integer-type field "marker"
    """
    
    print '#######################################################'
    print '     Modifying the boundary markers for grid in folder:'
    print '         %s'%suntanspath

    # Load the grid into an object
    grd = sunpy.Grid(suntanspath)
    
    # Find the edge points
    xe = np.mean(grd.xp[grd.edges],axis=1)
    ye = np.mean(grd.yp[grd.edges],axis=1)
    
    # Read the shapefile
    XY,newmarker = readShpPoly(bcfile,FIELDNAME='marker')
    if len(XY)<1:
        print 'Error - could not find any polygons with the field name "marker" in shapefile: %s'%bcfile
        return
    
    # Plot before updates
    #plt.figure()
    #grd.plotBC()
    #plt.plot(XY[0][:,0],XY[0][:,1],'m',linewidth=2)
    #plt.show()
    
    # Find the points inside each of the polygon and assign new bc
    for xpoly, bctype in zip(XY,newmarker):
        ind0 = grd.mark>0
        edges = np.asarray([xe[ind0],ye[ind0]])
        mark = grd.mark[ind0]
        ind1 = nxutils.points_inside_poly(edges.T,xpoly)
        mark[ind1]=bctype
        grd.mark[ind0]=mark
    
    # Save the new markers to edges.dat
    edgefile = suntanspath+'/edges.dat'
    grd.saveEdges(edgefile)
    print 'Updated markers written to: %s'%(edgefile)
    
    # Plot the markers
    plt.figure()
    grd.plotBC()
    plt.plot(XY[0][:,0],XY[0][:,1],'m',linewidth=2)
    figfile = suntanspath+'/BoundaryMarkerTypes.pdf'
    plt.savefig(figfile)
    
    print 'Marker plot saved to: %s'%(figfile)
    print 'Done.'
    print '#######################################################'
    

###################
# Testing stuff
