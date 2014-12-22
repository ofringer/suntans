"""
Cartesian grid utilities

Matt Rayson
Stanford University
October 2014
"""

import numpy as np
from scipy import sparse

class RegGrid(object):
    """
    Class for a regularly spaced cartesian grid
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

        i = np.floor( (x-self.ox)/self.dx )
        j = np.floor( (y-self.oy)/self.dy )

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

        return i.astype(int),j.astype(int)

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


