# -*- coding: utf-8 -*-
"""
Create a one cell quad grid

Created on Mon Nov 18 10:23:26 2013

@author: mrayson
"""

import numpy as np
from hybridgrid import HybridGrid
from gmsh import grd2suntans
from maptools import ll2utm

###
# Inputs
lon0 = -94.86367
lat0 = 29.29517
utmzone= 15
dx = 1000.
outpath = 'rundata'

periodic=True
###
xy = ll2utm(np.array((lon0,lat0)),utmzone)
x0 = xy[0,0]
y0 = xy[0,1]

# Create the cell outline
dx2 = dx/2.
xp = [x0-dx2,x0-dx2,x0+dx2,x0+dx2]
yp = [y0-dx2,y0+dx2,y0+dx2,y0+dx2]

cells = np.zeros((1,4),np.int)
cells[0,:] = [0,1,2,3]

nfaces=np.array([4])

# Convert to a suntans grid
grd = HybridGrid(xp,yp,cells,nfaces=nfaces)

# Manually set a couple of things here

#1) set the cell centers 
grd.xv[0]=x0
grd.yv[0]=y0

if periodic:
    #2) Set the neighbouring cells to itself. This should make the grid doubly
    #   periodic
    grd.neigh[:] = 0 

    #3) Set the edge type to be 0 and the edge connection to be the same cell
    grd.mark[:]=0
    grd.grad[:]=0

grd2suntans(grd,outpath)
