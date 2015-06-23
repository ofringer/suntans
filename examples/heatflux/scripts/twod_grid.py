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
outpath = 'rundata'

nx = 4
###
xy = ll2utm(np.array((lon0,lat0)),utmzone)
x0 = xy[0,0]
y0 = xy[0,1]

# Create the cell outline
dx = 1000.*nx
dx2 = dx/2.
xlims = [x0-dx2,x0+dx2]
ylims = [y0-dx2,y0+dx2]

xgrd = np.linspace(xlims[0],xlims[1],nx)
ygrd = np.linspace(ylims[0],ylims[1],nx)

cells=[]
xp = []
yp = []

def pntindx(j,i,nrows):
    return j*nrows + i

for jj in range(nx):
    for ii in range(nx):
        xp.append(xgrd[ii])
        yp.append(ygrd[jj])

for jj in range(nx-1):
    for ii in range(nx-1):
        cells.append([pntindx(jj,ii,nx), pntindx(jj+1,ii,nx),\
            pntindx(jj+1,ii+1,nx),pntindx(jj,ii+1,nx)])
        
Nc = (nx-1)**2
nfaces = 4*np.ones((Nc,),np.int)

cells = np.array(cells)
# Convert to a suntans grid
grd = HybridGrid(xp,yp,cells,nfaces=nfaces)

grd2suntans(grd,outpath)
