# -*- coding: utf-8 -*-
"""
Create a one cell quad grid

Created on Mon Nov 18 10:23:26 2013

@author: mrayson
"""

import numpy as np
import matplotlib.pyplot as plt
from hybridgrid import HybridGrid
from inpolygon import inpolygon

###
# Inputs
L = 5e4
W = 2.5e4
a =2e3
b =8e3
x0 = 0

# target grid scale
scale = 200.#
outpath = 'grids/quad/'
###

###
# Create the boundary outline
###
X = np.arange(x0-L/2,x0+L/2+scale,scale)


def gaussian(x,x0,a,b,pow=2.):
    return b*np.exp(-0.5*((x-x0)/a)**pow)

y_bot = gaussian(X,x0,a,b)
y_top = W*np.ones_like(y_bot)

y = np.hstack((y_bot,y_top[::-1],y_bot[0]))
x = np.hstack((X,X[::-1],X[0]))

xypoly =[(x[ii],y[ii]) for ii in range(x.shape[0])]


# Create the grid outline
xlims = [x.min(),x.max()]
ylims = [y.min(),y.max()]

xgrd = np.arange(xlims[0]-scale/2,xlims[1]+1.5*scale,scale)
ygrd = np.arange(ylims[0]-scale/2,ylims[1]+1.5*scale,scale)

nx = xgrd.shape[0]
ny = ygrd.shape[0]

# Create a mask polygon
X,Y = np.meshgrid(xgrd,ygrd)
XY = np.vstack((X.ravel(),Y.ravel())).T
mask = inpolygon(XY,xypoly)
mask = mask.reshape((ny,nx))

cells=[]
xp = []
yp = []

def pntindx(j,i,nrows):
    return j*nrows + i

for jj in range(ny):
    for ii in range(nx):
        #if mask[jj,ii]:
        xp.append(xgrd[ii])
        yp.append(ygrd[jj])

for jj in range(ny-1):
    for ii in range(nx-1):
        if mask[jj,ii] and mask[jj+1,ii] and mask[jj+1,ii+1] and mask[jj,ii+1]:
            cells.append([pntindx(jj,ii,nx), pntindx(jj+1,ii,nx),\
                pntindx(jj+1,ii+1,nx),pntindx(jj,ii+1,nx)])
        
Nc = len(cells)
nfaces = 4*np.ones((Nc,),np.int)

cells = np.array(cells)
# Convert to a suntans grid
grd = HybridGrid(xp,yp,cells,nfaces=nfaces)

grd.write2suntans(outpath)
