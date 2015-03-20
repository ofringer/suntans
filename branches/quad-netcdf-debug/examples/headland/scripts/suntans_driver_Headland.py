# -*- coding: utf-8 -*-
"""
Driver script for generating suntans input files

Created on Fri Oct 05 11:24:10 2012

@author: mrayson
"""

from sunpy import Grid
from sundepths import DepthDriver
from sunboundary import modifyBCmarker, Boundary, InitialCond
from metfile import SunMet

import numpy as np 
import matplotlib.pyplot as plt
import operator

PI = np.pi

######################################################################
#   General input variables
######################################################################

suntanspath = 'rundata'

#bcpolygonfile = 'grids/boundarytype2_headland.shp'

# Options to generate each type of file
makebathy=True # Interpolate raw depths on grid
makebnd=True # Boundary condition file
makewinds=False # Meteorological input file
makeinitial=True # Initial condition file

# Time parameters
starttime = '20130101.0000' # Start time of BC and wind file 'yyyymmdd.HHMM'
endtime = '20130201.0000' # End time of BC and wind file 'yyyymmdd.HHMM'
dt = 1800.0 # time step of BC and wind file (seconds)

u0 = 0.5
omega = 1.41e-4

# Grid parameters
L = 5e4
W = 2.5e4
a =2e3
b =8e3
x0 = 0

# Depth parameters
hmin = 1.0
hmax = 20.
w = 3000.

def gaussian(x,x0,a,b,pow=2.):
    return b*np.exp(-0.5*((x-x0)/a)**pow)

def calc_depth(x,y,hmin,hmax,w):
    y_bot = gaussian(x,x0,a,b)
    y_hat = y - y_bot # Distance from point to line

    ind = y_hat < w
    h = np.zeros_like(y)
    h[ind] = hmin+(hmax-hmin)/w*y_hat[ind]
    h[~ind] = hmax

    h[h<hmin]=hmin

    return h

def calc_depth2(y_hat,hmin,hmax,w):

    ind = y_hat < w
    h = np.zeros_like(y_hat)
    h[ind] = hmin+(hmax-hmin)/w*y_hat[ind]
    h[~ind] = hmax

    h[h<hmin]=hmin

    return h


####
# Verical grid options
####
Nkmax = 1 # number of layers
r = 1.30 # vertical stretching parameter

###
# Free surface function
###
def calc_fs(xv):
    a = 0.0
    Lx = 10000.
    return -xv*a/Lx

#################################
# Load the model grid and bathymetry into the object 'grd' and initialise the vertical spacing
##################################
grd = Grid(suntanspath)


##########
# Modify the boundary markers and create the boundary condition file
##########
# This changes the edge types based on the polygons in the shapefile
#modifyBCmarker(suntanspath,bcpolygonfile)

# Modify the left and right edges and convert to type 2
hgrd = grd.convert2hybrid()

grd.mark[grd.mark>0]=1 # reset all edges to land

# convert edges +/- half a grid cell from the edge
dx = hgrd.dg.max()
xmin = hgrd.xe.min()+dx/4.0
xmax = hgrd.xe.max()-dx/4.0

indleft = operator.and_(hgrd.mark>0, hgrd.xe < xmin) # all boundaries
grd.mark[indleft]=2
indright = operator.and_(hgrd.mark>0, hgrd.xe > xmax) # all boundaries
grd.mark[indright]=2

edgefile = suntanspath+'/edges.dat'
grd.saveEdges(edgefile)
print 'Updated markers written to: %s'%(edgefile)

#########################################
# Interpolate the depths onto the grid
#########################################
y_hat,idx = grd.find_nearest_boundary(markertype=1)
y_hat[grd.yv>b+w]=1e15
grd.dv = calc_depth2(y_hat,hmin,hmax,w)
#grd.dv = calc_depth(grd.xv,grd.yv,hmin,hmax,w)
grd.saveBathy('%s/depths.dat-voro'%suntanspath)




#Load the boundary object from the grid
#   Note that this zeros all of the boundary arrays
bnd = Boundary(suntanspath,(starttime,endtime,dt))

bnd.setDepth(grd.dv)

# Try out a linear interpolation 
#y0 = bnd.ye.min()
#y1 = bnd.ye.max()
#du = 0.0
#dy = y1-y0
#Utst = (bnd.ye.ravel()-y0)/dy*du

# Normal components
nx = hgrd.n1[grd.mark==2]
ny = hgrd.n2[grd.mark==2]

t = bnd.tsec-bnd.tsec[0]
for k in range(bnd.Nk):
   for ii in range(bnd.N2):

       umag = u0
       u = umag * np.sin(omega*t)
       sfac=1.0*bnd.de[ii]/hmax
       if indright[bnd.edgep[ii]]:
           sfac *= -1.0
       bnd.boundary_u[:,k,ii] = sfac*u*nx[ii]
       bnd.boundary_v[:,k,ii] = sfac*u*ny[ii]

# Set the u-velocity along
#bnd.boundary_u[:,:,:] = flowspeed

# type-3 
#h = calc_fs(bnd.xv)

#bnd.h[:] = np.ones_like(bnd.h) * h



# Write the boundary file
bcfile = 'Headland_BC.nc'
bnd.write2NC(suntanspath+'/'+bcfile)


##############################################
# Generate an initial condition netcdf file
##############################################
if makeinitial:
    ####
    # Initial condition options
    ####

    icfile = 'Headland_IC.nc' # Initial condition netcdf file name

    ####
    # End options
    ####

    # Initialise the class
    IC = InitialCond(suntanspath,starttime)

    
    # the initial condition arrays are stored in the following fields in the class:
    #   h, uc, vc, T, S
    # 
    # uc, vc, T and S have dimensions [1, Nk, Nc]
    # h has dimensions [1, Nc]
    
    # We want to set h to constant (other fields are zero by default)
    IC.h[:] = calc_fs(IC.xv) 

    # Set T as a circle
    x0 = -0
    y0 = b+w
    R = 500. # Radius

    dist = np.sqrt( (IC.xv-x0)**2 + (IC.yv-y0)**2.)

    IC.T[...,dist<=R] = 1.


    
    # Write the initial condition file
    IC.writeNC(suntanspath+'/'+icfile,dv=grd.dv)


