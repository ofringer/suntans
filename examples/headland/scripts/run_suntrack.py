# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 14:52:43 2013

@author: mrayson
"""

from suntrack import SunTrack, GridParticles
import numpy as np
from maptools import readShpPoly

########
# Inputs
ncfile = 'data/Headland_0*.nc'

outfile = 'data/Headland_particles.nc'

timeinfo = ('20130101.0200','20130102.0820',120.0)

dtout = 600.

dx = 100.
dy = 100.

# Grid parameters
W = 2.5e4
a =2e3
b =8e3
x0 = 0

########

# Initialse the particle locations using a polygon
# With a shapefile
#XY,newmarker = readShpPoly(polyfile,FIELDNAME=None )

# Polygon to fill with particles
XY = np.array([
    [x0-W/2, 0],\
    [x0-W/2, 2*b],\
    [x0+W/2, 2*b],\
    [x0+W/2, 0],\
    [x0-W/2, 0],\
    ])

x,y,z = GridParticles(ncfile,dx,dy,1,xypoly=XY)


# Initialise the particle tracking object
sun = SunTrack(ncfile, interp_method='mesh', interp_meshmethod='nearest',
is3D=False)

## Initialise the age polygon
#XY,newmarker = readShpPoly(agepolyfile,FIELDNAME='marker')
#agepoly=XY[0]
agepoly=None


# Run the model without animating the data
sun(x,y,z,timeinfo,agepoly=agepoly,outfile=outfile,dtout=dtout)

