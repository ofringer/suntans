# -*- coding: utf-8 -*-
"""
Create a western boundary current domain

Created on Mon Nov 18 10:23:26 2013

@author: mrayson
"""

import os
from shutil import copyfile
import numpy as np
import matplotlib.pyplot as plt
import operator

from soda.dataio.ugrid.ugridgen import cartesian_ugrid_gen
from soda.dataio.suntans.sunboundary import modifyBCmarker, Boundary, InitialCond
from soda.dataio.suntans.sunpy import Grid
from soda.dataio.suntans.metfile import SunMet

PI = np.pi
GRAV=9.81

####################################################
# Inputs

# Depth [m]
H = 5000.0

# Size of domain
dx = 50e3 # m
nz = 1

L = 1000e3 # m

suntanspath = 'data'

starttime = '20000101.000000'
endtime = '20020101.000000'
dt = 7*86400.

icfile = 'WestBC_IC.nc'

makewinds = True
metfile = 'WestBC_met.nc'
####################################################

#####
# Create the grid
if not os.path.isdir(suntanspath):
    print 'Creating new directory: %s'%suntanspath
    os.mkdir(suntanspath)
    copyfile('rundata/suntans.dat','%s/suntans.dat'%suntanspath)

xlims = [0,L]
ylims = [0,L]


# Create the grid
grd= cartesian_ugrid_gen(xlims, ylims, dx, suntanspath=suntanspath)

# Load the grid
grd = Grid(suntanspath)

grd.dv = H*np.ones_like(grd.xv)
grd.saveBathy('%s/depth.dat-voro'%suntanspath)

grd.dz = grd.calcVertSpace(nz, 1.0, H)
grd.saveVertspace('%s/vertspace.dat'%suntanspath)
grd.setDepth(grd.dz)


#########
# Create the initial conditions file
#########
IC = InitialCond(suntanspath,starttime)

###
#Update the latitude variable based on the prescribed Coriolis
omega=2*np.pi/86400.
beta = 2e-11
cor_f = 1e-4 + beta*grd.yv
rad2deg = 180./np.pi

IC.latv[:] = np.arcsin(cor_f/(2*omega)) *rad2deg
IC.lonv[:] = 0 # not used

print IC.latv

# Write the initial condition file
IC.writeNC(suntanspath+'/'+icfile,dv=grd.dv)

##########
# Make the wind file
##########
if makewinds:
    # Create a meteorological input file
    dxw = 100e3 # 100 km resolution
    x = np.arange(0,L+dxw,dxw)
    X,Y = np.meshgrid(x,x)
    xpt = X.ravel()
    ypt = Y.ravel()
    zpt = 10.0
    met = SunMet(xpt,ypt,zpt,(starttime,endtime,dt))

    # Create a constant velocity wind field
    tau0 = 1e-1 # N m-2
    rho_a = 1.
    Cd = 0.01
    U0 = np.sqrt(tau0/(rho_a*Cd))
    Uwind = -U0*np.cos(np.pi*ypt/L)
    
       
    met.Uwind.data[:,:] = Uwind[np.newaxis,...]
    met.Vwind.data[:] = 0.0
    met.RH.data[:] = 0.0
    met.cloud.data[:] = 0.0
    met.rain.data[:] = 0.0

    met.write2NC(suntanspath+'/'+metfile)


#grd.plotmesh()
#plt.show()

