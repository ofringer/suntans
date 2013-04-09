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

PI = 3.141592653589793

######################################################################
#   Input variables
######################################################################

suntanspath = 'rundata'

makebathy=True
makebnd=True
makewinds=True
makeinitial=True

# Time parameters
starttime = '20000601.0000' # Start time of BC and wind file 'yyyymmdd.HHMM'
endtime = '20000801.0000' # End time of BC and wind file 'yyyymmdd.HHMM'
dt = 3600.0 # time step of BC and wind file (seconds)


####
# Bathymetry interpolation options
####
depthfile = 'GIS/depths.shp'

depthmax=-0.1

# Interpolation method
interpmethod='idw' # 'nn', 'idw', 'kriging', 'griddata'

# Type of plot
plottype='mpl' # 'mpl', 'vtk2' or 'vtk3'

# Interpolation options
NNear=6
p = 1.0 #  power for inverse distance weighting
# kriging options
varmodel = 'spherical'
nugget = 0.1
sill = 0.8
vrange = 250.0

# Projection conversion info for input data
convert2utm=False
CS='NAD83'
utmzone=15
isnorth=True
vdatum = 'MSL'

# Smoothing options
smooth=True
smoothmethod='kriging' # USe kriging or idw for smoothing
smoothnear=4 # No. of points to use for smoothing

####
# Verical grid options
####
Nkmax = 24 # number of layers
r = 1.00 # vertical stretching parameter

####
# Open boundary options
####
bcpolygonfile = 'GIS/boundary_polygons.shp'
#bcpolygonfile = '../GIS/boundary_polygons_closed.shp'

bcfile = 'Estuary_BC.nc'

h0 = 0.25 # tidal range
omega = 2*PI/(24.0 * 3600.0) # Tidal frequency
T0 = 30 # Open boundary and initial condition background temperature
S0 = 32 # Open boundary and initial condition background salinity

edge_id = 1234 # Unique identifier of flux edge
Q = 100.0 # Volume flux rate boundary [m^3/s]

####
# Atmospheric input options
####

metfile = 'Estuary_MetForcing.nc'

Uwind = 0.0
Vwind = 5.0
RH = 70.0
Tair_mean = 30.0
Tair_amp = 5.0


####
# Initial condition options
####

icfile = 'Estuary_IC.nc'

# Sets initial temperature and salinity to T0 and S0, respectively

######################################################################
# End of input variables - do not change anything below here.
######################################################################

# Interpolate the depths onto the grid
if makebathy:
    D = DepthDriver(depthfile,interpmethod=interpmethod,plottype=plottype,NNear=NNear,\
    p=p,varmodel=varmodel,nugget=nugget,sill=sill,vrange=vrange,convert2utm=convert2utm,\
    utmzone=utmzone,isnorth=isnorth,vdatum=vdatum,smooth=smooth,smoothmethod=smoothmethod,\
    smoothnear=smoothnear)
    
    D(suntanspath,depthmax=depthmax)


# Load the model grid and bathymetry into the object 'grd' and initialise the vertical spacing
grd = Grid(suntanspath)

# Load the depth data into the grid object
grd.loadBathy(suntanspath+'/depths.dat-voro')
zmax = np.abs(grd.dv.min())

# Set up the vertical coordinates
dz = grd.calcVertSpace(Nkmax,r,zmax)
grd.setDepth(dz)

# Save vertspace.dat
grd.saveVertspace(suntanspath+'/vertspace.dat')



# Modify the boundary markers and create the boundary condition file
if makebnd:
    # This changes the edge types based on the polygons in the shapefile
    modifyBCmarker(suntanspath,bcpolygonfile)
    
    #Load the boundary object from the grid
    #   Note that this zeros all of the boundary arrays
    bnd = Boundary(suntanspath,(starttime,endtime,dt))
    
    #Fill the variables in the boundary object with your data here...
    #
    #Type 3 variables are at the cell-centre (xv,yv) and are named:
    #            uv, vc, wc, h, T, S
    #            
    #            Dimensions: [Nt,Nk,N3]
    #            
    #Type 2 variables are at the cell edges (xe, ye) and are named:
    #            boundary_u
    #            boundary_v
    #            boundary_w
    #            boundary_T
    #            boundary_S
    #            (no h)
    #            
    #            Dimensions: [Nt, Nk, N2]
    
    t = bnd.ncTime()

    # Set an oscillatory water level boundary along type 3 cells (T and S constant)
    for ii in range(bnd.N3):
        bnd.h[:,ii] = h0*np.sin(omega*t)
        for k in range(bnd.Nk):
            bnd.T[:,k,ii] = T0
            bnd.S[:,k,ii] = S0
            
    # Set a constant discharge along the type 2 boundary edge
    segmentID = [edge_id] # segment identifiers
    fluxrate = [Q] # constant flux rate at each of the segments
    for ID,flux in zip(segmentID,fluxrate):
        ind = np.argwhere(bnd.segp==ID)
        bnd.boundary_Q[:,ind]=flux
        
    # Set all type-2 boundary T and S constant
    bnd.boundary_S[:,:,:] = 0.0
    bnd.boundary_T[:,:,:] = T0
    
    # Write the boundary file
    bnd.write2NC(suntanspath+'/'+bcfile)
    
if makeinitial:
    # Generates an initial condition netcdf file
    
    # Initialise the class
    IC = InitialCond(suntanspath,starttime)
    
    # the initial condition arrays are stored in the following fields in the class:
    #   h, uc, vc, T, S
    
    # We just want to set T and S to constant (other fields are zero by default)
    IC.T[:] = T0
    IC.S[:] = S0
    
    # Write the initial condition file
    IC.writeNC(suntanspath+'/'+icfile)

if makewinds:
    # Create a meteorological input file
    xpt = grd.xv.mean()
    ypt = grd.yv.mean()
    zpt = 10.0
    met = SunMet(xpt,ypt,zpt,(starttime,endtime,dt))
    
    # Set a cyclic air temperature and all other variables constant (rain and cloud =0)
    omegaT = 2*PI/(24.0 * 3600.0) # diurnal frequency
    for ii in range(met.Tair.Npt):
        met.Tair.data[:,ii] = Tair_mean + Tair_amp * np.sin(omegaT * met.nctime)
        
    met.Uwind.data[:] = Uwind
    met.Vwind.data[:] = Vwind
    met.RH.data[:] = RH = RH
    met.cloud.data[:] = 0.0
    met.rain.data[:] = 0.0

    met.write2NC(suntanspath+'/'+metfile)
