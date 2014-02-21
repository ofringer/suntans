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

PI = np.pi

######################################################################
#   General input variables
######################################################################

suntanspath = 'data'

# Options to generate each type of file
makebathy=False # Interpolate raw depths on grid
makebnd=True # Boundary condition file
makewinds=True # Meteorological input file
makeinitial=True # Initial condition file

# Time parameters
starttime = '20130101.0000' # Start time of BC and wind file 'yyyymmdd.HHMM'
endtime = '20130201.0000' # End time of BC and wind file 'yyyymmdd.HHMM'
dt = 3600.0 # time step of BC and wind file (seconds)


####
# Verical grid options
####
Nkmax = 30 # number of layers
r = 1.30 # vertical stretching parameter


#########################################
# Interpolate the depths onto the grid
#########################################
if makebathy:
    ####
    # Bathymetry interpolation options
    ####

    # Some file containing depth data
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
    utmzone=10
    isnorth=True
    vdatum = 'MSL'

    # Smoothing options
    smooth=True
    smoothmethod='kriging' # USe kriging or idw for smoothing
    smoothnear=4 # No. of points to use for smoothing
    ###
    # End options
    ###

    D = DepthDriver(depthfile,interpmethod=interpmethod,plottype=plottype,NNear=NNear,\
    p=p,varmodel=varmodel,nugget=nugget,sill=sill,vrange=vrange,convert2utm=convert2utm,\
    utmzone=utmzone,isnorth=isnorth,vdatum=vdatum,smooth=smooth,smoothmethod=smoothmethod,\
    smoothnear=smoothnear)
    
    D(suntanspath,depthmax=depthmax)

#################################
# Load the model grid and bathymetry into the object 'grd' and initialise the vertical spacing
##################################
grd = Grid(suntanspath)

# Load the depth data into the grid object
grd.loadBathy(suntanspath+'/depth.dat-voro')
zmax = np.abs(grd.dv.max())
print zmax

# Set up the vertical coordinates
dz = grd.calcVertSpace(Nkmax,r,zmax) # Default version
#def calcvertspace(Nkmax,r,depth):
#    ktop=15
#    dztop=2.0


grd.setDepth(dz)

# Save vertspace.dat
grd.saveVertspace(suntanspath+'/vertspace.dat')

# Save cells.dat to ensure that is in the right format
grd.saveCells(suntanspath+'/cells.dat')

###################################################################
# Modify the boundary markers and create the boundary condition file
###################################################################
if makebnd:
    ####
    # Open boundary options
    ####

    # Shapefile with boundary condition types (set to None to skip)
    bcpolygonfile = 'gis/sfbay_boundary_poly.shp'

    bcfile = 'SFBay3D_BC.nc' # output boundary netcdf file name

    hamp = 0.5 # tidal range
    h0 = -5.0
    omega = 2*PI/(12.42 * 3600.0) # Tidal frequency
    T0 = 0 # Open boundary and initial condition background temperature
    S0 = 32 # Open boundary and initial condition background salinity

    # These must correspond to the edge_id field in the boundary polygon shapefile
    # For this case 101 : Sacramento Rv, 102 : San Joaquin Rv
    edge_id = [101,102] # Unique identifier of each flux edge (list)
    Q = [0.0,0.0] # Volume flux rate of each segment boundary [m^3/s]

    ####
    # End of options
    ####


    # This changes the edge types based on the polygons in the shapefile
    if not bcpolygonfile==None:
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
        bnd.h[:,ii] = h0 + hamp*np.cos(omega*t)
        for k in range(bnd.Nk):
            bnd.T[:,k,ii] = T0
            bnd.S[:,k,ii] = S0
            
    # Set a constant discharge along the type 2 boundary edge
    for ID,flux in zip(edge_id,Q):
        ind = np.argwhere(bnd.segp==ID)
        bnd.boundary_Q[:,ind]=flux
        
    # Set all type-2 (river) boundary T and S constant
    bnd.boundary_S[:,:,:] = 0.0
    bnd.boundary_T[:,:,:] = T0
    
    # Write the boundary file
    bnd.write2NC(suntanspath+'/'+bcfile)

##############################################
# Generate an initial condition netcdf file
##############################################
if makeinitial:
    ####
    # Initial condition options
    ####

    icfile = 'SFBay3D_IC.nc' # Initial condition netcdf file name

    h0 = -5.0 # initial free surface elevation
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
    IC.h[:] = h0

    # Now set T and S based on analytical functions
    def returnSalinity(x,y):
        x0 = 5.7e5
        S = 32.0*0.5*(1.0+np.tanh(-(x-x0)/1e4))
        S[y<4.2e6] = 32.
        return S

    def returnTemperature(x,y):
        x0 = 545691.
        y0 = 4185552.
        R = 1000.
        T = np.zeros_like(x)
        ind = (x-x0)**2 + (y-y0)**2 < R*R
        T[ind] = 1.0
        # ELSE 0
        return T

    # Horizontal coordinates are stored in the xv/yv attributes
    S = returnSalinity(IC.xv,IC.yv) # returns a 1-D array length Nc
    T = returnTemperature(IC.xv,IC.yv)

    # Now fill in the values in the array
    for k in range(IC.Nkmax):
        IC.T[:,k,:] = T
        IC.S[:,k,:] = S

    
    # Write the initial condition file
    IC.writeNC(suntanspath+'/'+icfile)

#######################################
# Create a meteorological input file
#######################################
if makewinds:
    ####
    # Atmospheric input options
    ####

    metfile = 'SFBay_MetForcing.nc'

    Uwind = 0.0
    Vwind = -5.0
    RH = 70.0
    Tair_mean = 20.0
    Tair_amp = 5.0
    ####
    # End of options
    ####


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
