# -*- coding: utf-8 -*-
"""
Create the fully coupled Galveston Input files

Created on Mon Apr 01 18:00:32 2013

@author: mrayson
"""

from sundriver import sundriver, dumpinputs

###
# Variable inputs for the class
starttime = '20100401.000000'
endtime = '20100501.000000'
dt = 3600.0
sunpath = 'rundata'

plotdir = 'plots/initial'
###

# Initiate the driver class
sun = sundriver()

# Switches to generate bathymetry, boundary, meteorology and initial condition input files
sun.makebathy=False
sun.makebnd=True
sun.makewinds=False
sun.makeinitial=False

sun.modifyedges = True# Option to modify the boundary edges

###
# General options
###
# Grid projection variables
sun.convert2utm=False
sun.CS='NAD83'
sun.utmzone=15
sun.isnorth=True
sun.vdatum = 'MSL'

# Verical grid options
sun.Nkmax = 10 # number of layers
sun.r = 1.0 # vertical stretching parameter
sun.setconstantdepth=True
sun.H0=10.0

###
# Bathymetry interpolation options
###
sun.depthfile = '../DATA/NOAA_25m_UTM_DEM.nc'
sun.depthmax=-0.1
sun.interpmethod='idw' # Interpolation method:  'nn', 'idw', 'kriging', 'griddata'
sun.plottype='mpl' # Type of plot: 'mpl', 'vtk2' or 'vtk3'

# Interpolation options
sun.NNear=3

# IF interpmethod = 'idw' 
sun.p = 1.0 #  power for inverse distance weighting

# IF interpmethod = 'kriging' 
sun.varmodel = 'spherical'
sun.nugget = 0.1
sun.sill = 0.8
sun.vrange = 250.0

# Smoothing options
sun.smooth=False
sun.smoothmethod='kriging' # USe kriging or idw for smoothing
sun.smoothnear=4 # No. of points to use for smoothing


####
# Open boundary options
####
sun.opt_bcseg = 'constant' # Segment boundary condition option: 'constant' or 'file'
sun.opt_bctype2 = 'constant' # Type 2 boundary condition option: 'constant'
sun.opt_bctype3 = 'OTIS' # Type 3 boundary condition option: 'constant', ,'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE'

sun.bcpolygonfile = 'janet/BoundaryPolygon.shp' # Shape file with fields 'marker' and 'edge_id'
sun.bcfile = 'Plume_BC.nc' # Input boundary condition file
    
# IF opt_bcseg = 'consant
sun.Q0 = 1260.0 # m3/s

# IF opt_bctype2/opt_bctype3 = 'constant'
sun.T0 = 30.0 # Open boundary and initial condition background temperature
sun.S0 = 30.0 # Open boundary and initial condition background salinity

# IF opt_bctype2 = 'file' or 'ROMSFILE'
sun.waterlevelstationID = 8772447 # USCG Freeport gage

# IF opt_bytype2 = 'filt'
sun.TairstationID = '722436-12906' # Houston-Ellington weather station

sun.filttype='low'
sun.cutoff = 7200.0

# Option to zero ROMS uv and eta at boundaries and initial conditions
sun.useROMSuv=False
sun.useROMSeta=False
####
# Initial condition options
####
sun.opt_ic = 'constant' #'constant', 'depth_profile', 'ROMS'

sun.icfile = 'PlumeQuad_IC.nc'

###
# Metfile
### 

# For plotting only
sun.metfile = 'Galveston_NARR_2010.nc'
###
# Input file names
### 
sun.romsfile = ['../DATA/txla_subset_HIS_2010.nc']
sun.otisfile = '/home/mrayson/GalvestonBay/DATA/Tides/Model_Mex'
sun.dbasefile = '../DATA/GalvestonObs.db'

######################
# Now call the class...
######################

sun(sunpath,starttime,endtime,dt)


###
# Dump figures of the input data
#dump = dumpinputs(suntanspath=sunpath,icfile=sun.icfile,bcfile=sun.bcfile,metfile=sun.metfile)

#0dump(plotdir)
