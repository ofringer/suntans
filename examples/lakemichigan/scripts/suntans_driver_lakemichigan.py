# -*- coding: utf-8 -*-
"""
Create the fully coupled Galveston Input files

Created on Mon Apr 01 18:00:32 2013

@author: mrayson
"""

from sundriver import sundriver, dumpinputs

###
# Variable inputs for the class
starttime = '20120301.000000'
endtime = '20120401.000000'
dt = 3600.0
sunpath = 'data'

plotdir = 'plots/inputs'
###

# Initiate the driver class
sun = sundriver()

# Switches to generate bathymetry, boundary, meteorology and initial condition input files
sun.makebathy=False
sun.makebnd=False
sun.makewinds=False
sun.makeinitial=True

sun.modifyedges = False # Option to modify the boundary edges

###
# General options
###
# Grid projection variables
sun.convert2utm=True
sun.CS='NAD83'
sun.utmzone=16
sun.isnorth=True
sun.vdatum = 'MSL'
sun.shapefieldname='depth'

# Verical grid options
sun.Nkmax = 40 # number of layers
sun.r = 1.05 # vertical stretching parameter

###
# Bathymetry interpolation options
###
sun.depthfile = 'gis/bathymetry_lake_michigan.shp'
sun.depthmax=0.1
sun.interpmethod='kriging' # Interpolation method:  'nn', 'idw', 'kriging', 'griddata'
sun.plottype='mpl' # Type of plot: 'mpl', 'vtk2' or 'vtk3'

# Interpolation options
sun.NNear=6

# Interpolate to nodes then take maximum depth for cell
sun.interpnodes=False

# IF interpmethod = 'idw' 
sun.p = 1.0 #  power for inverse distance weighting

# IF interpmethod = 'kriging' 
sun.varmodel = 'spherical'
sun.nugget = 0.1
sun.sill = 0.8
sun.vrange = 2500.00

# Smoothing options
sun.smooth=False
sun.smoothmethod='kriging' # USe kriging or idw for smoothing
sun.smoothnear=4 # No. of points to use for smoothing


# option to adjust channel depths using a shapefile
sun.adjust_depths=False
sun.channel_shpfile=''


####
# Open boundary options
####
sun.opt_bcseg = 'file' # Segment boundary condition option: 'constant' or 'file'
sun.opt_bctype2 = 'file' # Type 2 boundary condition option: 'constant'
sun.opt_bctype3 = 'ROMSOTISFILE' # Type 3 boundary condition option: 'constant',
#,'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'ROMSOTISFILE'

sun.bcpolygonfile = '' # Shape file with fields 'marker' and 'edge_id'
sun.bcfile = '' # Input boundary condition file
    
# IF opt_bcseg = 'consant
sun.Q0 = 0.0 # m3/s

# IF opt_bctype2/opt_bctype3 = 'constant'
sun.T0 = 4.0 # Open boundary and background temperature
sun.S0 = 0.0 # Open boundary and background salinity

# IF opt_bctype2 = 'file' or 'ROMSFILE'
#sun.waterlevelstationID = 8772447 # USCG Freeport gage
sun.waterlevelstationID = 8771450 # USCG Freeport gage

# IF opt_bytype2 = 'filt'
sun.TairstationID = '722436-12906' # Houston-Ellington weather station

sun.filttype='low'
sun.cutoff = 7200.0

# Option to zero ROMS uv and eta at boundaries and initial conditions
sun.useROMSuv=False
sun.useROMSeta=False

# Option to use OTS velocities along bounaries
sun.useOTISuv=False

####
# Initial condition options
####
sun.opt_ic = 'constant' #'constant', 'depth_profile', 'ROMS'

sun.icfile = 'LakeMichigan_IC_2012.nc'

sun.T0ic = 4.0 # Initial condition and background temperature
sun.S0ic = 0.0 # Initial condition and background salinity



sun.icfilterdx = None # spatial filter length scale for smoothing ic's

sun.agesourcepoly = None

###
# Metfile
### 

# For plotting only
sun.metfile = 'Galveston_NARR_2010.nc'
###
# Input file names
### 
sun.romsfile = ['../DATA/txla_subset_HIS_2009.nc']
sun.otisfile = '../DATA/Tides/Model_Mex'
sun.dbasefile = '../DATA/GalvestonObs.db'

######################
# Now call the class...
######################

sun(sunpath,starttime,endtime,dt)


###
# Dump figures of the input data
#dump = dumpinputs(suntanspath=sunpath,bcfile=sun.bcfile)
##
#dump(plotdir)
