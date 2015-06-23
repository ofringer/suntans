# -*- coding: utf-8 -*-
"""
Example of how to get National Weather Service observations, convert them into a
searchable format, and finally convert to a SUNTANS  meteorological netcdf format

Created on Mon Nov 18 13:50:17 2013

@author: mrayson
"""

import getNOAAWeatherStation as noaa
import netcdfio
from createSunMetFile import interpWeatherStations, write2NC
import os
from datetime import datetime,timedelta

###################################################
# Input variables

# Position of point of interest +/- dx degrees
lon0 = -94.86367
lat0 = 29.29517
utmzone = 15

# Time period to download
timestart = '20090101'
timeend = '20120101'
dt = 1.0/24.0

# Local directory where to store raw NWS-tar.gz data files
localdir = 'data-files'

# Netcdf file and shapefile containing the data and its location
ncfile = 'data-files/NCDCNWS_AirObs_20102011.nc'
shpfile ='data-files/NCDCNWS_AirObs_20102011.shp'

# Database file that stores station metadata
dbfile = 'data-files/WeatherObs.db'

# Name of the suntans netcdf met file
metncfile = 'rundata/Galveston_Met.nc'
###################################################

###
# Step 1: Download the NWS weather files and convert to a netcdf-4 format
###
dx = 0.1
# Set the date range slightly larger just in case
tstart = datetime.strptime(timestart,'%Y%m%d')-timedelta(days=30)
tend = datetime.strptime(timeend,'%Y%m%d')+timedelta(days=30)

latlon = [lon0-dx,lon0+dx,lat0-dx,lat0+dx]

data = noaa.noaaish2nc(latlon,[tstart.year,tend.year],localdir,ncfile,shpfile)

###
# Step 2: Create a database and insert the netcdf metadata into it
###
print 'Creating database: %s'%dbfile
    
if os.path.exists(dbfile):
    print 'Overwriting old file...'
    os.unlink(dbfile)
    
netcdfio.createObsDB(dbfile)
    
ncfiles=[ncfile]
for nc in ncfiles:
    print 'Inserting metadata from: %s'%nc    
    netcdfio.netcdfObs2DB(nc,dbfile)

###
# Step 3: Create a SUNTANS format meteorological input file
###
[coords, output, nctime] = interpWeatherStations(latlon,timestart,timeend,dt*24.0,utmzone,dbfile,maxgap=100,showplot=False)
# Write to the suntans meteorological netcdf format
write2NC(metncfile,coords,output,nctime)
