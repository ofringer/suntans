# -*- coding: utf-8 -*-
"""
Update a database with point observation netcdf files
"""

import netcdfio

####
# Inputs
create = False

dbfile = 'C:/Projects/GOMGalveston/DATA/GalvestonObs.db'

#ncfiles = ['C:/Projects/GOMGalveston/DATA/Ocean/USIOOS_OceanObs_20102011.nc',\
#    'C:/Projects/GOMGalveston/DATA/Winds/NCDCNWS_AirObs_20102011.nc',\
#    'C:/Projects/GOMGalveston/DATA/Ocean/IOOS_ADCP_20112012.nc',\
#    'C:/Projects/GOMGalveston/DATA/Ocean/IOOS_Tides_20112012.nc']
    
ncfiles = ['C:/Projects/GOMGalveston/DATA/WaterQuality/HGAC_WaterSampling.nc']
####

if create:
    print 'Creating database: %s'%dbfile
    netcdfio.createObsDB(dbfile)
    
for nc in ncfiles:
    print 'Inserting metadata from: %s'%nc    
    netcdfio.netcdfObs2DB(nc,dbfile)

print 'Done.'