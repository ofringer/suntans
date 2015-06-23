"""
Download NARR model data
"""
from createSunMetFile import narr2suntans

bbox = [-88.2, -84.5, 41.5,46.2]
utmzone = 16
tstart = '20120301'
tend = '20120308'
ncfile = 'rundata/LakeMichigan_NARR_2012.nc'

narr2suntans(ncfile,tstart,tend,bbox,utmzone)
