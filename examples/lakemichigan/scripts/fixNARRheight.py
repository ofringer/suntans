"""
Fix the height variables in the NARR files
"""

from netCDF4 import Dataset
import numpy as np

def fixHeight(infile):

    nc = Dataset(infile,mode='a')
    nc.variables['z_Tair'][:] = 2.0
    nc.variables['z_Uwind'][:] = 2.0
    nc.variables['z_Vwind'][:] = 2.0
    nc.variables['z_RH'][:] = 2.0
    nc.close()
    print 'Fixed file: %s.'%infile

narrfiles = [\
	#'../DATA/Galveston_NARR_2007.nc',\
	#'../DATA/Galveston_NARR_2008.nc',\
	#'../DATA/Galveston_NARR_2009.nc',\
	#'../DATA/Galveston_NARR_2010.nc',\
	#'../DATA/Galveston_NARR_2011.nc'\
	#'../DATA/Galveston_NARR_2014.nc'\
    'rundata/LakeMichigan_NARR_2012.nc',
	]

for ff in narrfiles:
    fixHeight(ff)
