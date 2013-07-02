# -*- coding: utf-8 -*-
"""
Suntans temporal filtering code

Created on Tue May 14 09:49:38 2013

@author: mrayson
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

from sunpy import Spatial, unsurf
from timeseries import timeseries
from suntans_ugrid import ugrid
import othertime


class sunfilter(Spatial):
    """
    Class for temporally filtering suntans model output
    """
    
    ftype='low'
    order=3
    cutoff_dt = 34.0*3600.0 # Cutoff time period in hours
    
    def __init__(self,ncfile,**kwargs):
        """
        Initialise the suntides class 
        
        See sunpy.Spatial class for list of kwargs
        """
        
        self.__dict__.update(kwargs)
        
        Spatial.__init__(self,ncfile,**kwargs)
        
        self.Nt = len(self.time)
            
    def __call__(self,tstart,tend,varname=None):
        """
        Calls the filter class
        """
        if tstart == -1:
            self.tstep=np.arange(0,self.Nt,1)
        else:
            self.tstep=self.getTstep(tstart,tend)
        
        if not varname == None:
            self.variable = varname
        
        self.loadData()
        
        # Load the data into a time series object (this has a filter method)
        T = timeseries(self.time[self.tstep],np.swapaxes(self.data,0,-1)) # Time dimension needs to be last (filtfilt may have a bug)
                
        dataout = T.filt(self.cutoff_dt,btype=self.ftype,axis=-1,order=self.order)  
        
        return np.swapaxes(dataout,-1,0) # Return with the dimensions in the right order
        
    def filter2nc(self,outfile,tstart,tend,substep=12,varlist=None,**kwargs):
        """
        Filters the variables in the list, varlist, and outputs the results to netcdf
        
        """
        self.__dict__.update(kwargs)
        
        if tstart == -1:
            self.tstep=np.arange(0,self.Nt,1)
        else:
            self.tstep=self.getTstep(tstart,tend)
        
        if varlist == None:
            varlist = ['eta','uc','vc','w']
            
        # Create the output file 
        self.writeNC(outfile)
        
        # Create the output variables
        for vv in varlist:
            print 'Creating variable: %s'%vv
            self.create_nc_var(outfile, vv, ugrid[vv]['dimensions'], ugrid[vv]['attributes'],\
                dtype=ugrid[vv]['dtype'],zlib=ugrid[vv]['zlib'],complevel=ugrid[vv]['complevel'],fill_value=ugrid[vv]['fill_value'])
        
        self.create_nc_var(outfile,'time', ugrid['time']['dimensions'], ugrid['time']['attributes'])
        
        # Loop through and filter each variable (do layer by layer on 3D variables for sake of memory)        
        nc = Dataset(outfile,'a')
        
        # Create the time variable first
        nctime = othertime.SecondsSince(self.time[self.tstep[0]:self.tstep[-1]:substep])

        nc.variables['time'][:] = nctime
        
        for vv in varlist:
            print 'Filtering variable: %s'%vv
            
            if len(ugrid[vv]['dimensions']) == 2:
                dataf = self.__call__(tstart,tend,varname=vv)
                nc.variables[vv][:] = dataf[::substep,:].copy()
            elif len(ugrid[vv]['dimensions']) == 3:
                for kk in range(0,self.Nkmax):
                    print '   layer: %d'%kk
                    self.klayer = [kk]
                    dataf = self.__call__(tstart,tend,varname=vv)
                    nc.variables[vv][:,kk,:] = dataf[::substep,:].copy()


        nc.close()
        print '#####\nComplete - Filtered data written to: \n%s \n#####'%outfile

def usage():
    print "--------------------------------------------------------------"
    print "sunfilter.py   -h                 # show this help message      "
    print "python sunfilter.py 'ncfilename.nc' 'outputfile.nc' [-v 'var1 var2 ...] [-s 12]"
        
if __name__ == '__main__':
    """
    Command line call to sunfilter
    """        
    import getopt, sys
    
    # Defaults
    tstart=-1
    tend = -1
    varnames=['eta','uc','vc','salt','temp','rho']
    substep=12
    
    try:
        opts,rest = getopt.getopt(sys.argv[3:],'hv:s:')
    except getopt.GetoptError,e:
        print e
        print "-"*80
        usage()
        exit(1)

    for opt,val in opts:
        if opt == '-h':
            usage()
            exit(1)
        elif opt == '-s':
            frqnames=int(val)
        elif opt == '-v':
            varnames = val.split()
    
    try:
	ncfile = sys.argv[1]
	outfile = sys.argv[2]
    except:
    	usage()
	exit(1)
    
    print ncfile, outfile, varnames, substep
    # Call the object
    sun=sunfilter(ncfile)
    sun.filter2nc(outfile,tstart,tend,substep=substep,varlist=varnames)