# -*- coding: utf-8 -*-
"""
Classes and tools for handling tidal harmonic data

Created on Mon Feb 11 10:01:04 2013

@author: mrayson
"""

import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import matplotlib.pyplot as plt

import netcdfio
from sunpy import Spatial, unsurf
import uspectra
from timeseries import timeseries, harmonic_fit
from suntans_ugrid import ugrid
import othertime

import pdb

class suntides(Spatial):
    """
    Class for calculating tidal harmonics from SUNTANS model output
    """    
    
    frqnames = None
    baseyear = 1990 # All phases are referenced to the 1st of the 1st of this year
    
    def __init__(self,ncfile,**kwargs):
        """
        Initialise the suntides class 
        
        See sunpy.Spatial class for list of kwargs
        """
        
        self.__dict__.update(kwargs)
        
        Spatial.__init__(self,ncfile,**kwargs)
        
        if self.hasDim('Ntide'):
            print 'Loading existing harmonic data...'
            self._loadVars()
            
        else:
            # Get the tidal fruequencies
            if self.frqnames == None:
    			# This returns the default frequencies from the uspectra class
                self.frq,self.frqnames = uspectra.getTideFreq(Fin=None)
            else:
                self.frq,self.frqnames = uspectra.getTideFreq(Fin=self.frqnames)
                
            self.Ntide = len(self.frqnames)
            
            self.reftime = datetime(self.baseyear,1,1)
                
            self.Nt = len(self.time)
        
    def __call__(self,tstart,tend,varnames=['eta','uc','vc']):
        """
        Actually does the harmonic calculation for the model time steps in tsteps
        (or at least calls the class that does the calculation)
        
        Set tstart = -1 to do all steps
        """
        if tstart == -1:
            self.tstep=np.arange(0,self.Nt,1)
        else:
            self.tstep=self.getTstep(tstart,tend)
        
        time = othertime.SecondsSince(self.time[self.tstep])
        
        self.varnames=varnames
        self._prepDict(varnames)
        
        for vv in varnames:
            ndim = self._returnDim(vv)
            self.variable=vv
            if ndim  == 2 or self.Nkmax==1:
                print 'Loading data...'
                self.loadData()
                print 'Performing harmonic fit on variable, %s...'%(self.variable)
                self.Amp[vv], self.Phs[vv], self.Mean[vv] = harmonic_fit(time,self.data,self.frq,phsbase=self.reftime)
                
            elif ndim == 3:
                for k in range(self.Nkmax):
                    self.klayer=[k]
                    print 'Loading data...'
                    self.loadData()
                    print 'Performing harmonic fit on variable, %s, layer = %d of %d...'%(self.variable,self.klayer[0],self.Nkmax)
                    self.Amp[vv][:,k,:], self.Phs[vv][:,k,:], self.Mean[vv][k,:] = harmonic_fit(time,self.data,self.frq,phsbase=self.reftime)
        
    def _prepDict(self,varnames):
        """
        Prepare the output dictionary
        """
        self.Amp = {}
        self.Phs = {}
        self.Mean = {}
        for vv in varnames:
            ndim = self._returnDim(vv)
            if ndim == 2:
                self.Amp.update({vv:np.zeros((self.Ntide,self.Nc))})
                self.Phs.update({vv:np.zeros((self.Ntide,self.Nc))})
                self.Mean.update({vv:np.zeros((self.Nc,))})
            elif ndim == 3:
                self.Amp.update({vv:np.zeros((self.Ntide,self.Nkmax,self.Nc))})
                self.Phs.update({vv:np.zeros((self.Ntide,self.Nkmax,self.Nc))})
                self.Mean.update({vv:np.zeros((self.Nkmax,self.Nc))})
    
    def _returnDim(self,varname):
        
        return self.nc.variables[varname].ndim 
        
    def _loadVars(self):
        """
        
        """
        self.frq = self.nc.variables['omega'][:]
        self.frqnames = self.nc.Constituent_Names.split()
        
        varlist = ['eta','uc','vc','w','temp','salt','rho']
        self.Amp={}
        self.Phs={}
        self.Mean={}
        for vv in varlist:
            name = vv+'_amp'
            if self.hasVar(name):
                print 'Loading %s harmonic data...'%(vv)
                self.Amp.update({vv:self.nc.variables[name][:]})
                
            name = vv+'_phs'
            if self.hasVar(name):
                self.Phs.update({vv:self.nc.variables[name][:]})
                
            name = vv+'_Mean'
            if self.hasVar(name):
                self.Mean.update({vv:self.nc.variables[name][:]})
                
        
    def plotAmp(self,vname='eta',k=0,con='M2',xlims=None,ylims=None,**kwargs):
        """
        Plots the amplitude of the constituent on a map
        """
        ii = findCon(con,self.frqnames)
        
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(self.Amp))
            self.clim.append(np.max(self.Amp))
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        if len(self.Amp[vname].shape)==2:
            zA = self.Amp[vname][ii,:]
        else:
            zA = self.Amp[vname][ii,k,:].ravel()
            
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,zA,xlim=xlims,ylim=ylims,\
                clim=self.clim,**kwargs)
        
        #titlestr='%s Amplitude\n%s [%s]'%(self.frqnames[ii],self.long_name,self.units)
        #plt.title(titlestr)
        
    def plotPhs(self,vname='eta',k=0,con='M2',phsunits='radians',xlims=None,ylims=None,**kwargs):
        """
        Plots the phase of the constituent on a map
        """
        ii = findCon(con,self.frqnames)
        
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        if len(self.Phs[vname].shape)==2:
            zP = self.Phs[vname][ii,:]
        else:
            zP = self.Phs[vname][ii,k,:].ravel()
            
        if phsunits=='radians':
            clim = [0,2*np.pi]
            self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,zP,xlim=xlims,ylim=ylims,\
                clim=clim,**kwargs)
                
        elif phsunits=='degrees':
            clim = [0,360.0]
            self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,zP*180/np.pi,xlim=xlims,ylim=ylims,\
                clim=clim,**kwargs)
        
        #titlestr='%s Phase\n%s [%s]'%(self.frqnames[ii],self.long_name,phsunits)
        #plt.title(titlestr)

    def plotMean(self,vname='eta',k=0,xlims=None,ylims=None,**kwargs):
        """
        Plots the amplitude of the constituent on a map
        """
        
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(self.Mean))
            self.clim.append(np.max(self.Mean))
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        if len(self.Amp[vname].shape)==2:
            zA = self.Mean[vname]
        else:
            zA = self.Mean[vname][k,:].ravel()
            
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,zA,xlim=xlims,ylim=ylims,\
                clim=self.clim,**kwargs)        
                
    def coRangePlot(self,vname ='eta', k=0, con='M2', clevs=20, phsint=1800.0,xlims=None,ylims=None, **kwargs):
        """
        Contour plot of amplitude and phase on a single plot
        
        phsint - phase contour step (in seconds)
        """
        from matplotlib import tri
        
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        ii = findCon(con,self.frqnames)
        
        if len(self.Amp[vname].shape)==2:
            zA = self.Amp[vname][ii,:]
            zP = self.Phs[vname][ii,:]
        else:
            zA = self.Amp[vname][ii,k,:].ravel()
            zP = self.Phs[vname][ii,k,:].ravel()
        
        # Create a matplotlib triangulation to contour the data       
        t =tri.Triangulation(self.xp,self.yp,self.cells)
        
        # Filled contour of amplitude
        fig = plt.gcf()
        ax = fig.gca()
        
        # Amplitude plot (note that the data must be on nodes for tricontourf)
        V = np.linspace(self.clim[0],self.clim[1],clevs)
        camp = plt.tricontourf(t, self.cell2node(zA), V, **kwargs)
        
        # Phase Plot
        Vphs = np.arange(0,2*np.pi,self.frq[ii]*phsint) # Phase contour lines
        cphs = plt.tricontour(t, self.cell2node(zP), Vphs,colors='k',linewidths=2.0)
                
        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        
        axcb = fig.colorbar(camp)
        
        #titlestr='%s Amplitude\n%s [%s]\nPhase contours interval: %3.1f hr'%(self.frqnames[ii],self.long_name,self.units,phsint/3600.)
        titlestr='%s Amplitude\nPhase contour interval: %3.1f hr'%(self.frqnames[ii],phsint/3600.)
        plt.title(titlestr)
    
    def tides2nc(self,outfile):
        """
        Saves the tidal harmonic data to netcdf
        """
        
        # Write the grid variables
        self.writeNC(outfile)

        nc = Dataset(outfile,'a')
        
        nc.Title = 'SUNTANS harmonic output'
        nc.Constituent_Names = ' '.join(self.frqnames)
        reftime = datetime.strftime(self.reftime,'%Y-%m-%d %H:%M:%S')
        nc.ReferenceDate = reftime
        nc.SimulationTime = '%s - %s'%(datetime.strftime(self.time[self.tstep[0]],'%Y-%m-%d %H:%M:%S'),datetime.strftime(self.time[self.tstep[-1]],'%Y-%m-%d %H:%M:%S'))
        # Add another dimension
        nc.createDimension('Ntide', self.Ntide)
        
        
        nc.close()
        
        # Create the output variables
        for vv in self.varnames:
            print 'Creating variable: %s'%vv

            ndim = self._returnDim(vv)
            if ndim == 2:
                dims = ('Ntide','Nc')
            elif ndim == 3:
                dims = ('Ntide','Nk','Nc')
                
            name = vv+'_amp'
            longname = '%s - harmonic amplitude'%vv
            self.create_nc_var(outfile, name, dims, {'long_name':longname,'units':self.nc.variables[vv].units},\
                dtype='f8',zlib=1,complevel=1,fill_value=999999.0)
                
            name = vv+'_phs'
            longname = '%s - harmonic phase'%vv
            self.create_nc_var(outfile, name, dims, {'long_name':longname,'units':'radians','reference_time':reftime},\
                dtype='f8',zlib=1,complevel=1,fill_value=999999.0)
                
            ndim = self._returnDim(vv)
            if ndim == 2:
                dims = ('Nc')
            elif ndim == 3:
                dims = ('Nk','Nc')
                
            name = vv+'_Mean'
            longname = '%s - Temporal mean'%vv
            self.create_nc_var(outfile, name, dims, {'long_name':longname,'units':self.nc.variables[vv].units},\
                dtype='f8',zlib=1,complevel=1,fill_value=999999.0)
        
        self.create_nc_var(outfile,'omega', ('Ntide',), {'long_name':'frequency','units':'rad s-1'})
        
        nc = Dataset(outfile,'a')
        nc.variables['omega'][:]=self.frq
        
        for vv in self.varnames:
            name = vv+'_amp'
            nc.variables[name][:]=self.Amp[vv]
            name = vv+'_phs'
            nc.variables[name][:]=self.Phs[vv]
            name = vv+'_Mean'
            nc.variables[name][:]=self.Mean[vv]
        nc.close()        
        
        print 'Completed writing harmonic output to:\n   %s'%outfile

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
            
    def __call__(self,tstart,tend,varname=None):
        """
        Calls the filter class
        """
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
        tstep=self.getTstep(tstart,tend)
        nctime = othertime.SecondsSince(self.time[tstep[0]:tstep[-1]:substep])

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
        
    
def findCon(name,conList):
    """
    Returns the index of a constituent from a list
    """    
    return (i for i, j in enumerate(conList) if j == name).next()


def QueryNC(dbfile,staname=None,yearrange=None,cons=None):
    """
    Query the tidal station data
    
    """
    outvar = ['NetCDF_Filename','NetCDF_GroupID','StationName','StationID']
    tablename = 'observations'
    #condition = "Variable_Name = '%s' and (StationName = '%s' or StationName = '%s' or StationName = '%s')" % (varname,staname1,staname2,staname3 )
    
    # Create the query
    varname1 = 'ssh_amp'
    varname2 = 'ssh_phs'
    if staname == None and yearrange == None:
        condition = "(Variable_Name = '%s' or Variable_Name = '%s')"%(varname1,varname2)
        ydim = None
    elif not staname == None and yearrange == None:
        condition = "(Variable_Name = '%s' or Variable_Name = '%s') and StationName = '%s'"%(varname1,varname2, staname)
        
        ydim = 'year'
    elif not staname == None and not yearrange == None:
        t1 = '%s-01-01 00:00:00'%yearrange[0]
        t2 = '%4d-01-01 00:00:00'%yearrange[1]
        condition = "(Variable_Name = '%s' or Variable_Name = '%s') and StationName = '%s' and time_start > %s and time_start < %s"%(varname1,varname2,staname,t1,t2)
        
        ydim = None
    elif staname == None and not yearrange == None:
        t1 = '%s-01-01 00:00:00'%yearrange[0]
        t2 = '%4d-01-01 00:00:00'%yearrange[1]
        condition = "(Variable_Name = '%s' or Variable_Name = '%s') and time_start > '%s' and time_start < '%s'"%(varname1,varname2,t1,t2)
        ydim = 'station'
        
    data, query = netcdfio.queryNC(dbfile,outvar,tablename,condition,fastmode=True)
    
    # Read the constituents from the netcdf file
    ncfile = query['NetCDF_Filename'][0]
    nc = Dataset(ncfile,'r')
    #names =  nc.__dict__.keys()
    names = nc.Tidal_Constituents.split(', ')
    #print nc.['Tidal Constituents']
    nc.close()
    
    # Find the constituent indices
    if cons == None:
        cons=names
    ind = []
    for nn in cons:
        ind.append((i for i, j in enumerate(names) if j == nn).next())                
     
    # Output the query data into a nicer format 
    #amp = [if dd.has_key('ssh_amp'): dd['ssh_amp'][[1,3,8]].ravel() for dd in data]
    amp = []
    phs = []
    time = []
    lon = []
    lat = []
    StationName=[]
    for ii,dd in enumerate(data):
        if dd.has_key('ssh_amp'): 
            amp.append(dd['ssh_amp'][ind].ravel())
            
            if ydim == 'year':
                time.append(dd['time'][0])
            elif ydim == 'station':
                lon.append(dd['longitude'])
                lat.append(dd['latitude'])
                StationName.append(query['StationName'][ii])
            else:
                lon.append(dd['longitude'])
                lat.append(dd['latitude'])
                StationName.append(query['StationName'][ii])
                time.append(dd['time'][0])

                
    
        if dd.has_key('ssh_phs'): 
            phs.append(dd['ssh_phs'][ind].ravel())
        
    amp = np.array(amp)
    phs = np.array(phs)
    
    # Get the time and station coordinates
    if ydim == 'station':
        time = data[0]['time'][0]
        lon = np.array(lon)
        lat = np.array(lat)
        
    elif ydim == 'year':
        lon = data[0]['longitude']
        lat = data[0]['latitude']
        StationName = query['StationName'][0]
        time = np.array(time)   
        
        
    return amp, phs, time, lon, lat, StationName, cons