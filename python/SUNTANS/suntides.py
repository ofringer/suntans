# -*- coding: utf-8 -*-
"""
Classes and tools for handling tidal harmonic data

Created on Mon Feb 11 10:01:04 2013

@author: mrayson
"""

import numpy as np
import netcdfio
from netCDF4 import Dataset
from sunpy import Spatial, unsurf
import uspectra
from datetime import datetime
import matplotlib.pyplot as plt

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
        
        # Get the tidal fruequencies
        if self.frqnames == None:
			# This returns the default frequencies from the uspectra class
            self.frq,self.frqnames = uspectra.getTideFreq(Fin=None)
        else:
            self.frq,self.frqnames = uspectra.getTideFreq(Fin=self.frqnames)
            
        self.Ntide = len(self.frqnames)
        
        self.reftime = datetime(self.baseyear,1,1)
        
        Spatial.__init__(self,ncfile,**kwargs)
        
        self.Nt = len(self.time)
        
    def __call__(self,tsteps):
        """
        Actually does the harmonic calculation for the model time steps in tsteps
        (or at least calls the class that does the calculation)
        
        Set tsteps = -1 to do all steps
        """
        if tsteps == -1:
            self.tstep=np.arange(0,self.Nt,1)
        else:
            self.tstep=tsteps
        
        # Load all of the time steps into memory
        self.loadData()
        
        # Initialize the amplitude and phase arrays
        self.Amp = np.zeros((self.Nc,self.Ntide))
        self.Phs = np.zeros((self.Nc,self.Ntide))
        
        # Initialise the harmonic object. This object does all of the work...
        t=self.time[self.tstep]
        U = uspectra.uspectra(t,self.data[:,0],frq=self.frq,method='lsqfast')
        
        # Loop through all grid cells
        print
        printstep = 5 
        printstep0 = 0
        print 'Calculating tidal harmonics for variable: "%s"\nTimerange: %s to %s (%d time points)'%\
        (self.variable,datetime.strftime(t[0],'%Y-%m-%d'),datetime.strftime(t[-1],'%Y-%m-%d'),len(self.tstep))
        for ii in range(self.Nc):
            perccomplete = float(ii)/float(self.Nc)*100.0
            if perccomplete > printstep0:
                print '%d %% complete...'%(int(perccomplete))
                printstep0+=printstep
				
            # Update the 'y' data in the grid class. This invokes a harmonic fit    
            U['y'] = self.data[:,ii]
			
			# Calculate the phase and amplitude. Note that the phase is offset to self.reftime
            self.Amp[ii,:],self.Phs[ii,:] = U.phsamp(phsbase=self.reftime)
			
        print '100% complete.'


    def plotAmp(self,con='M2',xlims=None,ylims=None,**kwargs):
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
            
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,self.Amp[:,ii],xlim=xlims,ylim=ylims,\
                clim=self.clim,**kwargs)
        
        titlestr='%s Amplitude\n%s [%s]'%(self.frqnames[ii],self.long_name,self.units)
        plt.title(titlestr)
        
    def plotPhs(self,con='M2',phsunits='radians',xlims=None,ylims=None,**kwargs):
        """
        Plots the phase of the constituent on a map
        """
        ii = findCon(con,self.frqnames)
        
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        if phsunits=='radians':
            clim = [0,2*np.pi]
            self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,self.Phs[:,ii],xlim=xlims,ylim=ylims,\
                clim=clim,**kwargs)
                
        elif phsunits=='degrees':
            clim = [0,360.0]
            self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,self.Phs[:,ii]*180/np.pi,xlim=xlims,ylim=ylims,\
                clim=clim,**kwargs)
        
        titlestr='%s Phase\n%s [%s]'%(self.frqnames[ii],self.long_name,phsunits)
        plt.title(titlestr)
        
    def coRangePlot(self, con='M2', clevs=20, phsint=1800.0,xlims=None,ylims=None, **kwargs):
        """
        Contour plot of amplitude and phase on a single plot
        
        phsint - phase contour step (in seconds)
        """
        from matplotlib import tri
        
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        ii = findCon(con,self.frqnames)
        
        # Create a matplotlib triangulation to contour the data       
        t =tri.Triangulation(self.xp,self.yp,self.cells)
        
        # Filled contour of amplitude
        fig = plt.gcf()
        ax = fig.gca()
        
        # Amplitude plot (note that the data must be on nodes for tricontourf)
        V = np.linspace(self.clim[0],self.clim[1],clevs)
        camp = plt.tricontourf(t, self.cell2node(self.Amp[:,ii]), V, **kwargs)
        
        # Phase Plot
        Vphs = np.arange(0,2*np.pi,self.frq[ii]*phsint) # Phase contour lines
        cphs = plt.tricontour(t, self.cell2node(self.Phs[:,ii]), Vphs,colors='k',linewidths=2.0)
                
        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        
        axcb = fig.colorbar(camp)
        
        #titlestr='%s Amplitude\n%s [%s]\nPhase contours interval: %3.1f hr'%(self.frqnames[ii],self.long_name,self.units,phsint/3600.)
        titlestr='%s Amplitude\nPhase contour interval: %3.1f hr'%(self.frqnames[ii],phsint/3600.)
        plt.title(titlestr)
        
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