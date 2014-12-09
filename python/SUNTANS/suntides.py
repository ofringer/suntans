#!/usr/bin/python
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
from timeseries import timeseries, harmonic_fit, ap2ep
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
        
        if self.hasVar('eta_amp'):
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
            if vv in ['ubar','vbar']:
                ndim=2
                if vv=='ubar': self.variable='uc'
                if vv=='vbar': self.variable='vc'
            else:
                ndim = self._returnDim(vv)
                self.variable=vv
                
            if ndim  == 2 or self.Nkmax==1:
                print 'Loading data from %s...'%vv
                if vv in ['ubar','vbar']:
                    data=self.loadDataBar()
                else:
                    data=self.loadData()
                print 'Performing harmonic fit on variable, %s...'%(self.variable)
                self.Amp[vv], self.Phs[vv], self.Mean[vv] = harmonic_fit(time,data,self.frq,phsbase=self.reftime)
                
            elif ndim == 3:
                for k in range(self.Nkmax):
                    self.klayer=[k]
                    print 'Loading data...'
                    data = self.loadData()
                    print 'Performing harmonic fit on variable, %s, layer = %d of %d...'%(self.variable,self.klayer[0],self.Nkmax)
                    self.Amp[vv][:,k,:], self.Phs[vv][:,k,:], self.Mean[vv][k,:] = harmonic_fit(time,data,self.frq,phsbase=self.reftime)
        
    def _prepDict(self,varnames):
        """
        Prepare the output dictionary
        """
        self.Amp = {}
        self.Phs = {}
        self.Mean = {}
        for vv in varnames:
            if vv in ['ubar','vbar']:
                ndim=2
            else:
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
        
        if varname in ['ubar','vbar']:
            return 2
        else:
            return self.nc.variables[varname].ndim 
        
    def _loadVars(self):
        """
        
        """
        self.frq = self.nc.variables['omega'][:]
        self.frqnames = self.nc.Constituent_Names.split()
        
        varlist = ['eta','uc','vc','w','T','S','rho','ubar','vbar']
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
                
        
    def plotAmp(self,vname='eta',k=0,con='M2',xlims=None,ylims=None,\
        barotropic=False,cbarpos=None,**kwargs):
        """
        Plots the amplitude of the constituent on a map

        vname: one of 'eta','ubar','vbar','uc','vc',...

        Set vname = 'ellmaj','ellmin','ellang','ellphs' to plot ellipse
        properties
        """
        ii = findCon(con,self.frqnames)
        
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(self.Amp))
            self.clim.append(np.max(self.Amp))
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims
            
        if vname in ['ellmaj','ellmin','ellang','ellphs']:
            ell=self.getEllipse(barotropic=barotropic,k=k,con=con)
            if vname=='ellmaj':
                zA=ell[0]
            elif vname=='ellmin':
                zA=ell[1]
            elif vname=='ellang':
                zA=ell[2]
            elif vname=='ellphs':
                zA=ell[3]

        elif len(self.Amp[vname].shape)==2:
            zA = self.Amp[vname][ii,:]
        else:
            zA = self.Amp[vname][ii,k,:].ravel()
            
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,zA,xlim=xlims,ylim=ylims,\
                clim=self.clim,colorbar=False,**kwargs)
        
        # Add a decent looking colorbar
        if not cbarpos==None:
            cbaxes = self.fig.add_axes(cbarpos) 
            cb = self.fig.colorbar(self.patches,cax = cbaxes,orientation='vertical')  
            #cb.ax.set_title('[$m s^{-1}$]')
    
        plt.sca(self.ax)

        return self.patches
 
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
                
    def getEllipse(self,barotropic=False,k=0,con='M2'):
        """
        Returns the ellipse parameters
        """
        iicon = findCon(con,self.frqnames)
        # Load the data
        if barotropic:
            uA = self.Amp['ubar'][iicon,:]
            uP = self.Phs['ubar'][iicon,:]
            vA = self.Amp['vbar'][iicon,:]
            vP = self.Phs['vbar'][iicon,:]
        else:
            uA = self.Amp['uc'][iicon,k,:]
            uP = self.Phs['uc'][iicon,k,:]
            vA = self.Amp['vc'][iicon,k,:]
            vP = self.Phs['vc'][iicon,k,:]
            
        # Calculate the ellipse parameters
        ell = ap2ep(uA,uP,vA,vP)
        # ell = (major, minor, angles, phase)
        return ell
 

    def plotEllipse(self,barotropic=False,k=0,con='M2',scale=1e4,subsample=4,\
            xlims=None,ylims=None,cbarpos=[0.15, 0.15, 0.03, 0.3],**kwargs):
        """
        Plots tidal ellipses on a map
        """
        from matplotlib.collections import EllipseCollection
        
        plt.ioff()
        fig = plt.gcf()
        ax = fig.gca()
        
        iicon = findCon(con,self.frqnames)
        
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(self.Amp))
            self.clim.append(np.max(self.Amp))
        if xlims==None or ylims==None:
            xlims=self.xlims 
            ylims=self.ylims

        ell = self.getEllipse(barotropic=barotropic,k=k,con=con)
            
        # Create the ellipse collection
        indices = range(0,self.Nc,subsample)
        widths = [ell[0][ii]*scale for ii in indices]
        heights = [ell[1][ii]*scale for ii in indices]
        angles = [ell[2][ii]*180.0/np.pi for ii in indices]
        #angles = [ell[2][ii] for ii in indices]
        offsets = [(self.xv[ii],self.yv[ii]) for ii in indices]

        
        collection = EllipseCollection(widths,heights,angles,units='xy',\
            offsets=offsets, transOffset=ax.transData,**kwargs)
        
        z=ell[0][indices]
        collection.set_array(np.array(z))
        collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        collection.set_edgecolors(collection.to_rgba(np.array(z))) 
        
        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        titlestr='%s - Semi-major Ellipse Amplitude'%(self.frqnames[iicon])
        plt.title(titlestr)
        
        ax.add_collection(collection)
        # Add a decent looking colorbar
        if not cbarpos==None:
            cbaxes = fig.add_axes(cbarpos) 
            cb = fig.colorbar(collection,cax = cbaxes,orientation='vertical')  
            cb.ax.set_title('[m s$^{-1}$]')
    
        plt.sca(ax)
        
        #axcb = fig.colorbar(collection)
        
        return collection
   
        
    def coRangePlot(self,vname ='eta', k=0, con='M2', clevs=20, phsint=1800.0,\
        cbarpos=[0.15, 0.15, 0.03, 0.3],xlims=None,ylims=None, **kwargs):
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
        #t =tri.Triangulation(self.xp,self.yp,self.cells)
        
        # Filled contour of amplitude
        fig = plt.gcf()
        ax = fig.gca()

        # Need to do this to avoid sending different keys to contourf
        if kwargs.has_key('titlestr'):
            titlestr=kwargs.pop('titlestr')
        else:
            titlestr=''
        if kwargs.has_key('colorbar'):
            colorbar=kwargs.pop('colorbar')
        else:
            colorbar=None
        
        # Amplitude plot (note that the data must be on nodes for tricontourf)
        V = np.linspace(self.clim[0],self.clim[1],clevs)
        #camp = plt.tricontourf(t, self.cell2node(zA), V, **kwargs)
        camp = self.contourf(z=zA,clevs=V,colorbar=colorbar,titlestr=titlestr,**kwargs)
        
        # Phase Plot
        Vphs = np.arange(0,2*np.pi,self.frq[ii]*phsint) # Phase contour lines
        #cphs = plt.tricontour(t, self.cell2node(zP), Vphs,colors='k',linewidths=2.0)
        cphs = self.contourf(z=zP,clevs=Vphs,filled=False,\
            colorbar=None,titlestr='',
            colors='k',linewidths=2.0)
                
        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        
        #axcb = fig.colorbar(camp)
        # Add a decent looking colorbar
        if not cbarpos==None:
            cbaxes = fig.add_axes(cbarpos) 
            cb = fig.colorbar(camp,cax = cbaxes,orientation='vertical')  
            cb.ax.set_title('[m]')
    
        plt.sca(ax)
 
        
        #titlestr='%s Amplitude\n%s [%s]\nPhase contours interval: %3.1f hr'%(self.frqnames[ii],self.long_name,self.units,phsint/3600.)
        titlestr='%s Amplitude\nPhase contour interval: %3.1f hr'%(self.frqnames[ii],phsint/3600.)
        plt.title(titlestr)
    
    def ustokes(self,barotropic=False,k=0,con='M2'):
        """
        Calculate the Stokes' velocity
        """
        iicon = findCon(con,self.frqnames)
             
        # Load the data
        if barotropic:
            uA = self.Amp['ubar'][iicon,:]
            uP = self.Phs['ubar'][iicon,:]
            vA = self.Amp['vbar'][iicon,:]
            vP = self.Phs['vbar'][iicon,:]
        else:
            uA = self.Amp['uc'][iicon,k,:]
            uP = self.Phs['uc'][iicon,k,:]
            vA = self.Amp['vc'][iicon,k,:]
            vP = self.Phs['vc'][iicon,k,:]
            
        # Calculate the phse gradient
        phiu_x,phiu_y = self.gradH(uP,k=k)
        phiv_x,phiv_y = self.gradH(vP,k=k)
        
        cff = 1.0/(2.0*self.frq[iicon])
        us = phiu_x*uA**2*cff
        vs = phiv_y*vA**2*cff
        
        return us, vs
        
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
                coords = 'omega xv yv'
            elif ndim == 3:
                dims = ('Ntide','Nk','Nc')
                coords = 'omega z_r xv yv'	

            if vv in ['ubar','vbar']:
                units='m s-1'
            else:
                units = self.nc.variables[vv].units

            name = vv+'_amp'
            longname = '%s - harmonic amplitude'%vv
            self.create_nc_var(outfile, name, dims,\
                {'long_name':longname,'units':units,'coordinates':coords},\
                dtype='f8',zlib=1,complevel=1,fill_value=999999.0)
                
            name = vv+'_phs'
            longname = '%s - harmonic phase'%vv
            self.create_nc_var(outfile, name, dims,\
                {'long_name':longname,'units':'radians','coordinates':coords,'reference_time':reftime},\
                dtype='f8',zlib=1,complevel=1,fill_value=999999.0)
                
            ndim = self._returnDim(vv)
            if ndim == 2:
                dims = ('Nc')
                coords = 'omega xv yv'
            elif ndim == 3:
                dims = ('Nk','Nc')
                coords = 'omega z_r xv yv'
                
            name = vv+'_Mean'
            longname = '%s - Temporal mean'%vv
            self.create_nc_var(outfile, name, dims,\
                {'long_name':longname,'units':units,'coordinates':coords},\
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

def usage():
    print "--------------------------------------------------------------"
    print "suntides.py   -h                 # show this help message      "
    print "python suntides.py 'ncfilename.nc' 'outputfile.nc' [-v 'var1 var2 ...] [-f 'K1 O1 M2...']"
    
if __name__ == '__main__':
    """
    Command line call to suntides
    """        
    import getopt, sys
    
    # Defaults
    tstart=-1
    tend = -1
    frqnames=None
    varnames=['eta','ubar','vbar']
    
    try:
        opts,rest = getopt.getopt(sys.argv[3:],'hv:f:')
    except getopt.GetoptError,e:
        print e
        print "-"*80
        usage()
        exit(1)

    for opt,val in opts:
        if opt == '-h':
            usage()
            exit(1)
        elif opt == '-f':
            frqnames=val.split()
        elif opt == '-v':
            varnames = val.split()
    
    try:
	ncfile = sys.argv[1]
	outfile = sys.argv[2]
    except:
    	usage()
	exit(1)
    
    print ncfile, outfile, varnames, frqnames
    # Call the object
    sun=suntides(ncfile,frqnames=frqnames)
    sun(tstart,tend,varnames=varnames)
    sun.tides2nc(outfile)
