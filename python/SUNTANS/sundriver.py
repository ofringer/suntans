# -*- coding: utf-8 -*-
"""
Driver classes for generating suntans input files

Matt Rayson
Stanford University
April 2013
"""

from sunpy import Grid, Spatial
from sundepths import DepthDriver
from sunboundary import modifyBCmarker, Boundary, InitialCond
from metfile import SunMet, metfile
import timeseries

import numpy as np 
import scipy.io as io
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import pdb

PI = 3.141592653589793

class sundriver(object):
    """
    Driver class for generating SUNTANS input files
    """
    
    # Switches to generate bathymetry, boundary, meteorology and initial condition input files
    makebathy=False
    makebnd=False
    makewinds=False
    makeinitial=False
    
    ###
    # General options
    ###
    # Grid projection variables
    convert2utm=False
    CS='NAD83'
    utmzone=51
    isnorth=False
    vdatum = 'MSL'
    
    # Verical grid options
    Nkmax = 1 # number of layers
    r = 1.01 # vertical stretching parameter
    setconstantdepth=False # Option to set constant depth
    H0 = 10.0 # Constant depth
    
    ###
    # Bathymetry interpolation options
    ###
    depthfile = None
    depthmax=0.1
    interpmethod='idw' # Interpolation method:  'nn', 'idw', 'kriging', 'griddata'
    plottype='mpl' # Type of plot: 'mpl', 'vtk2' or 'vtk3'
    
    # Interpolation options
    NNear=3
    
    # IF interpmethod = 'idw' 
    p = 1.0 #  power for inverse distance weighting
    
    # IF interpmethod = 'kriging' 
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0
    
    # Smoothing options
    smooth=True
    smoothmethod='kriging' # USe kriging or idw for smoothing
    smoothnear=4 # No. of points to use for smoothing


    ####
    # Open boundary options
    ####
    opt_bcseg = 'constant' # Segment boundary condition option: 'constant' or 'file'
    opt_bctype2 = 'constant' # Type 2 boundary condition option: 'constant' or 'file'
    opt_bctype3 = 'constant' # Type 3 boundary condition option: 'constant', ,'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'OTISFILE', 'ROMSOTISFILE'

    modifyedges = False # Option to modify the boundary edges
    bcpolygonfile = None # Shape file with fields 'marker' and 'edge_id'
    bcfile = 'SUNTANS_BC.nc' # Input boundary condition file
        
    # IF opt_bcseg = 'consant
    Q0 = 100.0 # m3/s
    
    # IF opt_bctype2/opt_bctype3 = 'constant'
    T0 = 30 # Open boundary and initial condition background temperature
    S0 = 34 # Open boundary and initial condition background salinity
    
    # IF opt_bctype3 = 'file' or 'ROMSFILE'
    waterlevelstationID = None
    
    # IF opt_bctype3 = 'harmonic'
    amp = 0.25
    omega = 2*PI/(24.0*3600.0)

    # IF opt_bctype2 = 'file'
    TairstatationID = None

    # Filter type (waterlevel)   
    filttype='low'
    cutoff=3600.0
    
    # Air temp cuttoff (bctype2 = file)
    tairfilttype = 'low'
    taircutoff = 14.0*24.0*3600.0

    ####
    # Atmospheric input options
    ####
    opt_met = 'constant' # Met file creation options: 'constant'
    metfile = 'SUNTANS_MetForcing.nc'
    
    # IF opt_met = 'consant'
    Uwind = 0.0
    Vwind = 5.0
    RH = 50.0
    Tair = 30.0
    Pair = 1010.0
    rain = 0.0
    cloud = 0.0
    
    ####
    # Initial condition options
    ####
    opt_ic = 'constant', 'depth_profile', 'ROMS'

    icfile = 'SUNTANS_IC.nc'

    icfilterdx = 0.0 # Filtering length scale
    
    ###
    # Input file names
    ### 
    romsfile = None
    otisfile = None
    dbasefile = None

    # Use ROMS u,v and eta
    useROMSuv=False
    useROMSeta=False

    # Use OTIS u & v
    useOTISuv=False
    
    ############################

    def __init__(self,**kwargs):
        
        self.__dict__.update(kwargs)
        
    def __call__(self,suntanspath,starttime,endtime,dt):

        self.suntanspath = suntanspath
        self.starttime = starttime
        self.endtime = endtime
        self.dt = dt
        
        # Step through and perform each step
        self._makebathy()
        
        if self.makebnd:
            self._makebnd()
        
        if self.makeinitial:
            self._makeinitial()
            
        if self.makewinds:
            self._makewinds()
            
        print '###########\n Completed generating input files. Check text for errors!!!!\n##########'
        
    def _makebathy(self):
        """
        Loads the grid object and interpolates the depths
        """
        # Interpolate the depths onto the grid
        if self.makebathy:
            if self.depthfile == None:
                raise Exception, 'need to set "depthfile" parameter'
            else:
                print 'Interpolation depths onto grid from file:\n%s'%self.depthfile
                
            D = DepthDriver(self.depthfile,interpmethod=self.interpmethod,\
            plottype=self.plottype,NNear=self.NNear,\
            p=self.p,varmodel=self.varmodel,nugget=self.nugget,sill=self.sill,vrange=self.vrange,\
            convert2utm=self.convert2utm,utmzone=self.utmzone,isnorth=self.isnorth,vdatum=self.vdatum,\
            smooth=self.smooth,smoothmethod=self.smoothmethod,smoothnear=self.smoothnear)
            
            D(self.suntanspath,depthmax=self.depthmax)
            
            self.grd = D.grd
            
            print 'SUNTANS depths saved to: %s'%(self.suntanspath+'/depths.dat-voro')
        
        elif self.setconstantdepth:
            print 'Using constant depth (%6.2f m)...'%self.H0
            self.grd = Grid(self.suntanspath)
            self.grd.dv = np.zeros_like(self.grd.xv)
            self.grd.dv[:] = self.H0
        else:
            print 'Loading grid from folder:\n%s'%self.suntanspath
            # Load the grid
            self.grd = Grid(self.suntanspath)
            # Load the depth data into the grid object
            self.grd.loadBathy(self.suntanspath+'/depths.dat-voro')


        zmax = np.abs(self.grd.dv.max())
        
        print 'Calculating vertical grid spacing for Nk = %d, r = %1.3f, %6.3f...'%(self.Nkmax,self.r,zmax)

        # Set up the vertical coordinates
        dz = self.grd.calcVertSpace(self.Nkmax,self.r,zmax)
        self.grd.setDepth(dz)
        
        # Save vertspace.dat
        self.grd.saveVertspace(self.suntanspath+'/vertspace.dat')

        # Write cells.dat and edges.dat to ensure they are in the right format
        print 'Overwriting cells.dat and edges.dat to ensure format consistency.'
        self.grd.saveCells(self.suntanspath+'/cells.dat')
        self.grd.saveEdges(self.suntanspath+'/edges.dat')

    def _makebnd(self):
        """
        Generate boundary condition files
        """
        
        if self.modifyedges:        
            modifyBCmarker(self.suntanspath,self.bcpolygonfile)
    
        #Load the boundary object from the grid
        bnd = Boundary(self.suntanspath,(self.starttime,self.endtime,self.dt))
        
        ###
        # Segment (flux) boundaries
        ###
        if self.opt_bcseg == 'constant':
            print 'Setting %d boundary segments to discharge of %6.3f m3/s'%(bnd.Nseg,self.Q0)
            bnd.boundary_Q[:]=self.Q0
            
        elif self.opt_bcseg == 'file':
            print 'Loading river segment data from file...\n'
            for ii, ID in enumerate(bnd.segp):
                print 'Loading discahrge data for boundary segment (%d) StationID: %d...'%(ii,ID)
                
                ts = timeseries.loadDBstation(self.dbasefile,ID,'discharge',timeinfo=(self.starttime,self.endtime,self.dt),\
                    filttype=self.filttype,cutoff=self.cutoff)
                    
                bnd.boundary_Q[:,ii]=ts.y.copy()
            
        else:
            print 'Unknown option: opt_bcseg = %s. Not setting boundary segment data.'%self.opt_bcseg
            
        ###
        # Type-3 boundaries
        ### 
        self.useROMS = False
        self.useOTIS = False
        self.useFILE = False
        self.useOTISFILE = False

        if self.opt_bctype3=='constant':
            print 'Setting constant type-3 boundary conditions...'  
            print 'Setting salinity = %f, temperature = %f'%(self.S0,self.T0)
            bnd.S[:]=self.S0
            bnd.T[:]=self.T0
            
        elif self.opt_bctype3=='depth_profile':
            print 'Setting type-3 boundary T/S from profile...'  
            
            self.loadTSprofile()
            for ii in range(0,bnd.N3):
                bnd.T[0,:,ii] = self.Tz
                bnd.S[0,:,ii] = self.Sz
        
        elif self.opt_bctype3 in ('ROMS'):
            self.useROMS = True

        elif self.opt_bctype3 in ('OTIS'):
            self.useOTIS = True
            
        elif self.opt_bctype3 in ('file'):
            self.useFILE = True
            
        elif self.opt_bctype3 in ('ROMSOTIS'):
            self.useROMS = True
            self.useOTIS = True
        
        elif self.opt_bctype3 in ('ROMSFILE'):
            self.useROMS = True
            self.useFILE = True
        
        elif self.opt_bctype3 in ('OTISFILE'):
            self.useOTISFILE = True

        elif self.opt_bctype3 in ('ROMSOTISFILE'):
            self.useOTISFILE = True
            self.useROMS = True

        else:
            print 'Unknown option: opt_bctype3 = %s. Not setting type-3 boundaries.'%self.opt_bctype3

            
        if self.useROMS:
            bnd.roms2boundary(self.romsfile,setUV=self.useROMSuv,seth=self.useROMSeta)
            
        if self.useOTIS:
            bnd.otis2boundary(self.otisfile,setUV=self.useOTISuv)

        if self.useOTISFILE:
            bnd.otisfile2boundary(self.otisfile,self.dbasefile,self.waterlevelstationID,setUV=self.useOTISuv)
            
        if self.useFILE:
            ID = self.waterlevelstationID
            print 'Loading waterlevel onto all type-3 points from stationID: %d...'%(ID)
            ts = timeseries.loadDBstation(self.dbasefile,ID,'waterlevel',timeinfo=(self.starttime,self.endtime,self.dt),\
                    filttype=self.filttype,cutoff=self.cutoff)
                    
            for ii in range(bnd.N3):
                bnd.h[:,ii] += ts.y.copy()
            
        ###
        # Type-2 boundaries
        ###
        self.useFILE2 = False

        if self.opt_bctype2 == 'constant':
            print 'Setting constant type-2 boundary conditions...'  
            print 'Setting salinity = %f, temperature = %f'%(self.S0,self.T0)
            bnd.boundary_S[:]=self.S0
            bnd.boundary_T[:]=self.T0
        elif self.opt_bctype2 == 'file':
            print 'Using file for type-2 boundary condition (temperature only)'
            print 'Setting salinity = %f'%(self.S0)
            bnd.boundary_S[:]=self.S0
            self.useFILE2 = True
        else:
            print 'Unknown option: opt_bctype2 = %s. Not setting type-2 boundaries.'%self.opt_bctype3
            
            
        if self.useFILE2:
            ID = self.TairstationID
            print 'Loading air temperature onto all type-2 points from stationID: %s...'%(ID)
            ts = timeseries.loadDBstation(self.dbasefile,ID,'Tair',timeinfo=(self.starttime,self.endtime,self.dt),\
            filttype=self.tairfilttype,cutoff=self.taircutoff)
                    
            for ii in range(bnd.N2):
                for kk in range(bnd.Nk):
                    bnd.boundary_T[:,kk,ii] += ts.y.copy()
            
        # Write to netcdf
        bnd.write2NC(self.suntanspath+'/'+self.bcfile)
            
    def _makeinitial(self):
        """
        Generate initial conditions
        """
        
        # Initialise the class
        IC = InitialCond(self.suntanspath,self.starttime)
    
        if self.opt_ic=='constant':
            print 'Setting constant initial conditions...'  
            print 'Setting salinity = %f, temperature = %f'%(self.S0,self.T0)
            IC.T[:]=self.T0
            IC.S[:]=self.S0
            
        elif self.opt_ic=='depth_profile':
            print 'Setting depth-varying initial conditions...'  
            
            self.loadTSprofile()
            for ii in range(0,IC.Nc):
                IC.T[0,:,ii] = self.Tz
                IC.S[0,:,ii] = self.Sz
                
        elif self.opt_ic=='ROMS':
            print 'Setting initial conditions from ROMS model output...'  
            IC.roms2ic(self.romsfile,setUV=self.useROMSuv,seth=self.useROMSeta,interpmethod='idw',NNear=5,p=2)
                #interpmethod=self.interpmethod,NNear=self.NNear,p=self.p,\
                #varmodel=self.varmodel,nugget=self.nugget,sill=self.sill,\
                #vrange=self.vrange)
            
        else:
            print 'Unknown option: opt_ic = %s. Not setting initial conditions.'%self.opt_ic

    
        # Filter the variables in space
        if self.icfilterdx>0:
            IC.filteric(self.icfilterdx)
        # Write the initial condition file
        IC.writeNC(self.suntanspath+'/'+self.icfile,dv=self.grd.dv)
        
    def _makewinds(self):
        """
        Generate a metfile
        """
        if self.opt_met=='constant':
            print 'Setting constant met forcing with: Uwind = %6.2f\nVwind = %6.2f\nTair = %6.2f\nPair = %6.2f\nRH = %6.2f\nrain = %6.2f\ncloud = %6.2f\n'\
                %(self.Uwind,self.Vwind,self.Tair,self.Pair,self.RH,self.cloud,self.rain)
            xpt = self.grd.xv.mean()
            ypt = self.grd.yv.mean()
            zpt = 10.0
            # Load the met object
            met = SunMet(xpt,ypt,zpt,(self.starttime,self.endtime,self.dt))
                
            # Fill the arrays
            met.Uwind.data[:] = self.Uwind
            met.Vwind.data[:] = self.Vwind
            met.Tair.data[:] = self.Tair
            met.Pair.data[:] = self.Pair
            met.RH.data[:] = self.RH 
            met.cloud.data[:] = self.cloud
            met.rain.data[:] = self.rain
        
            met.write2NC(self.suntanspath+'/'+self.metfile)
        else:
            print 'Unknown option: opt_met=%s. Failed to generate a metfile'%self.opt_met
            
    def loadTSprofile(self):
        """
        Loads a temperature salinity profile from a matfile and interpolates onto the grid depths
        """            
        try:
            data=io.loadmat(self.TSprofilefile)
            z = data['depth'].ravel()
            temp = data['temp_profile'].ravel()
            salt = data['salt_profile'].ravel()
            
            FT = interp1d(z,temp,kind=3)
            self.Tz=FT(self.grd.z_r)
            
            FS = interp1d(z,salt,kind=3)
            self.Sz=FS(self.grd.z_r)
            
        except:
            print 'Warning could not find variables: "depth", "temp_profile" or "salt_profile" in file:\n%s'%self.TSprofilefile
            self.Tz = self.T0
            self.Sz = self.S0
                    
            
class dumpinputs(object):
    """
    Dumps the input files
    """
    
    suntanspath = None
    icfile = None
    bcfile = None
    metfile = None
    
    def __init__(self,**kwargs):
        
        self.__dict__.update(kwargs)
   
    def __call__(self,plotdir):
        
        self.plotdir = plotdir
        
        if not self.bcfile == None:
            print 'Dumping boundary condition data to figures...'
            self._dumpbc()
        if not self.icfile == None:
            print 'Dumping initial condition data to figures...'
            self._dumpic()
            
        if not self.metfile == None:
            print 'Dumping metfile data to figures...'
            self._dumpmet()
            
    def _dumpbc(self):
        """
        Dump boundary condition plots
        
        For each variable (uc,vc,T,S,eta):
            - time series at one point (range of all other points)
            - scatter plots of all points
        """
        
        varnames = ['uc','vc','T','S','h','boundary_u','boundary_v','boundary_T','boundary_S']
        
	print self.suntanspath, self.bcfile
        bnd = Boundary(self.suntanspath+'/'+self.bcfile,0)
        
        self.tmin = bnd.time[0]
        self.tmax = bnd.time[-1]
        
        for vv in varnames:
            plt.figure()
            bnd.plot(varname=vv)
            outfile = '%s/BC_%s_TS.png'%(self.plotdir,vv)
            bnd.savefig(outfile)   
            
        # Boundary segments
        for ii in range(bnd.Nseg):
            fig=plt.figure()
            bnd.plot(varname='boundary_Q',j=ii,rangeplot=False,linewidth=2)
            plt.fill_between(bnd.time,bnd.boundary_Q[:,ii],alpha=0.5)
            outfile = '%s/BC_Flow_%d.png'%(self.plotdir,bnd.segp[ii])
            bnd.savefig(outfile)
            del fig
            
        # Scatter Plots
        for vv in varnames:
            fig=plt.figure()
            bnd.scatter(varname=vv)
            outfile = '%s/BC_%s_scatter.png'%(self.plotdir,vv)
            bnd.savefig(outfile) 
            del fig
            
    def _dumpic(self):
        """
        Dump initial condition plots
        
        For each variable (uc,vc,T,S,eta):
            - surface plot
            - seabed plot
        """
        
        varnames = ['uc','vc','temp','salt','eta']
        
        sun = Spatial(self.suntanspath+'/'+self.icfile,klayer=[0])
        
        # Plot the surface variables
        for vv in varnames:
            h=plt.figure()
            sun.variable=vv
            sun.loadData()
            sun.clim = [sun.data.min(),sun.data.max()]
            sun.plot()            
            outfile = '%s/IC_%s_surface.png'%(self.plotdir,vv)
            sun.savefig(outfile)
            del h
            
        # Plot the seabed variables
        sun.klayer=[-1]
        for vv in varnames:
            h=plt.figure()
            sun.variable=vv
            sun.loadData()
            sun.clim = [sun.data.min(),sun.data.max()]
            sun.plot()
            outfile = '%s/IC_%s_seabed.png'%(self.plotdir,vv)
            sun.savefig(outfile)  
            del h
    
    def _dumpmet(self):
        """
        Dump plots of meteorological data to figures
        """
        print self.suntanspath+'/'+self.metfile
        met = metfile(self.suntanspath+'/'+self.metfile)
        
        varnames = ['Uwind','Vwind','Tair','Pair','RH','rain','cloud']
        
        for vv in varnames:
            plt.figure()
            eval('met.%s.plot()'%vv)
            # Set the x limits to the boundary axis
            if self.__dict__.has_key('tmin'):
                plt.xlim((self.tmin,self.tmax))
                
            outfile = '%s/MET_%s_TS.png'%(self.plotdir,vv)
            plt.savefig(outfile)
            print 'Meteorological data image saved to file:%s'%outfile
            
            
        
        
        
        
