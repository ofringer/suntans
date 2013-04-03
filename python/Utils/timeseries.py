# -*- coding: utf-8 -*-
"""
Collection of tools for plotting and analysis of time series data

Created on Wed Feb 06 12:37:17 2013
Author: Matt Rayson
Stanford University
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from scipy import signal, interpolate
import othertime
from datetime import datetime, timedelta
from uspectra import uspectra, getTideFreq

import pdb

class timeseries(object):
    """
    Class for handling time series data
    
    Methods include:
        - Power spectral density (plot)
        - Filtering
        - Interpolation
        - Plotting
    """
    
    basetime = datetime(1900,1,1)
    VERBOSE=False
    
    def __init__(self,t,y,**kwargs):
        
        self.__dict__.update(kwargs)        
        self.t = t # independent variable (t,x, etc)
        self.y = y # dependent variable
        
        self.tsec = othertime.SecondsSince(self.t,basetime=self.basetime)
        
        self.ny = np.size(self.y)
        
        self._checkDT()
        
    def psd(self, plot=True,nbandavg=1,**kwargs):
        """
        Power spectral density
        
        nbandavg = Number of fft bins to average
        """
        
        if self.isequal==False and self.VERBOSE:
            print 'Warning - time series is unequally spaced consider using lomb-scargle fft'
        

        NFFT = int(2**(np.floor(np.log2(self.ny/nbandavg)))) # Nearest power of 2 to length of data
            
        Pyy,frq = mlab.psd(self.y-self.y.mean(),Fs=2*np.pi/self.dt,NFFT=NFFT,window=mlab.window_hanning,scale_by_freq=True)
        
        if plot:
            plt.loglog(frq,Pyy,**kwargs)
            plt.xlabel('Freq. [$cycles s^{-1}$]')
            plt.ylabel('PSD')
            plt.grid(b=1)
        
        return Pyy, frq
        
    def specgram(self, NFFT=256,noverlap=128,plot=True,vv=29, **kwargs):
        """
        Spectrogram plot
        """
        from matplotlib.colors import LogNorm
        
        Pyy,frq,tmid = mlab.specgram(self.y-self.y.mean(),Fs=2*np.pi/self.dt,window=mlab.window_hanning,\
            scale_by_freq=True,noverlap=noverlap,NFFT=NFFT)
            
        if plot==True:
            ax3=plt.gca()
            plt.contourf(tmid*2*np.pi,frq,Pyy,vv,norm=LogNorm())
            ax3.set_yscale('log')
            #ax3.set_xlim((1e-4,1e-2))
            #ax3.set_ylim([tmid.min(),tmid.max()])
            ax3.set_ylabel("$\\omega [rad.s^{-1}]$")
            ax3.set_xticklabels([])
    
        return Pyy,frq,tmid
        
    def filt(self,cutoff_dt, btype='low',order=3,axis=-1):
        """
        Butterworth filter the time series
        
        Inputs:
            cutoff_dt - cuttoff period [seconds]
            btype - 'low' or 'high'
        """
        
        if self.isequal==False and self.VERBOSE:
            print 'Warning - time series is unequally spaced. Use self.interp to interpolate onto an equal grid'
        
        Wn = self.dt/cutoff_dt
        (b, a) = signal.butter(order, Wn, btype=btype, analog=0, output='ba')
        
        return signal.filtfilt(b, a, self.y, axis=axis)
#        return signal.lfilter(b, a, self.y, axis=axis)

        
    def interp(self,tstart,tend,dt,method='linear',timeformat='%Y%m%d.%H%M%S'):
        """
        Interpolate the data onto an equally spaced vector
        
        tstart and tend - string with format 'yyyymmdd.HHMMSS'
        
        method - method passed to interp1d
        """
        
        # Create the time vector
        tnew = othertime.TimeVector(tstart,tend,dt,timeformat =timeformat)
        t = othertime.SecondsSince(tnew,basetime = self.basetime)
        
        F = interpolate.interp1d(self.tsec,self.y,kind=method)
        #F = interpolate.UnivariateSpline(self.tsec,self.y,k=method)
        
        return tnew, F(t)
        
        
    def tidefit(self,frqnames=None,basetime=None):
        """
        Perform a tidal harmonic fit to the data
        
        Returns the amp, phase, frequencies and fitted time series
        """
        
        # Get the tidal fruequencies
        if frqnames == None:
			# This returns the default frequencies from the uspectra class
            frq,frqnames = getTideFreq(Fin=None)
        else:
            frq,frqnames = getTideFreq(Fin=self.frqnames)
            
        # Call the uspectra method
        U = uspectra(self.tsec,self.y,frq=frq,method='lsqfast')
        
        amp,phs = U.phsamp(phsbase=basetime)
        
        return amp, phs, frq, U.invfft()   
        
    def plot(self,angle=17,**kwargs):
        """
        Plot
        
        Rotates date lables
        """
        
        plt.plot(self.t,self.y,**kwargs)
        plt.xticks(rotation=angle)
    
    def save2txt(self,txtfile):
        f = open(txtfile,'w')
        
        for t,v in zip(self.tsec,self.y):
            f.write('%10.6f\t%10.6f\n'%(t,v))
            
        f.close()
        
    def _checkDT(self):
        """
        Check that the time series is equally spaced
        """
        dt = np.diff(self.tsec)
        
        dt_unique = np.unique(dt)
        
        if np.size(dt_unique) == 1:
            self.isequal = True
        else:
            self.isequal = False
            
        self.dt = dt[1]


def harmonic_fit(t,X,frq,axis=0):
    """
    Least-squares harmonic fit on an array, X, with frequencies, frq. 
    
    X - vector [Nt] or array [Nt, (size)]
    t - vector [Nt]
    frq - vector [Ncon]
    
    where, dimension with Nt should correspond to axis = axis.
    """

    t = np.asarray(t)
    
    # Reshape the array sizes
    X = X.swapaxes(0, axis)
    sz = X.shape
    lenX = np.prod(sz[1:])
    
    if not len(t) == sz[0]:
        raise 'length of t (%d) must equal dimension of X (%s)'%(len(t),sz[0])
    
    X = np.reshape(X,(sz[0],lenX))
    
    # Resize the frequency array (this may need to go)
    frq=np.asarray(frq) 
    Nfrq = frq.shape[axis]
    sz0 = frq.shape
    if len(sz0)>1:
        frq = frq.swapaxes(0, axis)
        szf = frq.shape
        lenf = np.prod(szf[1:])
        if not lenf == lenX:
            raise 'the size of the non-varying array X must equal the size of frq'
        
        frq = np.reshape(frq,(szf[0],lenf))
    else:
        frq = np.repeat(frq,lenX)
        frq = np.reshape(frq,(Nfrq,lenX))
     
    # Initialize the amplitude and phase arrays
    Amp = np.zeros((Nfrq,lenX))
    Phs = np.zeros((Nfrq,lenX))
    
    # Initialise the harmonic object. This object does all of the work...
    U = uspectra(t,X[:,0],frq=frq[:,0],method='lsqfast')
    
    for ii in range(0,lenX):
        # Update the 'y' data in the grid class. This invokes a harmonic fit 
        U.frq = frq[:,ii]
        U['y'] = X[:,ii]
			
	   # Calculate the phase and amplitude
        Amp[:,ii], Phs[:,ii] = U.phsamp()
        
    # reshape the array
    Amp = np.reshape(Amp,(Nfrq,)+sz[1:])
    Phs = np.reshape(Phs,(Nfrq,)+sz[1:])
    
    # Output back along the original axis
    return Amp.swapaxes(axis,0), Phs.swapaxes(axis,0)
    
def loadDBstation(dbfile,stationID,varname,timeinfo=None,filttype=None,cutoff=3600.0):
    """
    Load station data from a database file
    
    Inputs:
        dbfile - location of database file
        stationID - Station ID in database
        varname - variable name e.g. 'waterlevel', 'discharge', 'salinity'
        
        timeinfo (optional) - tuple with (starttime,endtime,dt). Format 'yyyymmdd.HHMMSS'
            Use this to interpolate onto a constant time vector
        filttype (optional) - 'low' or 'high' 
            Set this to filter data
            
    Returns:
        timeseries object
        -1 on error
            
    """
    from netcdfio import queryNC
    
    outvar = ['NetCDF_Filename','NetCDF_GroupID','StationName']
    tablename = 'observations'
    #condition = 'Variable_Name = "%s" and StationID = "%s"' % (varname,stationID)
    condition = 'Variable_Name = "%s" and StationID LIKE "%%%s"' % (varname,stationID)
    
    print 'Querying database...'
    print condition
    data, query = queryNC(dbfile,outvar,tablename,condition)  
    
    if len(data)==0:
        print '!!! Warning - Did not find any stations matching query. Returning -1 !!!'
        return -1
    else:
        ts = timeseries(data[0]['time'],data[0][varname].ravel())
        
    if not timeinfo==None:
        print 'Interpolating station data between %s and %s\n'%(timeinfo[0],timeinfo[1])
        tnew,ynew = ts.interp(timeinfo[0],timeinfo[1],timeinfo[2])
        ts = timeseries(tnew,ynew)
        ts.dt = timeinfo[2] # This needs updating
        
    if not filttype==None:
        print '%s-pass filtering output data. Cutoff period = %f [s].'%(filttype,cutoff)
        yfilt = ts.filt(cutoff,btype=filttype)
        ts.y = yfilt.copy()
        
    return ts
    
def monthlyhist(t,y,ylim=0.1,xlabel='',ylabel='',title='',**kwargs):
    """
    Plots 12 histograms on a 6x2 matrix of variable, y, grouped by calendar month
    
    Inputs:
        y - vector of data
        t - vector of datetime objects
        kwargs - keyword arguments for numpy.hist
            
    """
    month = othertime.getMonth(t)
    fig=plt.gcf()
    
    for m in range(1,13):
        
        # Find the values
        ind = np.argwhere(month==m)
        data=y[ind]
        
        ax=plt.subplot(6,2,m)
        if len(data)>0:
            plt.hist(data,**kwargs)
        
        mon=datetime.strftime(datetime(1900,m,1),'%B')
        plt.title(mon)
        plt.ylim([0,ylim]) 
        
    
        if m not in (11,12):         
            ax.set_xticklabels([])
        else:
            plt.xlabel(xlabel)
            
            
        if m not in (1,3,5,7,9,11):
            ax.set_yticklabels([])
        else:
            plt.ylabel(ylabel)
            
        
        #Calc some stats
        textstr = 'Mean: %6.1f\nStd. Dev.: %6.1f\n'%(np.mean(data),np.std(data))
        plt.text(0.5,0.5,textstr,transform=ax.transAxes)
        
        # plot a title
        plt.figtext(0.5,0.95,title,fontsize=14,horizontalalignment='center')
        
    return fig
    
def window_index(serieslength,windowsize,overlap):
    """
    Determines the indices for start and end points of a time series window
    
    Inputs:
        serieslength - length of the vector [int]
        windowsize - length of the window [int]
        overlap - number of overlap points [int]
        
    Returns: pt1,pt2 the start and end indices of each window
    """

    p1=0
    p2=p1 + windowsize
    pt1=[p1]
    pt2=[p2]
    while p2 < serieslength:
        p1 = p2 - overlap
        p2 = min((p1 + windowsize, serieslength))
        pt1.append(p1)
        pt2.append(p2)
        
    return pt1, pt2
    
def window_index_time(t,windowsize,overlap):
    """
    Determines the indices for window start and end points of a time vector
    
    The window does not need to be evenly spaced
    
    Inputs:
        t - list or array of datetime objects
        windowsize - length of the window [seconds]
        overlap - number of overlap points [seconds]
        
    Returns: pt1,pt2 the start and end indices of each window
    """
    
    try:
        t=t.tolist()
    except:
        t=t
        
    t1=t[0]
    t2=t1 + timedelta(seconds=windowsize)
    pt1=[0]
    pt2=[othertime.findNearest(t2,t)]
    while t2 < t[-1]:
        t1 = t2 - timedelta(seconds=overlap)
        t2 = t1 + timedelta(seconds=windowsize)

        pt1.append(othertime.findNearest(t1,t))
        pt2.append(othertime.findNearest(t2,t))
        
    return pt1, pt2