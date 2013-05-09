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
        
        self.shape = self.y.shape
        self.ndim = len(self.shape)
        
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
            plt.contourf(tmid*2*np.pi,frq,Pyy,vv,norm=LogNorm(),**kwargs)
            ax3.set_yscale('log')
            #ax3.set_xlim((1e-4,1e-2))
            #ax3.set_ylim([tmid.min(),tmid.max()])
            ax3.set_ylabel("$\\omega [rad.s^{-1}]$")
            plt.colorbar()
            #ax3.set_xticklabels([])
    
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

        
    def interp(self,timein,method='linear',timeformat='%Y%m%d.%H%M%S'):
        """
        Interpolate the data onto an equally spaced vector
        
        timein is either:
            (3x1 tuple) - (tstart,tend,dt)
                tstart and tend - string with format 'yyyymmdd.HHMMSS'
        or
            datetime vector
        
        method - method passed to interp1d
        """
        
        # Create the time vector
        try:
            tstart=timein[0]
            tend=timein[1]
            dt=timein[2]
            tnew = othertime.TimeVector(tstart,tend,dt,timeformat =timeformat)
        except:
            tnew=timein
            
        t = othertime.SecondsSince(tnew,basetime = self.basetime)
        
        # Don't include nan points
        
        if self.ndim > 1:
            # Interpolate multidimensional arrays without a mask
            F = interpolate.interp1d(self.tsec,self.y.T,kind=method,axis=-1)
        else:
            mask = np.isnan(self.y) == False
            F = interpolate.interp1d(self.tsec[mask],self.y[mask],kind=method,axis=0)
        
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
        
        return amp, phs, frq, frqnames, U.invfft()   
        
    def running_harmonic(self,omega,windowlength=3*86400.0,overlap=12*3600.0, plot=True):
        """
        Running harmonic fit of the time series at frequency, omega. 
        
        windowlength - length of each time window [seconds]
        overlap - overlap between windows [seconds]
        """
        
        # Make sure that omega is a list
        try:
            len(omega)
        except:
            omega=[omega]
        
        pt1,pt2 = window_index_time(self.t,windowlength,overlap)
        npt = len(pt1)
        tmid = []
        amp = np.zeros((npt,))
        phs = np.zeros((npt,))
        ymean = np.zeros((npt,))
        ii=-1
        for t1,t2 in zip(pt1,pt2):
            ii+=1
            # Perform the harmonic fit on the segment
            U = uspectra(self.t[t1:t2],self.y[t1:t2],frq=omega,method='lsqfast')
            # Reference the phase to the start of the time series
            amp[ii],phs[ii] = U.phsamp(phsbase = self.t[0])
            
            # Return the mid time point
            ind = np.floor(t1 + (t2-t1)/2)
            tmid.append(self.t[ind])
            
            # Return the fitted time series
            ymean[ii] = self.y[t1:t2].mean()
            #yout[t1:t2] += ymean + U.invfft()

            
        tmid = np.asarray(tmid)
        

        
        if plot:
            plt.subplot(211)
            self.plot()
            plt.plot(tmid,ymean,'r')
            plt.fill_between(tmid,ymean-amp,y2=ymean+amp,color=[0.5,0.5,0.5],alpha=0.5)
            plt.legend(('Original Signal','Harmonic reconstruction'))
            
            ax=plt.subplot(212)
            plt.fill_between(tmid,amp,alpha=0.5)
            plt.xticks(rotation=17)
            plt.ylabel('Amplitude')
            ax.set_xlim([self.t[0],self.t[-1]])
            
            
        return tmid, amp, phs
            
    def running_mean(self,windowlength=3*86400.0,overlap=12*3600.0, plot=True):
        """
        Running mean and RMS of the time series
        
        windowlength - length of each time window [seconds]
        overlap - overlap between windows [seconds]
        """
        
        # Make sure that omega is a list

        pt1,pt2 = window_index_time(self.t,windowlength,overlap)
        npt = len(pt1)
        tmid = []
        ymean = np.zeros((npt,))
        yrms = np.zeros((npt,))

        ii=-1
        for t1,t2 in zip(pt1,pt2):
            ii+=1

            ymean[ii] = self.y[t1:t2].mean()
            yrms[ii] = crms(self.tsec[t1:t2],self.y[t1:t2]-ymean[ii])
            
            # Return the mid time point
            ind = np.floor(t1 + (t2-t1)/2)
            tmid.append(self.t[ind])
   
        tmid = np.asarray(tmid)
        
        if plot:
            plt.subplot(211)
            self.plot()
            plt.plot(tmid,ymean,'r')
            plt.fill_between(tmid,ymean-yrms,y2=ymean+yrms,color=[0.5,0.5,0.5],alpha=0.5)
            plt.legend(('Original Signal','Mean','$\pm$ RMS'))
            
            ax=plt.subplot(212)
            plt.fill_between(tmid,yrms,alpha=0.7)
            plt.xticks(rotation=17)
            plt.ylabel('RMS')
            ax.set_xlim([self.t[0],self.t[-1]])
            
            
        return tmid, ymean, yrms           
        
    def plot(self,angle=17,**kwargs):
        """
        Plot
        
        Rotates date lables
        """
        
        h1=plt.plot(self.t,self.y,**kwargs)
        plt.xticks(rotation=angle)
        
        return h1 
        
    def fillplot(self,angle=17,alpha=0.7,**kwargs):
        """
        
        """
        h1=plt.fill_between(self.t,self.y,alpha=alpha,**kwargs)
        plt.xticks(rotation=angle)
        
        return h1 
        
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
        
        try:
            self.dt = dt[1]
        except:
            self.dt = 0.0


def harmonic_fit(t,X,frq,mask=None,axis=0,phsbase=None):
    """
    Least-squares harmonic fit on an array, X, with frequencies, frq. 
    
    X - vector [Nt] or array [Nt, (size)]
    t - vector [Nt]
    frq - vector [Ncon]
    mask - array [(size non-time X)]
    phsbase - phase offset
    
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
    
    if not mask == None:
        mask = np.reshape(mask,(lenX,))
    else:
        mask = np.ones((lenX,))
    
    frq = np.array(frq)
    Nfrq = frq.shape[0]
    

    
    def buildA(t,frq):
        """
        Construct matrix A
        """
        nt=t.shape[0]
        nf=frq.shape[0]
        nff=nf*2+1
        A=np.ones((nt,nff))
        for ff in range(0,nf):
            A[:,ff*2+1]=np.cos(frq[ff]*t)
            A[:,ff*2+2]=np.sin(frq[ff]*t)
            
        return A
    
    def lstsqnumpy(A,y):    
        """    
        Solve the least square problem
        
        Return:
            the complex amplitude 
            the mean
        """
        N=A.shape[1]
        b = np.linalg.lstsq(A,y)
        A = b[0][1::2]
        B = b[0][2::2]
        
        return A+1j*B, b[0][0::N]
    
    def phsamp(C):
        return np.abs(C), np.angle(C)
        
    # Least-squares matrix approach
    A = buildA(t,frq)
    C, C0 = lstsqnumpy(A,X) # This works on all columns of X!!
    Amp, Phs= phsamp(C)

    # Reference the phase to some time
    if not phsbase == None:
        base = othertime.SecondsSince(phsbase)
	phsoff = phase_offset(frq,t[0],base)
	phsoff = np.repeat(phsoff.reshape((phsoff.shape[0],1)),lenX,axis=1)
	phs = np.mod(Phs+phsoff,2*np.pi)
    
    # Non-vectorized method (~20x slower)
#    Amp = np.zeros((Nfrq,lenX))
#    Phs = np.zeros((Nfrq,lenX))
#    for ii in range(0,lenX):    
#        if mask[ii]==True: 
#            C = lstsqnumpy(A,X[:,ii])
#            # Calculate the phase and amplitude
#            am, ph= phsamp(C)
#            Amp[:,ii] = am; Phs[:,ii] = ph
            
    
    # reshape the array
    Amp = np.reshape(Amp,(Nfrq,)+sz[1:])
    Phs = np.reshape(Phs,(Nfrq,)+sz[1:])
    C0 = np.reshape(C0,sz[1:])
    
    # Output back along the original axis
    return Amp.swapaxes(axis,0), Phs.swapaxes(axis,0), C0.swapaxes(axis,0)
    
def phase_offset(frq,start,base):
        """
        Compute a phase offset for a given fruequency
        """
        
        if type(start)==datetime:
            dx = start - base
            dx = dx.total_seconds()
        else:
            dx = start -base
        
        return np.mod(dx*np.array(frq),2*np.pi)
 
def loadDBstation(dbfile,stationID,varname,timeinfo=None,filttype=None,cutoff=3600.0,output_meta=False):
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
        ts = timeseries(data[0]['time'],data[0][varname].squeeze())
        
        
    if not timeinfo==None:
        print 'Interpolating station data between %s and %s\n'%(timeinfo[0],timeinfo[1])
        tnew,ynew = ts.interp((timeinfo[0],timeinfo[1],timeinfo[2]))
        ts = timeseries(tnew,ynew)
        ts.dt = timeinfo[2] # This needs updating
        
    if not filttype==None:
        print '%s-pass filtering output data. Cutoff period = %f [s].'%(filttype,cutoff)
        yfilt = ts.filt(cutoff,btype=filttype,axis=-1)
        ts.y = yfilt.copy()
    
    if output_meta:
        if data[0].has_key('elevation'):
            ele = data[0]['elevation']
        else:
            ele = 0.0
        meta = {'longitude':data[0]['longitude'],'latitude':data[0]['latitude'],'elevation':ele,'StationName':query['StationName'][0]}
        return ts, meta        
    else:
        return ts

def SpeedDirPlot(t,u,v,convention='current',units='m s^{-1}',color1='b',color2='r'):
    """
    Plots speed and direction on the same axes
    
    Inputs:
        t - time vector
        u,v - velocity cartesian components
        
    Returns:
        ax - list of axes  handles
        h - list of plot handles
    
    convention = 'current' or 'wind'
    
    See this example:
        http://matplotlib.org/examples/api/two_scales.html
    """
    import airsea
    
    Dir, Spd = airsea.convertUV2SpeedDirn(u,v,convention=convention)
    
    
    ax = range(2)
    h = range(2)
    fig = plt.gcf()
    ax[0] = fig.gca()
    
    
    # Left axes
    h[0] = ax[0].fill_between(t, Spd, color=color1,alpha=0.7)
    # Make the y-axis label and tick labels match the line color.
    ax[0].set_ylabel('Speed [$%s$]'%units, color=color1)
    for tl in ax[0].get_yticklabels():
        tl.set_color(color1)

    #Right axes
    ax[1] = ax[0].twinx() # This sets up the second axes
    ax[1].plot(t, Dir, '.',color=color2)
    ax[1].set_ylabel("Dir'n [$\circ$]", color=color2)
    ax[1].set_ylim([0,360])
    ax[1].set_yticks([0,90,180,270])
    ax[1].set_yticklabels(['N','E','S','W'])
    for tl in ax[1].get_yticklabels():
        tl.set_color(color2)
        
    plt.setp( ax[0].xaxis.get_majorticklabels(), rotation=17 )
        
    return ax, h

def ProfilePlot(t,y,z,scale=86400, axis=0,color=[0.5,0.5,0.5]):
    """
    Plot a series of vertical profiles as a time series
    
    scale - Sets 1 unit = scale (seconds)
    
    See this page on formatting:
        http://matplotlib.org/examples/pylab_examples/date_index_formatter.html
    """
    from matplotlib import collections
    from matplotlib.ticker import Formatter

    class MyFormatter(Formatter):
        def __init__(self, dates, fmt='%b %d %Y'):
            self.fmt = fmt
            self.dates = dates

        def __call__(self, x, pos=0):
            'Return the label for time x s'
            return datetime.strftime(datetime(1990,1,1)+timedelta(seconds=x),self.fmt)

    tsec = othertime.SecondsSince(t)
    formatter = MyFormatter(tsec)
    
    y = np.swapaxes(y,0,axis)
    
    lines=[]
    line2 =[]
    for ii, tt in enumerate(tsec):
        #xplot = set_scale(y[:,ii],tt)
        xplot = tt + y[:,ii]*scale
        lines.append(np.array((xplot,z)).T)
        line2.append(np.array([[tt,tt],[z[0],z[-1]]]).T)
        
    
    LC1 = collections.LineCollection(lines,colors=color,linewidths=1.5)
    LC2 = collections.LineCollection(line2,colors='k',linestyles='dashed') # Zero axis
    
    ax=plt.gca()
    ax.add_collection(LC1)
    ax.add_collection(LC2)
    ax.set_ylim((z.min(),z.max()))
    ax.xaxis.set_major_formatter(formatter)
    ax.set_xlim((tsec[0],tsec[-1]))
    plt.xticks(rotation=17)       
    
    return ax
    
    
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
    
def pol2cart(th,rho):
    """Convert polar coordinates to cartesian"""
    x = rho * np.cos(th)
    y = rho * np.sin(th)

    return x, y
    
def cart2pol(x,y):
    """
    Convert cartesian to polar coordinates
    """
    th = np.angle(x+1j*y)
    rho = np.abs(x+1j*y)
    
    return th, rho
    
def ap2ep(uamp,uphs,vamp,vphs):
    """
    Convert u/v amplitude phase information to tidal ellipses
    
    All angles are in radians
    
    Returns:
        SEMA, SEMI, INC, PHS, ECC
    
    Based on the MATLAB ap2ep function:
        https://www.mathworks.com/matlabcentral/fileexchange/347-tidalellipse/content/ap2ep.m
    """
    # Make complex amplitudes for u and v
    u = uamp*np.exp(-1j*uphs)
    v = vamp*np.exp(-1j*vphs)

    #Calculate complex radius of anticlockwise and clockwise circles:
    wp = (u+1j*v)/2.0     # for anticlockwise circles
    wm = np.conj(u-1j*v)/2.0  # for clockwise circles
    # and their amplitudes and angles
    Wp = np.abs(wp)
    Wm = np.abs(wm)
    THETAp = np.angle(wp)
    THETAm = np.angle(wm)
   
    # calculate ep-parameters (ellipse parameters)
    SEMA = Wp+Wm              # Semi  Major Axis, or maximum speed
    SEMI = Wp-Wm            # Semin Minor Axis, or minimum speed
    ECC = SEMI/SEMA          # Eccentricity

    PHA = (THETAm-THETAp)/2.0   # Phase angle, the time (in angle) when 
                               # the velocity reaches the maximum
    INC = (THETAm+THETAp)/2.0   # Inclination, the angle between the 
                               # semi major axis and x-axis (or u-axis).
                               
    return SEMA, SEMI, INC, PHA, ECC
    
    
    
def rms(x):
    """
    root mean squared
    """
    
    return np.sqrt(1.0/np.size(x) * np.sum(x**2))
    
def crms(t,y):
    """
    Continous function rms
    """
    fac = 1.0/(t[-1]-t[0])
    return np.sqrt(fac*np.trapz(y**2,x=t))
    
def tidalrmse(Ao,Am,Go,Gm):
    """
    Tidal harmonic RMSE
    Ao, Am - observed and modeled amplitude
    Go, Gm - observed and modeled phase (radians)
    """
    return np.sqrt( 0.5*(Ao**2 + Am**2) - Ao*Am*np.cos(Go-Gm) )
