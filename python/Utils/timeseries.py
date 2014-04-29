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
from matplotlib.lines import Line2D
from scipy import signal, interpolate
import othertime
from datetime import datetime, timedelta
from uspectra import uspectra, getTideFreq
import operator

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
    
    units=''
    long_name=''
    stationid = ''
    varname = ''
    Z=0.0
    
    def __init__(self,t,y,**kwargs):
        
        self.__dict__.update(kwargs)        
        self.t = t # independent variable (t,x, etc)
        self.y = y # dependent variable
        
        self.shape = self.y.shape
        self.ndim = len(self.shape)
        
        self.tsec = othertime.SecondsSince(self.t,basetime=self.basetime)
        
        self.ny = np.size(self.y)

        # make sure the original data is a masked array
        mask = ~np.isfinite(self.y)
        self.y = np.ma.MaskedArray(self.y,mask=mask)
        
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

    def autocorr(self,normalize=False,axis=-1):
        """
        Autocorrelation calculation
        """

        assert self.isequal,\
             'Data must be equally spaced to perform this function.'

        N = self.ny
        M = int(N)/10

        ymean = self.y.mean(axis=axis)
        y = self.y - ymean
        k = range(1,M)
        tau = np.asarray(k,dtype=np.float)*self.dt

        Cyy = [1./(N-kk) * np.sum(y[...,0:-kk]*y[...,kk::],axis=axis) for kk in k ]

        if normalize:
            return Cyy/y.var(), tau
        else:
            return Cyy ,tau
            
        
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
        
        if not btype == 'band':
            Wn = self.dt/cutoff_dt
        else:
            Wn = [self.dt/co for co in cutoff_dt]
            
        (b, a) = signal.butter(order, Wn, btype=btype, analog=0, output='ba')
        
        # filtfilt only likes to operate along the last axis
        ytmp = np.swapaxes(self.y,-1,axis)
        ytmp = signal.filtfilt(b, a, ytmp, axis=-1)
        return np.swapaxes(ytmp,axis,-1)
        #return signal.filtfilt(b, a, self.y, axis=axis)

    def godinfilt(self,filtwidths=[24,25]):
        """
        Apply the successive Godin-type filter to the data set
        
        The filter widths are defined by filtwidths (hours).
        
        For a 24-25-24 hour filter set filtwidths=[24,25] (default).
        """
        
        if self.isequal==False or self.dt != 3600.:
            # Puts the data onto an hourly matrix
            self._evenly_dist_data(3600.)
            
        ymean = self.running_mean(windowlength=24*3600.)
        self.y = ymean
        ymean = self.running_mean(windowlength=25*3600.)
        self.y = ymean
        ymean = self.running_mean(windowlength=24*3600.)
        self.y = ymean
        
            
#        # Calculate the weights
#        dt = self.dt/3600. # time step in hours
#        nbands = [int(fw/dt) for fw in filtwidths]
#        nb = map(float,nbands)
#        window = nbands[0]*2+nbands[1]-2 # number of points in the window
#        
#        weights = np.zeros((window,))
#        weights[0:nbands[0]]=1/nb[0]
#        weights[nbands[0]-1:nbands[1]+nbands[0]]+=1/nb[1]
#        weights[nbands[0]+nbands[1]-1::]+=1/nb[0]
#        weights/=weights.sum()
#        
#        # Emery and Thompson weights (assumes hourly data)
##        weights = np.zeros((36,))
##        for k in range(12):
##            weights[k] = 1./28800 * (1200 - (12-k)*(13-k) - (12+k)*(13+k) )
##            
##        for k in range(12,36):
##            weights[k] = 1./28800 * (36-k)*(37-k)
##        weights = np.hstack((weights[::-1],weights[1::]))
#        
#        
#        # masked values are included in the sum so we
#        # will just zero them
#        mask = self.y.mask.copy()
#        self.y[self.y.mask]=0.
#        
#        # Apply the weights along the last axis of a strided array
#        ytmp = self.y.copy()
#        ytmp = self._window_matrix(ytmp,window)
#        ytmp *=weights
#        ytmp = ytmp.sum(axis=-1)       
#        
#        # This result needs to be normalized to account for missing data,
#        # this is the same as calculating different weights for each section
#        ntmp= np.ones_like(self.y)
#        ntmp[mask] = 0.
#        norm = self._window_matrix(ntmp,window)
#        norm*=weights
#        norm = norm.sum(axis=-1)
#        
#        ytmp/=norm
#        
#        # Update the y variable
#        self.y = self._update_windowed_data(ytmp,window)
#        self.y.mask=mask
            
    def interp(self,timein,method='linear',timeformat='%Y%m%d.%H%M%S',axis=0):
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
            F = interpolate.interp1d(self.tsec,self.y,kind=method,axis=axis,\
                bounds_error=False,fill_value=0)
        else:
            #mask = np.isnan(self.y) == False
            mask = ~self.y.mask
            F = interpolate.interp1d(self.tsec[mask],self.y[mask],kind=method,axis=axis,\
            bounds_error=False,fill_value=0)
        
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
            frq,frqnames = getTideFreq(Fin=frqnames)
            
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
            
#    def running_mean(self,windowlength=3*86400.0,overlap=12*3600.0,\
#        plot=True,calcrms=True):
    def running_mean(self,windowlength=3*86400.0):
        """
        Running mean and RMS of the time series
        
        windowlength - length of each time window [seconds]
        overlap - overlap between windows [seconds]
        """
        mask = self.y.mask.copy()
        self.y[self.y.mask]=0.
        self.y.mask=mask
        
        windowsize = np.floor(windowlength/self.dt)
        ytmp = self.y.copy()
        ytmp = self._window_matrix(ytmp,windowsize)
        
        weights = 1./windowsize * np.ones((windowsize,))
        ytmp2 = np.sum(ytmp*weights,axis=-1)
        
        # This result needs to be normalized to account for missing data,
        # this is the same as calculating different weights for each section
        ntmp= np.ones_like(self.y)
        ntmp[mask] = 0.
        norm = self._window_matrix(ntmp,windowsize)
        #norm*= weights
        norm = norm.sum(axis=-1)
        norm /= windowsize
        
        ytmp2/=norm
                
        return self._update_windowed_data(ytmp2,windowsize)
        


#        pt1,pt2 = window_index_time(self.t,windowlength,overlap)
#        npt = len(pt1)
#        tmid = []
#        ymean = np.zeros((npt,))
#        yrms = np.zeros((npt,))
#
#        ii=-1
#        for t1,t2 in zip(pt1,pt2):
#            ii+=1
#
#            ymean[ii] = self.y[t1:t2].mean()
#            if calcrms:
#                yrms[ii] = crms(self.tsec[t1:t2],self.y[t1:t2]-ymean[ii])
#            
#            # Return the mid time point
#            ind = int(np.floor(t1 + (t2-t1)/2))
#            tmid.append(self.t[ind])
#   
#        tmid = np.asarray(tmid)
#        
#        if plot:
#            plt.subplot(211)
#            self.plot()
#            plt.plot(tmid,ymean,'r')
#            plt.fill_between(tmid,ymean-yrms,y2=ymean+yrms,color=[0.5,0.5,0.5],alpha=0.5)
#            plt.legend(('Original Signal','Mean','$\pm$ RMS'))
#            
#            ax=plt.subplot(212)
#            plt.fill_between(tmid,yrms,alpha=0.7)
#            plt.xticks(rotation=17)
#            plt.ylabel('RMS')
#            ax.set_xlim([self.t[0],self.t[-1]])
#            
#            
#        return tmid, ymean, yrms     
        
    def despike(self,nstd=4.,windowlength=3*86400.0,overlap=12*3600.0,\
        upper=np.inf,lower=-np.inf,maxdiff=np.inf,fillval=0.):
        """
        Despike time series by replacing any values greater than nstd*std.dev with the
        median of a running window. 
        
        nstd - number of standard deviations outside to replace        
        windowlength - length of each time window [seconds]
        overlap - overlap between windows [seconds]
        lower - lower value bound
        upper - upper value bound
        maxdiff - maximum difference between two points
        """
#        if self.isequal==False:
#            self._evenly_dist_data()
            
        nbad = 0
        
        # Now check the maximum difference
        ydiff = np.zeros_like(self.y)
        ydiff[1::] = np.abs(self.y[1::]-self.y[0:-1])
        ind = ydiff>=maxdiff        
        #self.y[ind]=fillval
        self.y.mask[ind] = True # Mask needs to be set after values are prescribed
        nbad += np.sum(ind)
        
        # First mask any NaN and values outside of bounds
        ind = operator.or_(self.y<=lower,self.y>=upper)
        #self.y[ind]=fillval
        self.y.mask[ind] = True
        nbad += np.sum(ind)
        
        ind =  np.isnan(self.y)
        #self.y[ind]=fillval
        self.y.mask[ind] = True
        
        nbad += np.sum(ind)
        
        # Now calculate the moving median and standard deviation
        windowsize = np.floor(windowlength/self.dt)
        ytmp = self.y.copy()
        ytmp = self._window_matrix(ytmp,windowsize)
        
        ytmp2 = np.mean(ytmp,axis=-1)
        ymean = self._update_windowed_data(ytmp2,windowsize)
        
        #ytmp2= np.std(ytmp,axis=-1)
        ytmp2 = np.apply_along_axis(np.std,-1,ytmp2)
        ystd = self._update_windowed_data(ytmp2,windowsize)
        
        # Mask values outsize of the
        ind = operator.or_(self.y >= ymean + nstd*ystd,\
                self.y <= ymean - nstd*ystd)
        
        #self.y[ind] = ymedian[ind]
        self.y.mask[ind] = True
        
        nbad += np.sum(ind)
        
        if self.VERBOSE:
            print 'Despiked %d points'%nbad
        
        
#        pt1,pt2 = window_index_time(self.t,windowlength,overlap)
#
#        for t1,t2 in zip(pt1,pt2):
#            
#            if t1==t2:
#                continue
#            
#            ytmp = self.y[t1:t2]
#           
#            # Mask any nan's and values outside of the bounds
#            ind = operator.or_(ytmp<=lower,ytmp>=upper)
#            ytmp.mask[ind] = True
#            ind2 =  np.isnan(ytmp)
#            ytmp.mask[ind2] = True
#
#
#            # Now mask any values outside of nstd away from the median
#            ymedian = np.mean(ytmp)
#            ystd = np.std(ytmp)
#            #print ymedian, ytmp.mean(), ystd, type(ymedian)
#            #if type(ymedian) == np.ma.core.MaskedConstant:
#            #    pdb.set_trace()
#            if ~np.any(ytmp.mask==False): # all points are masked
#                self.y[t1:t2]=ytmp.copy()               
#                continue
#
#            ind3 = operator.or_(ytmp >= ymedian + nstd*ystd,\
#                ytmp <= ymedian - nstd*ystd)
#            ytmp.mask[ind3] = True
#            
#            self.y[t1:t2]=ytmp.copy()

        
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
        
    def savetxt(self,txtfile):
        f = open(txtfile,'w')
        
        #for t,v in zip(self.tsec,self.y):
        #    f.write('%10.6f\t%10.6f\n'%(t,v))
        for ii in range(self.y.shape[0]):
            f.write('%s, %10.6f\n'%(datetime.strftime(self.t[ii],'%Y-%m-%d %H:%M:%S'),self.y[ii]))
            
        f.close()

    def copy(self):
        """
        Make a copy of the time-series object in memory
        """
        from copy import deepcopy
        return deepcopy(self)
        
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
            
    def _evenly_dist_data(self,dt):
        """
        Distribute the data onto an evenly spaced array
        
        No interpolation is performed
        """
        if self.VERBOSE:
            print 'inserting the data into an equally-spaced time vector (dt = %f s).'%self.dt
    
        t = self.tsec - self.tsec[0]
        # Put the data onto an evenly spaced, masked array
        tout = np.arange(t[0],t[-1]+dt,dt)
        
        tind = np.searchsorted(tout,t)
        
        yout = np.ma.MaskedArray(np.zeros_like(tout),mask=True)

        yout[tind] = self.y
        
        def updatetime(tsec):
            return timedelta(seconds=tsec) + self.t[0]
            
        self.t = map(updatetime,tout)
        self.y = yout
        
        self.tsec = tout
        
        self.ny = np.size(self.y)
        
        self.isequal = True
        
        self.dt = dt
        
    def _window_matrix(self,y,windowsize):
        """
        Returns the matrix as a strided array so that 'rolling' operations can
        be performed along the last axis
        """
        shape = y.shape[:-1] + (y.shape[-1] - windowsize + 1, windowsize)
        strides = y.strides + (y.strides[-1],)
        
        # The masked values get 
        return np.lib.stride_tricks.as_strided(y, shape=shape, strides=strides)
    
    def _update_windowed_data(self,ytmp,windowsize):
        """
        Re-inserts data that has been windowed into an array
        that is the same size as the original time series
        """
        y = np.zeros_like(self.y)
        indent = (windowsize-np.mod(windowsize,2))/2
        
        if np.mod(windowsize,2)==1:
            y[indent:-indent]=ytmp
        else:
            y[indent-1:-indent]=ytmp
        
        y = np.ma.MaskedArray(y,mask=self.y.mask)
        y.mask[0:indent]=True
        y.mask[-indent:]=True
        
        return y

class ModVsObs(object):
    """
    Class for handling and comparing two time series i.e. model vs observations
    """

    units=' '
    long_name=' '
    stationid = ' '
    varname = ' '
    Z=0.0

    def __init__(self,tmod,ymod,tobs,yobs,**kwargs):
        """
        Inputs:
            tmod,tobs - vector of datetime object
            ymod,yobs - vector of values

        Keywords:
            long_name: string containing variable's name (used for plotting)
            units: string containing variable's units (used for plotting)

        Note that tmod and tobs don't need to be the same length. yobs is
        linearly interpolated onto tmod.
        """
        self.__dict__.update(kwargs)


        # Set the range inclusive of both observation and model result
        time0 = max(tmod[0],tobs[0])
        time1 = min(tmod[-1],tobs[-1])

        if time1 < time0:
            print 'Error - the two datasets have no overlapping period.'
            return None
        
        # Clip both the model and observation to this daterange

        t0 = othertime.findNearest(time0,tmod)
        t1 = othertime.findNearest(time1,tmod)
        TSmod = timeseries(tmod[t0:t1],ymod[t0:t1])

        t0 = othertime.findNearest(time0,tobs)
        t1 = othertime.findNearest(time1,tobs)
        self.TSobs = timeseries(tobs[t0:t1],yobs[t0:t1])

        # Interpolate the observed value onto the model step
        #tobs_i, yobs_i = TSobs.interp(tmod[t0:t1],axis=0)
        #self.TSobs = timeseries(tobs_i, yobs_i)

        # Interpolate the modeled value onto the observation time step
        tmod_i, ymod_i = TSmod.interp(tobs[t0:t1],axis=0)
        self.TSmod = timeseries(tmod_i,ymod_i)

        self.N = self.TSmod.t.shape[0]

        self.calcStats()

    def plot(self,colormod='r',colorobs='b',legend=True,**kwargs):
        """
        Time-series plots of both data sets with labels
        """


        h1 = self.TSmod.plot(color=colormod,**kwargs)

        h2 = plt.plot(self.TSobs.t,self.TSobs.y,color=colorobs,**kwargs)

        plt.ylabel(r'%s [$%s$]'%(self.long_name,self.units)) # Latex formatting

        plt.grid(b=True)

        plt.title('StationID: %s'%self.stationid)

        if legend:
            plt.legend(('Model','Observation'),loc='lower right')

        return h1, h2
        
    def stackplot(self,colormod='r',colorobs='b',scale=None,ax=None,fig=None,**kwargs):
        """
        Stack plot of several time series
        """
        labels = ['z = %1.1f m'%z for z in self.Z.tolist()]
        
        fig,ax,ll = stackplot(self.TSobs.t,self.TSobs.y,ax=ax,fig=fig,\
            scale=scale,units=self.units,labels=labels,color=colorobs,*kwargs)
            
        fig,ax,ll = stackplot(self.TSmod.t,self.TSmod.y,ax=ax,fig=fig,\
            scale=scale,units=self.units,labels=labels,color=colormod,*kwargs)

    def calcStats(self):
        """
        Calculates statistics including:
            moments, RMS, CC, skill, ...
        """
        self.meanObs = np.mean(self.TSobs.y,axis=0)
        self.meanMod = np.mean(self.TSmod.y,axis=0)
        self.stdObs = np.std(self.TSobs.y,axis=0)
        self.stdMod = np.std(self.TSmod.y,axis=0)

        # RMSE
        self.rmse = rms(self.TSobs.y-self.TSmod.y,axis=0)

        # skill
        self.skill = 1.0 - np.sum( (self.TSobs.y-self.TSmod.y)**2.,axis=0) / \
            np.sum( (self.TSobs.y - self.meanObs)**2.,axis=0)

        # Correlation coefficient
        self.cc = 1.0/float(self.N) * np.sum( (self.TSobs.y-self.meanObs) * \
            (self.TSmod.y - self.meanMod) ,axis=0) / (self.stdObs * self.stdMod) 

    def printStats(self,f=None,header=True):
        """
        Prints the statistics to a markdown language style table
        """
        outstr=''

        if header:
            outstr += "|      | Mean Model | Mean Obs. | Std. Mod. | Std Obs | RMSE |   CC   | skill |\n"
            
            outstr += "|------| ---------- | --------- | --------- | ------- | --- | ----- | ------| \n"

        outstr += "| %s [%s] | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | \n"%\
            (self.stationid,self.units, self.meanMod, self.meanObs, self.stdMod, self.stdObs,\
            self.rmse,self.cc,self.skill)

        if f == None:
            print outstr
        else:
            f.write(outstr)


    def printStatsZ(self,f=None,header=True):
        """
        Prints the statistics to a markdown language style table
        """
        outstr=''

        if header:
            outstr += "| Depth | Mean Model | Mean Obs. | Std. Mod. | Std Obs | RMSE |   CC   | skill |\n"
            
            outstr += "|------| ---------- | --------- | --------- | ------- | --- | ----- | ------| \n"

        for ii,zz in enumerate(self.Z.tolist()):

            outstr += "| %3.1f [m] | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | \n"%\
            (zz, self.meanMod[ii], self.meanObs[ii], self.stdMod[ii],\
                self.stdObs[ii], self.rmse[ii],self.cc[ii],self.skill[ii])

        if f == None:
            print outstr
        else:
            f.write(outstr)

    def crosscorr(self,normalize=False,axis=-1):
        """
        Crosscorrelation calculation
        """

        assert self.TSobs.isequal,\
             'Data must be equally spaced to perform this function.'

        N = self.TSobs.ny
        M = int(N)/10

        ymean = self.TSobs.y.mean(axis=axis)
        y = self.TSobs.y - ymean
        xmean = self.TSmod.y.mean(axis=axis)
        x = self.TSmod.y - xmean

        k = range(1,M)
        tau = np.asarray(k,dtype=np.float)*self.TSobs.dt

        Cxy = [1./(N-kk) * np.sum(y[...,0:-kk]*x[...,kk::],axis=axis) for kk in k ]

        if normalize:
            return Cxy/(y.std()*x.std()), tau
        else:
            return Cxy ,tau
 
    def csd(self, plot=True,nbandavg=1,**kwargs):
        """
        Cross spectral density
        
        nbandavg = Number of fft bins to average
        """
        
        if self.isequal==False and self.VERBOSE:
            print 'Warning - time series is unequally spaced consider using lomb-scargle fft'
        

        NFFT = int(2**(np.floor(np.log2(self.ny/nbandavg)))) # Nearest power of 2 to length of data
            
        Pyy,frq = mlab.csd(self.TSobs.y-self.TSobs.y.mean(),Fs=2*np.pi/self.dt,NFFT=NFFT,window=mlab.window_hanning,scale_by_freq=True)
        
        if plot:
            plt.loglog(frq,Pyy,**kwargs)
            plt.xlabel('Freq. [$cycles s^{-1}$]')
            plt.ylabel('PSD')
            plt.grid(b=1)
        
        return Pyy, frq


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

    yout = data[0][varname].squeeze()
    # Zero nan
    yout[np.isnan(yout)] = 0.0
    
    if len(data)==0:
        print '!!! Warning - Did not find any stations matching query. Returning -1 !!!'
        return -1
    else:
        ts = timeseries(data[0]['time'],yout)
        
        
    if not timeinfo==None:
        print 'Interpolating station data between %s and %s\n'%(timeinfo[0],timeinfo[1])
        tnew,ynew = ts.interp((timeinfo[0],timeinfo[1],timeinfo[2]))
        ts = timeseries(tnew,ynew)
        ts.dt = timeinfo[2] # This needs updating
        
    if not filttype==None:
        print '%s-pass filtering output data. Cutoff period = %f [s].'%(filttype,cutoff)
        yfilt = ts.filt(cutoff,btype=filttype,axis=0)
        ts.y = yfilt.copy()
    
    if output_meta:
        if data[0].has_key('elevation'):
            ele = data[0]['elevation']
        else:
            ele = np.array([0.0])
        meta = {'longitude':data[0]['longitude'],'latitude':data[0]['latitude'],'elevation':ele,'StationName':query['StationName'][0]}
        return ts, meta        
    else:
        return ts


    

def stackplot(t,y,scale=None,gap=0.2,ax=None,fig=None,units='',labels=None,**kwargs):
    """
    Vertically stacked time series plot.
    
    Puts all of the time-series into one axes by working out a suitable spacing. 
    
    Inputs:
        y - 2d array [nt,ny] where ny is the number of time series
        t - datetime vector
        
    Returns: 
        fig, ax : figure and axes handles
        ll : plot handles to each line plot [list]
    """
    # Determine the scale factors and the heights of all of the axes
    ny = y.shape[1]
        
    if scale==None:
        scale = np.abs(y).max()
    
    if not labels == None:
        assert len(labels)==ny, ' number of labels (%d) must equal number of layers (%d)'%(len(labels),ny)
        
    # Height of each axes in normalized coordinates
    yheight = 1.0 / (ny + (ny+1.0)*gap)
    
    # Create a new figure
    if fig==None:
        fig=plt.figure()
    else:
        fig = plt.gcf()
    
    if ax == None:
        ax = fig.add_subplot(111,frame_on=False,ylim=[0,1.0],yticks=[])
        
    # Now add each line to the figure
    ll = [] # List of line objects
    
    def fakeaxes(yval,dy):
        cc=[0.5,0.5,0.5]
        ax.add_line(Line2D([0,1],[yval,yval],linewidth=0.5,color=cc,transform=ax.transAxes,linestyle='--'))
        yp = yval + dy/2.
        ym = yval - dy/2.
        ax.add_line(Line2D([0,0],[yp,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([1,1],[yp,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        #Little caps
        ax.add_line(Line2D([0,0.01],[yp,yp],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([0,0.01],[ym,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([0.99,1],[yp,yp],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([0.99,1],[ym,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        
    for N in range(1,ny+1):
        yoffset = N*(gap*yheight) + 0.5*yheight + (N-1)*yheight     
        # scaling factor
        vscale = yheight / (scale+yoffset)
        l = ax.plot(t,vscale*y[:,N-1]+yoffset,**kwargs)
        ll.append(l)
        #Adds an axes
        fakeaxes(yoffset,yheight)
        
        if not labels==None:
            plt.text(0.8,yoffset+0.5*yheight-0.02,labels[N-1],transform=ax.transAxes,fontstyle='italic')
              
    # Add a few extra features    
    ax.add_line(Line2D([0,1],[0.01,0.01],linewidth=0.5,color='k',transform=ax.transAxes))
    ax.add_line(Line2D([0,1],[1,1],linewidth=0.5,color='k',transform=ax.transAxes))
    plt.xticks(rotation=17)
    plt.ylabel('$Scale\ =\  %2.1f\  %s$'%(scale,units))
    
    return fig,ax,ll
    
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
    
def window_index_time_slow(t,windowsize,overlap):
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
    
    tsec = othertime.SecondsSince(t)
        
    t1=tsec[0]
    t2=t1 + windowsize
    pt1=[0]
    pt2=[np.searchsorted(tsec,t2)]
    while t2 < tsec[-1]:
        t1 = t2 - overlap
        t2 = t1 + windowsize

        pt1.append(np.searchsorted(tsec,t1))
        pt2.append(np.searchsorted(tsec,t2))
        
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
    
    
    
def rms(x,axis=None):
    """
    root mean squared
    """
    
    return np.sqrt(1.0/np.size(x) * np.sum(x**2,axis=axis))
    
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

def eofsvd(M):
    """
    Compute empirical orthogonal function using singular
    value decomposition technique
    
    Inputs:
        - M : matrix with time along first axis and observation points along
          second
    Returns:
        - PC : The principal component amplitude
        - s : the eigenvalues
        - E : the eigenvectors in each column (EOFs)
    """
    # Remove the mean from the columns
    M = M - M.mean(axis=0)

    # 
    U,s,V = np.linalg.svd(M,full_matrices=False)

    # The principal components are U*s
    PC = U*s

    # Each row of V are the eigenvectors (EOFs) so transpose them so that
    # columns are
    E = V.T

    # Note that the values of s from the svd are the eigenvalues ^0.5

    return PC,s*s,E


def loadtxt(txtfile):
    """
    Loads a text file with two columns
    Column 1 is time: seconds since 1990-1-1
    Column 2 is the data.
    """
    f = open(txtfile,'r')
    
    t=[]
    y=[]
    for line in f.readlines():
        line = line.strip()
        ll = line.split(',')
        t.append(datetime.strptime(ll[0],'%Y-%m-%d %H:%M:%S'))
        y.append(float(ll[1]))
        
    f.close()
    return timeseries(t,np.array(y))

