# -*- coding: utf-8 -*-
"""
Collection of tools for plotting and analysis of time series data

Created on Wed Feb 06 12:37:17 2013
Author: Matt Rayson
Stanford University
"""

import numpy as np
import matplotlib.pyplot as plt
import othertime
from datetime import datetime, timedelta
from uspectra import uspectra

import pdb

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