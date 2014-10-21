# -*- coding: utf-8 -*-
"""
"Other" time functions that are not part of the datetime module

Created on Fri Dec 07 15:39:15 2012

@author: mrayson
"""

from datetime import datetime,timedelta
import numpy as np

import pdb

def TimeVector(starttime,endtime,dt,istimestr=True,timeformat='%Y%m%d.%H%M%S'):
    """
    Create a vector of datetime objects
    """
    # build a list of timesteps
    if istimestr:
        t1 = datetime.strptime(starttime,timeformat)
        t2 = datetime.strptime(endtime,timeformat)
    else:
        t1 = starttime
        t2 = endtime
    
    time = []
    t0=t1
    while t0 <= t2:
        time.append(t0)
        t0 += timedelta(seconds=dt)

    # Make sure the last time step is on the end
    if time[-1]<t2:
        time.append(t2)

    return np.asarray(time)
    
    
def SecondsSince(timein,basetime = datetime(1990,1,1)):
    """
    Converts a list or array of datetime object into an array of seconds since "basetime"
    
    Useful for interpolation and storing in netcdf format
    """
    timeout=[]
    try:
        timein = timein.tolist()
    except:
        timein = timein
    
    try:
        for t in timein:
            dt = t - basetime
            timeout.append(dt.total_seconds())
    except:
        dt = timein - basetime
        timeout.append(dt.total_seconds()) 
        
    return np.asarray(timeout)
    
def MinutesSince(timein,basetime = datetime(1970,1,1)):
    """
    Converts a list or array of datetime object into an array of seconds since "basetime"
    
    Useful for interpolation and storing in netcdf format
    """
    return SecondsSince(timein, basetime = basetime)/60.0

def DaysSince(timein,basetime = datetime(1970,1,1)):
    """
    Converts a list or array of datetime object into an array of seconds since "basetime"
    
    Useful for interpolation and storing in netcdf format
    """
    return SecondsSince(timein, basetime = basetime)/86400.0


def datenum2datetime(timein):
    """
    Convert a matlab datenum float type to a python datetime object
    """
    try:
        timein = timein.tolist()
    except:
        timein = timein
     
    timeout=[]
    for datenum in timein:
        datenum-=366.0 # There is a leap year difference in the reference time
        day = np.floor(datenum)
        hms = datenum - day
        tday = datetime.fromordinal(int(day))
        
        timeout.append(tday+timedelta(days=hms))
        
    return timeout
    
def YearDay(timein):
    """
    Converts a list or array of datetime object into an array of yeardays (floats)
    
    Useful for plotting
    """
    basetime =datetime(timein[0].year,1,1)
    timeout=[]
    for t in timein:
        dt = t - basetime
        timeout.append(dt.total_seconds()/86400.0)
        
    return np.asarray(timeout)
    
def getMonth(timein):
    """
    Converts a list or array of datetime object into an array of months
    
    Useful for calculating monthly distributions
    """

    month=[]
    for t in timein:
        month.append(t.month)
        
    return np.asarray(month)
    
def getYear(timein):
    """
    Converts a list or array of datetime object into an array of years
    
    Useful for calculating monthly distributions
    """

    year=[]
    for t in timein:
        year.append(t.year)
        
    return np.asarray(year)
 
def findNearest(t,timevec):
    """
    Return the index from timevec the nearest time point to time, t. 
    
    """
    tnow = SecondsSince(t)
    tvec = SecondsSince(timevec)
    
    #tdist = np.abs(tnow - tvec)
    
    #idx = np.argwhere(tdist == tdist.min())
    
    #return int(idx[0])
    return np.searchsorted(tvec,tnow)[0]
    
def findGreater(t,timevec):
    """
    Return the index from timevec the first time step greater than t. 
    """
    tnow = SecondsSince(t)
    tvec = SecondsSince(timevec)
    
    idx = np.where(tvec > tnow)
    if len(idx)>0:
        return idx[0][0]
    else:
        return None

def monthlyVector(startyr,endyr,startmth,endmth):
    """
    Returns two datetime vectors with start and end dates at the start and end of months
    """
    month=[]
    year=[]

    m0 = startmth
    y0 = startyr
    while m0  <= endmth or y0 <=  endyr:
        month.append(m0)
        year.append(y0)
        if m0 == 12:
            m0 = 1
            y0 += 1
        else:
            m0 += 1
            
    time=[]
    for mm,yy in zip(month,year):
        time.append(datetime(yy,mm,1))
               
    
    return time[0:-1],time[1:]

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
 
