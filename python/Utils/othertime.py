# -*- coding: utf-8 -*-
"""
"Other" time functions that are not part of the datetime module

Created on Fri Dec 07 15:39:15 2012

@author: mrayson
"""

from datetime import datetime,timedelta
import numpy as np

def SecondsSince(timein,basetime = datetime(1990,1,1)):
    """
    Converts a list or array of datetime object into an array of seconds since "basetime"
    
    Useful for interpolation and storing in netcdf format
    """
    timeout=[]
    for t in timein:
        dt = t - basetime
        timeout.append(dt.total_seconds())
        
    return np.asarray(timeout)

def datenum2datetime(datenum):
    """
    Convert a matlab datenum float type to a python datetime object
    """
    datenum-=366.0 # There is a leap year difference in the reference time
    day = np.floor(datenum)
    hms = datenum - day
    tday = datetime.fromordinal(int(day))
    
    return tday+timedelta(days=hms)
    
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