# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 10:18:47 2012

@author: mrayson
"""
import time
from datetime import datetime, timedelta

def parse(tstr):
    """ Parse a date represented by a string to a decimal number.
    The time/datetime functions do not seem to support the format spat out by 
    database
    """
    # convert the string to a list so that it can be modified
    tlst = list(tstr)
    
    
    day = tstr[4:6]
    hour = tstr[12:14]
    dd = int(day)
    hh = int(hour)
    # Replace days
    if dd < 10:
        #tout.replace(day,"0"+day[1],1)
        tlst[4:6]="0"+day[1]
    
    # replace hours
    if hh < 10:
        tlst[12:14]="0"+hour[1]
    
    # Combine back into a string
    ttmp = "".join(tlst)

    # Convert to a time struct format
    t = time.strptime(ttmp.upper(), '%b %d %Y %I:%M%p')
    
    # Convert to a datetime format
    t2=datetime.fromtimestamp(time.mktime(t))
    
    # Return as seconds since 1970,1,1
#    tout = (t2.toordinal()-datetime(1970,1,1))*86400.0
    
    # Convert to the matlab time format
    tout = datetime2matlabdn(t2)
    
    return tout
    
def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + timedelta(days = 366)
   frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac