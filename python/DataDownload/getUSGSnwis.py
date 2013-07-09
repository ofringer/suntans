# -*- coding: utf-8 -*-
"""
Download and convert USGS stream gauge data

Created on Thu Nov 15 19:14:03 2012

@author: mrayson
"""

import numpy as np
from datetime import datetime
import othertime
import urllib2
import netcdfio


def getUSGSnwis(stationids,starttime,endtime,ncfile):
    """
    Main function for grabbing the data and saving it to a netcdf file
    """
    meta={}
    meta.update({'Station ID':[]})
    meta.update({'StationName':[]})
    meta.update({'Latitude':[]})
    meta.update({'Longitude':[]})
    
    time=[]
    discharge=[]
    
    for sid in stationids:
        meta = readUSGSmeta(sid,meta=meta)
        tt,dd = readUSGStxt(sid,starttime,endtime)
        time.append(tt)
        discharge.append(dd)
        
    USGS2netcdf(ncfile,meta,time,discharge)
    

def readUSGStxt(stationid,starttime,endtime):
    """
    Read the daily station data from a web-service
    
    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """
    
    scale_fac = 0.0283168 # Cubic feet to cubic metres conversion

    target_url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=%s&startDT=%s&endDT=%s&parameterCd=00060'%(stationid,starttime,endtime)            

    try:
        print 'Opening: %s'%target_url
        f = urllib2.urlopen(target_url)
    except:
        raise Exception, 'cannot open url:\n%s'%target_url  
        
        
    StationID=[]
    time=[]
    discharge=[]
    # Read and convert the raw data
    for s in f:
        line = s.split()
        if line[0]=='USGS':
            StationID.append(line[1])
            time.append(datetime.strptime(line[2],'%Y-%m-%d'))
            discharge.append(float(line[3])*scale_fac)
            
    f.close()
        
    return np.asarray(time), np.asarray(discharge)
    

def readUSGSmeta(stationid,meta=None):
    """
    Read the station meta data from a web-service
    
    See here:
        http://waterdata.usgs.gov/nwis/news/?automated_retrieval_info#Examples
    """

    target_url = 'http://waterservices.usgs.gov/nwis/site/?format=rdb&sites=%s'%stationid
    
    try:
        f = urllib2.urlopen(target_url)
    except:
        raise Exception, 'cannot open url:\n%s'%target_url    
    
    if meta==None:
        meta={}
        meta.update({'Station ID':[]})
        meta.update({'StationName':[]})
        meta.update({'Latitude':[]})
        meta.update({'Longitude':[]})
    
    for s in f:
        line = s.split('\t')
        if line[0]=='USGS':
            meta['Station ID'].append(line[1])
            meta['StationName'].append(line[2])
            meta['Latitude'].append(float(line[4]))
            meta['Longitude'].append(float(line[5]))
            
    return meta
    
def USGS2netcdf(ncfile,meta,time,discharge):
    """
    Convert the USGS files to netcdf4 format
    """

    
    shpfile = ncfile[:-2]+'shp'
    
    varname = 'discharge'
    longname = 'Stream Discharge Rate'
    units = 'm3 s-1'
    
    ncdict=[]
    ii=-1
    for tt,dd in zip(time,discharge):
        ii+=1
        timeout = othertime.MinutesSince(tt,basetime=datetime(1970,1,1))
        ncdict=netcdfio.createObsDict(varname,longname,units,[dd],[timeout],\
                      [meta['Latitude'][ii]],[meta['Longitude'][ii]],[0.0],[meta['Station ID'][ii]],[meta['StationName'][ii]],ncdict=ncdict )
                      
    ## Global atts
    globalatts={'Title':'USGS stream gage discharge data'}
    # Write to a netcdf file
    netcdfio.writePointData2Netcdf(ncfile,ncdict,globalatts)
    # Write to a shape file
    netcdfio.pointNC2shp(ncfile,shpfile)


#####
## Example call
#
####
## Input variables
#stationids = ['08066500',\
#            '08078000',\
#            '08067500',\
#            '08076000',\
#            '08075000',\
#            '08074500',\
#            '08076500',\
#            '08073600',\
#            '08075770']
#            
#starttime = '2000-01-01'
#endtime = '2012-01-01'
#ncfile = 'C:/Projects/GOMGalveston/DATA/River/USGS_Rivers_20002012.nc'
####
#
#getUSGSnwis(stationids,starttime,endtime,ncfile)