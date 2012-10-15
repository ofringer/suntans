# -*- coding: utf-8 -*-
r"""
Suite of tools for extracting time series observations from the US-IOOS opendap server:

http://opendap.co-ops.nos.noaa.gov/dods/IOOS/

Created on Mon Jul 02 10:18:47 2012
@author: mrayson

Example
--------
>>> startt = "201001"
>>> endt = "201201"
>>> latlon = [-95.60,-93.7,28.8,30]
>>> varlist = ['waterlevel', 'conductivity', 'watertemp', 'airtemp', 'airpressure',\
>>> 'windspeed','winddirn']
>>> ncfile = 'C:/Projects/GOMGalveston/DATA/Ocean/USIOOS_OceanObs_20102011.nc'
>>> shpfile = 'C:/Projects/GOMGalveston/DATA/Ocean/USIOOS_OceanObs_20102011.shp'
>>> 
>>> # Extract the data
>>> data = extractIOOS(varlist,startt,endt,latlon)
>>> 
>>> # Write to a netcdf file
>>> globalatts = {'title':'US-IOOS oceanographic observation data',\
>>> 'history':'Created on '+datetime.ctime(datetime.now()),\
>>> 'source':'http://opendap.co-ops.nos.noaa.gov/dods/IOOS/'}
>>> netcdfio.writePointData2Netcdf(ncfile,data,globalatts)
>>> 
>>> # Write the metadata to a shapefile
>>> netcdfio.pointNC2shp(ncfile,shpfile)
"""

from pydap.client import open_url
import numpy as np
import scipy.io as io
import time
from datetime import datetime, timedelta
import urllib2
from xml.dom import minidom
import netcdfio

from netCDF4 import Dataset
import shapefile


def extractIOOS(varlist,startt,endt,latlon):
    """ Main function for extracting data from the IOOS server"""
    # Get the station ID's to retrieve
    x,y,stationID,name = queryStations(latlon)
    
    # Build up a list of monthly start and end dates
    starttime,endtime=getStartEndDates(startt,endt)
    
    # Build up the output data as a list of dictionaries
    meta=[]
    for vv in varlist:
        for lon,lat,ID,nn in zip(x,y,stationID,name):
            coords = [{'Name':'longitude','Value':lon,'units':'degrees East'},\
            {'Name':'latitude','Value':lat,'units':'degrees North'},\
            {'Name':'time','Value':[],'units':'days since 0000-00-00 00:00:00'}]
            attribs = {'StationID':str(ID),'StationName':nn,'Data':[],'coordinates':'time, longitude, latitude','coords':coords} 
            
            meta.append(attribs)
    
    # Main for loop for extracting data
    k=-1
    data=[]
    for vv in varlist:    
        if k == 0:
            break
        for ID in stationID: 
                k=k+1
                tmp = {vv:meta[k]}
                ctr=0
                for t1,t2 in zip(starttime,endtime):
                    output,t,atts = getData(str(ID),t1,t2,vv)
                    if ctr==0:                
                        tmp[vv].update(atts)
                    
                    # Append the data to the list array
                    tmp[vv]['Data'] += output
                    # Append the time data
                    ctr=-1
                    for cc in tmp[vv]['coords']:
                        ctr+=1
                        if cc['Name']=='time':
                            tmp[vv]['coords'][ctr]['Value'] += t 
                if np.size(tmp[vv]['Data']) > 0:
                    data.append(tmp)    
    
    return data
    
def extractIOOSold(varlist,startt,endt,latlon,outfile):
    """ (Deprecated) Main function for extracting data from the IOOS server"""
        
    # Get the station ID's to retrieve
    x,y,stationID,name = queryStations(latlon)
    
    # Build up a list of monthly start and end dates
    starttime,endtime=getStartEndDates(startt,endt)
    
    # Build up the output data as a list of dictionaries
    data=[]
    for lon,lat,ID,nn in zip(x,y,stationID,name):
        tmpdict = {'ID':str(ID),'name':nn,'lon':lon,'lat':lat}
        for vv in varlist:
            tmpdict[vv]=np.zeros([0,1])
            tmpdict[vv+'_time']=np.zeros([0,1])
        
        data.append(tmpdict)
    
    # Main for loop for extracting data
    k=-1
    for ID in stationID:
        k=k+1
        for vv in varlist:
            for t1,t2 in zip(starttime,endtime):
                output,t,atts = getData(str(ID),t1,t2,vv)
                # Append the data to the list array
                data[k][vv] = np.concatenate((data[k][vv],output))
                data[k][vv+'_time']= np.concatenate((data[k][vv+'_time'],t))
    
    
    print '######################## Finsished Download #######################'            
    output={'data':data}
    io.savemat(outfile,output,appendmat=False)
    print 'Data saved to: %s' % outfile
    return output
    
def getData(station_id,starttime,endtime,outvar):
    """ 
    Access data from the NOAA opendap database.pointNC2shp(ncfile,shpfile)
    Usage example:
        data,t = getData("8639348","20120501","20120512","conductivity")
    
    Input variables:
        station_id: string
        starttime: string "yyyymmdd"
        endtime: string "yyyymmdd"
        outvar: string with variable name
            one of: 'waterlevel', 'conductivity', 'watertemp', 'airtemp', 'airpressure',
                'windspeed' or 'winddirn'
    Output Variables:
        data: vector of floats with data
        t: time vector (matlab datenum format)
   
    Note that this is not a netcdf file dataserver.
    See this website for guidance:
   https://oceana.mbari.org/confluence/display/OneStopShopping/Examples+using+pydap+from+Python+to+access+BOG+data+via+DRDS
    
    """
    # Items unique to each data site
    baseurl = "http://opendap.co-ops.nos.noaa.gov/dods/IOOS/"
    if outvar == 'waterlevel':
        # Six minute verified waterlevel data
        url = baseurl + "SixMin_Verified_Water_Level"
        seqname='WATERLEVEL_6MIN_VFD_PX'
        varname = 'WL_VALUE'
        attribs = {'long_name':'Water surface elevation','units':'m'}
    elif outvar == 'conductivity':
        # Conductivity (millisiemens/cm)
        url = baseurl + "Conductivity"
        seqname='CONDUCTIVITY_PX'
        varname = 'CONDUCTIVITY' 
        attribs = {'long_name':'Water conductivity','units':'millisiemens/cm'}
    elif outvar == 'watertemp':
        # Water temperature (degC)
        url = baseurl + "Water_Temperature"
        seqname='WATER_TEMPERATURE_PX'
        varname = 'WaterTemp'  
        attribs = {'long_name':'Water temperature','units':'degreesC'}
    elif outvar == 'airtemp':
        # Air Temperature (degC)
        url = baseurl + "Air_Temperature"
        seqname='AIR_TEMPERATURE_PX'
        varname = 'AirTemp'
        attribs = {'long_name':' Air temperature','units':'degreesC'}
    elif outvar == 'airpressure':
        # Air Presure (millibars)
        url = baseurl + "Barometric_Pressure"
        seqname='BAROMETRIC_PRESSURE_PX'
        varname = 'BP'
        attribs = {'long_name':'Air pressure','units':'mb'}
    elif outvar == 'windspeed':
        # Wind Speed (m/s)
        url = baseurl + "Wind"
        seqname = 'WIND_PX'
        varname = 'Wind_Speed'
        attribs = {'long_name':'Wind speed','units':'m/s'}
    elif outvar == 'winddirn':
        # Wind Direction (degrees)
        url = baseurl + "Wind"
        seqname='WIND_PX'
        varname = 'Wind_Direction'
        attribs = {'long_name':'Wind direction','units':'m/s'}
        
        
    # Open the database
    nc = open_url(url)
    
    #my_station = nc.WATERLEVEL_6MIN_VFD_PX[(nc.WATERLEVEL_6MIN_VFD_PX._STATION_ID == station_id) & \
    #    (nc.WATERLEVEL_6MIN_VFD_PX._DATUM == "MSL") & \
    #    (nc.WATERLEVEL_6MIN_VFD_PX._BEGIN_DATE ==starttime) & \
    #    (nc.WATERLEVEL_6MIN_VFD_PX._END_DATE==endtime)]
    
    print 'Retrieving data '+outvar+' @ site # '+station_id+' for date range: '+\
        starttime+' to '+endtime+'...'
    try:
        # Build a query with the server
        if outvar == 'waterlevel':
            # Water level requires a datum in the query
            my_station = nc[seqname][(nc[seqname]._STATION_ID == station_id) & \
                (nc[seqname]._DATUM == "MSL") & \
                (nc[seqname]._BEGIN_DATE ==starttime) & \
                (nc[seqname]._END_DATE==endtime)]
        else:
            my_station = nc[seqname][(nc[seqname]._STATION_ID == station_id) & \
                (nc[seqname]._BEGIN_DATE ==starttime) & \
                (nc[seqname]._END_DATE==endtime)] 
        
        #print "Query ok - downloading data..."        
        # Get the data
        #data = np.zeros((len(my_station['DATE_TIME']),1))
        #t = np.zeros((len(my_station['DATE_TIME']),1))
        k=0
        data=[]
        t=[]
        for dt, d in zip(my_station['DATE_TIME'], my_station[varname]):
            #data[k,0]=np.float(d)
            data.append(d)
            t.append(parseDate(dt))
            k=k+1
    except:
        print "The date range and/or the variable: " + varname + " are not available from station #: "+station_id
        #data = np.zeros([0,1])
        #t=np.zeros([0,1])
        data=[]
        t=[]
     
    return data, t, attribs

def getStations():
    """Returns a list of dictionaries with NOAA station ID data"""
    
    xmlfile = 'http://opendap.co-ops.nos.noaa.gov/stations/stationsXML.jsp'
    doc = minidom.parse(urllib2.urlopen(xmlfile))
    # Load the data into a list of dictionaries
    stations = []
    for node in doc.getElementsByTagName('station'):
        name = node.getAttribute('name')
        ID = node.getAttribute('ID')
        alist = node.getElementsByTagName('lat')
        lat=alist[0].childNodes[0].nodeValue    
        alist = node.getElementsByTagName('long')
        lon=alist[0].childNodes[0].nodeValue
        alist = node.getElementsByTagName('date_established')
        date=alist[0].childNodes[0].nodeValue
        stations.append([{'name':str(name),'ID':int(ID),'lat':float(lat),'lon':float(lon),'date':str(date)}])
    
    return stations

def queryStations(latlon):
    """ Build a query of stations inside a certain domain"""
    stations = getStations()
    stationID=[]
    x=[]
    y=[]
    name=[]
    for a in stations:
        if (a[0]['lon'] >= latlon[0]) & (a[0]['lon'] <= latlon[1]) & (a[0]['lat'] \
        >= latlon[2]) & (a[0]['lat'] <= latlon[3]):
            stationID.append(a[0]["ID"])
            x.append(a[0]["lon"])
            y.append(a[0]["lat"])
            name.append(a[0]["name"])
            
    return x, y, stationID, name
    
def getStartEndDates(startt,endt):
    """
    Create a list of strings with the first day and last day of a month between two dates
    """
    # convert to datetime objects
    t1 = datetime.strptime(startt,'%Y%m')
    t2 = datetime.strptime(endt,'%Y%m')
        
    tnow = t1
    tnow2 = tnow
    starttime=[]
    endtime=[]
    while tnow < t2:
        mo = tnow.month
        yr = tnow.year
        if mo == 12:
            tnow2=tnow2.replace(month=1)
            tnow2=tnow2.replace(year=yr+1)
        else:
            tnow2=tnow2.replace(month=mo+1)
        # Get the last day of the month here
        t3=datetime.fromordinal(tnow2.toordinal()-1) 
        starttime.append(tnow.strftime('%Y%m%d'))
        endtime.append(t3.strftime('%Y%m%d'))
        # Update tnow    
        tnow = tnow2
    
    return starttime, endtime

    
def parseDate(tstr):
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

   

