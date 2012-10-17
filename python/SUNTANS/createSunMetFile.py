# -*- coding: utf-8 -*-
"""
Create a SUNTANS meteorological input netcdf file

TODO:
-----
 Use the spectral interpolation class for filling in gaps
 
Created on Fri Jul 27 17:13:04 2012
@author: mrayson
"""

import numpy as np
from datetime import datetime, timedelta
from scipy import interpolate
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import shapefile
import netcdfio
from maptools import ll2utm

import pdb

def interpWeatherStations(latlon,timestart,timeend,dt,utmzone,dbfile, showplot=False):
    """ Temporally interpolate weather station data onto a specified time grid"""

    # Create the time variable
    timeList = []
    tnow=timestart
    while tnow<timeend:
        timeList.append(tnow)
        tnow += timedelta(hours=dt)
    
    nctime = convertTime(timeList)
    ntime = len(timeList)
        
    varnames = ['Tair','Pair','Uwind','Vwind','RH','rain','cloud']
    
    coords={}
    output = {}
    # Read in the semi-processed data
    for vv in varnames:
        print 'Interpolating variable %s...'%vv
        outvar = ['NetCDF_Filename','NetCDF_GroupID','StationName']
        tablename = 'observations'
        condition = 'Variable_Name = "%s"' % vv + \
            'and time_start <= "%s"'% datetime.strftime(timestart,'%Y-%m-%d %H:%M:%S') + \
            'and time_end >= "%s"'% datetime.strftime(timeend,'%Y-%m-%d %H:%M:%S') + \
            'and lon_start >= %3.6f '%latlon[0] + 'and lon_end <= %3.6f '%latlon[1] + \
            'and lat_start >= %3.6f '%latlon[2] + 'and lat_end <= %3.6f '%latlon[3]
        
        data, query = netcdfio.queryNC(dbfile,outvar,tablename,condition)
           
        ii=0   
        for dd in data:
            ind = np.isfinite(np.ravel(dd[vv]))
            timenow = convertTime(dd['time'])
            timegood = timenow[ind]  
            if nctime[0] <= timegood[0] or nctime[-1] >= timegood[-1]:
                data.pop(ii) 
            else:
                ii+=1    

        # Remove points that have large gaps
        ii=0  
        maxgap=40
        for dd in data:
            ind = np.isfinite(dd[vv])
            timenow = convertTime(dd['time'])
            
            i=-1
            for t in timenow:
                i+=1
                if t < nctime[0]:
                    t1=i
                   
            i=-1
            for t in timenow:
                i+=1
                if t < nctime[-1]:
                    t2=i
                      
            # Find the maximum gap size between the two time limits
            gapsize = 0
            gap=0
            for gg in ind[t1:t2]:
                if ~gg:
                    gap+=1
                    if gap > gapsize:
                        gapsize=gap
                else:
                    gap = 0   
            #print t1,t2,len(timenow),gapsize
            if gapsize > maxgap:
                print 'Removing data point - gap size %d is > %d'%(gapsize,maxgap)
                data.pop(ii)
            else:
                ii+=1
                
            #print gapsize,percgood
           

        # Work out the number of spatial points of each variable based on quality control
        coords['x_'+vv] = []
        coords['y_'+vv] = []
        coords['z_'+vv] = []              
        
        for dd in data:
            # Convert to utm
            ll = np.hstack((dd['longitude'],dd['latitude']))
            xy = ll2utm(ll,utmzone)
            coords['x_'+vv].append(xy[0][0])
            coords['y_'+vv].append(xy[0][1])
            #coords['x_'+vv].append(dd['longitude'])
            #coords['y_'+vv].append(dd['latitude'])
            coords['z_'+vv].append(dd['elevation'])
        
        varlen = len(data)
        
        # Initialize the output arrays
        output[vv] = {'Data':np.zeros((ntime,varlen))}
        
        # Loop trough and interpolate each variables onto the time array
        ctr=0
        for dd in data:
            # Interpolate the data 
            tmp = np.array(np.ravel(dd[vv]))
            timenow = convertTime(dd['time'])
            ind = np.isfinite(tmp)
            F = interpolate.interp1d(timenow[ind],tmp[ind],kind='linear')
            varinterp = F(nctime)
            #F = interpolate.InterpolatedUnivariateSpline(timenow[ind],tmp[ind])
            #varinterp = F(nctime)
            #F = interpolate.splrep(timenow[ind],tmp[ind],s=0)
            #varinterp = interpolate.splev(nctime,F,der=0)
            #F = interpolate.Rbf(timenow[ind],tmp[ind])
            #varinterp = F(nctime)
            output[vv]['Data'][:,ctr]=varinterp
            ctr+=1
            
            # Add the other info
            #output[vv].update({'long_name':dd[vv]['Longname'],'units':dd[vv]['Units']})
        
            if showplot:
                plt.figure()
                plt.hold('on')
                plt.plot(timenow,tmp)
                plt.plot(nctime,varinterp,'r')
                plt.title(dd['StationName']+' - '+vv)
                plt.show()
            
    # Return the data
    return coords, output, nctime
        
def interpWeatherStationsOld(latlon,timestart,timeend,dt,localdir,showplot=False):
    """ Temporally interpolate weather station data onto a specified time grid"""
    
    varnames = ['Tair','Pair','Uwind','Vwind','RH','rain','cloud']
    
    # Read in the semi-processed data
    data = readNOAAISH.readall(latlon,[timestart.year,timeend.year],localdir) 
    
    nctime,ntime = returnTime(timestart,timeend,dt)
        
    data = dataQC(data,nctime,varnames)
   
                                    
    # Work out the number of spatial points of each variable based on quality control
    coords={}
    for vv in varnames:
        coords['lon_'+vv] = []
        coords['lat_'+vv] = []              
    
    
    for dd in data:
        for vv in dd.keys():
            if vv in varnames:
                coords['lon_'+vv].append(dd['Longitude'])
                coords['lat_'+vv].append(dd['Latitude'])
    
    varlength={}
    ctr={}                    
    for vv in varnames:
        varlength[vv]=np.size(coords['lon_'+vv])
        ctr[vv]=0
        
    # Initialize the output arrays
    output = {}
    for vv in varnames:
        output[vv] = {'Data':np.zeros([ntime,varlength[vv]])}
    
    # Loop trough and interpolate each variables onto the time array
    for dd in data:
        # Interpolate the data 
        for vv in dd.keys():
            if vv in varnames:
                tmp = np.array(dd[vv]['Data'])
                timenow = np.double(dd[vv]['Time'])
                ind = np.isfinite(tmp)
                F = interpolate.interp1d(timenow[ind],tmp[ind],kind='linear')
                varinterp = F(nctime)
                output[vv]['Data'][:,ctr[vv]]=varinterp
                ctr[vv]+=1
                
                # Add the other info
                output[vv].update({'long_name':dd[vv]['Longname'],'units':dd[vv]['Units']})
    
                if showplot:
                    plt.figure()
                    plt.hold('on')
                    plt.plot(timenow,tmp)
                    plt.plot(nctime,varinterp,'r')
                    plt.title(dd['StationName']+' - '+vv)
                    plt.show()
        
    # Return the data
    return coords, output, nctime
    
def dataQC(data,nctime,varnames):
    """ Perform basic quality control on the raw data pass
     # Checks:
    # 1) Lat/lon isn't empty
    # 2) Variable isn't empty
    # 3) nctime is bounded by time in the file
    # 4) There is an adequate number of "good" data points
    """
    
    # 1) Lat/lon isn't empty
    ii=0
    for dd in data:
        if np.size(dd['Longitude'])<1:
            ii+=1
            data.pop(ii)
            
    # 2) Remove the variable if it is emptyrootgrp = Dataset('test.nc', 'w', format='NETCDF4')
    ii=-1       
    for dd in data:
        ii+=1
        for v in dd.keys():
            if v in varnames:        
                if np.size(dd[v]['Data']) < 1:
                    data[ii].pop(v)
    
    #3) Remove if there is an inadequate number (2) of good points 
    ii=-1
    for dd in data:
        ii+=1
        for v in dd.keys():
            if v in varnames:           
                ind = np.isfinite(dd[v]['Data'])
                if ind.sum() < 2:
                    data[ii].pop(v)                
                    
    # 3) Remove if it is outside of the time domain
    ii=-1   
    for dd in data:
        ii+=1
        for v in dd.keys():
            if v in varnames:
                ind = np.isfinite(dd[v]['Data'])
                timenow = np.array(dd[v]['Time'])
                timegood = timenow[ind]
                if nctime[0] <= timegood[0] or nctime[-1] >= timegood[-1]:
                    data[ii].pop(v)
                    
    return data
def returnTime(timestart,timeend,dt):
    """ Create the time variable"""
    basetime = datetime.datetime(1970,1,1)
    dt1 = timestart-basetime
    dt2 = timeend - basetime
    nctime = np.arange(dt1.total_seconds()/60.0,dt2.total_seconds()/60.0,dt*1440.0) # minutes since 1/1/1970
    ntime = np.size(nctime) 
    return nctime, ntime
    
def write2NC(ncfile,coords,output,nctime):

    """ Writes the data to a netcdf file"""
    print 'Writing to netcdf file: %s',ncfile
    # Create an output netcdf file
    nc = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
    
    # Define the dimensions
    nc.createDimension('nt',0)
    for vv in output.keys():
        dimname = 'N'+vv
        dimlength = np.size(coords['x_'+vv])
        nc.createDimension(dimname,dimlength)
        #print '%s, %d' % (dimname, dimlength)
        
    # Create the coordinate variables
    tmpvar=nc.createVariable('Time','f8',('nt',))
    # Write the data
    tmpvar[:] = nctime
    tmpvar.long_name = 'time'
    tmpvar.units = 'seconds since 1990-01-01 00:00:00'
    for vv in output.keys():
        dimname = 'N'+vv
        varx = 'x_'+vv
        vary = 'y_'+vv
        varz = 'z_'+vv
        #print dimname, varx, coords[varx]
        tmpvarx=nc.createVariable(varx,'f8',(dimname,))
        tmpvary=nc.createVariable(vary,'f8',(dimname,))
        tmpvarz=nc.createVariable(varz,'f8',(dimname,))
        tmpvarx[:] = coords[varx]
        tmpvary[:] = coords[vary]
        tmpvarz[:] = coords[varz]
        # Create the attributes
        tmpvarx.setncattr('long_name','Longitude at '+vv)
        tmpvarx.setncattr('units','degrees_north')
        tmpvary.setncattr('long_name','Latitude at '+vv)
        tmpvary.setncattr('units','degrees_east')
        tmpvarz.setncattr('long_name','Elevation at '+vv)
        tmpvarz.setncattr('units','m')
        
    # Create the main variables
    for vv in output.keys():
        dimname = 'N'+vv
        varx = 'x_'+vv
        vary = 'y_'+vv
        tmpvar = nc.createVariable(vv,'f8',('nt',dimname))
        # Write the data
        tmpvar[:] = output[vv]['Data']
        # Write the attributes
        for aa in output[vv].keys():
            if aa != 'Data':
                tmpvar.setncattr(aa,output[vv][aa])
        # Create the all important coordinate attribute
        tmpvar.setncattr('coordinates',varx+','+vary)
        
    nc.close()
    print 'Done.'
    return

def write2CSV(latlon,timestart,timeend,dt,localdir,csvfile):
    varnames = ['Tair','Pair','Uwind','Vwind','RH','rain','cloud']
     
     
    # Write to a CSV file   
    # Read in the semi-processed data
    data = readNOAAISH.readall(latlon,[timestart.year,timeend.year],localdir) 
    
    nctime,ntime = returnTime(timestart,timeend,dt)
        
    data = dataQC(data,nctime,varnames)
    
    f = open(csvfile, 'w')
    f.write('ID, Longitude, Latitude, Variable Name, Elevation, Site Name, Site ID\n')
    
    ID = 0
    for dd in data:
        for vv in dd:
            if vv in varnames:
                ID+=1
                lon = dd['Longitude']
                lat = dd['Latitude']
                varname = dd[vv]['Longname']
                if dd[vv].has_key('Height'):
                    ele = dd[vv]['Height']
                else:
                    ele = 0.0
                StationName = dd['StationName']
                StationID = dd['StationID']
                fstr = '%d, %10.6f, %10.6f, %s, %3.1f, %s, %s\n' % (ID,lon,lat,varname,ele,StationName,StationID)
                f.write(fstr)
    
    f.close()   
    return

def write2SHP(latlon,timestart,timeend,dt,localdir,shpfile):
    varnames = ['Tair','Pair','Uwind','Vwind','RH','rain','cloud']
         
    # Write to a SHP file   
    
    # Read in the semi-processed data
    data = readNOAAISH.readall(latlon,[timestart.year,timeend.year],localdir) 
    
    nctime,ntime = returnTime(timestart,timeend,dt)
        
    data = dataQC(data,nctime,varnames)
    
    ID = 0
    w = shapefile.Writer(shapefile.POINT)
    w.field('Long Name')
    w.field('Station Name')
    w.field('Station ID')
    for dd in data:
        for vv in dd:
            if vv in varnames:
                ID+=1
                lon = dd['Longitude']
                lat = dd['Latitude']
                varname = dd[vv]['Longname']
                if dd[vv].has_key('Height'):
                    ele = dd[vv]['Height']
                else:
                    ele = 0.0
                StationName = dd['StationName']
                StationID = dd['StationID']
                
                w.point(lon,lat)
                w.record(varname,StationName,StationID)
    #   
    w.save(shpfile)
    return
   
def convertTime(timein):
    """Converts a list of time object into an array of seconds since 1990-01-01"""
    basetime = datetime(1990,1,1)
    
    timeout=[]
    for t in timein:
        dt = t - basetime
        timeout.append(dt.total_seconds())
        
    return np.array(timeout)
    
    
#######
## Inputs
#latlon = [-95.60,-94.3,28.8,30]
##latlon = [-95.0,-94.3,28.8,29.5]
#timestart = datetime(2011,2,1,0,0,0)
#timeend = datetime(2011,3,1,0,0,0)
#dt = 1 # time step in hours
#utmzone = 15
#dbfile = 'C:/Projects/GOMGalveston/DATA/GalvestonObs.db'
#showplot=False
#ncfile='C:\Projects\GOMGalveston\MODELLING\WINDS\Galveston_Winds_Feb2011_UTM.nc'
#
#    
##
#[coords, output, nctime] = interpWeatherStations(latlon,timestart,timeend,dt,utmzone,showplot=False)
##
#
#write2NC(ncfile,coords,output,nctime)