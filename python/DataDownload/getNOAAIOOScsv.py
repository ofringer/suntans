"""
 Tools for extracting data through the NOAA GetObservation request function 
 
 Some data sets are only available through this interface e.g. the PORTS ADCP data
 
 http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/
"""

import urllib2  # the lib that handles the url stuff
from datetime import datetime, timedelta
import othertime
from xml.dom import minidom
import numpy as np
import airsea
import netcdfio

import pdb

def main(vartype,bbox,timestart,timeend,ncfile,dbfile=None):
    varlookup ={'CurrentsActive':'currents','WaterLevelActive':'water_surface_height_above_reference_datum',\
    'Salinity':'sea_water_salinity','Conductivity':'sea_water_electrical_conductivity'}
    
    
    # Find the stations
    staInfo = stationInfo(vartype)
    for vv in staInfo.keys():
        if staInfo[vv]['lon']>=bbox[0] and staInfo[vv]['lon']<=bbox[1] and staInfo[vv]['lat']>=bbox[2] and staInfo[vv]['lat']<=bbox[3]:
            print 'Station %s inside lat-lon range.' % staInfo[vv]['name']
        else:
            staInfo.pop(vv)
    
    ncdata=[]
    for ID in staInfo.keys():
        print 'Getting %s data from %s' % (vartype,staInfo[ID]['name'])
        
        # Grab the station data
        data = getAllTime(ID,varlookup[vartype],timestart,timeend)
        
        ###
        if len(data)>0:
            # Convert the output to the format required for the netcdf file
            if vartype == 'CurrentsActive':
                ncdatatmp=parseADCPdata(data,staInfo,ID)
            elif vartype == 'WaterLevelActive':
                ncdatatmp=parseWaterLev(data,staInfo,ID)
        
            ncdata+=ncdatatmp
        
    # Write to the output netcdf
    globalatts = {'title':'US-IOOS observation data',\
        'history':'Created on '+datetime.ctime(datetime.now()),\
        'source':'http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/index.jsp'}    
    netcdfio.writePointData2Netcdf(ncfile,ncdata,globalatts)
        
    # Write the metadata to a shapefile
    #shpfile = ncfile.split('.')[0]+'.shp'
    #netcdfio.pointNC2shp(ncfile,shpfile)
    
    # Update the database
    #createObsDB(dbfile)
    if not dbfile == None:
        print 'Updating database: %s'%dbfile
        netcdfio.netcdfObs2DB(ncfile,dbfile)                

def getAllTime(stationID,vartype,timestart,timeend): 
    
    """Wrapper function to grab data for extened time periods in chunks"""
    
    # Get the start and end times - can only extract 4-day chunks
    timeList=[]
    #
    #tnow = timestart
    #while tnow <= timeend:
    #    timeList.append(tnow)
    #    tnow += timedelta(days=4)
    timeList = othertime.TimeVector(timestart,timeend,4*86400.,istimestr=False)
    
    k=0
    for t1s,t2s in zip (timeList[0:-1],timeList[1:]):
        
        t1 = datetime.strftime(t1s,'%Y-%m-%dT%H:%M:%SZ')
        t2 = datetime.strftime(t2s,'%Y-%m-%dT%H:%M:%SZ')
        print 'Extracting dates %s to %s' % (t1,t2)
        # generate the url string
        target_url='http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/SOS?service=SOS&request=GetObservation&version=1.0.0&observedProperty='+vartype+'&offering=urn:ioos:station:NOAA.NOS.CO-OPS:'+stationID+'&responseFormat=text%2Fcsv&eventTime='+t1+'/'+t2
        
        # Get the data
        datatmp=getObsfromURL(target_url)
        
        # append the data
        if k == 0 and len(datatmp)>1:
            data=datatmp
            k+=1
        elif k > 0 and len(datatmp)>1:
            for vv in data.keys():
                data[vv]+=datatmp[vv]
    try:
        return data
    except:
        return []
    
def getObsfromURL(target_url):
    """ Reads the actual data from the CSV and returns a dictionary"""
    data={}
    #try:
    #print target_url
    csvfile = urllib2.urlopen(target_url) # it's a file like object and works just like a file
    k=-1
    
    for line in csvfile: # files are iterable
        k+=1
        if k == 0:
            # Get the variable names from the header
            varnames = line.split(',')
            for vv in varnames:
                data.update({vv:[]})
        else:
            # Append the data
            obs = line.split(',')
            for oo,vv in zip(obs,varnames):
                data[vv].append(parseData(oo,vv))
    #except:
    #    print 'No data exists for this period'                
       
    if len(data)==1:
        print 'No data exists for this period'
        data={}
    return data
    
    
def getType(vv):
    
    """ Gets the variable type from its name
        Types are:
            0 - str
            1 - int
            2 - float
            3 - time
    """
    
    vartype = {'station_id':0,'sensor_id':0,'"latitude (degree)"':2,'"longitude (degree)"':2,\
    'date_time':3,'"sensor_depth (m)"':2,'"direction_of_sea_water_velocity (degree)"':2,\
    '"sea_water_speed (cm/s)"':2,'"platform_orientation (degree)"':2,'"platform_pitch_angle (degree)"':2,\
    '"platform_roll_angle (degree)"':2,'"sea_water_temperature (C)"':2,'orientation':0,'"sampling_rate (Hz)"':2,\
    '"reporting_interval (s)"':1,'processing_level':0,'"bin_size (m)"':2,'"first_bin_center (m)"':2,\
    'number_of_bins':1,'"bin (count)"':1,'"bin_distance (m)"\n':2,\
    '"water_surface_height_above_reference_datum (m)"':2,'"vertical_position (m)"':2,'datum_id':0}
    
    if vartype.has_key(vv):
        return vartype[vv]
    #else:
        #print 'Unknown variable type: %s'%vv

def parseData(oo,vv):
    
    datatype = getType(vv)
    if datatype == 0:
        return oo
    elif datatype == 1:
        return int(oo)
    elif datatype == 2:
        return float(oo)
    elif datatype == 3:
        return datetime.strptime(oo,'%Y-%m-%dT%H:%M:%SZ')
    else:
        return oo

def maxList(listvec):
    mx = listvec[0]
    for l in listvec:
        if l > mx:
            mx = l
    return mx

def stationInfo(vartype='CurrentsActive'):
    
    """ Returns the information about the available stations from an xml file
    vartype must be one of: 
        'All','CurrentsActive','CurrentsSurvey','MetActive','WaterLevelActive'
    
    """
    
    xmlfile="http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/SOS?service=SOS&request=DescribeSensor&version=1.0.0&outputFormat=text/xml;subtype=\"sensorML/1.0.1\"&procedure=urn:ioos:network:NOAA.NOS.CO-OPS:"+vartype
    
    doc = minidom.parse(urllib2.urlopen(xmlfile))
    
    stations = {}
    for node in doc.getElementsByTagName('gml:Point'):
        # get the station in
        ID = node.getAttribute('gml:id')
        idstr = ID.split('-')
        # Get the coordinates
        alist = node.getElementsByTagName('gml:coordinates')
        coordstr=alist[0].childNodes[0].nodeValue
        lonlatstr = coordstr.split(' ')
        lat = float(lonlatstr[0])
        lon = float(lonlatstr[1])
        stations.update({idstr[1]:{'lon':lon,'lat':lat,'name':[]}})
    
    # Get the station names
    for node in doc.getElementsByTagName('sml:System'):
        ID = node.getAttribute('gml:id')
        idstr = ID.split('-')
         # Get the coordinates
        alist = node.getElementsByTagName('gml:description')
        namestr=alist[0].childNodes[0].nodeValue
        if stations.has_key(idstr[1]):
            stations[idstr[1]]['name']=namestr
        #print idstr[1], namestr
        #print node.attributes.keys()
        
    return stations

def parseADCPdata(data,staInfo,ID):
    basetime = datetime(1970,1,1)
    # Find the number of bins
    nz = maxList(data['"bin (count)"'])
    nt = len(data['"bin (count)"'])/nz
        
    # Get the depth data
    ele = -np.array(data['"bin_distance (m)"\n'][0:nz])
    # Load in the data
    time = np.zeros(nt)
    u = np.zeros((nt,nz))
    v = np.zeros((nt,nz))
    T = np.zeros(nt)
    ii=0
    for tt in range(0,nt):
        dt = data['date_time'][ii]-basetime
        time[tt] = dt.total_seconds()/60.0
        T[tt]=data['"sea_water_temperature (C)"'][ii]
        for zz in range(0,nz):
            uu,vv=airsea.convertSpeedDirn(data['"direction_of_sea_water_velocity (degree)"'][ii],data['"sea_water_speed (cm/s)"'][ii]/100)
            u[tt,zz]=uu
            v[tt,zz]=vv
            ii+=1
    # Insert each variable into a list
    ncdict = []
    # U
    coords = [{'Name':'longitude','Value':staInfo[ID]['lon'],'units':'degrees East'},\
                {'Name':'latitude','Value':staInfo[ID]['lat'],'units':'degrees North'},\
                {'Name':'elevation','Value':ele,'units':'metres','positive':'up'},\
                {'Name':'time','Value':time,'units':'minutes since 1970-01-01 00:00:00'}]
                
    attribs = {'StationID':ID,'StationName':staInfo[ID]['name'],'Data':u,\
    'coordinates':'time, elevation, latitude, longitude','long_name':'Eastward Water Velocity Component',\
    'units':'m/s','coords':coords} 
    ncdict.append({'water_u':attribs})
    # V
    coords = [{'Name':'longitude','Value':staInfo[ID]['lon'],'units':'degrees East'},\
                {'Name':'latitude','Value':staInfo[ID]['lat'],'units':'degrees North'},\
                {'Name':'elevation','Value':ele,'units':'metres','positive':'up'},\
                {'Name':'time','Value':time,'units':'minutes since 1970-01-01 00:00:00'}]
                
    attribs = {'StationID':ID,'StationName':staInfo[ID]['name'],'Data':v,\
    'coordinates':'time, elevation, latitude, longitude','long_name':'Northward Water Velocity Component',\
    'units':'m/s','coords':coords} 
    ncdict.append({'water_v':attribs})
    # T
    coords = [{'Name':'longitude','Value':staInfo[ID]['lon'],'units':'degrees East'},\
                {'Name':'latitude','Value':staInfo[ID]['lat'],'units':'degrees North'},\
                {'Name':'elevation','Value':ele[0],'units':'metres','positive':'up'},\
                {'Name':'time','Value':time,'units':'minutes since 1970-01-01 00:00:00'}]
                
    attribs = {'StationID':ID,'StationName':staInfo[ID]['name'],'Data':T,\
    'coordinates':'time, elevation, latitude, longitude','long_name':'Sea Water Temperature',\
    'units':'degrees C','coords':coords} 
    ncdict.append({'watertemp':attribs})
    
    return ncdict
    
def parseWaterLev(data,staInfo,ID):
    basetime = datetime(1970,1,1)
    # Find the number of bins
    nt = len(data['"water_surface_height_above_reference_datum (m)"'])
    
    # Get the depth data
    ele = 0.0
    # Load in the data
    time = np.zeros(nt)
    ssh = np.zeros(nt)
    ii=0
    for tt in range(0,nt):
        dt = data['date_time'][ii]-basetime
        time[tt] = dt.total_seconds()/60.0
        ssh[tt]=data['"water_surface_height_above_reference_datum (m)"'][ii]
        ii+=1
    # Insert each variable into a list
    ncdict = []
    # U
    coords = [{'Name':'longitude','Value':staInfo[ID]['lon'],'units':'degrees East'},\
                {'Name':'latitude','Value':staInfo[ID]['lat'],'units':'degrees North'},\
                {'Name':'elevation','Value':ele,'units':'metres','positive':'up'},\
                {'Name':'time','Value':time,'units':'minutes since 1970-01-01 00:00:00'}]
                
    attribs = {'StationID':ID,'StationName':staInfo[ID]['name'],'Data':ssh,\
    'coordinates':'time, elevation, latitude, longitude','long_name':'Water surface elevation',\
    'units':'m','coords':coords} 
    ncdict.append({'waterlevel':attribs})
    return ncdict

###

#if __name__=="__main__":
#    ###
#    # Inputs
#    bbox = [-95.40,-94.49,28.8,29.9]
#    vartype = 'WaterLevelActive'
#    vartype = 'Conductivity'
#    timestart = datetime(2000,1,1)
#    timeend = datetime(2000,1,6)
#    ncfile = 'C:\Projects\GOMGalveston\DATA\Ocean\IOOS_Salinity_20002012.nc'
#    dbfile = 'C:/Projects/GOMGalveston/DATA/GalvestonObs.db'
#    main(vartype,bbox,timestart,timeend,ncfile,dbfile)
