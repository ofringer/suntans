# -*- coding: utf-8 -*-
r"""
	Tools for downloading and converting NOAA/NWS weather station data.
	
	Data is downloaded from the NOAA ftp server:
		ftp://ftp.ncdc.noaa.gov/pub/data/ish/
	
	Example
    --------
	>>> from datetime import datetime
	>>> import getNOAAWeatherStation as noaa
	>>> 
	>>> latlon = [-95.60,-93.7,28.8,30]
	>>> timestart = datetime(2010,1,1,0,0,0)
	>>> timeend = datetime(2011,12,31,0,0,0)
	>>> dt = 1.0/24.0
	>>> localdir = 'C:/Projects/GOMGalveston/CODE/PYTHON/NOAAWeather/rawdata/'
	>>> showplot = False
	>>> ncfile = 'C:/Projects/GOMGalveston/DATA/Winds/NCDCNWS_AirObs_20102011.nc'
	>>> shpfile = 'C:/Projects/GOMGalveston/DATA/Winds/NCDCNWS_AirObs_20102011.shp'
	>>> 
	>>> data = noaa.noaaish2nc(latlon,[timestart.year,timeend.year],localdir,ncfile,shpfile)
"""
import csv 
import numpy as np
from ftplib import FTP
import gzip
from datetime import datetime
import os
import airsea
import netcdfio


def noaaish2nc(latlon,yearrange,localdir,ncfile,shpfile):
    
    """ Reads in the noaa data and spits it out to a netcdf file"""
    
    varnames = ['Tair','Pair','Uwind','Vwind','RH','rain','cloud']
    # Read in the semi-processed data
    data = readall(latlon,[timestart.year,timeend.year],localdir) 
    
    data = dataQC(data,varnames)
       
    # Create the dictionary format necessary to parse the data to a grouped netcdf file
    ncdict=[]
    for dd in data: 
        for vv in dd.keys():
            if vv in varnames:
                # Start building the dictionary
                if np.size(dd['Longitude'])>0 or np.size(dd['Latitude'])>0:
                    if dd[vv].has_key('Height'):
                        ele = dd[vv]['Height']
                    else:
                        ele=0.0
                    coords = [{'Name':'longitude','Value':dd['Longitude'],'units':'degrees East'},\
                    {'Name':'latitude','Value':dd['Latitude'],'units':'degrees North'},\
                    {'Name':'elevation','Value':ele,'units':'metres','positive':'up'},\
                    {'Name':'time','Value':dd[vv]['Time'],'units':'minutes since 1970-01-01 00:00:00'}]
                    
                    attribs = {'StationID':dd['StationID'],'StationName':dd['StationName'],'Data':dd[vv]['Data'],\
                    'coordinates':'time, elevation, longitude, latitude','long_name':dd[vv]['Longname'],\
                    'units':dd[vv]['Units'],'coords':coords} 
                    ncdict.append({vv:attribs})
                else:
                    print 'Station: %s missing coordinate info - skipping'%dd['StationName']
    
    globalatts = {'title':'NCDC/NWS integrated surface hourly observation data',\
    'history':'Created on '+datetime.ctime(datetime.now()),\
    'source':'ftp://ftp.ncdc.noaa.gov/pub/data/ish/'}            
    
    # Write the output to a netcdf file
    netcdfio.writePointData2Netcdf(ncfile,ncdict,globalatts)
    
    # Write the metadata to a shapefile
    netcdfio.pointNC2shp(ncfile,shpfile)
    
    return ncdict
    
    
def readall(latlon,yearrange,localdir):
    """ Main function for reading in the ISH file data"""
    
    # Connect to the ftp server
    ftpdir = '/pub/data/noaa/'
    ftpsrvr = 'ftp.ncdc.noaa.gov'
    ftp = FTP(ftpsrvr)
    ftp.login() 
    #ftp.cwd(ftpdir)
    #ftp.retrlines('LIST')
    
    # List the local files
    localfiles = os.listdir(localdir)
    
    # Get the names of the stations 
    stations = getFileNames(latlon,yearrange)
    
    # Loop through and check if the files are in the local directory, if not download
    data_all = []
    for dd in stations:
        kk=0
        for yy, ff in zip(dd['years'],dd['filenames']):
            kk+=1
            # Check if the file exists locally
            if ff in localfiles:
                print 'File found locally: %s' % ff
                gzfile = localdir+ff
            else:
                gzfile = localdir+ff
                ftpfile = ff
                ftp.cwd(ftpdir+str(yy))
                #ftp.retrlines('LIST')
                print 'Downloading file: %s...' % ftpfile
                ftp.retrbinary('RETR '+ftpfile,open(gzfile, 'wb').write)
                print 'Completed download.'
            
           
            if kk == 1:
                 # Now read the data
                 data = ishData2struct(gzfile,dd['station_id'],dd['station_name'])
            else:
                # Now read the data
                datanew = ishData2struct(gzfile,dd['station_id'],dd['station_name'])
                # Append the time series for each station
                data['Tair']['Data']=data['Tair']['Data']+datanew['Tair']['Data']
                data['Tair']['Time']=data['Tair']['Time']+datanew['Tair']['Time']
                data['Pair']['Data']=data['Pair']['Data']+datanew['Pair']['Data']
                data['Pair']['Time']=data['Pair']['Time']+datanew['Pair']['Time']
                data['Uwind']['Data']=data['Uwind']['Data']+datanew['Uwind']['Data']
                data['Uwind']['Time']=data['Uwind']['Time']+datanew['Uwind']['Time']
                data['Vwind']['Data']=data['Vwind']['Data']+datanew['Vwind']['Data']
                data['Vwind']['Time']=data['Vwind']['Time']+datanew['Vwind']['Time']
                data['RH']['Data']=data['RH']['Data']+datanew['RH']['Data']
                data['RH']['Time']+=datanew['RH']['Time']
                data['rain']['Data']+=datanew['rain']['Data']
                data['rain']['Time']+=datanew['rain']['Time']
                data['cloud']['Data']+=datanew['cloud']['Data']
                data['cloud']['Time']+=datanew['cloud']['Time']
            
        data_all.append(data)
    return data_all
    
             
def stationMeta(csvfile = "ish-history.csv"):
    """
    Read the NOAA station metadata csv file
    
    Created on Mon Jul 23 15:57:03 2012
    @author: mrayson
    """
    histf = csv.reader(open(csvfile,'r'))
    
    i=0
    missingval = -999.0
    station_id = []
    station = []
    lon = []
    lat = []
    elevation = []
    startyr = []
    endyr = []
    for row in histf:
        if i > 0:
            station_id.append(row[0]+'-'+row[1])
            station.append(row[2])
            
            if len(row[7])>0:
                lat.append(np.float(row[7])/1000)
            else:
                lat.append(missingval)
                
            if len(row[8])>0:
                lon.append(np.float(row[8])/1000)
            else:
                lon.append(missingval)
                
            if len(row[9])>0:
                elevation.append(np.float(row[9])/10)
            else:
                elevation.append(missingval)
                
            startyr.append(row[10][0:4])
            endyr.append(row[11][0:4])
        i += 1
        
    stationdata = {'station_id':station_id,'name':station, \
        'lon':lon,'lat':lat,'elevation':elevation,'start':startyr,'end':endyr}
    
    return stationdata
    
def parseAddData(line,recl):
    """
    Parses the additional data string at the end of the line
    """
    pos = 0
    data={}
    
    if 'RHUM'  in line:
        print '!!!!!!!!!!!!! Found Humidity !!!!!!!!!!!!!!!!!!!!!!!'
    for i in range(0,50):        
        fc = line[pos:pos+3] # field code
        pos=pos+3
        if fc in ('AA1', 'AA2', 'AA3', 'AA4'):
            #LIQUID-PRECIPITATION occurrence identifier
            data['LIQUID_PRECIPITATION_DURATION']=int(line[pos:pos+2]) # hrs
            data['LIQUID_PRECIPITATION_DURATION']=checkMissingVal(data['LIQUID_PRECIPITATION_DURATION'],99)
            pos+=2
            data['LIQUID_PRECIPITATION_VOLUME']=float(line[pos:pos+4])/10 # mm
            data['LIQUID_PRECIPITATION_VOLUME']==checkMissingVal(data['LIQUID_PRECIPITATION_VOLUME'],999.9)
            pos+=6          
        elif fc == 'AB1':
            #LIQUID-PRECIPITATION MONTHLY TOTAL (skip)
            pos+=7
        elif fc == 'AC1':
            #PRECIPITATION-OBSERVATION-HISTORY
            pos +=3
        elif fc == 'AD1':
            #LIQUID-PRECIPITATION GREATEST AMOUNT IN 24 HOURS, FOR THE MONTH
            pos +=19
        elif fc == 'AE1':
            #LIQUID-PRECIPITATION, NUMBER OF DAYS WITH SPECIFIC AMOUNTS, FOR THE MONTH                
            pos+=12
        elif fc == 'AG1':
            #PRECIPITATION-ESTIMATED-OBSERVATION
            pos+=4
        elif (fc == 'AH1') or (fc == 'AH2') or (fc == 'AH3') or (fc == 'AH4') or (fc=='AH5') or (fc=='AH6'):
            #LIQUID-PRECIPITATION MAXIMUM SHORT DURATION, FOR THE MONTH
            pos+=15
        elif (fc == 'AI1') or (fc == 'AI2') or (fc == 'AI3') or (fc == 'AI4')or (fc=='AI5') or (fc=='AI6'):
            #LIQUID-PRECIPITATION MAXIMUM SHORT DURATION, FOR THE MONTH
            pos+=15
        elif fc == 'AJ1':
            #SNOW-DEPTH
            pos+=14
        elif fc == 'AK1':
            #SNOW-DEPTH GREATEST DEPTH ON THE GROUND,                
            pos+=12
        elif (fc == 'AL1') or (fc == 'AL2') or (fc == 'AL3') or (fc == 'AL4'):
            # SNOW-ACCUMULATION
            pos+=7
        elif fc=='AM1':
            #SNOW-ACCUMULATION GREATEST AMOUNT IN 24 HOURS, FOR THE MONTH
            pos+=18
        elif fc=='AN1':
            pos+=9
        elif (fc == 'AO1') or (fc == 'AO2') or (fc == 'AO3') or (fc == 'AO4'):
            pos+=8
        elif (fc == 'AP1') or (fc == 'AP2') or (fc == 'AP3') or (fc == 'AP4'):
            pos+=6
        elif fc in ('AU1','AU2','AU3','AU4','AU5','AU6','AU7','AU8', 'AU9'):
            pos+=8
        elif (fc == 'AW1') or (fc == 'AW2') or (fc == 'AW3') or (fc == 'AW4'):
            pos+=3
        elif (fc == 'AX1') or (fc == 'AX2') or (fc == 'AX3') or (fc == 'AX4')or (fc=='AX5') or (fc=='AX6'):
            pos+=6
        elif (fc == 'AY1') or (fc == 'AY2'):
            pos+=5
        elif (fc == 'AZ1') or (fc == 'AZ2'):
            pos+=5       
        elif (fc == 'CB1') or (fc == 'CB2'):
            data['PRECIP_PERIOD']=int(line[pos:pos+2]) # minutes
            data['PRECIP_VOLUME']=float(line[pos:pos+2])/10 # mm
            pos+=2
        elif (fc == 'CF1') or (fc == 'CF2') or (fc == 'CF3'):   
             pos+=6
        elif (fc == 'CG1') or (fc == 'CG2') or (fc == 'CG3'):
             pos+=8
        elif (fc == 'CH1') or (fc == 'CH2'):
            # Relative humidity data
            pos+=9
            data['RH']=float(line[pos:pos+4])/10 # Percent
            print '!!!!!!!!!!!!! Found Humidity !!!!!!!!!!!!!!!!!!!!!!!'
            pos+=2
        elif (fc == 'CI1'):
            print '!!!!!!!!!!!!! Found Humidity !!!!!!!!!!!!!!!!!!!!!!!'
            pos+=28
        elif (fc == 'CN1'):
            pos+=18
        elif (fc == 'CN2'):
            pos+=18
        elif (fc == 'CN3'):
            pos+=16
        elif (fc == 'CN4'):
            pos+=16
        elif (fc == 'CR1'):
            pos+=7
        elif (fc == 'CT1')or(fc == 'CT2') or (fc == 'CT3'):
            pos+=7
        elif (fc == 'CU1')or(fc == 'CU2') or (fc == 'CU3'):
            pos+=13
        elif (fc == 'CV1')or(fc == 'CV2') or (fc == 'CV3'):
            pos+=26
        elif (fc == 'CW1'):
            pos+=14
        elif (fc == 'CX1')or(fc == 'CX2') or (fc == 'CX3'):
            pos+=26
        elif (fc == 'C01'):
            pos+=5
        elif fc in ('CO2','CO3','CO4','CO5','CO6','CO7','CO8','CO9'):
            pos+=8
        elif fc == 'ED1':
            pos+=8
        elif fc in ('GA1','GA2','GA3','GA4','GA5','GA6'):
            #00: None, SKC or CLR
            #01: One okta - 1/10 or less but not zero
            #02: Two oktas - 2/10 - 3/10, or FEW
            #03: Three oktas - 4/10
            #04: Four oktas - 5/10, or SCT
            #05: Five oktas - 6/10
            #06: Six oktas - 7/10 - 8/10
            #07: Seven oktas - 9/10 or more but not 10/10, or BKN
            #08: Eight oktas - 10/10, or OVC
            #09: Sky obscured, or cloud amount cannot be estimated
            #10: Partial obscuration
            #99: Missing
            data['SKY_COVER_LAYER'] = line[pos:pos+2]
            pos+=13
        elif fc in ('GD1','GD2','GD3','GD4','GD5','GD6'):
            #00: None, SKC or CLR
            #01: One okta - 1/10 or less but not zero
            #02: Two oktas - 2/10 - 3/10, or FEW
            #03: Three oktas - 4/10
            #04: Four oktas - 5/10, or SCT
            #05: Five oktas - 6/10
            #06: Six oktas - 7/10 - 8/10
            #07: Seven oktas - 9/10 or more but not 10/10, or BKN
            #08: Eight oktas - 10/10, or OVC
            #09: Sky obscured, or cloud amount cannot be estimated
            #10: Partial Obscuration
            #11: Thin Scattered
            #12: Scattered
            #13: Dark Scattered
            #14: Thin Broken
            #15: Broken
            #16: Dark Broken
            #17: Thin Overcast
            #18: Overcast
            #19: Dark overcast
            #99: Missing
            pos+=1
            data['SKY_COVER_SUMMARY'] = line[pos:pos+2]
            pos+=11
        elif fc == 'GF1':
            # Same codes as above
            data['SKY_COND_OBS'] = line[pos:pos+2]
            pos+=23
        elif fc in ('GG1','GG2','GG3','GG4','GG5','GG6'):
            pos+=15
        elif fc == 'GH1':
            # Hourly Solar Radiation Data
            data['SOLRAD'] = float(line[pos:pos+5])
            pos+=28
        elif fc == 'GJ1':
            pos+=5
        elif fc == 'GK1':
            pos+=4
        elif fc == 'GL1':
            pos+=6
        elif fc == 'GM1':
            # Solar irradiance
            pos+=4
            data['SOLIRRAD'] = float(line[pos:pos+4])
            pos+=26
        elif fc == 'GN1':
            # upwelling and downwell solar and longwave radiation
            pos+=28
        elif fc == 'GO1':
            # Net radiation
            pos+=4
            data['NET_SHORTWAVE'] = float(line[pos:pos+4])
            pos+=5
            data['NET_LONGWAVE'] = float(line[pos:pos+4])
            pos+=10
        elif fc == 'GP1':
            pos+=31
        elif fc == 'GQ1':
            pos+=13
        elif fc == 'GR1':
            pos+=14
            
        ########################################################################
        # Inserted these miscellaneuosly as they appear
        ########################################################################
        elif fc in ('KA1','KA2'):
            pos+=10
        elif fc in ('KB1','KB2','KB3'):
            pos+=10
        elif fc in ('KC1','KC2'):
            pos+=14
        elif fc in ('KD1','KD2'):
            pos+=9
        elif fc == 'KE1':
            pos+=12
        elif fc == 'MA1':
            pos+=12   
        elif fc == 'MD1':
            pos+=11
        elif fc == 'MG1':
            pos+=12
        elif fc == 'MH1':
            pos+=12
        elif fc == 'MK1':
            pos+=24
        elif fc == 'MV1':
            pos+=3
        elif fc in ('MW1','MW2','MW3','MW4','MW4','MW6','MW7'):
            pos+=3
        elif fc == 'OC1':
            pos+=5
        elif fc in ('OE1','OE2','OE3'):
            pos+=16
        elif fc in ('RH1','RH2','RH3'):
            # Relatuve humidity
            pos+=4
            data['RH1']=int(line[pos:pos+3])
            pos+=5
            print '!!!!!!!!!!!!! Found Humidity !!!!!!!!!!!!!!!!!!!!!!!'
        elif fc == 'SA1':
            # Sea surface temperature
            pos+=5
        elif fc == 'UA1':
            # Wave data
            pos+=10
        elif fc == 'UG1':
            # Wave data
            pos+=9
        elif fc in ('REM','SYN','MET'):
            # Remarks section
           return data
           break    
        elif len(fc)==0:
            return data 
            break
        else:
            
            print 'Unknown identifier: %s' % fc
            print line
            return data
            break
                
    return data
   
def checkMissingVal(data,missingval,tol=1e-3):
    """ Check the quality of the data"""
    if data >=(missingval-tol) and data <=(missingval+tol):
        return np.nan
    else:
        return data

# end of function
def readRawGZdata(gzfile):
    """
    Read the GZ raw format here
    """    
    f = gzip.open(gzfile, 'rb')
    ishdata=[] # Store all of the raw data as a list of dictionaries
    for line in f.readlines():
        
        data={}
        # Parse the line here
        tvc = int(line[0:4]) # total variable character
        recl = tvc-108 # Record length of the additional data
        data['yyyymmdd'] = line[15:23]
        data['hhmm'] = line[23:27]
        data['lat'] = float(line[28:34])/1000
        data['lon'] = float(line[34:41])/1000
        data['ele'] = float(line[46:51])
        data['winddir'] = float(line[60:63]) # degrees
        data['windspd'] = float(line[65:69])/10 # m/s
        data['airtemp'] = float(line[87:92])/10 # degC
        data['airdewpoint'] = float(line[93:98])/10 # degC
        data['airpres'] = float(line[99:104])/10 # hPa
        
        # Replace missing values with nans
        data['winddir'] = checkMissingVal(data['winddir'],999.0)
        data['windspd'] = checkMissingVal(data['windspd'],999.9)
        data['airtemp'] = checkMissingVal(data['airtemp'],999.9)
        data['airdewpoint'] = checkMissingVal(data['airdewpoint'],999.9)
        data['airpres'] = checkMissingVal(data['airpres'],9999.9)

        
        lineadd =  line[108:] # Additional variable data
        # Parse the additional data string
        addldata = parseAddData(lineadd,recl)
        
        # Append the additional data to the dictionary
        data.update(addldata)
        ishdata.append(data)
    
    f.close()
    return ishdata

 
def returnHumidity(dd):
    """ Returns humidity data if it exists in the dictionary"""
    rh = []
    if dd.has_key('RH'):
        rh = dd['RH']
    elif dd.has_key('RH1'):
        rh = dd['RH1']
    else:
        # Convert the dew point temperature to relative humidity
        Pmb = dd['airpres']/10 # hPa to mb
        rh = airsea.relHumFromTdew(dd['airtemp'],dd['airdewpoint'],Pmb)
    return rh        
        
def returnRainfall(dd):
    """Returns rainfall data in units kg/m2/s"""
    rho_fresh = 1000 # freshwater density
    rain = []
    if  dd.has_key('LIQUID_PRECIPITATION_VOLUME'):
        period = dd['LIQUID_PRECIPITATION_DURATION'] # hours
        volume = dd['LIQUID_PRECIPITATION_VOLUME'] # mm
        rain = volume/(1000*period*3600)*rho_fresh
    elif dd.has_key('PRECIP_VOLUME'):
        period = dd['PRECIP_PERIOD'] # minutes
        volume = dd['PRECIP_VOLUME'] # mm
        rain = volume/(1000*period*60)*rho_fresh
    return rain
    
def returnCloudCover(dd):
    """Returns cloud cover as a fraction"""
    
    cloud_codeA = {'00':0.0,'01':0.05,'02':0.25,'03':0.4,'04':0.5,'05':0.6,'06':0.75 \
    ,'07':0.95,'08':1.0,'09':np.nan,'10':np.nan,'99':np.nan}
    cloud = []
    if dd.has_key('SKY_COVER_LAYER'):
        cloud_code = dd['SKY_COVER_LAYER']
        if cloud_code in cloud_codeA:
            cloud = cloud_codeA[cloud_code] 
            
    elif dd.has_key('SKY_COVER_SUMMARY'):
        cloud_code = dd['SKY_COVER_SUMMARY']
        if cloud_code in cloud_codeA:
            cloud = cloud_codeA[cloud_code] 
            
    elif dd.has_key('SKY_COND_OBS'):
        cloud_code = dd['SKY_COND_OBS']
        if cloud_code in cloud_codeA:
            cloud = cloud_codeA[cloud_code] 
            
    return cloud
  
def ishData2struct(gzfile,station_id,station_name):
    """ 
    Convert the data into a useful structure array (dictionary)
    """      
    ### Input Variables
    #gzfile = './rawdata/722420-12923-2012.gz'
    #station_id = '722420-12923'
    #station_name = 'test'
    ###
    
    print 'Reading ISH gz file: %s...' % gzfile
       
    ishdata = readRawGZdata(gzfile)
    
    # Create the base structure (dictionary)
    station={}
    station['StationName']=station_name
    station['StationID']=station_id
    station['Latitude']=[]
    station['Longitude']=[]
    station['Tair'] = {'Data':[],'Time':[],'Units':'Celsius','Longname':'Air Temperature','Height':[],'TimeUnits':'minutes since 1970-01-01 00:00:00'}
    station['Pair'] = {'Data':[],'Time':[],'Units':'hPa','Longname':'Air Pressure','Height':[],'TimeUnits':'minutes since 1970-01-01 00:00:00'}
    station['Uwind'] = {'Data':[],'Time':[],'Units':'m s-1','Longname':'Eastward wind velocity component','Height':[],'TimeUnits':'minutes since 1970-01-01 00:00:00'}
    station['Vwind'] = {'Data':[],'Time':[],'Units':'m s-1','Longname':'Northward wind velocity component','Height':[],'TimeUnits':'minutes since 1970-01-01 00:00:00'}
    station['RH'] = {'Data':[],'Time':[],'Units':'percent','Longname':'Relative Humidity','Height':[],'TimeUnits':'minutes since 1970-01-01 00:00:00'}
    station['cloud'] = {'Data':[],'Time':[],'Units':'dimensionless','Longname':'Cloud cover fraction','TimeUnits':'minutes since 1970-01-01 00:00:00'}
    station['rain'] = {'Data':[],'Time':[],'Units':'kg m2 s-1','Longname':'rain fall rate','TimeUnits':'minutes since 1970-01-01 00:00:00'}
    
    # Loop through the ishdata
    ii=0
    basetime = datetime(1970,1,1)
    for dd in ishdata:
        if ii == 0:
            station['Latitude'] = dd['lat']
            station['Longitude'] = dd['lon']
        # Convert the time to a suitable time format
        t = datetime.strptime(dd['yyyymmdd']+dd['hhmm'],'%Y%m%d%H%M')
        dt = t-basetime
        tobs = dt.total_seconds()/60.0
        #       tobs = (t.toordinal()-basetime.toordinal())*1440.0
        
        # Check for missing values here
        
        # Append the required data fields
        station['Tair']['Data'].append(dd['airtemp'])
        station['Tair']['Time'].append(tobs)
        if ii == 0:
            station['Tair']['Height'] = dd['ele']
            
        station['Pair']['Data'].append(dd['airpres'])
        station['Pair']['Time'].append(tobs)
        if ii == 0:
            station['Pair']['Height'] = dd['ele']
        
        if np.size(dd['windspd']) == np.size(dd['winddir']):
            spd=dd['windspd']
            if spd >= 999: spd=np.NaN
            dirn=dd['winddir']
            if dirn >= 999: dirn=np.NaN
            # Note the flip of direction to go to cartesian vectors
            theta = np.mod(dirn-180,360)
            [Uwind,Vwind] = airsea.convertSpeedDirn(theta,spd)
            station['Uwind']['Data'].append(Uwind)
            station['Uwind']['Time'].append(tobs)
            if ii == 0:
                station['Uwind']['Height'] = dd['ele']
            
            station['Vwind']['Data'].append(Vwind)
            station['Vwind']['Time'].append(tobs)
            if ii == 0:
                station['Vwind']['Height'] = dd['ele']
            
        rh = returnHumidity(dd)     
        if np.size(rh) > 0:
           station['RH']['Data'].append(rh)
           station['RH']['Time'].append(tobs)
        if ii == 0:
            station['RH']['Height'] = dd['ele']
           
        rain = returnRainfall(dd)
        if np.size(rain) > 0:
           station['rain']['Data'].append(rain)
           station['rain']['Time'].append(tobs)
               
        cloud = returnCloudCover(dd)
        if np.size(cloud) > 0:
            station['cloud']['Data'].append(cloud)
            station['cloud']['Time'].append(tobs)
        
    # Print a summary
    print 'File Summary:'
    print   'No. of Tair observations - %d' % len(station['Tair']['Data'])       
    print   'No. of Pair observations - %d' % len(station['Pair']['Data']) 
    print   'No. of Uwind observations - %d' % len(station['Uwind']['Data']) 
    print   'No. of Vwind observations - %d' % len(station['Vwind']['Data']) 
    print   'No. of RH observations - %d' % len(station['RH']['Data']) 
    print   'No. of rain observations - %d' % len(station['rain']['Data']) 
    print   'No. of cloud observations - %d' % len(station['cloud']['Data']) 
    print '############# File read successfully################\n'
        
    return station
# End of function    
def getFileNames(latlon,yearrange):
    """ 
    Function to retrieve station names to download from the ftp site (see ftplib.retrievefile)
    """  
    # Read the station metadata CSV file
    csvfile = 'C:/Projects/GOMGalveston/CODE/PYTHON/NOAAWeather/ish-history.csv'
    data = stationMeta(csvfile)
    # Create a list of dictionaries each containing: years, filenames and station name
    stations = []
    for la, lo, sta_id, s, e, name  in zip(data['lat'],data['lon'],data['station_id'],data['start'],data['end'],data['name']):
       if (la >= latlon[2]) & (la <= latlon[3]) & (lo >=latlon[0]) & (lo <= latlon[1]):
           if len(e) > 0 and len(s) > 0:
               years = range(max(int(s),yearrange[0]),min(int(e),yearrange[1])+1)
               if len(years)>0:
                   filenames = []
                   for yy in years:
                       filenames.append(sta_id+'-'+str(yy)+'.gz')
                   dd = {'station_id':sta_id,'station_name':name,'years':years,'filenames':filenames}
                   stations.append(dd)
                   
    return stations
# End of function

def dataQC(data,varnames):
    """ Perform basic quality control on the raw data pass
     # Checks:
    # 1) Lat/lon isn't empty
    # 2) Variable isn't empty
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
                    
    return data
    
    
###


