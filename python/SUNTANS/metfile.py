# -*- coding: utf-8 -*-
"""

Tools for reading and writing a suntans meterological file

Created on Fri Oct 26 09:17:20 2012

@author: mrayson
"""

from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
import numpy as np
import othertime as otime

import pdb

class metfile(object):    
    """
    Class for handling a meteorological input file
    """
    
    mode = 'load'
    
    def __init__(self,infile=None,**kwargs):

        self.__dict__.update(kwargs)
        
        if self.mode == 'create':
            # Loads the variables into an object
            self.Uwind = metdata('Uwind',mode='create',longname='Eastward wind velocity component',units='m s-1')
            self.Vwind = metdata('Vwind',mode='create',longname='Northward wind velocity component',units='m s-1')
            self.Tair = metdata('Tair',mode='create',longname='Air temperature',units='degrees C')
            self.Pair = metdata('Pair',mode='create',longname='Air pressure',units='millibar')
            self.RH = metdata('RH',mode='create',longname='Relative humidity',units='Percent')
            self.rain = metdata('rain',mode='create',longname='Rain fall rate',units='kg m-2 s-1')
            self.cloud = metdata('cloud',mode='create',longname='Cloud cover fraction',units='Fraction (0-1)')
            
        elif self.mode == 'load':
            self.Uwind = metdata('Uwind',infile=infile)
            self.Vwind = metdata('Vwind',infile=infile)
            self.Tair = metdata('Tair',infile=infile)
            self.Pair = metdata('Pair',infile=infile)
            self.RH = metdata('RH',infile=infile)
            self.rain = metdata('rain',infile=infile)
            self.cloud = metdata('cloud',infile=infile)
            
    def __getitem__(self,y):
        x = self.__dict__.__getitem__(y)
        return x.data
    
class SunMet(metfile):    
    """
    Class for creating and writing a SUNTANS met file
    
    Assumes all variables are on the same grid
    """
    
    def __init__(self,x,y,z,timeinfo,**kwargs):
        
        metfile.__init__(self,mode='create')
        
        self.x = x
        self.y = y
        self.z = z        
        self.time = otime.TimeVector(timeinfo[0],timeinfo[1],timeinfo[2])
        self.nctime = otime.SecondsSince(self.time)
        
        self.varnames = ['Uwind','Vwind','Tair','Pair','RH','rain','cloud']
        
        # Update all of the metdata objects
        for vv in self.varnames:
            self.updateMetData(self[vv])     
        
    def updateMetData(self,metdataobject):
        """
        Updates the array in the metdata object
        """
        metdataobject.Nt = self.time.shape[0]
        metdataobject.Npt =np.size(self.x)
        
        metdataobject.data = np.zeros((metdataobject.Nt,metdataobject.Npt))
        metdataobject.x = np.asarray(self.x)
        metdataobject.y = np.asarray(self.y)
        metdataobject.z = np.asarray(self.z)
        metdataobject.time = self.time
    
    def write2NC(self,ncfile):

        """ Writes the data to a netcdf file"""
        print 'Writing to netcdf file: %s',ncfile
        # Create an output netcdf file
        nc = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
        
        # Define the dimensions
        nc.createDimension('nt',0)
        for vv in self.varnames:
            dimname = 'N'+vv
            dimlength = np.size(self.x)
            nc.createDimension(dimname,dimlength)
            #print '%s, %d' % (dimname, dimlength)
            
        # Create the coordinate variables
        tmpvar=nc.createVariable('Time','f8',('nt',))
        # Write the data
        tmpvar[:] = self.nctime
        tmpvar.long_name = 'time'
        tmpvar.units = 'seconds since 1990-01-01 00:00:00'
        for vv in self.varnames:
            dimname = 'N'+vv
            varx = 'x_'+vv
            vary = 'y_'+vv
            varz = 'z_'+vv
            #print dimname, varx, coords[varx]
            tmpvarx=nc.createVariable(varx,'f8',(dimname,))
            tmpvary=nc.createVariable(vary,'f8',(dimname,))
            tmpvarz=nc.createVariable(varz,'f8',(dimname,))
            tmpvarx[:] = self[vv].x
            tmpvary[:] = self[vv].y
            tmpvarz[:] = self[vv].z
            # Create the attributes
            tmpvarx.setncattr('long_name','Longitude at '+vv)
            tmpvarx.setncattr('units','degrees_north')
            tmpvary.setncattr('long_name','Latitude at '+vv)
            tmpvary.setncattr('units','degrees_east')
            tmpvarz.setncattr('long_name','Elevation at '+vv)
            tmpvarz.setncattr('units','m')
            
        # Create the main variables
        for vv in self.varnames:
            dimname = 'N'+vv
            varx = 'x_'+vv
            vary = 'y_'+vv
            tmpvar = nc.createVariable(vv,'f8',('nt',dimname))
            # Write the data
            tmpvar[:] = self[vv].data
            # Write the attributes

            tmpvarz.setncattr('long_name',self[vv].longname)
            tmpvarz.setncattr('units',self[vv].units)
            tmpvar.setncattr('coordinates',varx+','+vary)
            
        nc.close()
        print 'Complete. File written to:\n%s'%ncfile
        
    def __getitem__(self,y):
        x = self.__dict__.__getitem__(y)
        return x
        

    
class metdata(object):
    """
    Class for handling each data type
    
    To create an empty class with variable "Uwind" call:
        >>meto = metdata('Uwind',mode='create')
        
    To load data from a file into the class call:
        >>meto = metdata('Uwind','MetFile.nc')
    """
    
    mode = 'load' 
    longname = ''
    units = ''
    
    def __init__(self,varname,infile=None,**kwargs):
        
        self.__dict__.update(kwargs)
        self.infile=infile
        self.Name=varname
        
        if self.mode == 'create':
            # Create an empty object 
            self.data=[]
            self.x=[]
            self.y=[]
            self.z=[]
            self.time=[]
            self.longname=self.longname
            self.units=self.units
        elif self.mode == 'load':
            
            self.loadnc()
            
    def loadnc(self):
        """
        Load the variable from the netcdf file into the object
        """
        nc = Dataset(self.infile, 'r')
        
        # Load the coordinates
        vv = self.Name
        varx = 'x_'+vv
        vary = 'y_'+vv
        varz = 'z_'+vv
        self.x=nc.variables[varx][:]
        self.y=nc.variables[vary][:]
        self.z=nc.variables[varz][:]

        t = nc.variables['Time']
        self.time = num2date(t[:],t.units)
         
        # Load the data and attributes
        self.data=nc.variables[vv][:]
        try:
            self.longname=nc.variables[vv].long_name
            self.units=nc.variables[vv].units
        except:
            self.longname=''
            self.units=''                
        
        nc.close()
    
    def plot(self,site=0,rangeplot=True):
        """ 
        Time-series plotting method
        """        
        self.fig = plt.gcf()        
        plt.plot(self.time,self.data[:,site])
        if rangeplot:
            upper = np.max(self.data,1)
            lower = np.min(self.data,1)
            plt.fill_between(self.time,lower,y2=upper,color=[0.5, 0.5, 0.5],alpha=0.7)
            
        plt.ylabel('%s [%s]'%(self.Name,self.units))
        plt.xticks(rotation=17)
        plt.title('%s\n x: %10.2f, y: %10.2f, z: %10.2f'%(self.longname,self.x[site],self.y[site],self.z[site]))
        plt.grid(b=True)
        #plt.show()


        