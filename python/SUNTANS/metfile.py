# -*- coding: utf-8 -*-
"""

Tools for reading and writing a suntans meterological file

Created on Fri Oct 26 09:17:20 2012

@author: mrayson
"""

from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt

class metfile(object):    
    """
    Class for handling a meteorological input file
    """
    
    mode = 'load'
    
    def __init__(self,infile=None,**kwargs):

        if self.mode == 'create':
            
            self.Uwind = metdata('Uwind',mode='create')
            self.Vwind = metdata('Vwind',mode='create')
            self.Tair = metdata('Tair',mode='create')
            self.Pair = metdata('Pair',mode='create')
            self.RH = metdata('RH',mode='create')
            self.rain = metdata('rain',mode='create')
            self.cloud = metdata('cloud',mode='create')
            
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
       
class metdata(object):
    """
    Class for handling each data type
    
    To create an empty class with variable "Uwind" call:
        >>meto = metdata('Uwind',mode='create')
        
    To load data from a file into the class call:
        >>meto = metdata('Uwind','MetFile.nc')
    """
    
    mode = 'load' 
    
    def __init__(self,varname,infile=None,**kwargs):
        
        self.__dict__.update(kwargs)
        self.infile=infile
        self.Name=varname
        
        if self.mode == 'create':
            # Create an empty object 
            self.data=[]
            self.x=[]
            self.y=[]
            self.time=[]
            self.longname=''
            self.units=''
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
            self.long_name=nc.variables[vv].long_name
            self.units=nc.variables[vv].units
        except:
            self.long_name=''
            self.units=''                
        
        nc.close()
    
    def plot(self,site=0):
        """ 
        Time-series plotting method
        """        
        self.fig = plt.gcf()        
        plt.plot(self.time,self.data[:,site])
        plt.ylabel('%s [%s]'%(self.Name,self.units))
        plt.title('%s\n x: %10.2f, y: %10.2f, z: %10.2f'%(self.long_name,self.x[site],self.y[site],self.z[site]))
        plt.grid(b=True)
        plt.show()
        