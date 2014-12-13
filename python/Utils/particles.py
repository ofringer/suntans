"""
General module for interacting with particle tracking data

Matt Rayson
Stanford University
October 2014
"""

import os
import numpy as np
from netCDF4 import Dataset, num2date
from datetime import datetime
import othertime
from inpolygon import inpolygon

import pdb

class Particle(object):
    """
    General particle class
    """
    verbose = True
    has_age = False
    basetime = datetime(1990,1,1)

    def __init__(self,N,**kwargs):
        """
        Initialize the arrays
        """
        self.__dict__.update(kwargs)

        self.N = N
        self.X = np.zeros((N,))
        self.Y = np.zeros((N,))
        self.Z = np.zeros((N,))

        if self.has_age:
            self.age = np.zeros((N,))
            self.agemax = np.zeros((N,))

    def update_xyz(self,xnew,ynew,znew):
        """
        Update the particle location vectors

        Use this instead of doing it manually
        """
        self.X[:] = xnew[:]
        self.Y[:] = ynew[:]
        self.Z[:] = znew[:]

    def update_age(self,age):
        """
        Update the particle location vectors

        Use this instead of doing it manually
        """
        self.age[:] = age[:]
        self.agemax[:] = np.max([self.age,self.agemax],axis=0)

    def read_nc(self,tstep,ncfile=None):
        """
        Read the tstep from a netcdf file
        """
        # Check if the file is opened
        if not self.__dict__.has_key('_nc'):
            if ncfile == None:
                raise Exception, 'must set "ncfile" on call to read_nc().'
            else:
                # Open the file for reading
                self._nc = Dataset(ncfile,'r')

        # Read the time data 
        t = self.nc.variables['tp']
        self.time = num2date(t[tstep],t.units)

        self.X = self._nc.variables['xp'][:,tstep]
        self.Y = self._nc.variables['yp'][:,tstep]
        self.Z = self._nc.variables['zp'][:,tstep]

        if self.has_age:
            self.age = self._nc.variables['age'][:,tstep]
            self.agemax = self._nc.variables['agemax'][:,tstep]

    def write_nc(self,time,tstep,ncfile=None):
        """
        Writes the particle locations at the output time step, 'tstep'
        """
        if self.verbose:
            print '\tWriting netcdf output at tstep: %d...'%tstep

        # Check if the file is opened
        if not self.__dict__.has_key('_nc'):
            if ncfile == None:
                raise Exception, 'must set "ncfile" on call to write_nc().'
            else:
                self.create_nc(ncfile)
            
        t = othertime.SecondsSince(time,basetime=self.basetime)

        self._nc.variables['tp'][tstep]=t
        self._nc.variables['xp'][:,tstep]=self.X
        self._nc.variables['yp'][:,tstep]=self.Y
        self._nc.variables['zp'][:,tstep]=self.Z

        if self.has_age:
            self._nc.variables['age'][:,tstep]=self.age
            self._nc.variables['agemax'][:,tstep]=self.agemax

    def create_nc(self,ncfile):
        """
        Create the particle netcdf file

        NetCDF variable and dimension names are consistent with partrac data.
        """

        if self.verbose:
            print '\nInitialising particle netcdf file: %s...\n'%ncfile
            
        # Global Attributes
        nc = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
        nc.Description = 'Particle trajectory file'
        nc.Author = os.getenv('USER')
        nc.Created = datetime.now().isoformat()

        # Dimensions
        nc.createDimension('ntrac', self.N)
        nc.createDimension('nt', 0) # Unlimited
        
        # Create variables
        def create_nc_var( name, dimensions, attdict, dtype='f8'):
            tmp=nc.createVariable(name, dtype, dimensions)
            for aa in attdict.keys():
                tmp.setncattr(aa,attdict[aa])
          
        basetimestr = 'seconds since %s'%(datetime.strftime(self.basetime,\
            '%Y-%m-%d %H:%M:%S'))
        create_nc_var('tp',('nt'),{'units':basetimestr\
            ,'long_name':"time at drifter locations"},dtype='f8')
        create_nc_var('xp',('ntrac','nt'),{'units':'m',\
            'long_name':"Easting coordinate of drifter",'time':'tp'},dtype='f8')
        create_nc_var('yp',('ntrac','nt'),{'units':'m',\
            'long_name':"Northing coordinate of drifter",'time':'tp'},dtype='f8')
        create_nc_var('zp',('ntrac','nt'),{'units':'m',\
            'long_name':"vertical position of drifter (negative is downward from surface)",'time':'tp'},dtype='f8')
        if self.has_age:
            create_nc_var('age',('ntrac','nt'),{'units':'seconds',\
                'long_name':"Particle age",'time':'tp'},dtype='f8')
            create_nc_var('agemax',('ntrac','nt'),{'units':'seconds',\
                'long_name':"Maximum particle age",'time':'tp'},dtype='f8')

        # Keep the pointer to the open file as an attribute
        self._nc = nc
      
class ParticleAge(Particle):
    """
    Class for calculating particle age inside of a polygon
    """
    def __init__(self,agepoly,N):

        Particle.__init__(self,N,has_age=True)

        self.agepoly = agepoly

    def update_age(self,x,y,z,dt):
        """
        Find which particles are in the polygon and update their age
        """
        
        self.update_xyz(x,y,z)

        inpoly = self.inpolygon()
        
        age = self.age
        age[inpoly] += dt
        age[inpoly==False] = 0.

        Particle.update_age(self,age)

    def inpolygon(self):
        """
        Inpolygon logical test
        """
        xy = np.vstack((self.X,self.Y)).T

        return inpolygon(xy, self.agepoly)
        
