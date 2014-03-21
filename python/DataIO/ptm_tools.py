"""
Tools for reading and plotting output from Ed's PTM model
"""

import os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import pdb

class PtmBin(object):
    def __init__(self,fn,release_name=None):
        self.fn = fn

        if release_name is None:
            release_name = os.path.basename(fn)
            release_name = release_name.replace("_bin.out","")
        self.release = release_name

        self.fp = open(self.fn,'rb')
        self.fn_bytes = os.stat(self.fn).st_size

        self.read_bin_header()
        
        # File Format:
        #  int32: Nattr number of attributes
        #  Nattr*[int32 char80 char80]: attribute index, type, name
        #  Ntimesteps* {
        #    6*int32: year, month, day, hour, minute, Npart(t)
        #    Npart(t)* {
        #       int32, 3*float64, int32: id, xyz, active

        self.offsets = {} # map timestep => start of date header for that timestep
        self.offsets[0] = self.fp.tell()

        # Get the time information
        self.getTime()
        
    def read_bin_header(self):
        
        self.Nattr = int( np.fromstring(self.fp.read(4),np.int32) )

        # print "Nattr: ",self.Nattr

        atts = []
        for i in range(self.Nattr):
            idx = int( np.fromstring( self.fp.read(4), np.int32) )
            type_str = self.fp.read(80).strip()
            name_str = self.fp.read(80).strip()
            atts.append( (idx,type_str,name_str) )

    def scan_to_timestep(self,ts):
        """ Return true if successful, False if ts is beyond end of file.
        Set the file pointer to the beginning of the requested timestep.
        if the beginning of that timestep is at or beyond the end of the file
        return False, signifying that ts does not exist.
        """
        if not self.offsets.has_key(ts):
            for ts_scan in range(1,ts+1):
                if not self.offsets.has_key(ts_scan):
                    # if we don't have the offset of this step, go to the one
                    # before, and find out how big the previous frame was.
                    self.fp.seek( self.offsets[ts_scan-1])
                    tstep_header = np.fromstring( self.fp.read( 6*4 ), np.int32 )
                    Npart = tstep_header[5]
                    # print "Step %d has %d particles"%(ts_scan-1,Npart)
                    frame = 6*4 + Npart * (2*4 + 3*8)
                    self.offsets[ts_scan] = self.offsets[ts_scan-1] + frame
                    if self.offsets[ts_scan] >= self.fn_bytes:
                        #print "Hit end of file"
                        return False
        if self.offsets[ts] >= self.fn_bytes:
            return False

        self.fp.seek(self.offsets[ts])
        return True

    def count_timesteps(self):
        saved_pos = self.fp.tell()

        valid_ts = -1
        while 1:
            if self.scan_to_timestep(valid_ts+1):
                # next one is valid, keep going
                valid_ts += 1
            else:
                # valid_ts+1 doesn't exist, so valid_ts is the last valid timestep
                break
            
        self.fp.seek(saved_pos)
        # possible that this is 0!
        return valid_ts + 1

    def dt_seconds(self):
        dnum1,data = self.read_timestep(0)
        dnum2,data = self.read_timestep(1)
        return (dnum2-dnum1)*86400
            
    def read_timestep(self,ts=0):
        """ returns a datenum and the particle array
        """
        if not self.scan_to_timestep(ts):
            return None,None
        
        # Read the time
        dnum,Npart = self.readTime()

        part_dtype = [('id','i4'),
                      ('x','3f8'),
                      ('active','i4')]
        part_size = 2*4 + 3*8

        # print "reading %d particles"%Npart
        
        data = np.fromstring( self.fp.read( part_size * Npart), dtype=part_dtype)
        return dnum,data

    def readTime(self):
        """
        Reads the time header for one step and returns a datetime object
        """
        tstep_header = np.fromstring( self.fp.read( 6*4 ), np.int32 )

        year,month,day,hour,minute,Npart = tstep_header

        if minute == 60:
            hour += 1
            minute = 0
        if hour == 24:
            hour = 0
            day += 1
        
        return datetime(year,month,day,hour,minute),Npart


    def getTime(self):
        """
        Returns a list of datetime objects
        """
        self.nt = self.count_timesteps()
        self.time=[]
        for ts in range(self.nt):
            self.scan_to_timestep(ts)
            t,npart = self.readTime()
            self.time.append(t)

    def plot(self,ts,ax=None,xlims=None,ylims=None,fontcolor='k',\
        marker='.',color='m',**kwargs):
        """
        Plots the current time step
        """
        
        # Check for the plot handle
        if not self.__dict__.has_key('p_handle'):
            # Initialize the plot
            if ax==None:
                ax = plt.gca()
            h1 = ax.plot([],[],marker=marker,linestyle='None',color=color,**kwargs) 
            self.p_handle=h1[0]
            self.title = ax.set_title("",fontdict={'color':fontcolor})

        # Now just update the plot
        t,parts = self.read_timestep(ts=ts)
        x = parts['x'][:,0]
        y = parts['x'][:,1]
        if xlims == None:
            xlims = [x.min(), x.max()]
        if ylims == None:
            ylims = [y.min(), y.max()]

        self.p_handle.set_xdata(x)
        self.p_handle.set_ydata(y)
        self.title=ax.set_title('Particle positions at %s'%(datetime.strftime(t,'%Y-%m-%d %H:%M:%S')))
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_aspect('equal')


#hydrofile = '../InputFiles/untrim_hydro.nc'
#ptmfile = '../InputFiles/line_specify_bin.out'
#outfile = '../InputFiles/FISH_PTM.mov'
#
## Load the particle binary file
#pr = PtmBin(ptmfile)
#
#fig=plt.figure()
#ax=plt.gca()
#
## Load and plot the grid
#
#pr.plot(0,ax=ax)
#def updateLocation(ii):
#    pr.plot(ii,ax=ax)
#    return(pr.p_handle,pr.title)
#
#anim = animation.FuncAnimation(fig, updateLocation,\
#    frames=pr.nt, interval=2, blit=True)
#plt.show()
