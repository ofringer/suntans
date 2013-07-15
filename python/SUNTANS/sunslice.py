# -*- coding: utf-8 -*-
"""
Object for vertical slicing SUNTANS model output

Created on Wed Jul 03 12:29:14 2013

@author: mrayson
"""

from sunpy import Spatial
from trisearch import TriSearch
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.animation as animation

import pdb

class Slice(Spatial):
    """
    Suntans vertical slice class
    """
    
    def __init__(self,ncfile,xpt=None,ypt=None,Npt=100):
        
        Spatial.__init__(self,ncfile,klayer=[-99])
        
        # Calculate the horizontal coordinates of the slice
        self.Npt = Npt
        if xpt == None or ypt == None:
            self._getXYgraphically()
        else:
            self.xpt=xpt
            self.ypt=ypt
            self._getSliceCoords()
            
        # Initialise the slice interpolation object
        self._initInterp()
        
    def __call__(self,variable,tstep,method='nearest'):
        """
        Load the data and interpolate on the slice
        """
        self.tstep=tstep
        try:
            self.Ntslice = len(tstep)
        except:
            tstep=[tstep]
            self.Ntslice = 1
            
        
        self.data = self.interpSlice(variable,method=method)
        
        return self.data.squeeze()
        
        
    def pcolorslice(self,t=0,xaxis='xslice',titlestr=None,bathyoverlay=True,**kwargs):
        """
        Pcolor plot of the slice
        
        Returns a handle to the pcolor object and the colorbar
        """
        # Find the colorbar limits if unspecified

            
        a=self.data[t,:].squeeze()
        am = np.ma.array (a, mask=np.isnan(a))
        
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(am))
            self.clim.append(np.max(am))
        
        h1 = plt.pcolor(self[xaxis],self.zslice,am,vmin=self.clim[0],vmax=self.clim[1],**kwargs)
        
        #Overlay the bed
	if bathyoverlay:
	    self._overlayBathy(self[xaxis][0,:],facecolor=[0.5,0.5,0.5])
        
        # Set labels etc
        plt.xlabel(self._getXlabel(xaxis))
        plt.ylabel('Depth [m]')
        
        plt.xlim([self[xaxis].min(),self[xaxis].max()])
        plt.ylim([self.hslice.min(),0])
        
        axcb = plt.colorbar(h1)
        
        
        if titlestr==None:
            plt.title(self.__genTitle())
        else:
            plt.title(titlestr)
        
        return h1, axcb
            
    def contourslice(self,t=0,xaxis='xslice',clevs=20,titlestr=None,bathyoverlay=True,**kwargs):
        """
        Filled-contour plot of the slice
        
        Returns a handle to the pcolor object and the colorbar
        """
          
        a=self.data[t,:].squeeze()
        am = np.ma.array (a, mask=np.isnan(a))
        
        # Find the colorbar limits if unspecified
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(am))
            self.clim.append(np.max(am))
        
        V = np.linspace(self.clim[0],self.clim[1],clevs)
        h1 = plt.contourf(self[xaxis],self.zslice,am,V,vmin=self.clim[0],vmax=self.clim[1],**kwargs)
        
        #Overlay the bed
	if bathyoverlay:
	    self._overlayBathy(self[xaxis][0,:],facecolor=[0.5,0.5,0.5])
        
        # Set labels etc
        plt.xlabel(self._getXlabel(xaxis))
        plt.ylabel('Depth [m]')
        
        plt.xlim([self[xaxis].min(),self[xaxis].max()])
        plt.ylim([self.hslice.min(),0])
        
        axcb = plt.colorbar(h1)
        
        
        if titlestr==None:
            plt.title(self.__genTitle())
        else:
            plt.title(titlestr)
        
        return h1, axcb
    
    def xtplot(self,zlayer=0,xaxis='xslice',clevs=20,titlestr=None,**kwargs):    
        """
        x-t contour plot of the sliced variable along vertical layer, 'zlayer'. 
        
        zlayer can be:            
            [0 - Nkmax] - vertical layer number
            'seabed' - seabed value
            'diff' - top minus bottom difference)
                       
        """
        
	#kbed = np.max(self.Nk[self.cellind]-1,0)
	kbed = self.Nk[self.cellind]
        if zlayer == 'seabed':
            a= self.data[:,kbed,range(0,self.Npt)]
            zstring = 'seabed'
            
        elif zlayer == 'diff':
            atop = self.data[:,0,:]    
            abot = self.data[:,kbed,range(0,self.Npt)]
            a = atop - abot
            zstring = 'Surface value - seabed value'
            
        else:
            a = self.data[:,zlayer,:]
            zstring = '%3.1f [m]'%self.z_r[zlayer]
            
        am = np.ma.array (a, mask=np.isnan(a))
        
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(am))
            self.clim.append(np.max(am))
        
        V = np.linspace(self.clim[0],self.clim[1],clevs)
        
        h1 = plt.contourf(self[xaxis][0,:],self.time[self.tstep],am,V,vmin=self.clim[0],vmax=self.clim[1],**kwargs)
        
        plt.xlabel(self._getXlabel(xaxis))
        plt.ylabel('Time')
            
        axcb = plt.colorbar(h1)
        
        titlestr='%s [%s]\nLayer: %s'%(self.long_name,self.units,zstring)
        plt.title(titlestr)
        
        return h1, axcb
    
    def plotslice(self):
        """
        Plots the slice location on a map
        """
        self.contourf(z=-self.dv,clevs=30,titlestr='',cmap='gist_earth')
        plt.plot(self.xslice[0,:],self.yslice[0,:],'m--')
        
    def animate(self,**kwargs):
        """
        Animates the slice
        """
        
        # Initialise the plot object
        h1,cb = self.pcolorslice(**kwargs)
        
        fig = plt.gcf()
        ax = fig.gca()
        title=ax.set_title("")
        
        def updateScalar(ii):
            a=self.data[ii,:].squeeze()
            am = np.ma.array (a, mask=np.isnan(a))
            h1.set_array(am[am.mask==False])
            title.set_text(self.__genTitle(tt=self.tstep[ii]))
            
            return (title,h1)
            
        self.anim = animation.FuncAnimation(fig, updateScalar, frames=len(self.tstep), interval=50, blit=True)
                
    def interpSlice(self,variable,method='linear'):
        """
        Interpolates the data in raw data onto the slice array
        """
        tstep = self.tstep
        slicedata = np.zeros((self.Ntslice,self.Nkmax,self.Npt))

	#if method=='linear':
	#    cellind3d = np.repeat(self.cellind.reshape((1,self.Npt)),self.Nkmax,axis=0)
	#    k3d = np.arange(0,self.Nkmax)
	#    k3d = np.repeat(k3d.reshape((self.Nkmax,1)),self.Npt,axis=1)

        for tt in range(self.Ntslice):
            if self.Ntslice>1:
                print 'Slicing data at time-step: %d of %d...'%(tt,self.Ntslice)
            
            self.tstep=[tt]
            rawdata = self.loadData(variable=variable)

#	    if method == 'nearest':
#		for kk in range(self.Nkmax):
#                    slicedata[tt,kk,:] = rawdata[kk,self.cellind]
#	    elif method == 'linear':
#		slicedata[tt,:] = self.interpLinear(rawdata.ravel(),self.xslice.ravel(),self.yslice.ravel(),cellind3d.ravel(),k=k3d.ravel()).reshape((self.Nkmax,self.Npt))
#	    else:
#		raise Exception, ' unknown interpolation method: %s. Must be "nearest" or "linear"'%method
 
            for kk in range(self.Nkmax):
                if method == 'nearest':
                    slicedata[tt,kk,:] = rawdata[kk,self.cellind]
                elif method == 'linear':
                    slicedata[tt,kk,:] = self.interpLinear(rawdata[kk,:].squeeze(),self.xslice[0,:],self.yslice[0,:],self.cellind,k=kk)
                else:
                    raise Exception, ' unknown interpolation method: %s. Must be "nearest" or "linear"'%method
            
        mask = self.maskslice.reshape((1,self.Nkmax,self.Npt))
        mask = mask==False
        mask = mask.repeat(self.Ntslice,axis=0)
        slicedata[mask] = np.nan
        
        self.tstep=tstep
        
        return slicedata
        
        
    def _initInterp(self):
        """
        Initialise the interpolant
        
        Finds the horizontal indices of the slice points and 
        constructs the 3D mask array
        """
        
        # Find the cell index of each point along the slice
        self.Tri = TriSearch(self.xp,self.yp,self.cells)
        
        self.cellind = self.Tri(self.xslice,self.yslice)
        
        # Construct the 3D coordinate arrays
        self.xslice = np.repeat(self.xslice.reshape((1,self.Npt)),self.Nkmax,axis=0)
        self.yslice = np.repeat(self.yslice.reshape((1,self.Npt)),self.Nkmax,axis=0)
        self.distslice = np.repeat(self.distslice.reshape((1,self.Npt)),self.Nkmax,axis=0)
        self.zslice = np.repeat(-self.z_r.reshape((self.Nkmax,1)),self.Npt,axis=1)
        
        # Construct the mask array
        self.maskslice = np.zeros((self.Nkmax,self.Npt),dtype=np.bool)
        
        for kk in range(self.Nkmax):
            for ii in range(self.Npt):
                if kk <= self.Nk[self.cellind[ii]]:
                    self.maskslice[kk,ii]=True

        # Get the bathymetry along the slice
        self.hslice = -self.dv[self.cellind]
    
    def _getSliceCoords(self):
        """
        Fits a spline through the input slice points
        """
        n = self.xpt.shape[0]
        t = np.linspace(0,1,n)
        tnew = np.linspace(0,1,self.Npt)
        
        if n <= 3:
            kind='linear'
        else:
            kind=3
            
        Fx = interp1d(t,self.xpt,kind=kind)
        Fy = interp1d(t,self.ypt,kind=kind)
        
        self.xslice = Fx(tnew)
        self.yslice = Fy(tnew)
        
        # Calculate the distance along the slice
        self.distslice = np.zeros_like(self.xslice)
        self.distslice[1:] = np.sqrt( (self.xslice[1:]-self.xslice[:-1]) **2 + \
           (self.yslice[1:]-self.yslice[:-1]) **2 )
        
        self.distslice = self.distslice.cumsum()
        
    def _getXYgraphically(self):
        """
        Plot a map of the bathymetry and graphically select the slice points
        """
        self.contourf(z=-self.dv,clevs=30,titlestr='Select points for slice on map\nRight-click to finish; middle-click to remove last point.',cmap='gist_earth')
        x = plt.ginput(n=0,timeout=0,mouse_pop=2,mouse_stop=3)
        
        # Find the location of the slice
        xy = np.array(x)
        self.xpt=xy[:,0]
        self.ypt=xy[:,1]
        self._getSliceCoords()
        
        plt.plot(self.xslice,self.yslice,'m--')
        
        plt.title('Close figure to continue...')
        plt.show()
    
    def _getXlabel(self,xaxis):
        if xaxis == 'xslice':
            xlab = 'Easting [m]'
        elif xaxis == 'yslice':
            xlab = 'Northing [m]'
        elif xaxis == 'distslice':
            xlab = 'Distance along transect [m]'
        else:
            raise Exception, ' unknown "xaxis" value %d.\n Must be one of "xslice", "yslice" or "distslice".'%xaxis
            
        return xlab
        
    def _overlayBathy(self,xdata,**kwargs):
        """
        Pretty bathymetry overlay
        """
        
        plt.fill_between(xdata,self.hslice,y2=self.hslice.min(),**kwargs)
        
        
    def __genTitle(self,tt=None):
        
        if tt ==None:
            if type(self.tstep)==int:
                tt = self.tstep
            else:
                tt = self.tstep[0]
                
        titlestr='%s [%s]\nTime: %s'%(self.long_name,self.units,\
                datetime.strftime(self.time[tt],'%d-%b-%Y %H:%M:%S'))
                
        return titlestr
        
        
#####
## Testing data
#
####
## Inputs
#ncfile = 'C:/Projects/GOMGalveston/MODELLING/ForcingSensitivity/rundata/GalvCoarse_AprMay_2010_TWRH_007.nc'
#
#xpt = np.array([ 334640.1250722 ,  331097.21398258,  327837.73578013,
#        324861.69046485,  323444.526029  ,  321177.06293164,
#        319334.74916504,  318059.30117277,  319618.18205221,
#        322310.79448032,  326987.43711862,  329113.18377239,
#        332939.52774918])
#        
#ypt = np.array([ 3247017.02417632,  3247867.32283783,  3247867.32283783,
#        3250559.93526594,  3254386.27924273,  3259488.07121179,
#        3263456.13163216,  3267140.75916537,  3272525.98402159,
#        3275785.46222404,  3281170.68708027,  3287264.49415442,
#        3292649.71901064])
#
#outfile = 'TestSlice.mov'
####
#sun = Slice(ncfile,xpt=xpt,ypt=ypt,Npt=500)
#
#temp = sun('salt',range(23),method='linear')
#
## Test ploting
##sun.pcolors(t=1,xaxis='distslice')
##plt.show()
#
## Test animation
##sun.animate(xaxis='yslice')
##sun.saveanim(outfile)
##plt.show()
#
## Test the xtplot
##sun.xtplot(zlayer='seabed',xaxis='yslice')
##plt.show()
#
#sun.plotslice()
#plt.show()



#
#tempnodes = sun.cell2nodeki(temp[0,:],sun.cellind,0*np.ones_like(sun.cellind))

#cell_scalar = sun.loadData(variable='temp')
## Get the values at the nodes
#for kk in range(30):
#    print kk
#    tempnode = sun.cell2nodek(cell_scalar[kk,:],k=kk)
    


# Gradient routine
#dH_dx,dH_dy = sun.gradH(sun.dv,0)
##
#sun.plot(dH_dx,titlestr='')
#plt.show()

# Test out an interpolation tool
#temp = sun('salt',1)

#dT_dx,dT_dy = sun.gradH(cell_scalar[0,:],cellind=sun.cellind,k=0)

########################################################
# Deprecated functions
########################################################
#    def cell2nodek(self,cell_scalar,k=0):
#        """
#        ### NOT USED ###
#        
#        Map a cell-based scalar onto a node
#        
#        Uses sparse matrices to do the heavy lifting
#        """
#        if not self.__dict__.has_key('_datasparse'):
#            self._datasparse=[]
#            self._Asparse=[]
#            self._indsparse=[]
#            for kk in range(self.Nkmax):
#                self._datasparse.append(sparse.dok_matrix((self.Np,self.Nc),dtype=np.double))
#                self._Asparse.append(sparse.dok_matrix((self.Np,self.Nc),dtype=np.double))
#                self._indsparse.append(sparse.dok_matrix((self.Np,self.Nc),dtype=np.int))
#                
#                for i in range(self.Nc):
#                    for j in range(3):
#                        if kk <= self.Nk[i]:
#                            self._Asparse[kk][self.cells[i,j],i] = self.Ac[i]
#                
#                self._Asparse[kk]=self._Asparse[kk].tocoo()
#        
#        for i in range(self.Nc):
#            for j in range(3):
#                if k <= self.Nk[i]:
#                    self._datasparse[k][self.cells[i,j],i] = cell_scalar[i]
#                    #self._Asparse[k][self.cells[i,j],i] = self.Ac[i]
#                    #self._indsparse[k][self.cells[i,j],i] = i
#                
#        node_scalar = self._datasparse[k].tocoo().multiply(self._Asparse[k]).sum(axis=1) / self._Asparse[k].sum(axis=1)
#        
#        return np.array(node_scalar).squeeze()
#    def cell2nodekold(self,cell_scalar,k=0):
#        """
#        Map a cell-based scalar onto a node
#        
#        This is calculated via a mean of the cells connected to a node(point)
#        """
#        
#        # Area weighted interpolation
#        node_scalar = [np.sum(cell_scalar[self.pnt2cellsk(k,ii)]*self.Ac[self.pnt2cellsk(k,ii)])\
#            / np.sum( self.Ac[self.pnt2cellsk(k,ii)]) for ii in range(self.Np)]
#        return np.array(node_scalar)
#    
#    def cell2nodeki(self,cell_scalar,cellind,k):
#        """
#        NOT WORKING
#
#        ## Needs an additional lookup index to go from local cell index to the main 
#        grid cell index        
#        
#        Map a cell-based scalar onto a node
#        
#        This returns the nodes values of the nodes connected to cells, cellind at
#        vertical level, k
#        
#        This is useful for finding the nodal values for a limited number of points, i.e., 
#        during interpolation
#        """
#        ind = self.cells[cellind,:]
#        
#        node_scalar = np.zeros_like(ind)
#        
#        ni = cellind.shape[0]
#        
#        for ii in range(ni):
#            for jj in range(3):
#                node_scalar[ii,jj] = np.sum(cell_scalar[self.pnt2cellsk(k[ii],ind[ii,jj])]*self.Ac[self.pnt2cellsk(k[ii],ind[ii,jj])])\
#            / np.sum( self.Ac[self.pnt2cellsk(k[ii],ind[ii,jj])]) 
#                
#        
#        return node_scalar
#        
#    def pnt2cellsk(self, k, pnt_i):
#        """
#        Returns the cell indices for a point, pnt_i at level, k
#        
#        (Stolen from Rusty's TriGrid class)
#        """
#        if not self.__dict__.has_key('_pnt2cellsk'):
#            # build hash table for point->cell lookup
#            self._pnt2cellsk = []
#            for k in range(self.Nkmax):
#                self._pnt2cellsk.append({})
#                for i in range(self.Nc):
#                    for j in range(3):
#                        if not self._pnt2cellsk[k].has_key(self.cells[i,j]):
#                            self._pnt2cellsk[k][self.cells[i,j]] = []
#                        if k <= self.Nk[i]:
#                            self._pnt2cellsk[k][self.cells[i,j]].append(i)
#        return self._pnt2cellsk[k][pnt_i]
#        
#    def pnt2cellsparse(self,pnt_i):
#        
#        pdb.set_trace()
#        #self._pnt2cellsk = sparse.dok_matrix((self.Np,self.Nc),dtype=np.int)
#        testsparse = sparse.dok_matrix((self.Np,self.Nc),dtype=np.int)
#        
#        for i in range(self.Nc):
#            for j in range(3):
#                #self._pnt2cells[self.cells[i,j],i] = i
#                testsparse[self.cells[i,j],i] = i
#        
#        pdb.set_trace() 
