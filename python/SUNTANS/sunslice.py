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
from sunpy import unsurf
from hybridgrid import Line
from hybridgrid import Point as GPoint

from gridsearch import GridSearch
from shapely.geometry import LineString, Point

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
            
    def contourslice(self,z,t=0,xaxis='xslice',clevs=20,titlestr=None,bathyoverlay=True,\
        filled = True, outline = False,colorbar=True,**kwargs):
        """
        Filled-contour plot of the slice
        
        Returns a handle to the pcolor object and the colorbar
        """
          
        #a=self.data[t,:].squeeze()
        am = np.ma.array (z, mask=np.isnan(z))
        
        # Find the colorbar limits if unspecified
        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(am))
            self.clim.append(np.max(am))

        klayer,Nkmax = self.get_klayer()
        
        if type(clevs)==type(1): # is integer
            V = np.linspace(self.clim[0],self.clim[1],clevs)
        else:
            V = clevs
         
        if filled:
            h1 = plt.contourf(self[xaxis],-self.z_r[klayer],am,V,vmin=self.clim[0],vmax=self.clim[1],**kwargs)
        
        if outline:
            h2 = plt.contour(self[xaxis],-self.z_r[klayer],am,V,**kwargs)
            
        #Overlay the bed
        if bathyoverlay:
            self._overlayBathy(self[xaxis][:],facecolor=[0.5,0.5,0.5])
        
        # Set labels etc
        plt.xlabel(self._getXlabel(xaxis))
        plt.ylabel('Depth [m]')
        
        plt.xlim([self[xaxis].min(),self[xaxis].max()])
        plt.ylim([self.hslice.min(),0])
        
        if colorbar and filled:
            axcb = plt.colorbar(h1)
        
        
        if titlestr==None:
            plt.title(self.__genTitle())
        else:
            plt.title(titlestr)
        
        if filled and not outline:
            return h1, 
        elif not filled and outline:
            return h2, 
        elif filled and outline:
            return h1, h2, 
    
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
            
            self.tstep=[tstep[tt]]
            rawdata = self.loadData(variable=variable)
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
        
        
    def get_klayer(self):
        if self.klayer[0]==-99:
            klayer=range(self.Nkmax)
            Nkmax = self.Nkmax
        else:
            klayer=self.klayer
            Nkmax=len(klayer)
        return klayer,Nkmax


    def _initInterp(self):
        """
        Initialise the interpolant
        
        Finds the horizontal indices of the slice points and 
        constructs the 3D mask array
        """
        
        # Find the cell index of each point along the slice
        self.Tri = TriSearch(self.xp,self.yp,self.cells)
        
        self.cellind = self.Tri(self.xslice,self.yslice)

        klayer,Nkmax = self.get_klayer()
        
        # Construct the 3D coordinate arrays
        self.xslice = np.repeat(self.xslice.reshape((1,self.Npt)),self.Nkmax,axis=0)
        self.yslice = np.repeat(self.yslice.reshape((1,self.Npt)),self.Nkmax,axis=0)
        self.distslice = np.repeat(self.distslice.reshape((1,self.Npt)),self.Nkmax,axis=0)
        self.zslice = np.repeat(-self.z_r[klayer].reshape((self.Nkmax,1)),self.Npt,axis=1)
        
        # Construct the mask array
        self.calc_mask()

        # Get the bathymetry along the slice
        self.hslice = -self.dv[self.cellind]

    def calc_mask(self):
        """ Construct the mask array"""
        self.maskslice = np.zeros((self.Nkmax,self.Npt),dtype=np.bool)
        
        for kk in range(self.Nkmax):
            for ii in range(self.Npt):
                if kk <= self.Nk[self.cellind[ii]]:
                    self.maskslice[kk,ii]=True
    
    def _getSliceCoords(self,kind=3):
        """
        Fits a spline through the input slice points

        # Kind is the linear interpolation type
        """
        n = self.xpt.shape[0]
        t = np.linspace(0,1,n)
        tnew = np.linspace(0,1,self.Npt)
        
        if n <= 3:
            kind='linear' # Spline won't work with <= 3 points
        else:
            kind=kind
            
        Fx = interp1d(t,self.xpt,kind=kind)
        Fy = interp1d(t,self.ypt,kind=kind)
        
        self.xslice = Fx(tnew)
        self.yslice = Fy(tnew)

        self._getDistCoords()
        
    def _getDistCoords(self):
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
        
        plt.fill_between(xdata,self.hslice,y2=self.hslice.min(),zorder=1e6,**kwargs)
        
        
    def __genTitle(self,tt=None):
        
        if tt ==None:
            if type(self.tstep)==int:
                tt = self.tstep
            else:
                tt = self.tstep[0]
                
        titlestr='%s [%s]\nTime: %s'%(self.long_name,self.units,\
                datetime.strftime(self.time[tt],'%d-%b-%Y %H:%M:%S'))
                
        return titlestr
        
class SliceEdge(Slice):
    """
    Slice suntans edge-based data at all edges near a line

    Used for e.g. flux calculations along a profile
    """

    edgemethod=1
    def __init__(self,ncfile,xpt=None,ypt=None,Npt=100,klayer=[-99],**kwargs):
        
        self.Npt=Npt

        Spatial.__init__(self,ncfile,klayer=klayer,**kwargs)

        # Load the grid as a hybridgrid
        self.grd = GridSearch(self.xp,self.yp,self.cells,nfaces=self.nfaces,\
            edges=self.edges,mark=self.mark,grad=self.grad,neigh=self.neigh,\
                xv=self.xv,yv=self.yv)

        # Find the edge indices along the line
        self.update_xy(xpt,ypt)

    def update_xy(self,xpt,ypt):
        """
        Updates the x and y coordinate info in the object
        """
        if xpt == None or ypt == None:
            self._getXYgraphically()
        else:
            self.xpt=xpt
            self.ypt=ypt

        self._getSliceCoords(kind='linear')
        # List of the edge indices
        self.j,self.nodelist =\
            self.get_edgeindices(self.xslice,self.yslice,method=self.edgemethod)

        self.nslice = len(self.j)

        # Update the x and y axis of the slice
        self.xslice=self.xp[self.nodelist]
        self.yslice=self.yp[self.nodelist]

        self._getDistCoords()

        self.edgexy()

        # The x and y arrys need to be resized
        self.xslice = 0.5*(self.xslice[1:]+self.xslice[0:-1])
        self.yslice = 0.5*(self.yslice[1:]+self.yslice[0:-1])
        self.distslice = 0.5*(self.distslice[1:]+self.distslice[0:-1])

        # Get the mask
        self.calc_mask()

        # Calculate the area
        self.area = self.calc_area()

        # Calculate the normnal
        self.ne1, self.ne2, self.enormal = self.calc_normal(self.nodelist,self.j)

        # Get the bathymetry along the slice
        de = self.get_edgevar(self.dv)
        self.hslice = -de[self.j]


    def loadData(self,variable=None,setunits=True):
        """ 
        Load the specified suntans variable data as a vector

        Overloaded method for edge slicing - it is quicker to load time step by
        time step in a loop.
            
        """

        nc = self.nc

        if variable==None:
            variable=self.variable

        if setunits:
            try:
                self.long_name = nc.variables[variable].long_name
                self.units= nc.variables[variable].units
            except:
                self.long_name = ''
                self.units=''

        j=self.j
        # Check if cell-centered variable
        is3D=True
        isCell=False
        if self.hasVar(variable):
            if self.hasDim(variable,self.griddims['Ne']):
                isCell=False
            elif self.hasDim(variable,self.griddims['Nc']): 
                isCell=True
                nc1 = self.grad[j,0]
                nc2 = self.grad[j,1]
                # check for edges (use logical indexing)
                ind1 = nc1==-1
                nc1[ind1]=nc2[ind1]
                ind2 = nc2==-1
                nc2[ind2]=nc1[ind2]

            # Check if 3D
            if self.hasDim(variable,self.griddims['Nk']): # 3D
                is3D=True
            else:
                is3D=False

        klayer,Nkmax = self.get_klayer()

        def ncload(nc,variable,tt):
            if variable=='agemean':
                ac = nc.variables['agec'][tt,klayer,:]
                aa = nc.variables['agealpha'][tt,klayer,:]
                tmp = aa/ac
                tmp[ac<1e-12]=0.
                return tmp/86400.

            if variable=='area':
                eta = nc.variables['eta'][tt,:]
                dzf = self.getdzf(eta)
                dzf = Spatial.getdzf(self,eta)

                return self.df*dzf

            else:
                if self.hasDim(variable,self.griddims['Nk']): # 3D
                    return nc.variables[variable][tt,klayer,:]
                else:
                    return nc.variables[variable][tt,:]
                
        # For loop where the data is extracted 
        nt = len(self.tstep)
        ne = len(self.j)
        if is3D==True:
            self.data = np.zeros((nt,Nkmax,ne))
        else:
            self.data = np.zeros((nt,ne))

        for ii,tt in enumerate(self.tstep):
            #tmp=nc.variables[variable][tt,:,:]
            tmp = ncload(nc,variable,tt)
            # Return the mean for cell-based variables
            if isCell:
                self.data[ii,...] = 0.5*(tmp[...,nc1]+tmp[...,nc2])
            else:
                self.data[ii,...]=tmp[...,self.j]
            # Mask 3D data
            if is3D:
                maskval=0
                self.data[ii,self.maskslice]=maskval

        #fillval = 999999.0
        #self.mask = self.data==fillval
        #self.data[self.mask]=0.
        self.data[self.data==self._FillValue]=0.
        self.data = self.data.squeeze()
        
        return self.data

    def edgexy(self):
        """
        Nx2 vectors outlining each cell in the edge slice
        """
        def closePoly(xp,node,k):
            return  np.array([ [xp[node],\
                xp[node+1],xp[node+1], xp[node],xp[node]],\
                [-self.z_w[k],-self.z_w[k],-self.z_w[k+1],-self.z_w[k+1],-self.z_w[k]],\
                ]).T

        self.xye = [closePoly(self.distslice,jj,kk) for kk in range(self.Nkmax) \
            for jj in range(len(self.j)) ]

    def calc_normal(self,nodelist,j):
        """
        Calculate the edge normal
        """
        # Calculate the unit normal along the edge
        P1 = GPoint(self.xp[nodelist][0:-1],self.yp[nodelist][0:-1])
        P2 = GPoint(self.xp[nodelist][1:],self.yp[nodelist][1:])
        L = Line(P1,P2)
        ne1,ne2 = L.unitnormal()

        # Compute the unique normal of the dot product
        enormal = np.round(self.n1[j]*ne1 +\
            self.n2[j]*ne2)
        return ne1,ne2,enormal

    def mean(self,phi,axis='time'):
        """
        Calculate the mean of the sliced data along an axis

        axis: time, depth, area
            time : returns the time mean. size= (Nk, Nj)
            depth: returns the time and spatial mean. Size = (Nk)
            area: returns the area mean. Size = (Nt)

        """

        if axis=='time':
            return np.mean(phi,axis=0)
        elif axis=='area':
            area_norm = self.area / self.area.sum()
            return np.sum( np.sum(phi*area_norm,axis=-1),axis=-1)
        elif axis=='depth':
            dx = self.df[self.j]
            dx_norm = dx / dx.sum()
            return np.sum( self.mean(phi,axis='time')*dx_norm,axis=-1)

    def plot(self,z,titlestr=None,**kwargs):
        """
        Pcolor plot of the slice
        """

        if self.clim==None:
            self.clim=[]
            self.clim.append(np.min(z))
            self.clim.append(np.max(z))
        
        # Set the xy limits
        xlims=[self.distslice.min(),self.distslice.max()] 
        ylims=[-self.z_w.max(),-self.z_w.min()]
        
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xye,z.ravel(),xlim=xlims,ylim=ylims,\
            clim=self.clim,**kwargs)

        self.ax.set_aspect('auto')

    def plotedges(self,color='m',**kwargs):
        """
        plot for testing
        """
        self.plotmesh()
        #plt.plot(self.edgeline.xy,'r')
        for ee in self.j:
            plt.plot([self.xp[self.edges[ee,0]],self.xp[self.edges[ee,1]]],\
                [self.yp[self.edges[ee,0]],self.yp[self.edges[ee,1]]],color=color,\
                **kwargs)


    def calc_area(self,eta=None):
        """
        Calculate thee cross-sectional area of each face
        """
        if eta==None:
            eta = np.zeros((self.nslice,)) # Assumes the free-surface is zero

        dzf = self.getdzf(eta)

        area = dzf * self.df[self.j]

        area[self.maskslice]=0
        
        return area

    def getdzf(self,eta):
        """ Get the cell thickness along each edge of the slice"""
        dzf = Spatial.getdzf(self,eta,j=self.j)
        dzf[self.maskslice]=0
        return dzf

    def get_width(self):
        """
        Calculate the width of each edge as a 2d array

        Missing cells are masked
        """
        df = self.df[self.j]
        width =  np.ones((self.Nkmax,1)) * df[np.newaxis,:]
        width[self.maskslice]=0
        return width

    def calc_mask(self):
        """ Construct the mask array"""
        klayer,Nkmax=self.get_klayer()
        self.maskslice = np.zeros((Nkmax,len(self.j)),dtype=np.bool)
        
        for k,kk in enumerate(klayer):
            for ii,j in enumerate(self.j):
                if kk >= self.Nke[j]:
                    self.maskslice[k,ii]=True
 

    def get_edgeindices(self,xpt,ypt,method=1):
        """
        Return the indices of the edges (in order) along the line

        method - method for line finding algorithm
               0 - closest point to line
               1 - closest point without doing a u-turn  
        """
        # Load the line as a shapely object
        #edgeline = asLineString([self.xslice,self.yslice])
        Npt = xpt.shape[0]
        xyline = [(xpt[ii],ypt[ii]) for ii in range(Npt)]
        self.edgeline = LineString(xyline)

        # Find the nearest grid Node to the start and end of the line
        xy_1 = np.vstack((xpt[0],ypt[0])).T
        node0 = self.grd.findnearest(xy_1)
        xy_2 = np.vstack((xpt[-1],ypt[-1])).T
        endnode = self.grd.findnearest(xy_2)

        # This is the list containing all edge nodes
        nodelist = [node0[0]]

        def connecting_nodes(node,nodelist):
            """ finds the nodes connecting to the node"""
            edges = self.grd.pnt2edges(node)
            cnodes = []
            for ee in edges:
                for nn in self.grd.edges[ee]:
                    if nn not in nodelist:
                        cnodes.append(nn)
            return cnodes

        def min_dist(nodes,line):    
            """Returns the index of the node with the minimum distance
            to the line"""

            # Convert all nodes to a point object
            points = [Point((self.xp[nn],self.yp[nn])) for nn in nodes]

            # Calculate the distance
            dist = [line.distance(pp) for pp in points]
            for ii,dd in enumerate(dist):
                if dd == min(dist):
                    return nodes[ii]

        def min_dist_line(cnode,nodes,line):    
            """Returns the index of the node with the minimum distance
            to the line"""

            # Convert all nodes to a point object
            points = [Point((0.5*(self.xp[nn]+self.xp[cnode]),\
                0.5*(self.yp[nn]+self.yp[cnode]))) for nn in nodes]
            #lines = [LineString([(self.xp[cnode],self.yp[cnode]),\
            #    (self.xp[nn],self.yp[nn])]) for nn in nodes]

            # Calculate the distance
            dist = [line.distance(pp) for pp in points]
            for ii,dd in enumerate(dist):
                if dd == min(dist):
                    return nodes[ii]

        def min_dist_angle(cnode,nodes,line):    
            """Returns the index of the node with the minimum distance
            to the line"""

            # Convert all nodes to a point object
            points = [Point((0.5*(self.xp[nn]+self.xp[cnode]),\
                0.5*(self.yp[nn]+self.yp[cnode]))) for nn in nodes]

            # Calculate the distance
            dist = [line.distance(pp) for pp in points]
            dist = np.array(dist)

            # Calculate the angle along the line of the new coordinate
            def calc_ang(x1,x2,y1,y2):
                return np.arctan2( (y2-y1),(x2-x1) )

            angle1 = [calc_ang(self.xp[cnode],self.xp[nn],\
                self.yp[cnode],self.yp[nn]) for nn in nodes]

            # Calculate the heading of the line near the two points
            def calc_heading(P1,P2,L):
                d1 = L.project(P1)
                d2 = L.project(P2)
                if d1 <= d2:
                    P3 = L.interpolate(d1)
                    P4 = L.interpolate(d2)
                else:
                    P3 = L.interpolate(d2)
                    P4 = L.interpolate(d1)

                return calc_ang(P3.xy[0][0],P4.xy[0][0],P3.xy[1][0],P4.xy[1][0])

            P1 = Point((self.xp[cnode],self.yp[cnode]))
            angle2 = [calc_heading(P1,Point( (self.xp[nn],self.yp[nn]) ),line) \
                for nn in nodes]

            angdiff = np.array(angle2) - np.array(angle1)

            # Use the minimum distance unless the point is a u-turn
            rank = np.argsort(dist)
            for nn in range(dist.shape[0]):
                if np.abs(angdiff[rank[nn]]) <= np.pi/2:
                    return nodes[rank[nn]] 
            # if they all u-turn return the min dist
            return nodes[rank[0]]

        # Loop through and find all of the closest points to the line
        MAXITER=10000
        for ii in range(MAXITER):
            cnodes = connecting_nodes(nodelist[-1],nodelist)
            #if method==0:
            #    newnode = min_dist(cnodes,self.edgeline)
            if method==0:
                newnode = min_dist_line(nodelist[-1],cnodes,self.edgeline)
            elif method==1:
                newnode = min_dist_angle(nodelist[-1],cnodes,self.edgeline)
            #print 'Found new node: %d...'%newnode
            if newnode==None:
                break
            if ii>1:
                if self.mark[self.grd.find_edge([newnode,nodelist[-1]])] not in [0,5]:
                    print 'Warning: reached a boundary cell. Aborting edge finding routine'
                    break

            nodelist.append(newnode)
            if newnode == endnode:
                #print 'Reached end node.'
                break
                
        # Return the list of edges connecting all of the nodes
        return [self.grd.find_edge([nodelist[ii],nodelist[ii+1]]) for ii in\
            range(len(nodelist)-1)], nodelist

class MultiSliceEdge(SliceEdge):
    """
    Slice suntans edge-based data at all edges near a line

    Used for e.g. flux calculations along a profile
    """

    def __init__(self,ncfile,xpt=None,ypt=None,Npt=100,klayer=[-99],**kwargs):
        
        self.Npt=Npt

        Spatial.__init__(self,ncfile,klayer=klayer,**kwargs)

        # Load the grid as a hybridgrid
        self.grd = GridSearch(self.xp,self.yp,self.cells,nfaces=self.nfaces,\
            edges=self.edges,mark=self.mark,grad=self.grad,neigh=self.neigh,\
                xv=self.xv,yv=self.yv)

        # Find the edge indices along the line
        self.update_xy(xpt,ypt)

    def update_xy(self,xpts,ypts):
        """
        Updates the x and y coordinate info in the object
        """
        self.slices=[]
        for xpt,ypt in zip(xpts,ypts):
            self.xpt=xpt
            self.ypt=ypt
            self._getSliceCoords(kind='linear')
            # List of the edge indices
            j,nodelist = self.get_edgeindices(self.xslice,self.yslice,method=self.edgemethod)
            self.j = j # Need this to calculate other quantities
            self.nslice = len(self.j)

            # Update the x and y axis of the slice
            self.xslice=self.xp[nodelist]
            self.yslice=self.yp[nodelist]
            

            self._getDistCoords()

            # Get the mask
            self.calc_mask()

            # Get the area and the normal
            area = self.calc_area()
            ne1, ne2, enormal = self.calc_normal(nodelist,j)
            dx = self.df[j]

            # Store all of the info as a dictionary
            self.slices.append({'j':j,'nodelist':nodelist,\
                'xslice':self.xslice,'yslice':self.yslice, \
                'distslice':self.distslice,'area':area,'normal':enormal,'dx':dx})

        # Find the unique j values and sort them
        j=[]
        for ss in self.slices:
            j = j + ss['j']

        j = np.array(j)
        self.j = np.sort(np.unique(j))

        # Find the index of each slice into this j array
        for ii,ss in enumerate(self.slices):
            ind = np.searchsorted(self.j,np.array(ss['j'])) 
            self.slices[ii].update({'subind':ind})

        self.j = self.j.tolist()

    def loadData(self,**kwargs):
        """
        Overloaded method of MultiEdgeSlice

        loads the data and then inserts it into each slice dictionary

        returns a list of 3D arrays
        """
        SliceEdge.loadData(self,**kwargs)
        data=[]
        for ii,ss in enumerate(self.slices):
            data.append(self.data[...,ss['subind']])
        return data

    def mean(self,phi,axis='area'):
        """
        Calculate the mean of the sliced data along an axis

        axis: depth, area
            depth: returns the time and spatial mean. Size = (Nslice,Nk)
            area: returns the area mean. Size = (Nslice,Nt)

        """
        Nslice = len(self.slices)
        Nt = len(self.tstep)
        Nk = self.Nkmax

        if axis=='area':
            data = np.zeros((Nslice,Nt))
            for ii,ss in enumerate(self.slices):
                area_norm = ss['area']/ ss['area'].sum()
                data[ii,:] = np.sum( np.sum(phi[ii]*area_norm,axis=-1),axis=-1)
        elif axis=='depth':
            data = np.zeros((Nslice,Nk))
            for ii,ss in enumerate(self.slices):
                phimean = np.mean(phi[ii],axis=0) # time mean
                dx = self.df[ss['j']]
                ne = len(ss['j'])
                dx = np.repeat(dx[np.newaxis,:],Nk,axis=0)
                dx[phimean==0.]=0
                dx_sum = dx.sum(axis=-1)
                dx_sum = np.repeat(dx_sum[:,np.newaxis],ne,axis=-1)
                dx_norm = dx / dx_sum
                data[ii,:] =  np.sum( phimean*dx_norm,axis=-1)
        else:   
            raise Exception, 'axis = %s not supported'%axis

        return data

           
        
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
