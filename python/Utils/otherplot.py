"""
Other plotting routines outside of matplotlib
"""
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


class StreakPlot(object):
    """
    Class for generating a streak plot
    """

    # Width of the start and end of the tail
    widthmin = 0.1
    widthmax = 1.5
    # Plotting parameters
    colors = '0.9'
    alpha = 0.9

    xlim=None
    ylim=None

    def __init__(self,xp,yp,**kwargs):
        """
        StreakPlot:

        xp, yp - spatial location matrices with dimension Nt x Nxy

        Streaks are along the first dimension Nt
        """

        self.__dict__.update(**kwargs)

        self.xp = xp
        self.yp = yp
        self.nx,self.nt = np.shape(self.xp)

        self.create_xy()


    def create_xy(self):
        # Create the start and end points
        self.points = np.concatenate([self.xp[...,np.newaxis],\
            self.yp[...,np.newaxis]],axis=-1)

        # Create the linewidths
        lwidths = np.linspace(self.widthmin,self.widthmax,self.nt-1)
        lwidths = lwidths**2 # this thins out the tail a bit
        self.lwidths = lwidths/lwidths[-1] * self.widthmax

        # Create a list of segment arrays
        self.segments=[ np.concatenate([self.points[ii,:-1,np.newaxis,:],\
            self.points[ii,1:,np.newaxis,:]],axis=1)\
            for ii in range(self.nx) ]

        # Create a list of line collections
        self.lcs = [ LineCollection(segment, linewidths=self.lwidths,\
            colors=self.colors,alpha=self.alpha)\
            for segment in self.segments ]

    def update(self,xp,yp):
        """
        Update the line collections with new x,y points
        """
        self.xp=xp
        self.yp=yp
        # Create the start and end points
        self.points = np.concatenate([xp[...,np.newaxis],yp[...,np.newaxis]],axis=-1)

        # Create a list of segment arrays
        self.segments=[ np.concatenate([self.points[ii,:-1,np.newaxis,:],\
            self.points[ii,1:,np.newaxis,:]],axis=1)\
            for ii in range(self.nx) ]

        for lc, seg in zip(self.lcs,self.segments):
            lc.set_segments(seg)

    def plot(self, ax):
        """Inserts each line collection into current plot"""

        if self.xlim == None:
            self.xlim = [self.xp.min(),self.xp.max()]
        if self.ylim == None:
            self.ylim = [self.yp.min(),self.yp.max()]

        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        ax.set_aspect('equal')

        map(ax.add_collection,self.lcs)
        #for lc in self.lcs:
        #    ax.add_collection(lc)

        return ax

def streakplot(xp,yp,ax=None,**kwargs):
    """
    Functional call to StreakPlot class
    """
    if ax==None:
        ax=plt.gca()

    S = StreakPlot(xp,yp,**kwargs)
    S.plot(ax)
    return S


def axcolorbar(cbobj,pos=[0.7, 0.8, 0.2, 0.04],ax=None,fig=None,orientation='horizontal',**kwargs):
	"""
	Inserts a colorbar with a position relative to an axes and not a figure
	
	Inputs:
		cbobj - plot object for colorbar
		pos - position vector [x0, y0, width, height] in dimensionless coordinates
		ax - axes to insert colobar
		figure - figure 
		**kwargs - arguments for plt.colorbar
	
	Returns a colorbar object
	
	Derived from this post:
		http://stackoverflow.com/questions/22413211/cant-fix-position-of-colorbar-in-image-with-multiple-subplots
	"""
	if fig == None:
		fig=plt.gcf()
	if ax == None:
		ax=plt.gca()
		
	fig.tight_layout()  # You call fig.tight_layout BEFORE creating the colorbar

	# You input the POSITION AND DIMENSIONS RELATIVE TO THE AXES
	x0, y0, width, height = pos

	# and transform them after to get the ABSOLUTE POSITION AND DIMENSIONS
	Bbox = transforms.Bbox.from_bounds(x0, y0, width, height)
	trans = ax.transAxes + fig.transFigure.inverted()
	l, b, w, h = transforms.TransformedBbox(Bbox, trans).bounds

	# Now just create the axes and the colorbar
	cbaxes = fig.add_axes([l, b, w, h])
	cbar = plt.colorbar(cbobj, cax=cbaxes,orientation=orientation, **kwargs)
	cbar.ax.tick_params(labelsize=9)

	return cbar
