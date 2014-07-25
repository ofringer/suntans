"""
Other plotting routines outside of matplotlib
"""
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt

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
