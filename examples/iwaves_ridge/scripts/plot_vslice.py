"""
Plot a vertical slice of the data
"""

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from soda.dataio.suntans.sunslice import SliceEdge

import matplotlib
# Set some default parameters
matplotlib.rcParams['text.color']='white'
matplotlib.rcParams['savefig.facecolor']='black'
matplotlib.rcParams['savefig.edgecolor']='black'
matplotlib.rcParams['figure.facecolor']='black'
matplotlib.rcParams['figure.edgecolor']='black'
matplotlib.rcParams['axes.facecolor']='black'
matplotlib.rcParams['axes.edgecolor']='white'
matplotlib.rcParams['axes.labelcolor']='white'
matplotlib.rcParams['xtick.color']='white'
matplotlib.rcParams['ytick.color']='white'
matplotlib.rcParams['font.family']='serif'

def main(ncfile):
    #
    ###
    # Inputs

    #ncfile = 'rundata/IWave_nosponge_0000.nc'

    xpt = np.array([100., 2.8e5])
    ypt = np.array([100., 100.])

    #varname = 'uc'
    #clim = [-0.08, 0.08]

    varname = 'w'
    clim = [-0.001, 0.001]

    #varname = 'salt'
    #clim = None

    #varname = 'temp'
    #clim = None

    cmap = 'RdBu'

    t0 = 0
    ######

    # Load the slice object
    print ncfile
    sun = SliceEdge(ncfile, xpt=xpt, ypt=ypt)
    sun.tstep = [t0]

    if clim is not None:
        sun.clim = clim

    # This removes masked values
    cmap = getattr(matplotlib.cm, cmap)
    cmap.set_bad('k')

    #print sun.j

    # Load the data
    data = sun.loadData(variable = varname, method='max')

    # Build the figure
    fig = plt.figure(figsize = (12,8), num = 'SUNTANS Vertical Slice Tool')

    # Time slider axes
    axtime = plt.subplot2grid((7,3), (6,0), colspan=2, rowspan=1)

    ax = plt.subplot2grid((7,3), (0,0), colspan=3, rowspan=5)
    #sun.contourslice(data,t=0)
    h1, axcb, title = sun.pcolorslice(data, cmap = cmap, bathyoverlay=True, interpolation='none' )

    # Create the time slider
    valstr = ' of %d'%(sun.Nt-1)
    ts = Slider(axtime, 'Time', 0, sun.Nt-1, valinit=t0, valfmt='%d'+valstr,facecolor='0.5',)

    # On change of slider: load new data and set the plot objects accordingly
    def update_slider(val):
        t = int(np.floor(val))
        #print t
        if not sun.tstep[0] == t:
            sun.tstep = [t]
            sun.loadData(variable = varname, method='mean')
            h1.set_array(sun.data)
            title.set_text('%s [%s]\n%s'%(sun.long_name, sun.units,\
                datetime.strftime(sun.time[t], '%Y-%m-%d %H:%M:%S')))
            fig.canvas.draw_idle()

    ts.on_changed(update_slider)

    plt.show()

if __name__=='__main__':
    import sys
    ncfile = sys.argv[1]
    print ncfile
    main(ncfile)
