"""
Generates a series of animations
"""

from sunpy import Spatial
import numpy as np
import matplotlib.pyplot as plt

runname='Plume'
ncfile = '%s/%s.nc.0'%('rundata',runname)

plt.figure()
sun = Spatial(ncfile,klayer=[0],variable='salt')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.clim = [28.0,32.0]
sun.animate(vector_overlay=False)
sun.saveanim('plots/%s_salt.mov'%runname)

sun.clim=None
plt.figure()
sun = Spatial(ncfile,klayer=[0],variable='uc')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/%s_uc.mov'%runname)

sun.clim=None
plt.figure()
sun = Spatial(ncfile,klayer=[0],variable='eta')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/%s_eta.mov'%runname)

