"""
Generates a series of animations
"""

from sunpy import Spatial
import numpy as np
import matplotlib.pyplot as plt

runname='SFBay3D'
ncfile = '%s/%s_0*.nc'%('data',runname)
k = 1 # depth layer

plt.figure()
sun = Spatial(ncfile,klayer=[k],variable='salt')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.clim = [28.0,32.0]
sun.animate(vector_overlay=False)
sun.saveanim('plots/%s_salt.mov'%runname)

sun.clim=None
plt.figure()
sun = Spatial(ncfile,klayer=[k],variable='uc')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/%s_uc.mov'%runname)

sun.clim=None
plt.figure()
sun = Spatial(ncfile,klayer=[k],variable='eta')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/%s_eta.mov'%runname)

