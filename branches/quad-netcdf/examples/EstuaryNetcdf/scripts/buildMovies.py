"""
Generates a series of animations from the suntans output

"""

from sunpy import Spatial
import numpy as np
import matplotlib.pyplot as plt

ncfile = 'rundata/Estuary_SUNTANS_00*nc'

plt.figure()
sun = Spatial(ncfile,klayer=[1],variable='salt')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/salt.mov')

plt.figure()
sun = Spatial(ncfile,klayer=[1],variable='temp')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/temp.mov')

plt.figure()
sun = Spatial(ncfile,variable='tau_y')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/tau_y.mov')

plt.figure()
sun = Spatial(ncfile,klayer=[1],variable='vc')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate(vector_overlay=False)
sun.saveanim('plots/vc.mov')

plt.figure()
sun = Spatial(ncfile,variable='eta')
sun.tstep=np.arange(0,len(sun.time))
sun.loadData()
sun.animate()
sun.saveanim('plots/eta.mov')

