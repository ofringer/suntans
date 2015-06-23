"""
Plots a hybrid grid 
"""

from sunpy import Grid
import matplotlib.pyplot as plt

###
# 
#ncfile = 'rundata/Plume.nc.0'
ncfile = 'grids/hybrid'

sun = Grid(ncfile)

plt.figure()
#sun.plot()
sun.plotmesh(facecolors='none')
plt.show()
