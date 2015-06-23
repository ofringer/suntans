"""
Calculates the dual of a grid
"""

from sunpy import Grid
import matplotlib.pyplot as plt
import numpy as np
from hybridgrid import HybridGrid


        
###

ncfile='grid'

outpath='grid/hex'

sun = Grid(ncfile)
sun.neigh = sun.neigh.astype(int)

# Load the data into a hybrid grid class
grd = HybridGrid(sun.xp,sun.yp,sun.cells,nfaces=sun.nfaces)

grd.create_dual_grid(minfaces=5,outpath=outpath)


sun = Grid(outpath)

plt.figure()
sun.plotmesh(facecolors='none',linewidths=0.2)
plt.plot(sun.xv,sun.yv,'m.')


plt.figure()
sun.plothist()
plt.show()

