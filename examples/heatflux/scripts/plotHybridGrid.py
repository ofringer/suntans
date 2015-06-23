"""
Plots a hybrid grid 
"""

from sunpy import Grid
import matplotlib.pyplot as plt
import numpy as np
from hybridgrid import HybridGrid
from maptools import Polygon2GIS
###
# 
ncfile = 'rundata'
###


sun = Grid(ncfile)
sun.neigh = sun.neigh.astype(int)

# Load the data into a hybrid grid class
grd = HybridGrid(sun.xp,sun.yp,sun.cells,nfaces=sun.nfaces,
                 edges=sun.edges,mark=sun.mark,grad=sun.grad,neigh=sun.neigh)


# Update the suntans grid with the hybrid grid features
#sun.mark = grd.mark
#sun.grad = grd.grad
#sun.neigh = grd.neigh
#sun.saveEdges(ncfile+'/edges.dat')

### Test the rotation
##cell_ccw = grd.ensure_ccw()
#
#Ac = grd.calc_area()
#
#ind = 1
#plt.figure()
#for i in range(grd.nfaces[ind]):
#    plt.plot(grd.xp[grd.cells[ind,i]],grd.yp[grd.cells[ind,i]],'ro')
#    plt.text(grd.xp[grd.cells[ind,i]],grd.yp[grd.cells[ind,i]],'%d'%i)   
#    
##    plt.plot(grd.xp[cell_ccw[ind,i]],grd.yp[cell_ccw[ind,i]],'bo')
##    plt.text(grd.xp[cell_ccw[ind,i]],grd.yp[cell_ccw[ind,i]]-2,'%d'%i,color='b')   
#plt.show()

## plot a cell and it's neighbors
ind = 4
# Find and plot the edge locations
nodes = sun.cells[ind,:]
nfaces=grd.nfaces[ind]
e1 = grd.find_edge([nodes[0],nodes[1]])
e2 = grd.find_edge([nodes[1],nodes[2]])

if nfaces==4:
    e3 = grd.find_edge([nodes[2],nodes[3]])
    e4 = grd.find_edge([nodes[3],nodes[0]])
else:
    e3 = grd.find_edge([nodes[2],nodes[0]])

grd.edge_centers()

plt.figure()
#sun.plot()
sun.plotmesh(facecolors='none',linewidths=0.2)

#plt.plot(grd._xva,grd._yva,'b.')
#plt.plot(grd._xvb,grd._yvb,'b.')
plt.plot(grd.xv,grd.yv,'m.')
#plt.plot(grd.xc,grd.yc,'r.')
#plt.plot(grd.xca,grd.yca,'ro')

# Plot the neighbours of a selected cell
#plt.plot(grd.xv[ind],grd.yv[ind],'bo')
#plt.plot(grd.xv[sun.neigh[ind,0:nfaces]],grd.yv[sun.neigh[ind,0:nfaces]],'ro')
#plt.plot(grd.xp[sun.cells[ind,0:nfaces]],grd.yp[sun.cells[ind,0:nfaces]],'mo')
#plt.plot(grd.xe[e1],grd.ye[e1],'gs',markersize=10)
#plt.plot(grd.xe[e2],grd.ye[e2],'gs',markersize=10)
#plt.plot(grd.xe[e3],grd.ye[e3],'gs',markersize=10)
#if nfaces==4:
#    plt.plot(grd.xe[e4],grd.ye[e4],'gs',markersize=10)

plt.figure()
sun.plothist()
plt.show()

# Write to a shapefile
#Polygon2GIS(sun.xy,shpfile,15)
