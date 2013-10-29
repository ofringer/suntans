"""
Example of how to plot the grid and open boundaries in python
"""

from sunpy import Grid
import matplotlib.pyplot as plt
from maptools import Polygon2GIS

grdfolder = 'rundata/'
kmlfile = 'gis/sfbay_suntans_grid.kml'
shpfile = 'gis/sfbay_suntans_grid.shp'
utmzone = 10

# Load the grid into a sunpy object
sun = Grid(grdfolder)

# save the grid to a google earth file
#Polygon2GIS(sun.xy,kmlfile,utmzone) # the grid polygons are stored in the 'xy' attribute

# save the grid to a gis shapefile
#Polygon2GIS(sun.xy,shpfile,utmzone) 

# Plot the mesh
plt.figure()
sun.plotmesh()

# Plot the boundary conditions
sun.plotBC()
plt.show()
