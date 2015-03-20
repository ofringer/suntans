"""
Generate a constriction grid

Script that generates simple inputs from TOM and generates a triangle 
grid
"""

import numpy as np
from shapely import geometry
from maptools import shapely2shp
import os
from tom import Tom

import matplotlib.pyplot as plt

###
# Inputs
# Inputs
L = 5e4
W = 2.5e4
a =2e3
b =8e3
x0 = 0

# target grid scale
scale = 1000.#
scale_inner = 100.
outpath = 'grids/tri/'
#
##########
if not os.path.isdir(outpath):
    print 'Creating new directory: %s'%outpath
    os.mkdir(outpath)

###
# Create the boundary outline
###
###
# Create the boundary outline
###
X = np.arange(x0-L/2,x0+L/2+scale,scale_inner)


def gaussian(x,x0,a,b,pow=2.):
    return b*np.exp(-0.5*((x-x0)/a)**pow)

y_bot = gaussian(X,x0,a,b)
y_top = W*np.ones_like(y_bot)

y = np.hstack((y_bot,y_top[::-1],y_bot[0]))
x = np.hstack((X,X[::-1],X[0]))

xy =[(x[ii],y[ii]) for ii in range(x.shape[0])]

###
# Generate the tom input shapefiles
###
# Create the boundary shapefile
poly = geometry.Polygon(xy)
bndfile = '%s/boundary.shp'%outpath
shapely2shp(poly,bndfile)

# Create the scale shape file
line = geometry.LineString(xy)
lineint = geometry.LineString([[x0-a/2.,b],[x0+a/2.,b]])
scalefile = '%s/scale.shp'%outpath
atts = {'scale':[scale]}
shapely2shp([line],scalefile,atts=atts)

#atts = {'scale':[scale,scale_inner]}
#shapely2shp([line,lineint],scalefile,atts=atts)

# Create an interiors and telescoping shapefile
lineint = geometry.LineString([[x0-a/2.,b],[x0+a/2.,b]])
atts = {'scale':[scale_inner]}
scalefile = '%s/interiors.shp'%outpath
shapely2shp([lineint],scalefile,atts=atts)


####
# Run tom
####
os.chdir(outpath)
tom = Tom()
#tom.run(['interactive','-s','scale.shp','-b','boundary.shp',])
tom.run(['interactive','-s','scale.shp','-b','boundary.shp',\
    #'-i','interiors.shp',\
    '-a','interiors.shp',\
    '-t','1.03'])
    #])

#plt.plot(X,y_bot)
#plt.plot(X,y_top)
#plt.plot(x,y,'k')
#plt.show()
