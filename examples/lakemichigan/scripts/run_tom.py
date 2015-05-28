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

outpath = 'grid'
scalefile = '../gis/scale.shp'
bdyfile = '../gis/Lake_Michigan_wIslands_UTM.shp'
interiorfile = '../gis/interiors.shp'
####
# Run tom
####
os.chdir(outpath)
tom = Tom()
#tom.run(['interactive','-s','scale.shp','-b','boundary.shp',])
tom.run(['interactive','-s',scalefile,'-b',bdyfile,\
    #'-i','interiors.shp',\
    #'-a',interiorfile,\
    #'-t','1.03'])
    ])

