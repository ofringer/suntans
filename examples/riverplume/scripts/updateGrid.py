"""
Updates suntans grid text files to the quad format

Edges.dat has an extra column for the 'edge_id' variable
"""

from sunpy import Grid
import os

#
grdfolder = 'grids/quad/'
outfolder = 'rundata/'
maxfaces=4

grd = Grid(grdfolder,MAXFACES=maxfaces)

try:
    os.mkdir(outfolder)
except:
    pass
os.system('cp %s/*.dat %s'%(grdfolder,outfolder))

grd.saveEdges(outfolder+'/edges.dat')
grd.saveCells(outfolder+'/cells.dat')

print 'Saved new grid in: %s'%outfolder

