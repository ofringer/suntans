
from sunboundary import modifyBCmarker
import trigrid
import numpy as np
import matplotlib.pyplot as plt

inpath = 'C:/Projects/GOMGalveston/MODELLING/GalvestonCoarseFinal/rundata/'
outpath = 'C:/Projects/GOMGalveston/MODELLING/GalvestonCoarseFinal/grid2'

g = trigrid.TriGrid(suntans_path=inpath)
for j in np.nonzero(g.edges[:,2]==3)[0]:
    g.delete_edge(j)
g.renumber()
g.plot() # matplotlib plot of the edges to see if it did the right thing
g.write_suntans(outpath)
