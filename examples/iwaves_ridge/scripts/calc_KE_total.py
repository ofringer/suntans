"""
Compute the total volume integrated KE for a domain
"""

import numpy as np
import matplotlib.pyplot as plt
import xray

from soda.dataio.suntans.sunpy import Spatial


def plotfile(ncfile, color):
    sun = Spatial(ncfile)

    # Calculate the volume
    V = sun.dz[:,np.newaxis]*sun.Ac[np.newaxis,:]
    Vall = V.sum()

    # Load the data as a xray object
    ds = xray.open_dataset(ncfile)

    KE = 0.5*(ds['uc']**2 + ds['vc']**2)

    int_KE = KE*V[np.newaxis,...]
    int_KE = int_KE.values
    int_KE[np.isnan(int_KE)] = 0

    KEt = int_KE[:,20:80,:].sum(axis=1).sum(axis=1)

    plt.plot(KEt, color)


plotfile('rundata/IWave_sponge_0000.nc', 'r')
plotfile('rundata/IWave_nosponge_0000.nc', 'b')
plt.show()
