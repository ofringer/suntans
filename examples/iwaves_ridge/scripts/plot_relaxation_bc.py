"""
Plot the exponential function
"""

import numpy as np
import matplotlib.pyplot as plt

def sponge(x, R):
    return np.exp(-8 * (x/R)**1.0)


sponge_dist = 4e4
dist = np.linspace(0,2*sponge_dist, 200)

plt.plot(dist, sponge(dist, sponge_dist) )
plt.show()
