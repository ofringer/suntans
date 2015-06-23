"""
Try plotting particles as streaks using matplotlib
"""
from suntrack import PtmNC
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation
from otherplot import streakplot
from maptools import plotmap

import pdb

###########
# Inputs
ncfile = 'data/Headland_particles.nc'
outfile = 'plots/Headland_particle_streak.mov'

taillength = 18 # Number of time steps for tail
subset = 60 # Only plot every nth particle

# Plot specific stuff
#shpfile= '../../DATA/GalvestonBasemapCoarse.shplocations
W = 2.5e4
a =2e3
b =8e3
x0=0

xlims = [x0-W/2, x0+W/2]
ylims = [0, 2*b]
###########


# Itialize the particle class##
sun = PtmNC(ncfile)

tstart = range(0,sun.nt-taillength)
ts = range(tstart[0],tstart[0]+taillength)

# Read the particle locations at the initial step
xp = sun.read_step(ts,'xp')[::subset,:]
yp = sun.read_step(ts,'yp')[::subset,:]
zp = sun.read_step(ts,'zp')[::subset,:]


# Plot specific stuff
fig,ax = plt.subplots()
ax.set_axis_bgcolor('#001933')
ax.set_xticklabels([])
ax.set_yticklabels([])

# This plots a map
#plotmap(shpfile)

# Initialize the streak plot
S=streakplot(xp,yp,ax=ax,xlim=xlims,ylim=ylims)

####
# Animation code
####

def updateposition(i):
    print i
    ts = range(tstart[i],tstart[i]+taillength)
    xp = sun.read_step(ts,'xp')[::subset,:]
    yp = sun.read_step(ts,'yp')[::subset,:]

    S.update(xp,yp)

    return S.lcs

anim  = animation.FuncAnimation(fig,updateposition,frames=len(tstart),\
    interval=50,blit=True)

print 'Building animation...'
anim.save(outfile,fps=6,bitrate=3600)
print 'Done'

#plt.show()

