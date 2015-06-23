"""
Contour plot (z-t) of the heat flux model
"""

from sunpy import TimeSeries
import matplotlib.pyplot as plt
import numpy as np

ncfile = 'rundata/LakeMichigan_2012_0*.nc'
XY = np.array([5.0e5, 4.9e6])

# Load the suntans time series object
TS = TimeSeries(ncfile,XY)


# Plot temperature
plt.figure()

TS['variable']='temp'
ax1=plt.subplot(311)
TS.contourf()
plt.colorbar()
plt.ylabel('Depth [m]')

# Plot salinity
TS['variable']='salt'
ax2=plt.subplot(312)
TS.contourf()
plt.colorbar()
plt.ylabel('Depth [m]')

# Plot the heat flux value
ax3=plt.subplot(313)
TS['variable']='Hsw'
Hsw = TS.data.copy()
p0 = plt.plot(TS.time,Hsw,color='r')

TS['variable']='Hlw'
Hlw = TS.data.copy()
p1 = plt.plot(TS.time,Hlw,color='b')

TS['variable']='Hl'
Hl = TS.data.copy()
p2 = plt.plot(TS.time,Hl,color='m')

TS['variable']='Hs'
Hs = TS.data.copy()
p3 = plt.plot(TS.time,Hs,color='g')
plt.xticks(rotation=17)
plt.ylabel('Heat Flux $[W \ m^{-2}]$')

## Plot the heat flux value
#ax3=plt.subplot(313)
#TS['variable']='Hsw'
#Heat0 = np.zeros_like(TS.data)
#Heat = Heat0 + TS.data
#p0 = plt.fill_between(TS.time,Heat,y2=Heat0,color='r',alpha=0.7)
#
#TS['variable']='Hlw'
#Heat0 = Heat.copy()
#Heat = Heat0 + TS.data
#p1= plt.fill_between(TS.time,Heat,y2=Heat0,color='b',alpha=0.7)
#
#TS['variable']='Hl'
#Heat0 = Heat.copy()
#Heat = Heat0 + TS.data
#p2 = plt.fill_between(TS.time,Heat,y2=Heat0,color='m',alpha=0.7)
#
#TS['variable']='Hs'
#Heat0 = Heat.copy()
#Heat = Heat0 + TS.data
#p3 = plt.fill_between(TS.time,Heat,y2=Heat0,color='g',alpha=0.7)
#plt.xticks(rotation=17)
#plt.ylabel('Heat Flux $[W \ m^{-2}]$')


# legend
plt.legend((p0[0],p1[0],p2[0],p3[0]),('$H_{sw}$','$H_{lw}$','$H_L$','$H_S$'))

# Set the axis width
bbox2 = ax2.get_position().bounds
bbox3 = list(ax3.get_position().bounds)
bbox3[2]=bbox2[2]
ax3.set_position(bbox3)

ax1.set_xticklabels([])
ax2.set_xticklabels([])

plt.ylim([-1e3,1e3])

# Compare Tair to Hsw to test timing
plt.figure()
ax1 = plt.gca()
TS['variable']='Tair'
Tair = TS.data.copy()
ax2 = ax1.twinx()
ax1.plot(TS.time,Hsw,'r')
ax2.plot(TS.time,Tair,'b')

plt.show()
