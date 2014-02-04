"""
Calculates the heat and salt budget using average output
"""

from sunpy import Spatial
from maptools import maskPoly
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
from scipy.integrate import cumtrapz

import pdb

###
# Inputs
ncfile = 'data/Heatflux_AVG.nc.0'
xypoly = np.array(([500,9500,9500,500,500],[-10,-10,100,100,-10]))
#cellindex=range(0,9)
cellindex = [5]
t0 = 1
###

sun = Spatial(ncfile,klayer=[-99])
sun.tstep = range(t0,sun.Nt)
time = sun.time[sun.tstep]

# Constants
dt = sun.globalatts['dt']*sun.globalatts['ntaverage']
RHO0 = 1000.0
Cp = 4186.0

fac = (RHO0*Cp)

area = sun.Ac[cellindex]
sumarea = np.sum(area)

# Load the depth-average and flux quantities
s_dz = sun.loadDataRaw(variable='s_dz')
T_dz = sun.loadDataRaw(variable='T_dz')

Sflux = sun.loadDataRaw(variable='s_F')
Tflux = sun.loadDataRaw(variable='T_F')
Mflux = sun.loadDataRaw(variable='U_F') #[m3 s-1] * [s] = [m3]

# Load the total volume
eta = sun.loadDataRaw(variable='eta')


# Get the flux values at the cell
face = sun.face[cellindex,:]
normal = 1.0*sun.normal[cellindex,:]
eta = eta[:,cellindex]
s_dz = s_dz[:,cellindex]
T_dz = T_dz[:,cellindex]

Mflux = Mflux[:,:,face]
Mflux*=normal
Sflux = Sflux[:,:,face]
Sflux*=normal
Tflux = Tflux[:,:,face]
Tflux*=normal

if sun.hasVar('Hs'):
    # Load the surface flux quantities
    Hs = sun.loadDataRaw(variable='Hs')
    Hsw = sun.loadDataRaw(variable='Hsw')
    Hl = sun.loadDataRaw(variable='Hl')
    Hlw = sun.loadDataRaw(variable='Hlw')

    Qs = (Hs+Hl+Hlw+Hsw)/fac

    # Surface flux contribution
    T_surf =    (Qs[:,cellindex]*sun.Ac[cellindex])
    # Time-integrate
    #T_surf= T_surf[1:]
    #T_surf[0] = 0
    #T_surf = np.cumsum(T_surf*dt)
    #T_surf *= dt
    #T_surf = cumtrapz(T_surf,dx=dt)
else:
    T_surf=np.zeros_like(eta)

# Compute the total mass/tracer flux in/out of the cell
if len(Mflux.shape)==4:
    Mass = np.sum(np.sum(Mflux,axis=-1),axis=-2)
    Salt = np.sum(np.sum(Sflux,axis=-1),axis=-2)
    Temp = np.sum(np.sum(Tflux,axis=-1),axis=-2)
else:
    print 'Run is 2D!'
    Mass = np.sum(Mflux,axis=-1)
    Salt = np.sum(Sflux,axis=-1)
    Temp = np.sum(Tflux,axis=-1)


# Convert the volume to a height
Height = Mass/area
Height = np.mean(Height,axis=-1).squeeze()
# Calculate the total salt and heat content

# Area integrate all terms
Mass = np.sum(Mass,axis=-1)
Mass = Mass.squeeze()
Salt = np.sum(Salt,axis=-1)
Salt = Salt.squeeze()
Temp = np.sum(Temp,axis=-1)
Temp = Temp.squeeze()
T_surf = np.sum(T_surf,axis=-1)
T_surf = T_surf.squeeze()

eta = np.sum(eta*area,axis=-1)
eta = eta.squeeze()
s_dz = np.sum(s_dz*area,axis=-1)
s_dz = s_dz.squeeze()
T_dz = np.sum(T_dz*area,axis=-1)
T_dz = T_dz.squeeze()

# Calculate the LHS and RHS terms of the mass and scalar budgets
LHS_mass = eta[:]
LHS_salt = s_dz[:]
LHS_temp = T_dz[:]

Height[0]=0
Salt[0]=0
Temp[0]=0
T_surf[0]=0
RHS_mass =  eta[0] - np.cumsum(Height[:]*dt) 
RHS_salt =  s_dz[0] - np.cumsum(Salt[:]*dt) 
RHS_temp =  T_dz[0] - np.cumsum(Temp[:]*dt) + np.cumsum(T_surf[:]*dt)

timeflux = time[:]
timelhs = time[:]

# Error plot
#plt.figure(figsize=(6,12))
#plt.subplot(311)
#plt.plot(timeflux,LHS_mass-RHS_mass,linewidth=2)
#plt.title('Mass budget')
#
#plt.subplot(312)
#plt.plot(timeflux,LHS_salt-RHS_salt,linewidth=2)
#plt.title('Salt budget')
#
#plt.subplot(313)
#plt.plot(timeflux,LHS_temp-RHS_temp,linewidth=2)
#plt.show()


##
plt.figure(figsize=(6,12))
plt.subplot(311)
plt.plot(timelhs,LHS_mass,linewidth=2)
plt.plot(timeflux,RHS_mass,'k')
plt.title('Mass budget')

plt.subplot(312)
plt.plot(timelhs,LHS_salt,linewidth=2)
plt.plot(timeflux,RHS_salt,'k')
plt.title('Salt budget')
plt.legend(('$\int\ \phi\ dz$','$\sum Fluxes$'))

plt.subplot(313)
plt.plot(timelhs,LHS_temp,linewidth=2)
plt.plot(timeflux,RHS_temp,'k',linewidth=1.5)
plt.title('Heat budget')
plt.show()

