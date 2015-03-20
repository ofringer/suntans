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
ncfile = 'data/Heatflux_AVG.0'
#cellindex=range(0,9)
cellindex=[4]
t0 = 0
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

depth = sun.dv[cellindex]
volume = area*depth # Cell volume
sumvolume = np.sum(volume)

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
    T_surf =    Qs[:,cellindex]*area # units [C m3 s-1]
    # Time-integrate
    #T_surf= T_surf[1:]
    #T_surf[0] = 0
    #T_surf = np.cumsum(T_surf*dt)
    #T_surf *= dt
    #T_surf = cumtrapz(T_surf,dx=dt)
else:
    T_surf=np.zeros_like(eta)

if sun.hasVar('EP'):
    EPs0 = sun.loadDataRaw(variable='EP')
    s_surf = EPs0[:,cellindex]*area  # units [psu m3 s-1]
else:
    s_surf = np.zeros_like(eta)

# Compute the total mass/tracer flux in/out of each cell
# (Sums by nsides the nlayers)
if len(Mflux.shape)==4:
    Mass = np.sum(np.sum(Mflux,axis=-1),axis=-2)
    Salt = np.sum(np.sum(Sflux,axis=-1),axis=-2)
    Temp = np.sum(np.sum(Tflux,axis=-1),axis=-2)
else:
    print 'Run is 2D!'
    Mass = np.sum(Mflux,axis=-1)
    Salt = np.sum(Sflux,axis=-1)
    Temp = np.sum(Tflux,axis=-1)


# Calculate the total salt and heat content

# Area average all total
#s_dz = np.sum(s_dz*area,axis=-1).squeeze()/sumvolume # units S
#T_dz = np.sum(T_dz*area,axis=-1).squeeze()/sumvolume # units T
#eta = np.sum(eta*area,axis=-1).squeeze()/sumarea
s_dz = np.sum(s_dz*area,axis=-1).squeeze() # m3 S
T_dz = np.sum(T_dz*area,axis=-1).squeeze() # m3 C
eta = np.sum(eta*area,axis=-1).squeeze() # m3 [volume]


# Sum the fluxes in and out of all cells
Mass = np.sum(Mass,axis=-1)
Mass = Mass.squeeze()
Salt = np.sum(Salt,axis=-1)
Salt = Salt.squeeze()
Temp = np.sum(Temp,axis=-1)
Temp = Temp.squeeze()
T_surf = np.sum(T_surf,axis=-1)
T_surf = T_surf.squeeze()
s_surf = np.sum(s_surf,axis=-1)
s_surf = s_surf.squeeze()

##########
# units are:
##########
# s_dz [ppt m3]
# T_dz [C m3]
# eta [m3]

# Mass [m3 s-1]
# Salt [ppt m3 s-1]
# Temp [C m3 s-1]

# Compute each of the terms in the budget

# Tendency
Tend_eta = (eta[:-1]-eta[1:]).squeeze()/dt # m3 s-1
Tend_s = (s_dz[:-1]-s_dz[1:]).squeeze()/dt # psu m3 s-1
Tend_T = (T_dz[:-1]-T_dz[1:]).squeeze()/dt # C m3 s-1

# Advective fluxes
Adv_mass = Mass[1:]# m3 s-1
Adv_s = Salt[1:]# psu s-1
Adv_T = Temp[1:]# C s-1

# Surface fluxes (note change of sign)
Sflux_T = -T_surf[1:]# C m3 s-1
Sflux_s = s_surf[1:]# psu m3 s-1

# Compute the error in each budget(relative)
error_eta =(Tend_eta - Adv_mass) / Tend_eta.mean() * 100
error_T = (Tend_T - Adv_T - Sflux_T) / Tend_T.mean() * 100
error_s = (Tend_s - Adv_s - Sflux_s) / Tend_s.mean() * 100

# Output time
time = time[1:]

# Free-surface
fig1=plt.figure()
f1ax1=fig1.add_subplot(2,1,1)
plt.title('Volume budget')
plt.plot(time,Tend_eta,'b',linewidth=2)
plt.plot(time,Adv_mass,'r')
plt.ylabel('$m^3 \ s^{-1}$')
plt.legend(('Tendency','Advection'))

ax2=fig1.add_subplot(2,1,2,sharex=f1ax1)
plt.plot(time,error_eta)
plt.ylabel('error')

fig2=plt.figure()
f2ax1=fig2.add_subplot(2,1,1)
plt.title('Temperature budget')
plt.plot(time,Tend_T,'b',linewidth=2)
plt.plot(time,Adv_T,'r')
plt.plot(time,Adv_T + Sflux_T,'k')
plt.ylabel(r'$^\circ C \ m^3 \ s^{-1}$')
plt.legend(('Tendency','Advection','Adv. + Sflux'))

f2ax1=fig2.add_subplot(2,1,2,sharex=f2ax1)
plt.title('Temperature budget')
plt.plot(time,error_T)
plt.ylabel('error')

fig3=plt.figure()
f3ax1=fig3.add_subplot(2,1,1)
plt.title('Salt budget')
plt.plot(time,Tend_s,'b',linewidth=2)
plt.plot(time,Adv_s,'r')
plt.plot(time,Adv_s + Sflux_s,'k')
plt.ylabel('$psu \ m^3 \ s^{-1}$')
plt.legend(('Tendency','Advection','Adv. + Sflux'))

f2ax1=fig3.add_subplot(2,1,2,sharex=f3ax1)
plt.plot(time,error_s)
plt.ylabel('error')

plt.figure()
sun.plotmesh()
plt.plot(sun.xv[cellindex],sun.yv[cellindex],'m.')
plt.show()

# Compute the error in the budget
#error_eta = (eta[:-1]-eta[1:]).squeeze() - Mass[1:]*dt/sumarea
#error_s = (s_dz[:-1]-s_dz[1:]).squeeze() - Salt[1:]*dt/sumarea
#error_T = (T_dz[:-1]-T_dz[1:]).squeeze() - (Temp[1:]+T_surf[1:])*dt/sumarea
#error_T2 = (T_dz[:-1]-T_dz[1:]).squeeze() - Temp[1:]*dt/sumarea

#print 'Mass error:'
#print error_eta
#print 'Salt error:'
#print error_s
#print 'Temp error:'
#print error_T



