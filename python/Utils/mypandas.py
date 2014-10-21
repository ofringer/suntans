"""
Extensions to the wonderful Pandas module
"""

import pandas as pd
from scipy import signal
import numpy as np

def godin(phi):
    """
    Godin type filter of a time series
    """
    dt = get_dt(phi)
    window24 = 24.*3600//dt
    window25 = 25.*3600//dt
    phi_f = pd.rolling_mean(phi,window24,center=True)
    phi_f = pd.rolling_mean(phi_f,window25,center=True)
    phi_f = pd.rolling_mean(phi_f,window24,center=True)
    return phi_f

def butterfilt(phi,cutoff_dt,btype='low',order=3):
    """
    Butter worth filter the time series
    """
    dt = get_dt(phi)

    if not btype == 'band':
        Wn = dt/cutoff_dt
    else: # Band-pass expects a list of cuttoff frequencies
        Wn = [dt/co for co in cutoff_dt]
        
    (b, a) = signal.butter(order, Wn, btype=btype, analog=0, output='ba')
    
    # filtfilt only likes to operate along the last axis
    ytmp = phi.values.swapaxes(-1,0)
    ytmp = signal.filtfilt(b, a, ytmp, axis=-1)

    # Return a pandas object like the input
    phi_filt = phi.copy()
    phi_filt[:] = ytmp.swapaxes(0,-1)

    return phi_filt
 
def rms(phi,cutoff_dt):
    """
    rolling rms
    """
    dt = get_dt(phi)
    window = cutoff_dt//dt
    phi_rms = pd.rolling_mean(phi**2,window,center=True)

    return np.sqrt(phi_rms)

def get_dt(phi):
    """Calculates the time step of a pandas object"""
    if type(phi) == pd.Panel:
        time = phi.items
    else:
        time = phi.index

    jdate = time.to_julian_date()*86400.
    return jdate[1]-jdate[0]

