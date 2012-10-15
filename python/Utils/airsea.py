# -*- coding: utf-8 -*-
"""
Collection of tools for calculating various air-sea quantities

Created on Fri Jul 27 13:56:06 2012
@author: mrayson
"""

import numpy as np

def satVapPres(Ta,P=1010):
    """ Calculates the saturation vapor pressure of air at temperature, T

    Inputs: T - temperature [C]
            P - air pressure [mb]
            
    """
    ew=np.power(10,((0.7859+0.03477*Ta)/(1+0.00412*Ta)))
    
    fw=1 + 1e-6*P*(4.5+0.0006*Ta**2)
    
    ew=fw*ew;
    
    return ew
    
def relHumFromTdew(T,Tdew,P=1010):
    """ Calculates the relative humidity (%) from Dew point temperature"""
    
    e_dew = satVapPres(Tdew,P)
    e_dry = satVapPres(T,P)
    
    rh = e_dew / e_dry *100 
    
    return rh
    
def convertSpeedDirn(theta,rho):
    
    """
    (modifed from MATLAB compass2cart function)
    %COMPASS2CART convert speed and direction data (degN) into
    % cartesian coordinates.
    %   COMPASS2CART(THETA,RHO) convert the vector rho (e.g. speed) with
    %      direction theta (degree North) into cartesian coordinates u and v.
    %      note: theta is in degrees and between 0 and 360.
    """
    
    if theta >= 0 and theta <90:
        theta=np.abs(theta-90)
    elif theta >= 90 and theta <= 360:
        theta=np.abs(450-theta)
    
    u,v = pol2cart(theta*np.pi/180,rho)
    
    return u, v
    
def pol2cart(th,rho):
    """Convert polar coordinates to cartesian"""
    x = rho * np.cos(th)
    y = rho * np.sin(th)

    return x, y