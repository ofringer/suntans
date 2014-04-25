# -*- coding: utf-8 -*-
"""
Collection of tools for calculating various air-sea quantities

Main Reference:
    Kantha and Clayson, 2000, "Small Scale Processes in Geophysical Fluid Flows",
    Academic Press
    
Created on Fri Jul 27 13:56:06 2012
Author: Matt Rayson
Stanford University

"""

import numpy as np
import operator
import pdb

def buoyancyFlux(SSS,SST,Q,EP,dz):
    """
    Calculate the surface buoyancy flux
    
    Inputs:
    	SSS - sea surface salinity (psu)
	SST - sea surface temperater (celsius)
	Q - net heat flux (W m-2)
	EP - evaporation minus precipitation (m s-1)
	dz - grid cell height [m]
    Returns:
    	B_heat - heat buoyancy flux [W m-2]
	B_salt - salt buoyancy flux [W m-2]
    Ref:
    	Gill, 1982
    """
    try:
    	import seawater # CSIRO seawater toolbox
    except:
    	raise Exception, ' need to install CSIRO seawater toolbox'

    # Constants
    Cpinv = 1./4200.0 # Specific heat capacity [J kg-1 K-1]
    g = 9.81
    RHO0 = 1000.0

    # Convert EP from [m s-1] -> [kg m-2 s-1]
    EP = EP*RHO0

    # Calculate thermal expansion and saline contraction coefficient
    alpha = seawater.alpha(SSS,SST,0*SST)
    beta = seawater.beta(SSS,SST,0*SST)

    # Buoyancy flux equation (see Gill eq 2.7.1, p36)
    #  Note that units are [kg m s-3] or [W m-1]
    B_heat = Cpinv * g * alpha * Q
    B_salt = g*beta*EP*SSS

    # returns the fluxes in units [W m-2]
    return B_heat/dz, B_salt/dz

def heatFluxes(Uwind,Vwind,Ta,Tw,Pa,RH,cloud):
    """ 
    Calculate the non-penetrative radiative and turbulent heat flux terms:
        H_lwu - upward longwave radiation
        H_lwd - downward longwave radiation
        H_l - latent heat flux
        H_s - sensible heat flux
        
    **Sign Convention**
        -ve -> out of water
        +ve -> into water
        
    Inputs:
        Uwind - eastward wind velocity [m/s]
        Vwind - northward wind velocity [m/s]
        Ta - air temp [C]
        Tw - water temp [C]
        Pa - air pressure [mb]
        RH - relative humidity [%]
        cloud - cloud cover fraction [0-1]
        
    """
    
    # Coefficients
    Ce = 1.10e-3 # Dalton number
    Ch = 1.10e-3 # Stanton number
    
    Le = 2.5e6
    cp = 4.2e3
    rhoa = 1.20
    
    r_LW = 0.04
    
    # Compute wind speed
    S = np.sqrt(Uwind**2+Vwind**2)
    
    # Latent heat flux
    q = qspec(Ta,RH,Pa)
    qs = 0.98*qsat(Tw,Pa) # Need to scale seawater by 0.98
    Hl = latentBulk(qs,q,S,rhoa,Ce,Le) 
    dq = qs-q
    
    # Sesible heat flux
    Hs = sensibleBulk(Tw,Ta,S,rhoa,Ch,cp)
    dT = Tw - Ta
    
    # longwave
    Hlwu = longwaveUp(Tw)
    
    Hlwd = longwaveDown(Ta,cloud,r_LW)
    
    
    return Hl, Hs, Hlwu, Hlwd, dq, dT, S
    
    
def latentBulk(qs, q, S, rhoa=1.2, Ce=1.5e-3, Le=2.5e6):
    """
    Latent heat flux from water using the bulk exchange formulation
    
    Inputs:
        qs - Saturation specific humidity [kg kg-1]
        q - Air specific humidity [kg kg-1]
        S - Wind speed magnitude [m s^-1]
        rhoa - air density [kg m^-3]
        Ce - Dalton number
        Le - Latent heat of evaporation [J kg^-1]
        
    To calculate given Ta, Tw, RH, Pa, U, V:
    ---------------------------------------- 
    q = qspec(Ta,RH,Pa)
    qs = 0.98*qsat(Tw,Pa) # Need to scale seawater by 0.98
    S = sqrt(U**2+V**2)
    Hl = latentBulk(qs,q,S)     
    
    """
    
    return -rhoa*Le*Ce*S*(qs-q)
    
def sensibleBulk(Tw,Ta,S,rhoa=1.2,Ch=1.5e-3,cpa=1004.67):
    """
    Sensible heat flux from water using the bulk exchange formulation
    
    Inputs:
        Tw - Water temp [C]
        Ta - Air temp [C]
        S - Wind speed magnitude [m s^-1]
        rhoa - air density [kg m^-3]
        Ch - Stanton number
        cpa - Specific heat of air [J kg-1 K-1]
    """    
    return -rhoa*cpa*Ch*S*(Tw-Ta)    

def stressBulk(u,S,rhoa=1.2,Cd=1.1e-3):
    """
    Calculate the wind stress component using the bulk exchange formulation
    
    Inputs:
        u - wind velocity x/y component [m s-1]
        S - wind speed magnitude [m s-1]
        Cd - drag coefficient
        rhoa - air density [kg m-3]
    """
    return rhoa*Cd*S*u
    
def qsat(T,Pa):
    """
    Compute the specific humidity at saturation for a given T
    Inputs:
        T - temperature [degC]
        P - air pressure [mb or hPa]
    """
    ew = 6.1121*(1.0007+3.46e-6*Pa)*np.exp((17.502*T)/(240.97+T)) # in mb
    return 0.62197*(ew/(Pa-0.378*ew))                    # mb -> kg/kg 
    
def qspec(T,RH,Pa):
    """
    Calculate the specific humidity of air
    Inputs
        T - air temperature [degC]
        RH - relative humidity [%]
        Pa - air pressure [mb or hPa]
    """
    
    return 0.01*RH*qsat(T,Pa)
    
    
def satVapPres(Ta,P=1010):
    """ Calculates the saturation vapor pressure of air at temperature, T

    Inputs: T - temperature [C]
            P - air pressure [mb]
            
    """
    ew=np.power(10,((0.7859+0.03477*Ta)/(1+0.00412*Ta)))
    
    fw=1 + 1e-6*P*(4.5+0.0006*Ta**2)
    
    ew=fw*ew;
    
    return ew
    
def longwaveUp(Tw, epsilonw=0.97, sigma=5.67e-8):
    """
    Calculate upward longwave radiation from water
    
    Inputs:
        Tw - water temp [C]
        epsilonw - emissivity of air
        sigma - Stefan-Boltzman constan [W m-2 K-4]
    """
    
    return -epsilonw*sigma*(Tw+273.16)**4
    
def longwaveDown2(Ta,C,Pa,RH,r_LW=0.03):
    """
    Calculate downward longwave radiation that reaches the ocean
    
    Inputs:
        Ta - air temp [C]
        C - cloud cover [fraction]
        Pa - air pressure [mb or hPa]
        RH - relative humidity [%]
        r_LW - reflected fraction
    """
    
    sigma=5.67051e-8
    
    pv = vaporPres(Pa,RH,Ta)
    LW = sigma*(Ta+273.16)**4*(0.746+6.6*pv) # clear sky radiation

    return LW*(1+0.1762*C**2)*(1-r_LW) # correction for clouds and reflection

def longwaveDown(Ta,C,r_LW=0.03):
    """
    Downward longwave radiation using Martin and McCutcheon Formula
    """    
    alpha0 = 0.937e-5
    sigma=5.67051e-8
    
    epsilona = alpha0*(1.0 + 0.17 * C**2)*(Ta+273.16)**2

    return epsilona*sigma*(1.0-r_LW)*(Ta+273.16)**4

def vaporPres(P,RH,Ta):
    """
    Calculate vapor pressure
    Equation B33 in Kantha and Clayson
    
    Inputs:
        P - air pressure [mb]
        RH - relative humidity [%]
        Ta - air temp [C]
    """
    
    epsilon=0.62197
    
    q = qspec(Ta,RH,P)
    r = q/(1-q)
    return P*(r/(r+epsilon))
    
    
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
    
    try:
        if theta >= 0 and theta <90:
            theta=np.abs(theta-90)
        elif theta >= 90 and theta <= 360:
            theta=np.abs(450-theta)
    except:
        idx = operator.and_(theta>=0.,theta<90.)
        theta[idx] = np.abs(theta[idx]-90.)

        idx = operator.and_(theta>=90.,theta<=360.)
        theta[idx] = np.abs(450.-theta)
    
    u,v = pol2cart(theta*np.pi/180,rho)
    
    return u, v
    
def convertUV2SpeedDirn(u,v,convention='current'):
    """
    Convert velocity components into speed and direction (in degrees north)
    
    Set 'convention' = 'wind' to flip directions to "from" instead of "to"
    
    Adapted from matlab code:
    function [theta,rho] = cart2compass(u,v)
    %CART2COMPASS convert cartesian coordinates into
    % speed and direction data (degN).
    %    [THETA,RHO] = CART2COMPASS convert the vectors u and v
    %      from a cartesian reference system into rho (e.g. speed) with
    %      direction theta (degree North).
    %
    %   See also CART2POL
    %
    
    % Author: Arnaud Laurent
    % Creation : March 20th 2009
    % MATLAB version: R2007b
    %
    
    [theta,rho] = cart2pol(u,v);
    
    theta = theta*180/pi;
    
    idx = find(theta<0);
    theta(idx) = 360 + theta(idx);
    
    idx = find(theta>=0&theta<90);
    theta_comp(idx,1) = abs(theta(idx) - 90);
    
    idx = find(theta>=90&theta<=360);
    theta_comp(idx,1) = abs(450 - theta(idx));
    
    theta = theta_comp;
    """
    pi=np.pi
    
    theta,rho = cart2pol(u,v)
    
    theta = theta*180.0/pi
        
    #idx = np.argwhere(theta<0.0)
    idx = np.where(theta<0.0)
    theta[idx] = theta[idx] + 360.0
    
    #idx = np.argwhere(theta>=0.0&theta<90.0)
    fltr=operator.and_(theta>=0.0, theta<90.0)
    #idx = np.argwhere(fltr)
    idx = np.where(fltr)
    theta[idx] = np.abs(theta[idx] - 90.0)
    
    #idx = np.argwhere(theta>=90.0&theta<=360.0)
    fltr=operator.and_(theta>=90.0, theta<=360.0)
    #idx = np.argwhere(fltr)
    idx = np.where(fltr)
    theta[idx] = np.abs(450.0 - theta[idx])
    
    # flip the direction    
    if convention=='wind':
        theta = np.mod(theta+180.0, 360.0)
        
    return theta, rho    
    
def pol2cart(th,rho):
    """Convert polar coordinates to cartesian"""
    x = rho * np.cos(th)
    y = rho * np.sin(th)

    return x, y

def cart2pol(x,y):
    """
    Convert cartesian to polar coordinates
    """
    th = np.angle(x+1j*y)
    rho = np.abs(x+1j*y)
    
    return th, rho
