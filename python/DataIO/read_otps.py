# -*- coding: utf-8 -*-
"""
Tools for handling the OSU tidal prediction software (OTPS) output data
(http://volkov.oce.orst.edu/tides/)

This software is based on the tide model driver (TMD) matlab code from here:
	http://polaris.esr.org/ptm_index.html

Matt Rayson
Stanford University
March 2013
"""

import os
import numpy as np

from interpXYZ import interpXYZ
import othertime
from datetime import datetime

import pdb

otis_constits = { 'M2':{'index':1,'omega':1.405189e-04,'v0u':1.731557546},\
     'S2':{'index':2,'omega':1.454441e-04,'v0u':0.000000000},\
     'N2':{'index':3,'omega':0.00013787970,'v0u':6.050721243},\
     'K2':{'index':4,'omega':0.0001458423,'v0u':3.487600001},\
     'K1':{'index':5,'omega':7.292117e-05,'v0u':0.173003674},\
     'O1':{'index':6,'omega':6.759774e-05,'v0u':1.558553872},\
     'P1':{'index':7,'omega':7.252295e-05,'v0u':6.110181633},\
     'Q1':{'index':8,'omega':6.495854e-05,'v0u':5.877717569}}    
    
def tide_pred(modfile,lon,lat,time,z=None,conlist=None):
    """
    Performs a tidal prediction at all points in [lon,lat] at times in vector [time]
    
    """
    
    # Read and interpolate the constituents
    u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist = extract_HC(modfile,lon,lat,z=z,conlist=conlist)
    
    # Initialise the output arrays
    sz = lon.shape
    nx = np.prod(sz)
    nt = time.shape[0]
    ncon = omega.shape[0]
    
    h_re = h_re.reshape((ncon,nx))
    h_im = h_im.reshape((ncon,nx))
    u_re = u_re.reshape((ncon,nx))
    u_im = u_im.reshape((ncon,nx))
    v_re = v_re.reshape((ncon,nx))
    v_im = v_im.reshape((ncon,nx))
        
    # Calculate nodal correction to amps and phases
    #baseyear = time[0].year
    #t1992 = othertime.SecondsSince(datetime(baseyear,1,1),basetime=datetime(1992,1,1))/86400.0
    t1992 = othertime.SecondsSince(time[0],basetime=datetime(1992,1,1))/86400.0
    
    pu,pf,v0u = nodal(t1992+48622.0,conlist)
       
    # Calculate the time series
    tsec = othertime.SecondsSince(time,basetime=datetime(1992,1,1)) # Needs to be referenced to 1992
    h=np.zeros((nt,nx))
    u=np.zeros((nt,nx))
    v=np.zeros((nt,nx))
    for nn,om in enumerate(omega):
        for ii in range(0,nx):
            h[:,ii] += pf[nn]*h_re[nn,ii] * np.cos(om*tsec + v0u[nn] + pu[nn]) - \
                pf[nn]*h_im[nn,ii] * np.sin(om*tsec + v0u[nn] + pu[nn])
                
            u[:,ii] += pf[nn]*u_re[nn,ii] * np.cos(om*tsec + v0u[nn] + pu[nn]) - \
                pf[nn]*u_im[nn,ii] * np.sin(om*tsec + v0u[nn] + pu[nn])
                
            v[:,ii] += pf[nn]*v_re[nn,ii] * np.cos(om*tsec + v0u[nn] + pu[nn]) - \
                pf[nn]*v_im[nn,ii] * np.sin(om*tsec + v0u[nn] + pu[nn])
    
    szo = (nt,)+sz
    return h.reshape(szo), u.reshape(szo), v.reshape(szo)
    

def tide_pred_old(modfile,lon,lat,time,z=None,conlist=None):
    """
	
	### UNUSED ###
    Performs a tidal prediction at all points in [lon,lat] at times in vector [time]
    
    """
    
    # Read and interpolate the constituents
    u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist = extract_HC(modfile,lon,lat,z=z,conlist=conlist)
    
    # Initialise the output arrays
    sz = lon.shape
    nx = np.prod(sz)
    nt = time.shape[0]
    ncon = omega.shape[0]
    
    h_re = h_re.reshape((ncon,nx))
    h_im = h_im.reshape((ncon,nx))
    u_re = u_re.reshape((ncon,nx))
    u_im = u_im.reshape((ncon,nx))
    v_re = v_re.reshape((ncon,nx))
    v_im = v_im.reshape((ncon,nx))
    
    # Nodal correction to amps and phases here...
    baseyear = time[0].year
    amp, phase = cart2pol(h_re, h_im)
    amp,phase = nodal_correction(baseyear,conlist, amp, phase)
    h_re, h_im = pol2cart(amp, phase)
    
    amp, phase = cart2pol(u_re, u_im)
    amp, phase = nodal_correction(baseyear,conlist, amp, phase)
    u_re, u_im = pol2cart(amp, phase)
        
    amp, phase = cart2pol(v_re, v_im)
    amp, phase = nodal_correction(baseyear,conlist, amp, phase)
    v_re, v_im = pol2cart(amp, phase)
    
    
    # Calculate the time series
    tsec = othertime.SecondsSince(time,basetime=datetime(baseyear,1,1))
    h=np.zeros((nt,nx))
    u=np.zeros((nt,nx))
    v=np.zeros((nt,nx))
    for nn,om in enumerate(omega):
        for ii in range(0,nx):
            h[:,ii] += h_re[nn,ii] * np.cos(om*tsec) + h_im[nn,ii] * np.sin(om*tsec)
            u[:,ii] += u_re[nn,ii] * np.cos(om*tsec) + u_im[nn,ii] * np.sin(om*tsec)
            v[:,ii] += v_re[nn,ii] * np.cos(om*tsec) + v_im[nn,ii] * np.sin(om*tsec)
    
    szo = (nt,)+sz
    return h.reshape(szo), u.reshape(szo), v.reshape(szo)
    
    
def extract_HC(modfile,lon,lat,z=None,conlist=None):
    """
    Extract harmonic constituents from OTIS binary output and interpolate onto points in lon,lat
    
    set "z" to specifiy depth for transport to velocity conversion
    
    set "constituents" in conlist
    
    Returns:
        u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist
    
    """
    
    ###
    # Make sure the longitude is between 0 and 360
    lon = np.mod(lon,360.0)
    
    ###
    # Read the filenames from the model file
    pathfile = os.path.split(modfile)
    path = pathfile[0]
    
    f = open(modfile,'r')
    hfile = path+'/' + f.readline().strip()
    uvfile = path+'/' + f.readline().strip()
    grdfile = path+'/' + f.readline().strip()
    f.close()
    
    ###
    # Read the grid file 
    X,Y,depth, mask = read_OTPS_grd(grdfile)
    #X[X>180.0] = 180.0 - X[X>180.0]
    mask = mask == 1
    
    # Create an interpolation object
    sz = lon.shape
    lon = lon.ravel()
    lat = lat.ravel()
    nx = lon.size
    F= interpXYZ(np.vstack((X[mask],Y[mask])).T,np.vstack((lon,lat)).T,method='idw',NNear=3,p=1.0)
    
    # Interpolate the model depths onto the points if z is none
    if z == None:
        z = F(depth[mask])
    else:
        z = np.abs(z) # make sure they are positive
        
    ###
    # Check that the constituents are in the file
    conOTIS = get_OTPS_constits(hfile)
    
    if conlist == None:
        conlist = conOTIS
        
    for vv in conlist:
        if not vv in conOTIS:
            print 'Warning: constituent name: %s not present in OTIS file.'%vv
            conlist.remove(vv)
            
    ###
    # Now go through and read the data for each 
    
    # Initialse the arrays
    ncon = len(conlist)
    u_re = np.zeros((ncon,nx))
    u_im = np.zeros((ncon,nx))
    v_re = np.zeros((ncon,nx))
    v_im = np.zeros((ncon,nx))
    h_re = np.zeros((ncon,nx))
    h_im = np.zeros((ncon,nx))
    omega = np.zeros((ncon,))
    
    for ii, vv in enumerate(conlist):
        idx = otis_constits[vv]['index']
        omega[ii] = otis_constits[vv]['omega']
        print 'Interpolating consituent: %s...'%vv
        
        # Read and interpolate h
        X ,Y, tmp_h_re, tmp_h_im = read_OTPS_h(hfile,idx)
        h_re[ii,:] = F(tmp_h_re[mask])
        h_im[ii,:] = F(tmp_h_im[mask])
                
        # Read and interpolate u and v - Note the conversion from transport to velocity
        X ,Y, tmp_u_re, tmp_u_im, tmp_v_re, tmp_v_im = read_OTPS_UV(uvfile,idx)
        u_re[ii,:] = F(tmp_u_re[mask]) / z
        u_im[ii,:] = F(tmp_u_im[mask]) / z
        v_re[ii,:] = F(tmp_v_re[mask]) / z
        v_im[ii,:] = F(tmp_v_im[mask]) / z
        
        
    # Return the arrays in their original shape
    szout = (ncon,) + sz
    return u_re.reshape(szout), u_im.reshape(szout), v_re.reshape(szout), \
        v_im.reshape(szout), h_re.reshape(szout), h_im.reshape(szout), omega, conlist


def nodal_correction(year,conlist,amp, phase):
        """
		### UNUSED ###
        Applies a lunar nodal correction to the amplitude and phase        
        
        Code modified from Rusty Holleman's GET_COMPONENTS code below...
        #
        # GET_COMPONENTS
        #   [UAMP,UPHASE,VAMP,VPHASE,HAMP,HPHASE]=GET_COMPONENTS(YEAR,OMEGAT,LUN_NODE,V0U,AG)
        #   calculates the tidal amplitudes and phases from the interpolated OTIS
        #   data in the AG matrix.
        #
        #   This code has been adapted from Brian Dushaw's matlab scripts
        #   obtained from http://909ers.apl.washington.edu/~dushaw/tidegui/tidegui.html
        #
        #function [uamp,uphase,vamp,vphase,hamp,hphase]=get_components(YEAR,omegat,lun_node,v0u,AG)
        """
        
        import tide_consts as tc

        #oneday=np.array( [335.62,  0,  322.55,  1.97,   334.63,  0.99,  -0.99,  321.57])
        oneday = {'M2':335.62, 'S2':0, 'N2':322.55, 'K2':1.97, 'O1':334.63, 'K1':0.99, 'P1':-0.99, 'Q1':321.57}
        #if year < 1970 or year > 2037:
        #    print 'Constants for prediction year are not available'
            #return None
        
        # Find the index
        JJ=[]
        od=np.zeros((len(conlist),1))
        for ii,vv in enumerate(conlist):
            
            jj=[item for item in range(len(tc.const_names)) if tc.const_names[item] == vv]
            if len(jj) > 0:
                JJ.append(jj)
                
            if oneday.has_key(vv):
                od[ii]=(np.pi/180)*oneday[vv]
                  
        I = int( np.where(year==tc.years)[0] )
        vou=tc.v0u[JJ,I]
        lunnod=tc.lun_nodes[JJ,I]
                
        vou=(np.pi/180)*vou
        
        
        #oneday=(np.pi/180)*oneday
        
        hamp = amp*lunnod
        #hphase = - oneday[JJ] + vou[JJ] - G
        hphase =  -od + vou - phase

        return hamp, hphase
    
def read_OTPS_UV(uvfile,ic):
    """
    Reads the tidal transport constituent data from an otis binary file
    
    ic = constituent number
    
    Returns: X, Y, h_re and h_im (Real and imaginary components)
    
    See this post on byte ordering
         http://stackoverflow.com/questions/1632673/python-file-slurp-w-endian-conversion
    """
    f = open(uvfile,'rb')
    #f = hfile
    # Try numpy
    ll = np.fromfile(f,dtype=np.int32,count=1)
    nm = np.fromfile(f,dtype=np.int32,count=3)
    th_lim = np.fromfile(f,dtype=np.float32,count=2)
    ph_lim = np.fromfile(f,dtype=np.float32,count=2)
    
    # Need to go from little endian to big endian
    ll.byteswap(True)
    nm.byteswap(True)
    th_lim.byteswap(True)
    ph_lim.byteswap(True)
    
    n = nm[0]
    m = nm[1]
    nc = nm[2]
    
    if ic < 1 or ic > nc:
        raise Exception,'ic must be > 1 and < %d'%ic
    
    # Read the actual data
    nskip = (ic-1)*(nm[0]*nm[1]*16+8) + 8 + ll - 28
    f.seek(nskip,1)
    
    htemp = np.fromfile(f,dtype=np.float32,count=4*n*m)
    htemp.byteswap(True)
    
    f.close()
    
    htemp = np.reshape(htemp,(m,4*n))
    U_re = htemp[:,0:4*n-3:4]
    U_im = htemp[:,1:4*n-2:4]
    V_re = htemp[:,2:4*n-1:4]
    V_im = htemp[:,3:4*n:4]
        
    X,Y = np.meshgrid(np.linspace(th_lim[0],th_lim[1],n),np.linspace(ph_lim[0],ph_lim[1],m))
    
    return X, Y, U_re, U_im, V_re, V_im

def read_OTPS_grd(grdfile):
    """
    Reads the grid data from an otis binary file

    Returns: X, Y, hz, mask
    
    See this post on byte ordering
         http://stackoverflow.com/questions/1632673/python-file-slurp-w-endian-conversion
    """
    f = open(grdfile,'rb')
    #
    ## Try numpy
    f.seek(4,0)
    n = np.fromfile(f,dtype=np.int32,count=1)
    m = np.fromfile(f,dtype=np.int32,count=1)
    lats = np.fromfile(f,dtype=np.float32,count=2)
    lons = np.fromfile(f,dtype=np.float32,count=2)
    dt = np.fromfile(f,dtype=np.float32,count=1)
    
    n.byteswap(True)
    m.byteswap(True)
    lats.byteswap(True)
    lons.byteswap(True)
    dt.byteswap(True)
    
    nob = np.fromfile(f,dtype=np.int32,count=1)
    nob.byteswap(True)
    if nob == 0: 
       f.seek(20,1)
       iob = []
    else:
       f.seek(8,1)
       iob = np.fromfile(f,dtype=np.int32,count=2*nob)
       iob.byteswap(True)
       iob = np.reshape(iob,(2,nob))
       f.seek(8,1)
    
    hz = np.fromfile(f,dtype=np.float32,count=n*m)
    f.seek(8,1)
    mask = np.fromfile(f,dtype=np.int32,count=n*m)
    
    hz.byteswap(True)
    mask.byteswap(True)
    
    hz = np.reshape(hz,(m,n))
    mask = np.reshape(mask,(m,n))
    
    f.close()
    
    X,Y = np.meshgrid(np.linspace(lons[0],lons[1],n),np.linspace(lats[0],lats[1],m))
    
    return X, Y ,hz, mask

def read_OTPS_h(hfile,ic):
    """
    Reads the elevation constituent data from an otis binary file
    
    ic = constituent number
    
    Returns: X, Y, h_re and h_im (Real and imaginary components)
    
    See this post on byte ordering
         http://stackoverflow.com/questions/1632673/python-file-slurp-w-endian-conversion
    """
    f = open(hfile,'rb')
    #f = hfile
    # Try numpy
    ll = np.fromfile(f,dtype=np.int32,count=1)
    nm = np.fromfile(f,dtype=np.int32,count=3)
    th_lim = np.fromfile(f,dtype=np.float32,count=2)
    ph_lim = np.fromfile(f,dtype=np.float32,count=2)
    
    # Need to go from little endian to big endian
    ll.byteswap(True)
    nm.byteswap(True)
    th_lim.byteswap(True)
    ph_lim.byteswap(True)
    
    n = nm[0]
    m = nm[1]
    nc = nm[2]
    
    if ic < 1 or ic > nc:
        raise Exception,'ic must be > 1 and < %d'%ic
        #return -1
    
    # Read the actual data
    nskip = (ic-1)*(nm[0]*nm[1]*8+8) + 8 + ll - 28
    f.seek(nskip,1)
    
    htemp = np.fromfile(f,dtype=np.float32,count=2*n*m)
    htemp.byteswap(True)
    
    #
    f.close()
    
    htemp = np.reshape(htemp,(m,2*n))
    h_re = htemp[:,0:2*n-1:2]
    h_im = htemp[:,1:2*n:2]
    
    X,Y = np.meshgrid(np.linspace(th_lim[0],th_lim[1],n),np.linspace(ph_lim[0],ph_lim[1],m))
    
    return X ,Y, h_re, h_im


def get_OTPS_constits(hfile):
    """
    Returns the list of constituents in the file
    """
    f = open(hfile,'rb')
    ll = np.fromfile(f,dtype=np.int32,count=1)
    nm = np.fromfile(f,dtype=np.int32,count=3)  
    
    ll.byteswap(True)
    nm.byteswap(True)
    
    f.close()
    
    ncon = nm[2]
    conList = []    
    for ii in range(1,ncon+1):
        for vv in otis_constits:
            if otis_constits[vv]['index']==ii:
                conList.append(vv)
                
    return conList
  
def cart2pol(re,im):
    
    amp = np.abs(re + 1j*im)
    phs = np.angle(re + 1j*im)
    
    return amp, phs
    
def pol2cart(amp,phs):
    
    re = amp * np.cos(phs)
    im = amp * np.sin(phs)
    
    return re, im
 
def astrol(time):
    """
    %function  [s,h,p,N]=astrol(time);
    %  Computes the basic astronomical mean longitudes  s, h, p, N.
    %  Note N is not N', i.e. N is decreasing with time.
    %  These formulae are for the period 1990 - 2010, and were derived
    %  by David Cartwright (personal comm., Nov. 1990).
    %  time is UTC in decimal MJD.
    %  All longitudes returned in degrees.
    %  R. D. Ray    Dec. 1990
    %  Non-vectorized version. Re-make for matlab by Lana Erofeeva, 2003
    % usage: [s,h,p,N]=astrol(time)
    %        time, MJD
    circle=360;
    T = time - 51544.4993;
    % mean longitude of moon
    % ----------------------
    s = 218.3164 + 13.17639648 * T;
    % mean longitude of sun
    % ---------------------
    h = 280.4661 +  0.98564736 * T;
    % mean longitude of lunar perigee
    % -------------------------------
    p =  83.3535 +  0.11140353 * T;
    % mean longitude of ascending lunar node
    % --------------------------------------
    N = 125.0445D0 -  0.05295377D0 * T;
    %
    s = mod(s,circle);
    h = mod(h,circle);
    p = mod(p,circle);
    N = mod(N,circle);
    
    """
    
    circle=360;
    T = time - 51544.4993;
    # mean longitude of moon
    # ----------------------
    s = 218.3164 + 13.17639648 * T;
    # mean longitude of sun
    # ---------------------
    h = 280.4661 +  0.98564736 * T;
    # mean longitude of lunar perigee
    # -------------------------------
    p =  83.3535 +  0.11140353 * T;
    # mean longitude of ascending lunar node
    # --------------------------------------
    N = 125.0445 -  0.05295377 * T;
    #
    s = np.mod(s,circle);
    h = np.mod(h,circle);
    p = np.mod(p,circle);
    N = np.mod(N,circle);
    
    return s,h,p,N
    
    
def nodal(time,con):
    """
    Nodal correction 
    
    Derived from the tide model driver matlab scipt: nodal.m
    """
    
    rad = np.pi/180.0
    
    s,h,p,omega=astrol(time)
#    
#    omega = 
   #
   #     determine nodal corrections f and u 
   #     -----------------------------------
    sinn = np.sin(omega*rad);
    cosn = np.cos(omega*rad);
    sin2n = np.sin(2*omega*rad);
    cos2n = np.cos(2*omega*rad);
    sin3n = np.sin(3*omega*rad);
    
    ndict={'M2':{'f':np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2),\
        'u':np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))/rad},\
        'S2':{'f':1.0, 'u':0.0},\
        'K1':{'f':np.sqrt((1.+.1158*cosn-.0029*cos2n)**2 + (.1554*sinn-.0029*sin2n)**2),\
            'u':np.arctan((-.1554*sinn+.0029*sin2n)/(1.+.1158*cosn-.0029*cos2n))/rad},\
        'O1':{'f':np.sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 + (0.189*sinn-0.0058*sin2n)**2),\
            'u':10.8*sinn - 1.3*sin2n + 0.2*sin3n},\
        'N2':{'f':np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2),\
            'u':np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))/rad},\
        'P1':{'f':1.0, 'u':0.0},\
        'K2':{'f':np.sqrt((1.+.2852*cosn+.0324*cos2n)**2 + (.3108*sinn+.0324*sin2n)**2),\
            'u':np.arctan(-(.3108*sinn+.0324*sin2n) /(1.+.2852*cosn+.0324*cos2n))/rad},\
        'Q1':{'f':np.sqrt((1.+.188*cosn)**2+(.188*sinn)**2),\
            'u':np.arctan(.189*sinn / (1.+.189*cosn))/rad} }
            
    # Prepare the output data
    ncon = len(con)
    pu = np.zeros((ncon,1))
    pf = np.ones((ncon,1))
    v0u = np.zeros((ncon,1))
    for ii,vv in enumerate(con):
        if ndict.has_key(vv):
            pu[ii,:] = ndict[vv]['u']*rad
            pf[ii,:] = ndict[vv]['f']
            
        if otis_constits.has_key(vv):
            v0u[ii,:] = otis_constits[vv]['v0u']
     
    return pu, pf, v0u
            
            
    
#%%
#f=zeros(nT,53);
#f(:,1) = 1;                                     % Sa
#f(:,2) = 1;                                     % Ssa
#f(:,3) = 1 - 0.130*cosn;                        % Mm
#f(:,4) = 1;                                     % MSf
#f(:,5) = 1.043 + 0.414*cosn;                    % Mf
#f(:,6) = sqrt((1+.203*cosn+.040*cos2n).^2 + ...
#              (.203*sinn+.040*sin2n).^2);        % Mt
#
#f(:,7) = 1;                                     % alpha1
#f(:,8) = sqrt((1.+.188*cosn).^2+(.188*sinn).^2);% 2Q1
#f(:,9) = f(:,8);                                % sigma1
#f(:,10) = f(:,8);                               % q1
#f(:,11) = f(:,8);                               % rho1
#f(:,12) = sqrt((1.0+0.189*cosn-0.0058*cos2n).^2 + ...
#                 (0.189*sinn-0.0058*sin2n).^2);% O1
#f(:, 13) = 1;                                   % tau1
#% tmp1  = 2.*cos(p*rad)+.4*cos((p-omega)*rad);
#% tmp2  = sin(p*rad)+.2*sin((p-omega)*rad);% Doodson's
#tmp1  = 1.36*cos(p*rad)+.267*cos((p-omega)*rad);% Ray's
#tmp2  = 0.64*sin(p*rad)+.135*sin((p-omega)*rad);
#f(:,14) = sqrt(tmp1.^2 + tmp2.^2);                % M1
#f(:,15) = sqrt((1.+.221*cosn).^2+(.221*sinn).^2);% chi1
#f(:,16) = 1;                                    % pi1
#f(:,17) = 1;                                    % P1
#f(:,18) = 1;                                    % S1
#f(:,19) = sqrt((1.+.1158*cosn-.0029*cos2n).^2 + ...
#                (.1554*sinn-.0029*sin2n).^2);  % K1
#f(:,20) = 1;                                    % psi1
#f(:,21) = 1;                                    % phi1
#f(:,22) = 1;                                    % theta1
#f(:,23) = sqrt((1.+.169*cosn).^2+(.227*sinn).^2); % J1
#f(:,24) = sqrt((1.0+0.640*cosn+0.134*cos2n).^2 + ...
#                (0.640*sinn+0.134*sin2n).^2 ); % OO1
#f(:,25) = sqrt((1.-.03731*cosn+.00052*cos2n).^2 + ...
#                (.03731*sinn-.00052*sin2n).^2);% 2N2
#f(:,26) = f(:,25);                                % mu2
#f(:,27) = f(:,25);                                % N2
#f(:,28) = f(:,25);                                % nu2
#f(:,29) = 1;                                    % M2a
#f(:,30) = f(:,25);                                % M2
#f(:,31) = 1;                                    % M2b
#f(:,32) = 1;                                    % lambda2
#temp1 = 1.-0.25*cos(2*p*rad)-0.11*cos((2*p-omega)*rad)-0.04*cosn;
#temp2 = 0.25*sin(2*p*rad)+0.11*sin((2*p-omega)*rad)+ 0.04*sinn;
#f(:,33) = sqrt(temp1.^2 + temp2.^2);              % L2
#f(:,34) = 1;                                    % T2
#f(:,35) = 1;                                    % S2
#f(:,36) = 1;                                    % R2
#f(:,37) = sqrt((1.+.2852*cosn+.0324*cos2n).^2 + ...
#                (.3108*sinn+.0324*sin2n).^2);  % K2
#f(:,38) = sqrt((1.+.436*cosn).^2+(.436*sinn).^2); % eta2
#f(:,39) = f(:,30).^2;                            % MNS2
#f(:,40) = f(:,30);                              % 2SM2
#f(:,41) = 1;   % wrong                          % M3
#f(:,42) = f(:,19).*f(:,30);                     % MK3
#f(:,43) = 1;                                    % S3
#f(:,44) = f(:,30).^2;                           % MN4
#f(:,45) = f(:,44);                              % M4
#f(:,46) = f(:,44);                              % MS4
#f(:,47) = f(:,30).*f(:,37);                     % MK4
#f(:,48) = 1;                                    % S4
#f(:,49) = 1;                                    % S5
#f(:,50) = f(:,30).^3;                           % M6
#f(:,51) = 1;                                    % S6
#f(:,52) = 1;                                    % S7
#f(:,53) = 1;                                    % S8
#%
#u=zeros(nT,53);
#u(:, 1) = 0;                                       % Sa
#u(:, 2) = 0;                                       % Ssa
#u(:, 3) = 0;                                       % Mm
#u(:, 4) = 0;                                       % MSf
#u(:, 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n;      % Mf
#u(:, 6) = atan(-(.203*sinn+.040*sin2n)./...
#             (1+.203*cosn+.040*cos2n))/rad;        % Mt
#u(:, 7) = 0;                                       % alpha1
#u(:, 8) = atan(.189*sinn./(1.+.189*cosn))/rad;     % 2Q1
#u(:, 9) = u(:,8);                                  % sigma1
#u(:,10) = u(:,8);                                  % q1
#u(:,11) = u(:,8);                                  % rho1
#u(:,12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n;       % O1
#u(:,13) = 0;                                       % tau1
#u(:,14) = atan2(tmp2,tmp1)/rad;                    % M1
#u(:,15) = atan(-.221*sinn./(1.+.221*cosn))/rad;    % chi1
#u(:,16) = 0;                                       % pi1
#u(:,17) = 0;                                       % P1
#u(:,18) = 0;                                       % S1
#u(:,19) = atan((-.1554*sinn+.0029*sin2n)./...
#           (1.+.1158*cosn-.0029*cos2n))/rad;       % K1
#u(:,20) = 0;                                       % psi1
#u(:,21) = 0;                                       % phi1
#u(:,22) = 0;                                       % theta1
#u(:,23) = atan(-.227*sinn./(1.+.169*cosn))/rad;    % J1
#u(:,24) = atan(-(.640*sinn+.134*sin2n)./...
#           (1.+.640*cosn+.134*cos2n))/rad;         % OO1
#u(:,25) = atan((-.03731*sinn+.00052*sin2n)./ ...
#           (1.-.03731*cosn+.00052*cos2n))/rad;     % 2N2
#u(:,26) = u(:,25);                                 % mu2
#u(:,27) = u(:,25);                                 % N2
#u(:,28) = u(:,25);                                 % nu2
#u(:,29) = 0;                                       % M2a
#u(:,30) = u(:,25);                                   % M2
#u(:,31) = 0;                                       % M2b
#u(:,32) = 0;                                       % lambda2
#u(:,33) = atan(-temp2./temp1)/rad ;                % L2
#u(:,34) = 0;                                       % T2
#u(:,35) = 0;                                       % S2
#u(:,36) = 0;                                       % R2
#u(:,37) = atan(-(.3108*sinn+.0324*sin2n)./ ...
#             (1.+.2852*cosn+.0324*cos2n))/rad;     % K2
#u(:,38) = atan(-.436*sinn./(1.+.436*cosn))/rad;    % eta2
#u(:,39) = u(:,30)*2;                               % MNS2
#u(:,40) = u(:,30);                                 % 2SM2
#u(:,41) = 1.5d0*u(:,30);                           % M3
#u(:,42) = u(:,30) + u(:,19);                       % MK3
#u(:,43) = 0;                                       % S3
#u(:,44) = u(:,30)*2;                               % MN4
#u(:,45) = u(:,44);                                 % M4
#u(:,46) = u(:,30);                                 % MS4
#u(:,47) = u(:,30)+u(:,37);                         % MK4
#u(:,48) = 0;                                       % S4
#u(:,49) = 0;                                       % S5
#u(:,50) = u(:,30)*3;                               % M6
#u(:,51) = 0;                                       % S6
#u(:,52) = 0;                                       % S7
#u(:,53) = 0;                                       % S8 
### 
# Testing
###


#grdfile = 'C:/Projects/GOMGalveston/DATA/Tides/DATA/grid_Mex'
#hfile = 'C:/Projects/GOMGalveston/DATA/Tides/DATA/h_Mex2010'
#uvfile = 'C:/Projects/GOMGalveston/DATA/Tides/DATA/UV_Mex2010'
#ic=6

#X, Y ,hz, mask = read_OTPS_grd(grdfile)
#plt.figure()
#plt.contourf(X,Y,hz,30)
#plt.colorbar()

#X,Y,h_re, h_im = read_OTPS_h(hfile,ic)
#plt.figure()
#plt.contourf(X,Y,h_re,30)
##plt.imshow(h_re)
#plt.colorbar()
#plt.show()

#X,Y,U_re, U_im, V_re, V_im = read_OTPS_UV(uvfile,ic)
#plt.figure()
#plt.contourf(X,Y,U_re,30)
##plt.imshow(h_re)
#plt.colorbar()
#plt.show()

# Inputs:
#modfile ='C:/Projects/GOMGalveston/DATA/Tides/Model_Mex'
#lon = np.array([-90.0,-91.0,-91.0])
#lat = np.array([27.0,28.0,28.5])
#z = None
#conlist = ['M2','S2','K1','O1','r2d2']
#time = othertime.TimeVector('20000101.0000','20000201.0000',3600.)
##
##tides = extract_HC(modfile,lon,lat)
#
#huv = tide_pred(modfile,lon,lat,time,z=None,conlist=None)