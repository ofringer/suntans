"""
Signal processing tools

Matt Rayson
Stanford University
"""
import numpy as np

import pdb

def powerspec2D(phi,dx=1.,dz=1.,window=np.hanning,quadrant=0):
    """
    Compute the 2D power spectrum of matrix phi

    Returns
        S_kz - 2d power spectra units [units**2 dx^{-1}]
        kx - x wave number (frequency)
        kz - z wave number (frequency)
        quadrant to return : 0 = upper right
                             1 = upper left
                             2 = lower left
                             3 = lower right
    """
    Mz,Mx = phi.shape

    # Remove the mean
    rowMean = phi.mean(axis=0)
    result = phi-rowMean
    #colMean = result.mean(axis=1)
    #result = result - colMean[...,np.newaxis]

    # Compute the window
    windowVal = window2d(Mx,Mz,windowfunc=window)
    result = result*windowVal

    # FFT part
    result = np.fft.fft2(result)

    # This is the power in units^2 / dx
    result = dx*dz*np.abs(result)**2

    # Get the frequencies(wavenumbers)
    kx = np.fft.fftfreq(int(Mx),d=dx/(2*np.pi))
    kz = np.fft.fftfreq(int(Mz),d=dz/(2*np.pi))

    if quadrant==0:
        # Only retain the positive frequencies
        result = result[1:Mz//2,1:Mx//2]
        kx = kx[1:Mx//2]
        kz = kz[1:Mz//2]
    elif quadrant==1:
        # positve z, negative x
        result = result[1:Mz//2,Mx//2::]
        kx = np.abs(kx[Mx//2::])
        kz = kz[1:Mz//2]
    elif quadrant==2:
        # Only retain the negative frequencies
        result = result[Mz//2::,Mx//2::]
        kx = np.abs(kx[Mx//2::])
        kz = np.abs(kz[Mz//2::])
    elif quadrant==3:
        # negative z, positive x
        result = result[Mz//2::,1:Mx//2]
        kx = kx[1:Mx//2]
        kz = np.abs(kz[Mz//2::])

    return result, kx, kz

def rotary_spectra(tsec,u,v,K=3,power=2.):
    """
    Calculates the rotary spectral kinetic energy from velocity.
    
    See Alford and Whitmont, 2007, JPO for details.
    """
    
    M = tsec.shape[0]
    dt = tsec[1]-tsec[0]
    t=np.arange(M,dtype=np.double)
    M_2 = np.floor(M/2)
    
    # Put the velocity in rotary form
    u_r = u+1j*v
    
    h_tk = window_sinetaper(M,K=K)
   
    # Weight the time-series and perform the fft
    u_r_t = u_r[...,np.newaxis,:]*h_tk
    S_k = np.fft.fft(u_r_t,axis=-1)
    S_k = dt *np.abs(S_k)**power
    S = np.mean(S_k,axis=-2)
        
    omega = np.fft.fftfreq(int(M),d=dt/(2*np.pi))
    
    domega = 1/(M*dt)
    
    # Extract the clockwise and counter-clockwise component
    omega_ccw = omega[0:M_2]
    omega_cw = omega[M_2::] # negative frequencies
    S_ccw = S[...,0:M_2]
    S_cw = S[...,M_2::]
    
    return omega_cw,omega_ccw,S_cw,S_ccw,domega

def integrate_rotspec(omega_cw,omega_ccw,S_cw,S_ccw,domega,omega_low=None,omega_high=None):
    """
    Integrate the kinetic energy of a rotary spectra beween two bands:
        omega_low and omega_high
    """

    #find the index limits for the integration range
    # note that ccw frequencies go from low to high
    #           cw frequences go from high to low
    if omega_low==None:
        t0_ccw=0
        t0_cw=0
    else:
        t0_cw = np.argwhere(omega_cw>=-omega_high)[0]
        t0_ccw = np.argwhere(omega_ccw>=omega_low)[0]
    if omega_high==None:
        t1_cw = omega_cw.shape[0]
        t1_ccw = omega_ccw.shape[0]
    else:
        t1_cw = np.argwhere(omega_cw<=-omega_low)[-1]
        t1_ccw = np.argwhere(omega_ccw<=omega_high)[-1]
    
    # Integrate under the spectrum to get the kinetic energy
    KE_ccw = 0.5*np.sum(S_ccw[...,t0_ccw:t1_ccw]*domega,axis=-1)
    KE_cw = 0.5*np.sum(S_cw[...,t0_cw:t1_cw]*domega,axis=-1)
    
    return KE_ccw, KE_cw
 
def eofsvd(M):
    """
    Compute empirical orthogonal function using singular
    value decomposition technique
    
    Inputs:
        - M : matrix with time along first axis and observation points along
          second
    Returns:
        - PC : The principal component amplitude
        - s : the eigenvalues
        - E : the eigenvectors in each column (EOFs)
    """
    # Remove the mean from the columns
    M = M - M.mean(axis=0)

    # 
    U,s,V = np.linalg.svd(M,full_matrices=False)

    # The principal components are U*s
    PC = U*s

    # Each row of V are the eigenvectors (EOFs) so transpose them so that
    # columns are
    E = V.T

    # Note that the values of s from the svd are the eigenvalues ^0.5

    return PC,s*s,E

def window2d(M,N,windowfunc=np.hanning,**kwargs):
    """
    2D window function
    """
    wc = windowfunc(N,**kwargs)
    wr = windowfunc(M,**kwargs)
    maskr,maskc = np.meshgrid(wr,wc)
    return maskr*maskc
    
def window_sinetaper(M,K=3):
    # Generate the time-domain taper for the FFT
    h_tk = np.zeros((K,M))
    for k in range(K):
        h_tk[k,:] = np.sqrt(2./(M+1.))*np.sin( (k+1)*np.pi*t / (M+1.) ) 

    return h_tk
 


