""" 
    Spectral analysis tools for unevenly spaced data
    
    M.Rayson
    Stanford University
    September 2012
"""

### TO DO LIST
# "Interpolation" function - needs to calculate a phase adjustment
#
# Spectral band smooting - moving average
#
# Band-pass spectral filter



import numpy as np
from datetime import datetime  
import matplotlib.pyplot as plt
import operator
from scipy.interpolate import interp1d
import pdb

# For testing 
#import netcdfio      
    
class uspectra(object):
    """
        Object for performing uneven spectra operations on a dataset
    """
    # Properties
    method = 'lsq' # 'lsq' - least squares method; 'lomb' - Lomb-Scargle method 
    nfft = None # Number of bands
    frq = None # 2*np.pi*np.array([10.0]) #None # Angular frequency bands (2*pi/f)    
    verbose = False
    nt = 2 # Determines the frequency resolution as nt*sampling_period
    window = None # Window function None or 'hanning'
    nbandavg = 1 # Number of frequency bands to average together
    
    # Units for plotting only
    tunits = 'rad s'
    yunits = 'm s$^{-1}$'
    

    def __init__(self,t,y,**kwargs):
        
        self.t = t # independent variable (t,x, etc)
        self.y = y # dependent variable
        
        self.__dict__.update(kwargs)
        
        if type(self.t[0]) == datetime:
            # convert to seconds since the start
            self.t0 = getT0time(self.t)
        else:
            # Convert to units since the start
            self.t0 = getT0(self.t)
        
        if self.method == 'lsq':
            # Least-squares method - need to find the frequencies
            if self.nfft == None and self.frq == None:
                # Work out the number of frequency bands and frequencies
                self.getTs()
                self.getk()            
                self.getNFFT()
                self.getfrq()
            elif self.frq == None:
                self.getTs()  
                self.getfrq()
                
            self.N = len(self.t)    
            self.getWindow()
            
            self.C = lstsqfftseq(self.t0,self.y.copy()*self.w_n,self.frq,self.verbose)
            self.C= np.ravel(self.C)
        elif self.method=='lomb':
            # Lomb scargle method
            self.N = len(self.t)    
            self.getWindow()
            
            self.C, self.frq, wk2,phs = lspr(self.t0,self.y.copy()*self.w_n,ofac=2.0,verbose=self.verbose)
            self.frq *= 2.0*np.pi
            
            # Convert the power amplitude to spectral amp
            #amp =np.sqrt(2.0*wk2*2.0*np.var(self.y)/self.N)
            #phs = phs
            
            amp =2.0*np.abs(self.C)/self.N
            phs = np.angle(self.C)
            self.C = amp*np.cos(phs) + 1j*amp*np.sin(phs)
            
        
        # Band average (not sure if this is correct)
        self.bandAvg()
        
    def getTs(self):
        """ Get the sampling rate from the data"""
        sample_rate =np.diff(self.t0)
        self.Ts = sample_rate.mean()
        #self.Ts = sample_rate.min()
        
    def getk(self):
        """ Frequency resolution """
        self.k = 1/(self.nt*self.Ts) # Nyquist frequency
    
    def getNFFT(self):
        """ Find the number of frequency bins """
        #N = len(self.t0)
        #logpow = np.floor(np.log2(N/2))
        #self.nfft = int(2**logpow) # number of frequencies
        self.nfft = np.floor(self.t0[-1]*self.k)
        
    def getWindow(self):
        """ Gets the window function w_n"""
        if self.window == None:
            self.w_n = np.ones((self.N,))
        elif self.window == 'hanning':
            n =  np.floor(self.N * self.t0/self.t0[-1])
            self.w_n = 0.5 * (1.0 - np.cos(2*np.pi*n / (self.N-1.0)))
            
    
    def getfrq(self):
        """ Get the fruequeny value"""
        #Tsmin = self.Ts*2
        #Tsmax = self.t0[-1]
        #self.frq=np.linspace(np.pi/Tsmax,np.pi/Tsmin,self.nfft) 
        #self.frq=np.linspace(0.5*np.pi/Tsmax,0.5*np.pi/Tsmin,self.nfft) 
        
        # This ensures orthogonality
        f = 2*np.pi*self.k/self.nfft
        self.frq = np.arange(0,self.nfft)*f
        
    def plot(self,**kwargs):
        """ Plots the spectrum"""
        plt.loglog(self.frq,abs(self.C),**kwargs)
        plt.xlabel('Freq. [%s$^{-1}$]'%self.tunits)
        plt.ylabel('Amp [%s]'%(self.yunits))
        plt.grid(b=True)
     
    def plotinterp(self,N=10,**kwargs):
        """ Plots the spectrum"""
        N*=len(self.t0)
        plt.plot(self.t0,self.invfft(),kwargs)
        tint = np.linspace(min(self.t0),max(self.t0),N)
        yint = self.interp(tint)
        plt.plot(tint, yint)
        plt.xlabel('Freq. [%s$^{-1}$]'%self.tunits)
        plt.ylabel('Amp [%s]'%(self.yunits))
        plt.grid(b=True)

    def invfft(self):
        """ 
        Reconstruct the signal from the spectral coefficients
        """
        nt=len(self.t0)
        y = np.zeros((nt,))
        
        ii=-1
        for f in self.frq:
            ii+=1
            y+=np.real(self.C[ii])*np.cos(f*self.t0) + np.imag(self.C[ii])*np.sin(f*self.t0)
        
        # Needs to be scaled by the window weights
        return y/self.w_n
    

    def idealfilter(self,klow, khigh):
        fltr = operator.and_(self.frq<khigh, self.frq > klow)
        #fltr = np.reshape(fltr,(len(self.frq),1))
        self.C *= fltr
        

    def butterfilter(self,klow, khigh, order=3):
        """
        Apply butterworth band-pass filter
        """
        #http://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter    
        from scipy.signal import butter, freqz
        
        def butter_bandpass(lowcut, highcut):
            nyq =  self.frq[-1]
            low =  lowcut / nyq
            high = highcut / nyq
            b, a = butter(order, [low, high], btype='band')
            return b, a
        
        def butter_gain(lowcut, highcut):
            b, a = butter_bandpass(lowcut, highcut)
            #scalefact = np.pi/max(self.frq)
            scalefact = np.pi/self.frq[-1]
            w, h = freqz(b, a, worN=scalefact*self.frq)
            return h

        #fltr = np.reshape(butter_gain(klow, khigh),(len(self.frq),1))
        #pdb.set_trace()
        #self.C *= fltr
        self.C *= butter_gain(klow, khigh)
        

    def interp(self,t,nbands=None):
        """ 
        Reconstruct the signal from the spectral coefficients
        """
        nt=len(t)
        y = np.zeros((nt,))
        
        # Output needs to be scaled by the inverse of the window weights
        # The weights need to be interpolated onto the new x variable
        f = interp1d(self.t0,self.w_n,kind='linear')
        w_n = 1.0/f(t)
            
        if nbands == None:
            ii=-1
            for f in self.frq:
                ii+=1
                y+=np.real(self.C[ii])*w_n*np.cos(f*t) + np.imag(self.C[ii])*w_n*np.sin(f*t)
        else:
            # Only use the highest nbands
            C,frq = self.rankBands(nbands)
            ii=-1
            for f in frq:
                ii+=1
                y+=np.real(C[ii])*np.cos(f*t) + np.imag(C[ii])*np.sin(f*t)

        return y
    
    def mse(self):
        """
        Compute the mean square error of the spectral time series
        """
 
        return 1/float(self.N) * np.sum( np.power(self.y - self.interp(self.t0),2.0) )

        
    def rankBands(self,N):
        """ 
        Ranks the frequencies from highest to lowest magnitude and returns the 
        greatest N
        """
        amp = np.abs(self.C)
        
        ind = np.argsort(np.ravel(amp))
        amp = np.sort(np.ravel(amp))
        frq = self.frq[ind]
        return amp[-N:], frq[-N:]
    def bandAvg(self):
        """
        Band average the spectral coefficients
        def movingaverage(interval, window_size):
        """
        window = np.ones(int(self.nbandavg))/float(self.nbandavg)
        self.C =  np.convolve(self.C, window, 'same')
        
        
def lstsqfftseq(t,y,frq,verbose=False):
    """
        Performs a sequential least-squares fft
        Slower but more memory efficient than doing it for all frequencies
    
    """
    # Construct matrix A
    nt=len(t)
    nf=len(frq)
    C = np.zeros((nf,1),dtype='complex')
    A=np.ones((nt,2))
    ii=-1
    printstep = 5 
    printstep0 = 0
    if verbose:
        print 'Computing uneven spectra for %d bands...' % nf
    for ff in range(0,nf):
        ii+=1
        perccomplete = float(ii)/float(nf)*100.0
        if perccomplete > printstep0:
            if verbose:
                print '%d %% complete...'%(int(perccomplete))
                printstep0+=printstep
        
        # Construct the harmonic signal
        A[:,0]=np.cos(frq[ff]*t)
        A[:,1]=np.sin(frq[ff]*t)
        
        # Solve the least square problem
        b = np.linalg.lstsq(A,y)
        C[ii] = b[0][0] + 1j *  b[0][1]
        
        # Remove the fitted signal        
        y2=b[0][0]*np.cos(ff*t) + b[0][1]*np.sin(ff*t)
        y-=y2
     
    return C
    
def lstsqfft(t,y,frq):
    # Do a least squares
    
    # Construct matrix A
    nt=len(t)
    nf=len(frq)
    nff=nf*2+1
    A=np.ones((nt,nff))
    for ff in range(0,nf):
        A[:,ff*2+1]=np.cos(frq[ff]*t)
        A[:,ff*2+2]=np.sin(frq[ff]*t)
        
    # Solve the least square problem
    b = np.linalg.lstsq(A,y)
    A = b[0][1::2]
    B = b[0][2::2]
    
    return A+1j*B




def getT0time(t):
    #print 'Convert the time to seconds since the start...'
    tsec = np.zeros((len(t)))
    ii=-1
    for tt in t:
        ii+=1
        dt = tt-t[0]
        tsec[ii]=dt.total_seconds()
    return tsec

def getT0(t):
    #print 'Convert the t to "units" since the start...'
    tsec = np.zeros((len(t)))
    ii=-1
    for tt in t:
        ii+=1
        dt = tt-t[0]
        tsec[ii]=dt
    return tsec

def getTideFreq(Fin=None):
    """
    Return a vector of frequencies of common tidal constituents
    """
    twopi= 2*np.pi
    tidedict = {'M2':twopi/(12.42*3600.0), 
                'S2':twopi/(12.00*3600.0), 
                'N2':twopi/(12.66*3600.0),  
                'K2':twopi/(11.97*3600.0), 
                'K1':twopi/(23.93*3600.0), 
                'O1':twopi/(25.85*3600.0), 
                'P1':twopi/(24.07*3600.0), 
                'Q1':twopi/(26.87*3600.0), 
                'MF':twopi/(327.90*3600.0), 
                'MM':twopi/(661.30*3600.0),
                'M4':twopi/(6.21*3600.0)
                }
    if Fin==None:
        Fin=tidedict.keys()
        
    frq = []
    for f in Fin:
        frq.append(tidedict[f])
        
    return np.sort(np.asarray(frq))

    
def lspr(x,y,ofac=4.0,hifac=1,verbose=True):
    """
    Lomb-Scargle periodogram modified from K. Hocke's matlab function
    
    Created on Mon Nov 26 20:02:32 2012
    
    @author: mrayson
    
    function [wk1,wk2,ph,vari,hifac,ofac,F,teta]=lspr(x,y,ofac)
    %///////////////////////////////////////////////////////////////
    % Lomb-Scargle periodogram (work version for reconstruction)  
    % procedure is based on the Numerical Recipe's programs 
    % period.f (please see there for comments and explanations; 
    % Numerical Recipes, 2nd edition, Chapter 13.8, by Press et al., Cambridge, 1992)
    % Here, the program code is adopted from the related IDL program lnp.pro
    % and is translated to Matlab. New features are the 
    % phase determination (Hocke, Ann. Geophys. 16, 356-358,1998) 
    % and the output of a complex Fourier spectrum F. 
    % This spectrum can be used for inverse FFT and reconstruction of an evenly 
    % spaced time series (Scargle, 1989).
    %    
    % ATTENTION: 
    % -> Because of the long story of program development and some open problems  
    % -> of phase definition and construction of the FFT spectrum, 
    % -> the program must be regarded as a working and discussion version 
    % -> without any warranty!  
    % -> Particularly the phase determination with the Lomb-Scargle
    % -> periodogram has been done in a heuristic manner. Switching between 
    % -> the phase reference systems may introduce errors which I am not aware yet.   
    % -> Scargle (1989) gives more informations on the details  of the problem. 
    %    (K. Hocke, Nov. 2007).
    %
    %  program call: 
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %  [wk1,wk2,ph,vari,hifac,ofac,F]=lspr(x,y);
    %  or [wk1,wk2,ph,vari,hifac,ofac,F]=lspr(x,y,ofac);
    %;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 
    %  input:
    %  x: e.g. time vector 
    %  y: observational data y(x)
    %  ofac: optional (oversampling factor , integer, default is 4)
    %
    %  output:
    %  wk1: frequency axis ( a vector of increasing linear frequencies)
    %  wk2: Lomb normalized power as function of wk1-vector
    %  ph:  phase vector as function of wk1-vector. The phase 
    %       is in radians and is defined to be the argument 
    %       of the cosine wave at the time x=0 ! 
    %  vari: sigma^2,  variance of y, necessary to derive 
    %       the amplitude from the normalized power wk2 
    %  F:   complex Pseudo-Fourier spectrum 
    %
    %  please check the phases and their signs before interpreting the phases! 
    %
    %  keywords:
    %  ofac: oversampling factor , integer  
    %        The default is 4
    %  hifac: integer, 1 for frequencies up to the Nyquist frequency 
    %         (2 for 2*Nyquist frequency)
    % 
    % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    """

    xstart=x[0]
    x=x-xstart  # simplifies the FFT construction (wx(1)=0 ) 

    twopi=6.2831853071795865
    n=len(x)
    #s=size(y); if s(1)==1;y=y';end  % transpose the series if necessary
    nout=int(0.5*ofac*hifac*n)
    nmax=nout

    wi=np.zeros((nmax,))
    wpi=np.zeros((nmax,))
    wpr=np.zeros((nmax,))
    wr=np.zeros((nmax,))
    wtemp=np.zeros((nmax,))
    px=np.zeros((nmax,))
    py=np.zeros((nmax,))
    ph=np.zeros((nmax,))
    ph1=np.zeros((nmax,))
    Fx=np.zeros((nmax,))
    Fy=np.zeros((nmax,))
    teta=np.zeros((nmax,))
    ave=np.mean(y)
    vari=np.var(y)
    
    
    xmax=np.max(x)
    xmin=np.min(x)
    xdif=xmax-xmin
    xave=0.5*(xmax+xmin)

    pymax=0.0
    pnow=1.0/(xdif*ofac)
    arg=twopi*((x-xave)*pnow)
    wpr=-2.0*np.sin(0.5*arg)**2
    wpi=np.sin(arg)
    wr=np.cos(arg)
    wi=wpi
    yy = y-ave
    
    printstep = 5 
    printstep0 = 0
    if verbose:
        print 'Computing uneven spectra using Lomb-Scargle algorithm...' 
        
    for i in range(nout):
        perccomplete = float(i)/nout*100.0
        if perccomplete > printstep0:
            if verbose:
                print '%d %% complete...'%(int(perccomplete))
                printstep0+=printstep

        px[i]=pnow	
        sumsh=np.sum(wr*wi)
        sumc=np.sum((wr-wi)*(wr+wi))
        wtau=0.5*np.arctan2(2.0*sumsh,sumc)  
        swtau=np.sin(wtau)      
        cwtau=np.cos(wtau)
        ss=wi*cwtau-wr*swtau
        cc=wr*cwtau+wi*swtau
        sums=np.sum(ss**2)
        sumc=np.sum(cc**2)
        sumsy=np.sum(yy*ss)
        sumcy=np.sum(yy*cc)
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        iy=sumsy/np.sqrt(sums) # imaginary part of Lomb-Scargle spectral component
        ry=sumcy/np.sqrt(sumc) # real part 
        py[i]=0.5*(ry**2+iy**2)/vari # power
        # here, the FFT phase is computed from the Lomb-Scargle Phase 
        # at each new frequency 'pnow' by adding the phase shift 'arg0'     
        phLS=np.arctan2(iy,ry)            # phase of Lomb-Scargle spectrum 
        arg0=twopi*(xave+xstart)*pnow +wtau  # phase shift with respect to 0
        arg1=twopi*xave*pnow +wtau   # phase shift for FFT reconstruction 
        ph[i]=np.mod(phLS+arg0, twopi)  # phase with respect to 0
        ph1[i]=np.mod(phLS+arg1, twopi) # phase for complex FFT spectrum	
        pnow=pnow + 1.0/(ofac*xdif)    # next frequency
        teta[i]=np.mod(arg1,twopi)

        


    dim=2*nout + 1.0     #dimension of FFT spectrum
    fac=np.sqrt(vari*dim/2.0) 
    a=fac*np.sqrt(py)# amplitude vector for FFT
    Fx=a*np.cos(ph1)# real part of FFT spectrum
    Fy=a*np.sin(ph1)# imaginary part of FFT spectrum 
    ph=np.mod(ph +5.0*twopi, twopi)# for value range 0,..., 2 pi	
    wk1=px 
    wk2=py
    return Fx + 1j*Fy, wk1, wk2,ph # fft, freqeuncy
    
    #% Fourier spectrum F: arrangement of the Lomb-Scargle periodogram   
    #% as a Fourier spectrum in the manner of Matlab:
    #% (it is not fully clear yet if and how the complex Fourier spectrum 
    #% can be exactly constructed from the Lomb-Scargle periodogram.  
    #% The present heuristic approach works well for the FFT back transformation   
    #% of F and reconstruction of an evenly spaced series 'yfit' in the time domain (after 
    #% multiplication by a constant, 'yfit' fits to 'y') 
      
    #Fxr=flipud(Fx); Fxr(1)=[];  
    #Fyr=flipud(Fy); Fyr(1)=[];
    #complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram: 
    #F=[complex(ave,0)' complex(Fx,Fy)' complex(Fxr,-Fyr)'];  