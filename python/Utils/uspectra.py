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
    method = 'lsq' # 'lsq' - least squares method; 'lomb' - Lomb-Scargle method; 'lsqfast' - faster the lsq but uses more memory
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
            self.istime = True
        else:
            # Convert to units since the start
            self.t0 = getT0(self.t)
            self.istime = False
        
        if self.method in ['lsq','lsqfast']:
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
            
            if self.method == 'lsq':
                self.C = lstsqfftseq(self.t0,self.y.copy()*self.w_n,self.frq,self.verbose)
            elif self.method =='lsqfast':
                self.C = lstsqfft(self.t0,self.y.copy()*self.w_n,self.frq)
                
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
        if not self.window==None:
            f = interp1d(self.t0,self.w_n,kind='linear')
            w_n = 1.0/f(t)
        else:
            w_n = 1.0
            
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
        
    def rmse(self):
        """
        Compute the root mean square error of the spectral time series
        """
        return np.sqrt(self.mse())
    
    def percvar(self):
        """
        Compute the percentage of variance captured by the fitted signal
        """
        return np.cov(self.interp(self.t0))/np.cov(self.y)*100.0

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
        
    def phsamp(self,phsbase=None):
        """
        Compute the phase and amplitude of the complex number
        
        The phase is offset relative to phsbase (if specified)
        
        phsbase can be datetime object or real
        """
        
        amp = np.abs(self.C)
        phs = np.angle(self.C)+np.pi # [0, 2*pi]
        
        if not phsbase == None:
            phs = np.mod(phs+phase_offset(self.frq,self.t[0],phsbase),2*np.pi)
        
        return amp, phs
    
    def __setitem__(self,key,value):
        """
        Performs another least-squares fft if the value of y is changed
        """
        if key == "y":
            self.y=value
            #print 'Updating C...'
            self.C = lstsqfft(self.t0,self.y.copy()*self.w_n,self.frq)
        else:
            self.__dict__[key]=value
        
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

def phase_offset(frq,start,base):
        """
        Compute a phase offset for a given fruequency
        """
        
        if type(start)==datetime:
            dx = start - base
            dx = dx.total_seconds()
        else:
            dx = start -base
        
        return np.mod(dx*np.array(frq),2*np.pi)
        
def getTideFreq(Fin=None):
    """
    Return a vector of frequency [rad s-1] of common tidal constituents
    
    """
    twopi= 2*np.pi
#    tidedict = {'M2':twopi/(12.42*3600.0), 
#                'S2':twopi/(12.00*3600.0), 
#                'N2':twopi/(12.66*3600.0),  
#                'K2':twopi/(11.97*3600.0), 
#                'K1':twopi/(23.93*3600.0), 
#                'O1':twopi/(25.85*3600.0), 
#                'P1':twopi/(24.07*3600.0), 
#                'Q1':twopi/(26.87*3600.0), 
#                'MF':twopi/(327.90*3600.0), 
#                'MM':twopi/(661.30*3600.0),
#                'M4':twopi/(6.21*3600.0)
#                }
                
    tidedict = {'J1':                           15.5854433,
    'K1':                           15.0410686,
    'K2':                           30.0821373,
    'L2':                           29.5284789,
    'M1':                           14.4966939,
    'M2':                           28.9841042,
    'M3':                           43.4761563,
    'M4':                           57.9682084,
    'M6':                           86.9523126,
    'M8':                          115.9364169,
    'N2':                           28.4397295,
    '2N2':                          27.8953548,
    'O1':                           13.9430356,
    'OO1':                          16.1391017,
    'P1':                           14.9589314,
    'Q1':                           13.3986609,
    '2Q1':                          12.8542862,
    'R2':                           30.0410667,
    'S1':                           15.0000000,
    'S2':                           30.0000000,
    'S4':                           60.0000000,
    'S6':                           90.0000000,
    'T2':                           29.9589333,
    'LAM2':                         29.4556253,
    'MU2':                          27.9682084,
    'NU2':                          28.5125831,
    'RHO':                         13.4715145,
    'MK3':                          44.0251729,
    '2MK3':                         42.9271398,
    'MN4':                          57.4238337,
    'MS4':                          58.9841042,
    '2SM2':                         31.0158958,
    'MF':                            1.0980331,
    'MSF':                           1.0158958,
    'MM':                            0.5443747,
    'SA':                            0.0410686,
    'SSA':                           0.0821373,
    'SA-IOS':                        0.0410667,
    'MF-IOS':                        1.0980331,
    'S1-IOS':                       15.0000020,
    'OO1-IOS':                      16.1391017,
    'R2-IOS':                       30.0410667,
    'A7':                            1.6424078,
    '2MK5':                         73.0092771,
    '2MK6':                         88.0503457,
    '2MN2':                         29.5284789,
    '2MN6':                         86.4079379,
    '2MS6':                         87.9682084,
    '2NM6':                         85.8635632,
    '2SK5':                         75.0410686,
    '2SM6':                         88.9841042,
    '3MK7':                        101.9933813,
    '3MN8':                        115.3920422,
    '3MS2':                         26.9523126,
    '3MS4':                         56.9523126,
    '3MS8':                        116.9523126,
    'ALP1':                         12.3827651,
    'BET1':                         14.4145567,
    'CHI1':                         14.5695476,
    'H1':                           28.9430375,
    'H2':                           29.0251709,
    'KJ2':                          30.6265120,
    'ETA2':                         30.6265120,
    'KQ1':                          16.6834764,
    'UPS1':                         16.6834764,
    'M10':                         144.9205211,
    'M12':                         173.9046253,
    'MK4':                          59.0662415,
    'MKS2':                         29.0662415,
    'MNS2':                         27.4238337,
    'EPS2':                         27.4238337,
    'MO3':                          42.9271398,
    'MP1':                          14.0251729,
    'TAU1':                         14.0251729,
    'MPS2':                         28.9430356,
    'MSK6':                         89.0662415,
    'MSM':                           0.4715211,
    'MSN2':                         30.5443747,
    'MSN6':                         87.4238337,
    'NLK2':                         27.8860711,
    'NO1':                          14.4966939,
    'OP2':                          28.9019669,
    'OQ2':                          27.3509801,
    'PHI1':                         15.1232059,
    'KP1':                          15.1232059,
    'PI1':                          14.9178647,
    'TK1':                          14.9178647,
    'PSI1':                         15.0821353,
    'RP1':                          15.0821353,
    'S3':                           45.0000000,
    'SIG1':                         12.9271398,
    'SK3':                          45.0410686,
    'SK4':                          60.0821373,
    'SN4':                          58.4397295,
    'SNK6':                         88.5218668,
    'SO1':                          16.0569644,
    'SO3':                          43.9430356,
    'THE1':                         15.5125897,
    '2PO1':                         15.9748271,
    '2NS2':                         26.8794590,
    'MLN2S2':                       26.9523126,
    '2ML2S2':                       27.4966873,
    'SKM2':                         31.0980331,
    '2MS2K2':                       27.8039339,
    'MKL2S2':                       28.5947204,
    'M2(KS)2':                      29.1483788,
    '2SN(MK)2':                     29.3734880,
    '2KM(SN)2':                     30.7086493,
    'NO3':                          42.3827651,
    '2MLS4':                        57.4966873,
    'ML4':                          58.5125831,
    'N4':                           56.8794590,
    'SL4':                          59.5284789,
    'MNO5':                         71.3668693,
    '2MO5':                         71.9112440,
    'MSK5':                         74.0251729,
    '3KM5':                         74.1073101,
    '2MP5':                         72.9271398,
    '3MP5':                         71.9933813,
    'MNK5':                         72.4649024,
    '2NMLS6':                       85.3920422,
    'MSL6':                         88.5125831,
    '2ML6':                         87.4966873,
    '2MNLS6':                       85.9364169,
    '3MLS6':                        86.4807916,
    '2MNO7':                       100.3509735,
    '2NMK7':                       100.9046319,
    '2MSO7':                       101.9112440,
    'MSKO7':                       103.0092771,
    '2MSN8':                       116.4079379,
    '2(MS)8':                      117.9682084,
    '2(MN)8':                      114.8476675,
    '2MSL8':                       117.4966873,
    '4MLS8':                       115.4648958,
    '3ML8':                        116.4807916,
    '3MK8':                        117.0344499,
    '2MSK8':                       118.0503457,
    '2M2NK9':                      129.8887361,
    '3MNK9':                       130.4331108,
    '4MK9':                        130.9774855,
    '3MSK9':                       131.9933813,
    '4MN10':                       144.3761464,
    '3MNS10':                      145.3920422,
    '4MS10':                       145.9364169,
    '3MSL10':                      146.4807916,
    '3M2S10':                      146.9523126,
    '4MSK11':                      160.9774855,
    '4MNS12':                      174.3761464,
    '5MS12':                       174.9205211,
    '4MSL12':                      175.4648958,
    '4M2S12':                      175.9364169,
    'M1C':                          14.4920521,
    '3MKS2':                        26.8701754,
    'OQ2-HORN':                     27.3416965,
    'MSK2':                         28.9019669,
    'MSP2':                         29.0251729,
    '2MP3':                         43.0092771,
    '4MS4':                         55.9364169,
    '2MNS4':                        56.4079379,
    '2MSK4':                        57.8860711,
    '3MN4':                         58.5125831,
    '2MSN4':                        59.5284789,
    '3MK5':                         71.9112440,
    '3MO5':                         73.0092771,
    '3MNS6':                        85.3920422,
    '4MS6':                         85.9364169,
    '2MNU6':                        86.4807916,
    '3MSK6':                        86.8701754,
    'MKNU6':                        87.5788246,
    '3MSN6':                        88.5125831,
    'M7':                          101.4490066,
    '2MNK8':                       116.4900752,
    '2(MS)N10':                    146.4079379,
    'MNUS2':                        27.4966873,
    '2MK2':                         27.8860711}
    
    # Set default constituents    
    if Fin==None:
        #Fin=tidedict.keys()
        Fin = ['M2','S2','N2','K2','K1','O1','P1','Q1','M4']
        
    frq = []
    Fout = Fin[:]
    for f in Fin:
        if tidedict.has_key(f):
            frq.append(twopi/(360.0/tidedict[f]*3600.0))
        else:
            'Warning: could not find constituent name: %s'%f
            Fout.remove(f)
        
    return frq, Fout

    
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