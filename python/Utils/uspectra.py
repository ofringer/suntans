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
import pdb

# For testing 
#import netcdfio      
    
class uspectra(object):
    """
        Object for performing uneven spectra operations on a dataset
    """
    # Properties
    nfft = None # Number of bands
    frq = None # 2*np.pi*np.array([10.0]) #None # Angular frequency bands (2*pi/f)    
    verbose = False
    
    # Units for plotting only
    tunits = 'rad s'
    yunits = 'm s$^{-1}$'
    

    def __init__(self,t,y,**kwargs):
        
        self.t = t # independent variable (t,x, etc)
        self.y = y # dependent variable
        self.nt = 2 # Determines the frequency resolution as nt*sampling_period
        self.window = None # Window function None or 'hanning'
        self.__dict__.update(kwargs)
        
        if type(self.t[0]) == datetime:
            # convert to seconds since the start
            self.t0 = getT0time(self.t)
        else:
            # Convert to units since the start
            self.t0 = getT0(self.t)
        
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
            self.w_n = 1.0
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
        plt.loglog(self.frq,abs(self.C))
        plt.xlabel('Freq. [%s$^{-1}$]'%self.tunits)
        plt.ylabel('Amp [%s]'%(self.yunits))
        plt.grid(b=True)
     
    def plotinterp(self,N=10,**kwargs):
        """ Plots the spectrum"""
        N*=len(self.t0)
        plt.plot(self.t0,self.invfft())
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
        
        return y
    

    def idealfilter(self,klow, khigh):
        fltr = operator.and_(khigh < self.frq, self.frq < klow)
        fltr = np.reshape(fltr,(len(self.frq),1))
        #pdb.set_trace()
        self.C *= fltr
        print self.C

    def butterfilter(self,klow, khigh, order=3):
        """
        Apply butterworth band-pass filter
        """
        #http://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter    
        from scipy.signal import butter, lfilter, freqz
        from scipy.signal import freqz
        
        def butter_bandpass(lowcut, highcut):
            #nyq = 0.5 * fs
            #nyq = 2.0/self.Ts
            nyq =  self.frq[-1]
            low =  lowcut / nyq
            high = highcut / nyq
            b, a = butter(order, [low, high], btype='band')
            return b, a
        
        def butter_gain(lowcut, highcut):
            b, a = butter_bandpass(lowcut, highcut)
            scalefact = np.pi/max(self.frq)
            w, h = freqz(b, a, worN=scalefact*self.frq)
            #plt.figure()
            #plt.plot(w,abs(h))

            return h

        fltr = np.reshape(butter_gain(klow, khigh),(len(self.frq),1))
        #plt.figure()
        #plt.plot(self.frq,abs(fltr))
        #plt.title('butter')
        #plt.figure()
        self.C *= fltr
        

    def interp(self,t,nbands=None):
        """ 
        Reconstruct the signal from the spectral coefficients
        """
        nt=len(t)
        y = np.zeros((nt,))
        
        if nbands == None:
            ii=-1
            for f in self.frq:
                ii+=1
                y+=np.real(self.C[ii])*np.cos(f*t) + np.imag(self.C[ii])*np.sin(f*t)
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
    
    


# Testing code    
#dbfile = 'C:/Projects/GOMGalveston/DATA/GalvestonObs.db'
#outvar = ['NetCDF_Filename','NetCDF_GroupID','StationName']
#tablename = 'observations'
#
#staname = 'Galveston Bay Entr Channel LB 11'
##staname = 'USCG Freeport, TX'
#varname = 'water_u'
#condition = 'Variable_Name = "%s" and StationName = "%s"' % (varname,staname)
#datau, query = netcdfio.queryNC(dbfile,outvar,tablename,condition)
#varname = 'water_v'
#condition = 'Variable_Name = "%s" and StationName = "%s"' % (varname,staname)
#datav, query = netcdfio.queryNC(dbfile,outvar,tablename,condition)
#
#
#
## Input variables
#t = datau[0]['time']
#y = np.ravel(datau[0]['water_u'][:,0,0,0])
#
#u = uspectra(t,y,nfft=1000)
#
#plt.figure()
#u.plot()
#plt.show()
