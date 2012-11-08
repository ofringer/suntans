# -*- coding: utf-8 -*-
"""
    Kriging interpolation library
"""
from scipy import spatial
import numpy as np

import pdb

class kriging(object):
    
    """ Class for kriging interpolation"""
    
    ### Properties ###
    maxdist = 1000
    NNear = 12
    
    # Variogram paramters
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0
    
    verbose = True
    
    def __init__(self,XYin,XYout,**kwargs):
        self.__dict__.update(kwargs)
        
        self.XYin = XYin
        self.XYout = XYout
        
        self._buildWeights()
        
    def __call__(self,Zin):
        """
        Calls the interpolation function with the scalar in Zin
        """
        self.Z = np.zeros((self.Nc,))
        for ii in range(0,self.Nc):
            self.Z[ii] = np.dot(self.W[:,ii],Zin[self.ind[ii,:]])
            
        return self.Z
                
    def _buildWeights(self):
        """ Calculates the kriging weights for all of the points in the grid"""
        # Compute the spatial tree
        kd = spatial.cKDTree(self.XYin)
        
        # Perform query on all of the points in the grid
        dist,self.ind=kd.query(self.XYout,distance_upper_bound=self.maxdist,k=self.NNear)
        
        self.Nc = np.size(self.ind,axis=0)
        print '%d interpolation points.'%self.Nc
        # Now loop through and get the weights for each point
        self.W = np.zeros((self.NNear,self.Nc))

        # Print percentages
        p0=0
        pstep=5
        for ii in range(0,self.Nc):
            
            if self.verbose:
                pfinish = float(ii)/float(self.Nc)*100.0
                if  pfinish> p0:
                    print '%3.1f %% complete...'%pfinish
                    p0+=pstep
                                
            W = self.getWeights(dist[ii,:],self.XYin[self.ind[ii,:],0],self.XYin[self.ind[ii,:],1])
            self.W[:,ii] = W.T 
                
   
#            if all(dist[ii,:]!=np.inf):
#                self.W[:,ii] = self.getWeights(dist[ii,:],self.XYin[ind[ii,:],0],self.XYin[ind[ii,:],1])
#                
#            else:
#                Z[ii] = np.nan

        
        
    def getWeights(self,dist,xin,yin):
        
        """ Calculates the kriging weights point by point"""
        
        Ns = len(dist)
        
        # Construct the LHS matrix C
        C=np.ones((Ns+1,Ns+1))
        for i in range(0,Ns):
            C[i,i]=0
            for j in range(i+1,Ns):
                D = np.sqrt((xin[i]-xin[j])**2+(yin[i]-yin[j])**2)
                C[i,j] = self.semivariogram(D)
                C[j,i] = C[i,j]

        C[Ns,Ns]=0

        # Calculate the inverse of C 
        Cinv = np.linalg.inv(C)
        
        # Loop through each model point and calculate the vector D
        gamma = np.ones((Ns+1,1))
        
        for j in range(0,Ns):
            gamma[j,0]= self.semivariogram( dist[j])
        # Solve the matrix to get the weights
        W = np.dot(Cinv,gamma)
        W = W[:-1,:]
        
        #print np.size(gamma,axis=0),np.size(gamma,axis=1)   
        return 1.0/float(Ns)*np.ones((Ns,1))
        
    def semivariogram(self,D):
        """ Semivariogram functions"""
        if self.varmodel == 'spherical':
            if D > self.vrange:
                F = self.sill
            else:
                tmp = D/self.vrange
                F = self.nugget + (self.sill-self.nugget)*(1.5*tmp - 0.5*tmp**3)
        return F
  
        
        