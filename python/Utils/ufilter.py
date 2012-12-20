# -*- coding: utf-8 -*-
"""
Tools for filtering unstructured grid data

Created on Tue Dec 18 13:19:39 2012

@author: mrayson
"""

import numpy as np
from scipy import spatial, sparse

import pdb

class ufilter(object):
    """
    Unstructured grid filter class
    """
    
    c = 4 # uses point c x p for filter
    
    def __init__(self,X,delta_f,**kwargs):
        
        self.X = X
        self.delta_f = delta_f
        
        self.__dict__.update(kwargs)
        
        s = np.shape(self.X)
        self.n = s[0] # Number of points
        if len(s)>1:
            self.m = s[1] # Number of dimensions
        else:
            self.m = 1
        
        self.GetP()
        
        self.BuildFilterMatrix()
        
    def __call__(self,y):
        """
        Performs the filtering operation on data in vector y
        """
        return self.G*y
        
    def BuildFilterMatrix(self):
        """
        Builds a sparse matrix, G, used to filter data in vector, y, via:
            y_filt = G x y
            
        """
        # Compute the spatial tree
        kd = spatial.cKDTree(self.X)
        eps=1e-6
        # Initialise the sparse matrix
        self.G = sparse.lil_matrix((self.n,self.n))
        
        for nn in range(self.n):
            #print nn
            # Find all of the points within c * p distance from point
            dx, i = kd.query(self.X[nn,:]+eps,k=self.n+1,distance_upper_bound=self.c*self.p)
            
            # Calculate the filter weights on these points
            ind = dx != np.inf
            Gtmp = self.Gaussian(dx[ind])
            
            # Insert the points into the sparse matrix
            I = i[ind]
            self.G[nn,I] = Gtmp
 
    def GetP(self):
        """
        Calculate the 'p' parameter
        """
        #self.p = self.delta_f**2/40.0
        self.p = self.delta_f
        
    def Gaussian(self,dx):
        """
        Calculate the Gaussian filter weights
        """      
        #Gtmp = 1.0 / (4.0*np.pi*self.p) * np.exp(-dx/np.sqrt(self.p))/dx
        #Gtmp = Gtmp[1:] # discard closest point (self)
        if self.m == 1:
            # 1D gaussian
            coef = 1.0/(np.sqrt(2.0*np.pi)*self.p)
        elif self.m == 2:
            # 2D gaussian
            coef = 1.0/(2.0*np.pi*self.p**2)
            
        Gtmp = coef * np.exp(- dx**2 / (2*self.p**2))
        
        return Gtmp / np.sum(Gtmp)