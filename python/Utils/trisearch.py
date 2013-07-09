# -*- coding: utf-8 -*-
"""
Class for performing searches in a triangulation using only matplotlib and scipy functions

Created on Thu May 09 15:17:58 2013

@author: mrayson
"""

from matplotlib.tri import Triangulation
from scipy.spatial import cKDTree
from matplotlib.nxutils import points_inside_poly
import numpy as np
import operator as op

import pdb

class TriSearch(Triangulation):
    
    verbose=False
    
    def __init__(self, x, y, cells,**kwargs):
        Triangulation.__init__(self,x,y,cells)
        
        self.Nc = cells.shape[0]
        
    def __call__(self,xin,yin):
        """
        Performs the search 
        """
        
        # Step 1) Find the nearest node
        if self.verbose:
            print 'Finding nearest nodes...'
        self.xpt = xin
        self.ypt = yin
        
        self.cellind=self.tsearch(self.xpt,self.ypt)
        
        return self.cellind
    
    def tsearchold(self,xin,yin):
        """
        DEPRECATED
        
        non-vectorised version of tsearch
        """
        
        xyin  = np.vstack((xin,yin)).T
        node =  self.findnearest(xyin)
        Np = xin.shape[0]
        
        # Step 2) Find the cells surrounding each node and check if the point is inside each one       
        cellind = -1*np.ones((Np,),dtype=np.int32)
        if self.verbose:
            print 'Finding cells...'
        for nn in range(Np):
            #print '\t%d of %d'%(nn,Np)
            cell =  self.pnt2cells(node[nn])
            for cc in cell:
                if self.inCell(cc,xyin[nn]):
                    cellind[nn]=cc
                    continue
                
        return cellind
        
    def tsearch(self,xin,yin,MAXNODES=8):
        """
        Vectorized version of tseach
        
        """
        xyin  = np.vstack((xin,yin)).T
        node =  self.findnearest(xyin)
        Np = xin.shape[0]

        cell = -1*np.ones((Np,MAXNODES),dtype=np.int32)
        for nn in range(Np):
            p2c = self.pnt2cells(node[nn])
            cell[nn,0:len(p2c)]=p2c
            
        cellind = -1*np.ones((Np,),dtype=np.int32)
        for ii in range(MAXNODES):
            ind = op.and_(cell[:,ii]!=-1,cellind==-1)
            if any(ind):
                ind2 = self.inCellVec(cell[ind,ii],xin[ind],yin[ind])
                
                ind3=np.where(ind)
                
                cellind[ind3[0][ind2]]= cell[ind3[0][ind2],ii]
 
        return cellind
        
    def pnt2cells(self,pnt_i):
        """
        Returns the cell indices for a point, pnt_i
        
        (Stolen from Rusty's TriGrid class)
        """
        if not self.__dict__.has_key('_pnt2cells'):
            # build hash table for point->cell lookup
            self._pnt2cells = {}
            for i in range(self.Nc):
                for j in range(3):
                    if not self._pnt2cells.has_key(self.triangles[i,j]):
                        #self._pnt2cells[self.cells[i,j]] = set()
                        self._pnt2cells[self.triangles[i,j]] = []
                    #self._pnt2cells[self.cells[i,j]].add(i)
                    self._pnt2cells[self.triangles[i,j]].append(i)
        return self._pnt2cells[pnt_i]
        

    def findnearest(self,xy,NNear=1):
        """
        Returns the node indices of the closest points to the nx2 array xy
        
        Uses the scipy KDTree routine
        """
        
        if not self.__dict__.has_key('kd'):
            self.kd = cKDTree(np.vstack((self.x,self.y)).T)
    
        # Perform query on all of the points in the grid
        dist,ind=self.kd.query(xy,k=NNear)
        
        return ind
        
    def inCell(self,cellind,xy):
        """
        Check whether a point is inside a cell

        Basically a wrapper for points_inside_polyo function        
        """
        
        xyverts=np.zeros((4,2))                
        xyverts[:,0]= [self.x[self.triangles[cellind,0]],\
                        self.x[self.triangles[cellind,1]],\
                        self.x[self.triangles[cellind,2]],\
                        self.x[self.triangles[cellind,0]]]
        xyverts[:,1]= [self.y[self.triangles[cellind,0]],\
                        self.y[self.triangles[cellind,1]],\
                        self.y[self.triangles[cellind,2]],\
                        self.y[self.triangles[cellind,0]] ]
                        
        return points_inside_poly(xy.reshape(1,2),xyverts)
        
    def inCellVec(self,cellinds,x,y):
        """
        Check whether a point is inside a cell

        Basically a wrapper for points_inside_poly function        
        """
        
        nx = x.shape[0]
        xyvertsall=np.zeros((nx,4,2))  
        xyvertsall[:,0,0]= self.x[self.triangles[cellinds,0]]
        xyvertsall[:,1,0]= self.x[self.triangles[cellinds,1]]
        xyvertsall[:,2,0]= self.x[self.triangles[cellinds,2]]
        xyvertsall[:,3,0]= self.x[self.triangles[cellinds,0]]
        xyvertsall[:,0,1]= self.y[self.triangles[cellinds,0]]
        xyvertsall[:,1,1]= self.y[self.triangles[cellinds,1]]
        xyvertsall[:,2,1]= self.y[self.triangles[cellinds,2]]
        xyvertsall[:,3,1]= self.y[self.triangles[cellinds,0]] 
        
        
        xyverts=np.zeros((4,2))    
        def _xyverts(ii):
            xyverts[:,:] = xyvertsall[ii,:,:]
            return xyverts            
                            
        xy = np.zeros((1,2))
        def _xy(ii):
            xy[0,0]=x[ii]
            xy[0,1]=y[ii]
            return xy
            
        nx = x.shape[0]          
        inpoly = np.zeros((x.shape),dtype=np.bool)
        for ii in range(nx):
            inpoly[ii] = points_inside_poly(_xy(ii),_xyverts(ii))
        
        return inpoly
        
    def updatexy(self,xnew,ynew):
        """
        Finds the triangle index when x and y are updated
        
        Attempt at being faster than raw search performed during __call__
        """
        
        # Check that size of the arrrays match
        assert xnew.size == self.xpt.size, ' size of xnew must be the same as xin'        
        
        self.Nx = xnew.size
        
        # Check if the particle has crossed any edges (ie is it in the same cell)
        innewcell = np.zeros(xnew.shape,dtype=np.bool)
        innewcell[:]=True
        newcell = self.cellind.copy()
        
        # Check if the particle has crossed any edges (ie is it in the same cell)
        changedcell, neigh = self.checkEdgeCrossingVec(self.cellind,xnew,ynew,self.xpt,self.ypt)
        newcell[changedcell] = self.neighbors[self.cellind[changedcell],neigh[changedcell]]
        
        # Check if particle is actually in new cell
        innewcell[changedcell] = self.inCellVec(newcell[changedcell],xnew[changedcell],ynew[changedcell])
            
        cellind2 = self.tsearch(xnew[innewcell==False],ynew[innewcell==False])
        newcell[innewcell==False] = cellind2
        
        # Update the class attributes
        self.xpt=xnew
        self.ypt=ynew
        self.cellind=newcell

        return newcell
            
        
    def checkEdgeCrossing(self,cell_i,xnew, ynew, xold, yold):
        """
        Check to see if a particle has crossed any edge of a cell
        """

        p1 = Point(xold,yold)
        p2 = Point(xnew,ynew)
        A = Point(self.x[self.triangles[cell_i,0]],self.y[self.triangles[cell_i,0]])
        B = Point(self.x[self.triangles[cell_i,1]],self.y[self.triangles[cell_i,1]])
        C = Point(self.x[self.triangles[cell_i,2]],self.y[self.triangles[cell_i,2]])
        
        if intersect(p1,p2,A,B):
            return True, 0
        elif intersect(p1,p2,B,C):
            return True, 1
        elif intersect(p1,p2,C,A):
            return True, 2
        else:
            return False, -1
            
    def checkEdgeCrossingVec(self,cell_i,xnew, ynew, xold, yold):
        """
        Check to see if a particle has crossed any edge of a cell
        
        Vectorised version
        """

        p1 = Point(xold,yold)
        p2 = Point(xnew,ynew)
        A = Point(self.x[self.triangles[cell_i,0]],self.y[self.triangles[cell_i,0]])
        B = Point(self.x[self.triangles[cell_i,1]],self.y[self.triangles[cell_i,1]])
        C = Point(self.x[self.triangles[cell_i,2]],self.y[self.triangles[cell_i,2]])
        
        changedcell = np.zeros(xnew.shape,dtype=np.bool)
        neigh = -1*np.ones(xnew.shape,dtype=np.int)
        
        leftcell = intersectvec(p1,p2,A,B)
        changedcell[leftcell] = True
        neigh[leftcell] = 0
        
        leftcell = intersectvec(p1,p2,B,C)
        changedcell[leftcell] = True
        neigh[leftcell] = 1
        
        leftcell = intersectvec(p1,p2,C,A)
        changedcell[leftcell] = True
        neigh[leftcell] = 2
        
        return changedcell, neigh
        
        
#####
# Line crossing code from
#
# http://www.bryceboe.com/wordpress/wp-content/uploads/2006/10/intersect.py
#####
class Point:
	def __init__(self,x,y):
		self.x = x
		self.y = y

def ccw(A,B,C):
	return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

def intersect(A,B,C,D):
	return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def intersectvec(A,B,C,D):
	return op.and_( op.ne(ccwvec(A,C,D),ccwvec(B,C,D)),op.ne(ccwvec(A,B,C),ccwvec(A,B,D)) )

def ccwvec(A,B,C):
    return op.gt( (C.y-A.y)*(B.x-A.x),(B.y-A.y)*(C.x-A.x) )

#a = Point(0,0)
#b = Point(0,1)
#c = Point(1,1)
#d = Point(1,0)
#
#
#print intersect(a,b,c,d)
#print intersect(a,c,b,d)
#print intersect(a,d,b,c)
