# -*- coding: utf-8 -*-
"""
Hybrid Grid class

Attempt at replicating Rusty's TriGrid

Created on Tue Oct 22 18:29:07 2013

@author: mrayson
"""

import numpy as np
from scipy import sparse

import matplotlib.pyplot as plt
import pdb

###
# HybridGrid class global definitions
###
# edge markers:
CUT_EDGE = 37 # the marker for a cut edge
OPEN_EDGE = 3
LAND_EDGE = 1
DELETED_EDGE = -1

# edge-cell markers ( the cell ids that reside in the edge array
BOUNDARY = -1 # cell marker for edge of domain
UNMESHED = -2 # cell marker for edges not yet meshed

class TriGridError(Exception):
    pass

class NoSuchEdgeError(TriGridError):
    pass

class NoSuchCellError(TriGridError):
    pass

class HybridGrid(object):
    """
    Class for dealing with grids that have cells with an 
    arbitrary number of sides.
    
    
    """
    
    _pnt2cells = None
    _pnt2edges = None
    
    def __init__(self,xp,yp,cells,nfaces=None,edges=None,\
        mark=None,grad=None,neigh=None,**kwargs):
        
        self.xp = xp
        self.yp = yp
        
        self.cells = cells
        self.Nc = len(cells)
        
        if nfaces==None:
            self.nfaces = 3*np.ones((self.Nc,),np.int)
            self.MAXFACES = 3
        else:
            self.nfaces = np.array(nfaces)
            self.MAXFACES = max(nfaces)
    
        if edges == None or mark == None or grad == None:
            #self.make_edges_from_cells()
            self.make_edges_from_cells_sparse()
        else:
            self.edges=edges
            self.mark=mark
            self.grad=grad
        
        if neigh == None:
            self.make_neigh_from_cells()
        else:
            self.neigh=neigh
        
        # Calculate the coordintes
        self.calc_centroids()
        #self.edge_centers()
        
        # Make sure the BCs are ok
        self.check_missing_bcs()
    
    def calc_centroids(self):
        """
        Calculate the centroid coordinates of each cell
        
        This needs to be done differently depending on geometry
         - Triangles : use circumcenter
         - Quads : calculate circumcenter of two triangles and take geometric mean
         - Other : mid-point (average coordinate)
        """ 
        xp = np.array(self.xp)
        yp = np.array(self.yp)

        self.xv = np.zeros((self.Nc,))
        self.yv = np.zeros((self.Nc,))
        
        for N in range(3,self.MAXFACES+1):
            ind = self.nfaces==N
            cells = self.cells[ind,0:N]
            
            #xtmp,ytmp = centroid(xp[cells],yp[cells],N)
            if N ==3:
                # Use the circumcenter for triangles
                xtmp,ytmp = circumcenter(xp[cells[:,0]],yp[cells[:,0]],\
                    xp[cells[:,1]],yp[cells[:,1]],xp[cells[:,2]],yp[cells[:,2]])
            elif N == 4:
                # Use the average circumcenters of the two interior triangles for quads
                xtmp1,ytmp1 = circumcenter(xp[cells[:,0]],yp[cells[:,0]],\
                    xp[cells[:,1]],yp[cells[:,1]],xp[cells[:,2]],yp[cells[:,2]])
                xtmp2,ytmp2 = circumcenter(xp[cells[:,0]],yp[cells[:,0]],\
                    xp[cells[:,2]],yp[cells[:,2]],xp[cells[:,3]],yp[cells[:,3]])
                
                # Arithmetic mean
                #xtmp = np.vstack((xtmp1,xtmp2)).mean(axis=0)
                #ytmp = np.vstack((ytmp1,ytmp2)).mean(axis=0)
                #geometric mean
                xtmp = np.sqrt(xtmp1*xtmp2)
                ytmp = np.sqrt(ytmp1*ytmp2)
                
            else:
                xtmp = xp[cells].mean(axis=-1)
                ytmp = yp[cells].mean(axis=-1)
                
            self.xv[ind] = xtmp
            self.yv[ind] = ytmp
        
    
    def edge_centers(self):
        
        self.xe = 0.5 * (self.xp[self.edges[:,0]] + self.xp[self.edges[:,1]])
        self.ye = 0.5 * (self.yp[self.edges[:,0]] + self.yp[self.edges[:,1]])
        
    def make_edges_from_cells_sparse(self):
        """
        Utilize sparse matrices to find the unique edges
        
        This is slower than make_edges_from_cells() but is more robust for 
        hybrid grid types. 
        """
        A = sparse.lil_matrix((self.Npoints(),self.Npoints()))
        Vp1 = sparse.lil_matrix((self.Npoints(),self.Npoints()))
        Vp2 = sparse.lil_matrix((self.Npoints(),self.Npoints()))
        
        for i in range(self.Ncells()):

            for j in range(self.nfaces[i]):

                p1 = self.cells[i,j]
                p2 = self.cells[i,(j+1)%self.nfaces[i]]
                
                if A[p1,p2]==0:
                    Vp1[p1,p2] = i
                    Vp1[p2,p1] = i
                else:
                    Vp2[p1,p2] = i
                    Vp2[p2,p1] = i
                
                # Graph is undirected so should be symmetric
                A[p1,p2] += 1.0
                A[p2,p1] += 1.0
    
        I,J,V = sparse.find(A)
        Ne = I.shape[0]
        #edges = [[I[k], J[k], A[I[k], J[k]]!=2, Vp1[I[k], J[k]]-1,Vp2[I[k], J[k]]-1] for k in range(Ne) if I[k]<=J[k]]
        edges = [[I[k], J[k], A[I[k], J[k]]!=2, Vp1[I[k], J[k]],Vp2[I[k], J[k]]] for k in range(Ne) if I[k]<=J[k]]
        
        edges = np.array(edges,dtype=int)
        Ne = edges.shape[0]
        self.edges = np.array([edges[ii,0:2] for ii in range(Ne)])
        self.mark = np.array([edges[ii,2] for ii in range(Ne)])
        self.grad = np.array([edges[ii,3:5] for ii in range(Ne)])
        
        # Now go back and set grad[2] = -1 for boundary cells
        for k in range(Ne):
            if self.mark[k]==1 and self.grad[k,1]==0:
                self.grad[k,1]=-1
        
            
    ######################
    # TriGrid functions  #
    # (with adjustments) #
    ######################
    def make_edges_from_cells(self):
        # iterate over cells, and for each cell, if it's index
        # is smaller than a neighbor or if no neighbor exists,
        # write an edge record
        edges = []
        default_marker = 0

        # this will get built on demand later.
        self._pnt2edges = None
        
        for i in range(self.Ncells()):
            # find the neighbors:
            # the first neighbor: need another cell that has
            # both self.cells[i,0] and self.cells[i,1] in its
            # list.
            my_set = set([i])
            #n = [-1,-1,-1]
            n = self.nfaces[i] * [-1]
            #for j in 0,1,2:
            #for j in range(self.MAXFACES):
            for j in range(self.nfaces[i]):
                pnt_a = self.cells[i,j]
                #pnt_b = self.cells[i,(j+1)%3]
                #pnt_b = self.cells[i,(j+1)%self.MAXFACES]
                pnt_b = self.cells[i,(j+1)%self.nfaces[j]]
                
                    
                adj1 = self.pnt2cells(pnt_a) # cells that use pnt_a
                adj2 = self.pnt2cells(pnt_b) # cells that use pnt_b

                # the intersection is us and our neighbor
                #  so difference out ourselves...
                neighbor = adj1.intersection(adj2).difference(my_set)
                # and maybe we ge a neighbor, maybe not (we're a boundary)
                if len(neighbor) == 1:
                    n = neighbor.pop()
                else:
                    n = -1
                    
                if n==-1 or i<n:
                    # we get to add the edge:
                    edges.append((pnt_a,
                                  pnt_b,
                                  default_marker,
                                  i,n))

        #self.edge = np.array(edges,np.int32)
        Ne = len(edges)
        edges = np.array(edges)
        self.edges = np.array([edges[ii,0:2] for ii in range(Ne)])
        self.mark = np.array([edges[ii,2] for ii in range(Ne)])
        self.grad = np.array([edges[ii,3:5] for ii in range(Ne)])
            
            
    def check_missing_bcs(self):
        """
        Check for missing BCs
        """        
        missing_bcs = (self.mark==0) & (self.grad[:,1]<0)
        n_missing = missing_bcs.sum()
        if any(missing_bcs):
            print "WARNING: %d edges are on the boundary but have marker==0"%n_missing
            print "Assuming they are closed boundaries!"

            self.mark[missing_bcs] = 1

        
    def make_neigh_from_cells(self):
        """
        Find the neighbouring cells
        """
        self.neigh = np.zeros((self.Ncells(),self.MAXFACES),np.int)
        for i in range(self.Ncells()):
            # find the neighbors:
            # the first neighbor: need another cell that has
            # both self.cells[i,0] and self.cells[i,1] in its
            # list.
            my_set = set([i])
            #n = self.MAXFACES * [-1]
            n = self.nfaces[i] * [-1]
            #for j in range(self.MAXFACES):
            for j in range(self.nfaces[i]):
                adj1 = self.pnt2cells(self.cells[i,j])
                adj2 = self.pnt2cells(self.cells[i,(j+1)%self.nfaces[i]])
                #adj2 = self.pnt2cells(self.cells[i,(j+1)%self.MAXFACES])
                neighbor = adj1.intersection(adj2).difference(my_set)
                if len(neighbor) == 1:
                    n[j] = neighbor.pop()
            
            self.neigh[i,0:self.nfaces[i]] = n
        
    def pnt2cells(self,pnt_i):
        if self._pnt2cells is None:
            # build hash table for point->cell lookup
            self._pnt2cells = {}
            for i in range(self.Ncells()):
                #for j in range(self.MAXFACES):
                for j in range(self.nfaces[i]):
                    if not self._pnt2cells.has_key(self.cells[i,j]):
                        self._pnt2cells[self.cells[i,j]] = set()
                    self._pnt2cells[self.cells[i,j]].add(i)
        return self._pnt2cells[pnt_i]
        
    def pnt2edges(self,pnt_i):
        if self._pnt2edges is None:
            # print "building pnt2edges"
            
            p2e = {}
            for e in range(self.Nedges()):
                # skip deleted edges
#                if self.edges[e,2] == DELETED_EDGE:
#                    continue
                if self.mark[e] == DELETED_EDGE:
                    continue
                
                for p in self.edges[e,:2]:
                    if not p2e.has_key(p):
                        p2e[p] = []
                    p2e[p].append(e)
            self._pnt2edges = p2e

        if self._pnt2edges.has_key(pnt_i):
            return self._pnt2edges[pnt_i]
        else:
            return []
    
    def find_edge(self,nodes):
    
        el0 = self.pnt2edges(nodes[0])
        el1 = self.pnt2edges(nodes[1])
        for e in el0:
            if e in el1:
                return e
        raise NoSuchEdgeError,str(nodes)
        
    def Nedges(self):
        return len(self.edges)
    def Ncells(self):
        return len(self.cells)
    def Npoints(self):
        return len(self.xp)
        
def polyarea(x,y,N):
    """
    Calculate the centroid of an arbitrary side polygon.
    
    Uses the formula here:
        http://www.seas.upenn.edu/~sys502/extra_materials/Polygon%20Area%20and%20Centroid.pdf
    """
    
    A = np.zeros(x.shape[:-1])
    for i in range(N-1):
        A += x[...,i]*y[...,i+1] - x[...,i+1]*y[...,i]
        
    return 0.5*A
    

def centroid(x,y,N):
    """
    **NOT USED (and not working...)**
    Calculate the centroid of an arbitrary side polygon.
    
    Uses the formula here:
        http://www.seas.upenn.edu/~sys502/extra_materials/Polygon%20Area%20and%20Centroid.pdf
    """
    
    A = polyarea(x,y,N)
    
    Cx = np.zeros(x.shape[:-1])
    Cy = np.zeros(x.shape[:-1])
    for i in range(N-1):
        tmp = x[...,i]*y[...,i+1] - x[...,i+1]*y[...,i]
        
        Cx += (x[...,i]+x[...,i+1])*tmp
        Cy += (y[...,i]+y[...,i+1])*tmp
        
    fac = 1. / (6. * A)
    
    return fac*Cx, fac*Cy
    
def circumcenter(p1x,p1y,p2x,p2y,p3x,p3y):
    refx = p1x.copy()
    refy = p1y.copy()
    
    p1x -= refx # ==0.0
    p1y -= refy # ==0.0
    p2x -= refx
    p2y -= refy
    p3x -= refx
    p3y -= refy

    vcx = np.zeros( p1x.shape, np.float64)
    vcy = np.zeros( p1y.shape, np.float64)
    
    # taken from TRANSFORMER_gang.f90
    dd=2.0*((p1x-p2x)*(p1y-p3y) -(p1x-p3x)*(p1y-p2y))
    b1=p1x**2+p1y**2-p2x**2-p2y**2
    b2=p1x**2+p1y**2-p3x**2-p3y**2 
    vcx=(b1*(p1y-p3y)-b2*(p1y-p2y))/dd + refx
    vcy=(b2*(p1x-p2x)-b1*(p1x-p3x))/dd + refy
    
    return vcx,vcy