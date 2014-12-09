# -*- coding: utf-8 -*-
"""
Hybrid Grid class

Attempt at replicating Rusty's TriGrid

Created on Tue Oct 22 18:29:07 2013

@author: mrayson
"""

import numpy as np
from scipy import sparse
import operator as op

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
        mark=None,grad=None,neigh=None,xv=None,yv=None,**kwargs):
        
        self.xp = np.array(xp)
        self.yp = np.array(yp)
        
        self.cells = cells
        self.Nc = len(cells)
    
        if nfaces==None:
            self.nfaces = 3*np.ones((self.Nc,),np.int)
            self.MAXFACES = 3
        else:
            self.nfaces = np.array(nfaces,dtype=np.int)
            self.MAXFACES = max(self.nfaces)
            
        # Make sure the nodes are rotated counter-clockwise
        self.cells = self.ensure_ccw()
        
        if edges == None or mark == None or grad == None:
            self.make_edges_from_cells()
            #self.make_edges_from_cells_sparse()
        else:
            self.edges=edges
            self.mark=mark
            self.grad=grad
        
        if neigh == None:
            self.make_neigh_from_cells()
        else:
            self.neigh=neigh


        # Face->edge connectivity
        self.face = self.cell_edge_map()
        
        # Calculate the coordintes
        self.edge_centers()
        if xv==None:
            self.calc_centroids()
        else:
            self.xv = xv
            self.yv = yv

        # Calculate distance and other metrics
        self.calc_dg()
        self.calc_def()
        self.calc_dfe()
        self.calc_df()
        self.calc_tangent()
        self.calc_unitnormal()
        self.calc_normal()

        
        # Make sure the BCs are ok
        self.check_missing_bcs()
     
    ###################################
    # Geometry functions
    ###################################
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

        self._xva = np.zeros((self.Nc,))
        self._yva = np.zeros((self.Nc,))
        self._xvb = np.zeros((self.Nc,))
        self._yvb = np.zeros((self.Nc,))
        self.xc = np.zeros((self.Nc,))
        self.yc = np.zeros((self.Nc,))
        
        self.xca = np.zeros((self.Nc,))
        self.yca = np.zeros((self.Nc,))
        
        for N in range(3,self.MAXFACES+1):
            ind = self.nfaces==N
            cells = self.cells[ind,0:N]
            
            #xtmp,ytmp = centroid(xp[cells],yp[cells],N)
            if N ==3:
                # Use the circumcenter for triangles
                xtmp,ytmp = circumcenter(xp[cells[:,0]],yp[cells[:,0]],\
                    xp[cells[:,1]],yp[cells[:,1]],xp[cells[:,2]],yp[cells[:,2]])
                    
                xtmpa = xtmp
                ytmpa = ytmp
                #self.xca[ind],self.yca[ind]  = centroid(xp[cells[:,0:N]],yp[cells[:,0:N]],N) 
                # Triangle just use circumcenter
                self.xca[ind] = xtmp
                self.yca[ind] = ytmp
                
            elif N == 4:
                xtmp1,ytmp1 = circumcenter(xp[cells[:,0]],yp[cells[:,0]],\
                    xp[cells[:,1]],yp[cells[:,1]],xp[cells[:,2]],yp[cells[:,2]])
                xtmp2,ytmp2 = circumcenter(xp[cells[:,0]],yp[cells[:,0]],\
                    xp[cells[:,2]],yp[cells[:,2]],xp[cells[:,3]],yp[cells[:,3]])
                    
                xtmp = np.sqrt(xtmp1*xtmp2)                
                ytmp = np.sqrt(ytmp1*ytmp2)
                
                xtmp1,ytmp1 = circumcenter(xp[cells[:,1]],yp[cells[:,1]],\
                    xp[cells[:,2]],yp[cells[:,2]],xp[cells[:,3]],yp[cells[:,3]])
                xtmp2,ytmp2 = circumcenter(xp[cells[:,3]],yp[cells[:,3]],\
                    xp[cells[:,0]],yp[cells[:,0]],xp[cells[:,1]],yp[cells[:,1]])
                    
                xtmpa = np.sqrt(xtmp1*xtmp2)                
                ytmpa = np.sqrt(ytmp1*ytmp2)
                
                self.xca[ind],self.yca[ind]  = centroid(xp[cells[:,0:N]],yp[cells[:,0:N]],N) 
                
            else:
                xtmp = xp[cells].mean(axis=-1)
                ytmp = yp[cells].mean(axis=-1)
                
                xtmpa = xtmp
                ytmpa = ytmp
                
                self.xca[ind],self.yca[ind]  = centroid(xp[cells[:,0:N]],yp[cells[:,0:N]],N)
                
                
            self._xva[ind] = xtmp
            self._yva[ind] = ytmp
            
            self._xvb[ind] = xtmpa
            self._yvb[ind] = ytmpa
            
            # Save centroids in xc,yc attributes
            self.xc[ind] = xp[cells].mean(axis=-1)
            self.yc[ind] = yp[cells].mean(axis=-1)
            
        ####
        # Option 1) use the most orthogonal mid-point
        ####
#        # Check the orthogonality of each of the alternative mid-points
#        ang1 = self.check_orthogonality(self._xva,self._yva)
#        ang2 = self.check_orthogonality(self._xvb,self._yvb)
#        
#        # Now go through and decide which point to use based on orthogonality
#        self.orthog = np.zeros((self.Nc,))
#        self.xv = np.zeros((self.Nc,))
#        self.yv = np.zeros((self.Nc,))
#        for i in range(self.Nc):
#            if ang1[i]<ang2[i]:
#                self.orthog[i]=ang1[i]
#                self.xv[i]=self._xva[i]
#                self.yv[i]=self._yva[i]
#            else:
#                self.orthog[i]=ang2[i]
#                self.xv[i]=self._xvb[i]
#                self.yv[i]=self._yvb[i]
                
        ####
        # Option 2) Use the area centroid (circumcenter for triangles)
        ### 
        self.xv = self.xca.copy()
        self.yv = self.yca.copy()
                

    def calc_area(self):
        """
        Calculates the area of each cell
        """
        xp = np.array(self.xp)
        yp = np.array(self.yp)
        
        Ac = np.zeros((self.Nc,))
        
        for N in range(3,self.MAXFACES+1):
            ind = self.nfaces==N
            cells = self.cells[ind,0:N]
            Ac[ind] = polyarea(xp[cells[:,0:N]],yp[cells[:,0:N]],N)            
            #Ac[ind] = signed_area(xp[cells[:,0:N]],yp[cells[:,0:N]],N)       
        
        return Ac
    
    def ensure_ccw(self):
        """
        Ensure that the nodes are rotated counter-clockwise
        """
        Ac = self.calc_area()
        cells_ccw = np.zeros_like(self.cells)
        for i in range(self.Nc):
            ctmp=self.cells[i,0:self.nfaces[i]].copy()
            if Ac[i] < 0:
                #print 'Cell %d is clock-wise - reversing.'%i
                cells_ccw[i,0:self.nfaces[i]] = ctmp[::-1] # reverse order
            else:
                cells_ccw[i,0:self.nfaces[i]] =  ctmp
        
        return cells_ccw
    
    def edge_centers(self):
        
        xp = np.array(self.xp)
        yp = np.array(self.yp)
        
        self.xe = 0.5 * (xp[self.edges[:,0]] + xp[self.edges[:,1]])
        self.ye = 0.5 * (yp[self.edges[:,0]] + yp[self.edges[:,1]])

    def calc_dg(self):
        """
        Manually calculate the distance between voronoi points, 'dg'
        """
        
        grad = self.grad.copy()
        Ne = self.Nedges()

        for ii in range(Ne):
            if grad[ii,0]==-1:
                grad[ii,0]=grad[ii,1]
            elif grad[ii,1]==-1:
                grad[ii,1]=grad[ii,0]
                
                
        x1 = self.xv[grad[:,0]]
        x2 = self.xv[grad[:,1]]
        y1 = self.yv[grad[:,0]]
        y2 = self.yv[grad[:,1]]
        
        dx=x1-x2
        dy=y1-y2
        
        self.dg = np.sqrt( dx*dx + dy*dy )

    def calc_def(self):
        """
        Calculate the edge to face(cell) distance
        
        dimensions: Nc x maxfaces
        """
        ne = np.array(self.face)

        cellmask = self.face==-1

        ne[cellmask]=0


        self.DEF = dist(self.xv,self.xe[ne].T,self.yv,self.ye[ne].T).T

        self.DEF = np.ma.masked_array(self.DEF,mask=cellmask)

    def calc_dfe(self):
        """
        Calculate the face(cell) to edge distance

        dimensions: Ne x 2
        """
        grad = self.grad.copy()
        mask = grad==-1
        grad[mask]=0

        Ne = self.Nedges()

        self.dfe = np.zeros((Ne,2))

        self.dfe[:,0] = dist(self.xv[grad[:,0]],self.xe,\
            self.yv[grad[:,0]],self.ye)
        self.dfe[:,1] = dist(self.xv[grad[:,1]],self.xe,\
            self.yv[grad[:,1]],self.ye)

        self.dfe[mask] = 0


    def calc_df(self):
        """
        Calculate the length of each edge segment
        """
        x = self.xp[self.edges]
        y = self.yp[self.edges]

        self.df = dist(x[:,0],x[:,1],y[:,0],y[:,1])

    def calc_tangent(self):
        """
        Calculate the tangential vector for the edges of each cell
        """
        dx = np.zeros(self.cells.shape)    
        dy = np.zeros(self.cells.shape)  

        dx[:,0:-1] = self.xp[self.cells[:,1::]] - self.xp[self.cells[:,0:-1]]               
        dy[:,0:-1] = self.yp[self.cells[:,1::]] - self.yp[self.cells[:,0:-1]]               

        for ii in range(self.Nc):
            dx[ii,self.nfaces[ii]-1] = self.xp[self.cells[ii,0]] -\
                self.xp[self.cells[ii,self.nfaces[ii]-1]]  
            dy[ii,self.nfaces[ii]-1] = self.yp[self.cells[ii,0]] -\
                self.yp[self.cells[ii,self.nfaces[ii]-1]]  

        mag = np.sqrt(dx*dx + dy*dy)
        
        self.tx = dx/mag
        self.ty = dy/mag

        #self.nx = -self.ty
        #self.ny = self.tx

        #return self._tx, self._ty, self._mag

    def calc_unitnormal(self):
        """
        Calculate the unit normal vector at each edge
        """

        dx = self.xp[self.edges[:,0]] - self.xp[self.edges[:,1]]
        dy = self.yp[self.edges[:,0]] - self.yp[self.edges[:,1]]

        mag = np.sqrt(dx*dx + dy*dy)
        
        self.nx = -dy/mag
        self.ny = dx/mag


    def calc_normal(self):
        """
        Create the normal array
        """
        Nc = self.Ncells()

        self.normal=np.zeros((Nc,self.MAXFACES))
        
        for ii in range(Nc):
            for nf in range(self.nfaces[ii]):
                if self.grad[self.face[ii,nf],1]==ii:
                    self.normal[ii,nf]=-1
                else:
                    self.normal[ii,nf]=1

        
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
    
    def check_orthogonality(self,xv,yv):
        """
        Checks the orthogonality of the grid cells with index 'cell' and mid-points
        'xv, yv'
        
        Returns the maximum deviation from 90 degrees of each line connecting the 
        edge point to the cell mid-point
        """
        print 'calculating orthogonality...'
        nc = xv.shape[0]
        orthoang = np.zeros((nc,))
        pi_on_2 = 0.5*np.pi
        
        for i in range(nc):
            
            maxangle = 0.0
            for j in range(self.nfaces[i]):
                # Find the two node values
                pnt_a = self.cells[i,j]
                pnt_b =self.cells[i,(j+1)%self.nfaces[j]]
                
                # Find the edge index
                edg = self.find_edge([pnt_a,pnt_b])
                
                # Create a two lines and find the angle between them
                P0 = Point(self.xe[edg],self.ye[edg])
                #P1a = Point(self.xp[pnt_a],self.yp[pnt_a])
                P1 = Point(self.xp[pnt_b],self.yp[pnt_b])
                P2 = Point(xv[i],yv[i])
                L1 = Line(P0,P1)
                L2 = Line(P0,P2)

#                plt.figure()
#                plt.plot(P0.x,P0.y,'ro')
#                plt.plot(P1a.x,P1a.y,'bo')
#                plt.plot(P1.x,P1.y,'bo')
#                plt.plot(P2.x,P2.y,'ko')
#                plt.show()
                
                ang = np.abs(pi_on_2 - L1.angle(L2))
                if ang > maxangle:
                    maxangle=ang
                    
            orthoang[i]=maxangle
        
        return orthoang
        
    ###########################
    # Input output functions
    ###########################
    def write2suntans(self,suntanspath):
        """
        Write to suntans grid format ascii files
        """
        ### Save cells.dat into hybrid grid format
        f = open(suntanspath+'/cells.dat','w')
    
        for ii in range(self.Ncells()):
            outstr = '%d %10.6f %10.6f '%(self.nfaces[ii],self.xv[ii],self.yv[ii])
            for nn in range(self.nfaces[ii]):
                outstr += '%d '%self.cells[ii,nn]
    
            for nn in range(self.nfaces[ii]):
                outstr += '%d '%self.neigh[ii,nn]
    
            outstr += '\n'
            f.write(outstr)
        
        f.close()
        
        # Save edges.dat
        f = open(suntanspath+'/edges.dat','w')
    
        for ee,m,gg in zip(self.edges,self.mark,self.grad):
            e1=ee[0]
            e2=ee[1]
            g1=gg[0]
            g2=gg[1]
            f.write('%d %d  %d  %d  %d  0\n'%(e1,e2,m,g1,g2))
    
        f.close()
        
        # Save to points.dat
        f = open(suntanspath+'/points.dat','w')
    
        for x,y in zip(self.xp,self.yp):
            f.write('%10.6f %10.6f  0\n'%(x,y))
    
        f.close()
        #print 'Completed gmsh to suntans conversion.'
        
    def create_dual_grid(self,minfaces=3,outpath=None):
        """
        Create a new grid using the dual of the current grid. 
        
        Returns a new hybridgrid object. Set outpath to save directly to a suntans
        grid.
        """
        # Locate the points forming each cell
        
        # The centre points are now the nodes
        xp = self.xv
        yp = self.yv
        
        # ...and the nodes are now the centre point
        xv = self.xp
        yv = self.yp
        
        # Find the number of faces of each cell
        Np = self.Npoints()
        nfaces = np.array([len(self.pnt2cells(ii)) for ii in range(Np)])
        
        maxfaces = nfaces.max()
        
        # Reorder the nodes into anti-clockwise order
        def reordercells(ii):
            cell = np.array(list(self.pnt2cells(ii)))
            # Find the order of the points that form a non-intersecting polygon
            xy = np.array([xp[cell],yp[cell]])
            # calculate the angles from the centre point and sort them
            ang = np.arctan2(xy[0,:]-xv[ii],xy[1,:]-yv[ii])
            
            order = np.argsort(ang)
            
            return cell[order]
            
        cells_list = map(reordercells,range(Np))
        
        cells = -1*np.ones((Np,maxfaces),np.int)
        for ii in range(Np):
            cells[ii,0:nfaces[ii]]=cells_list[ii]
            
        # Now remove cells with less than 'minfaces' points
        ind = nfaces>=minfaces
        
        cells = cells[ind,:]
        nfaces =  nfaces[ind]
        xv = xv[ind]
        yv = yv[ind]

        dualgrd = HybridGrid(xp,yp,cells,nfaces=nfaces,xv=xv,yv=yv)
        
        if not outpath == None:
            dualgrd.write2suntans(outpath)
            
        return dualgrd
             
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
                pnt_b = self.cells[i,(j+1)%self.nfaces[i]]
                
                    
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
        
    def cell2edges(self,cell_i):
        if self.cells[cell_i,0] == -1:
            raise "cell %i has been deleted"%cell_i
        
        # return indices to the three edges for this cell:
        pnts = self.cells[cell_i] # the three vertices

        # the k-th edge is opposite the k-th point, like in CGAL
        nf = self.nfaces[cell_i]
        edges = [ self.find_edge( (pnts[(i+1)%nf], pnts[(i+2)%nf]) ) for i in range(nf) ]
        return edges

    _cell_edge_map = None
    def cell_edge_map(self):
        """ cell2edges for the whole grid
        return an integer valued [Nc,3] array, where [i,k] is the edge index
        opposite point self.cells[i,k]

        N.B. this is not kept up to date when modifying the grid.
        """
        if self._cell_edge_map is None:
            cem = 999999*np.ones( (self.Ncells(),self.MAXFACES), np.int32)

            for i in xrange(self.Ncells()):
                cem[i,0:self.nfaces[i]] = self.cell2edges(i)
            self._cell_edge_map = cem
        return self._cell_edge_map

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
        #raise NoSuchEdgeError,str(nodes)
        return -1
        
    def Nedges(self):
        return len(self.edges)
    def Ncells(self):
        return len(self.cells)
    def Npoints(self):
        return len(self.xp)

#def signed_area(x,y,N):
#    i = np.arange(N)
#    ip1 = (i+1)%(N)
#    #return 0.5*(points[i,0]*points[ip1,1] - points[ip1,0]*points[i,1]).sum()
#    return 0.5*(x[...,i]*y[...,ip1] - x[...,ip1]*y[...,i])
    
        
def polyarea(x,y,N):
    """
    Calculate the area of an arbitrary side polygon.
    
    Uses the formula here:
        http://www.seas.upenn.edu/~sys502/extra_materials/Polygon%20Area%20and%20Centroid.pdf
    """
    
    A = np.zeros(x.shape[:-1])
    for i in range(N):
        ip1 = (i+1)%(N)
        A += x[...,i]*y[...,ip1] - x[...,ip1]*y[...,i]
        
    return 0.5*A
    

def centroid(x,y,N):
    """
    **NOT USED **
    Calculate the centroid of an arbitrary side polygon.
    
    Uses the formula here:
        http://www.seas.upenn.edu/~sys502/extra_materials/Polygon%20Area%20and%20Centroid.pdf
    """
    
    A = polyarea(x,y,N)
    
    Cx = np.zeros(x.shape[:-1])
    Cy = np.zeros(x.shape[:-1])
    for i in range(N):
        ip1 = (i+1)%(N)
        tmp = x[...,i]*y[...,ip1] - x[...,ip1]*y[...,i]
        
        Cx += (x[...,i]+x[...,ip1])*tmp
        Cy += (y[...,i]+y[...,ip1])*tmp
        
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
    
class Point:
    def __init__(self,x,y):
        self.x = x
        self.y = y
  
class Line:
    def __init__(self,P1,P2):
        self.x = P2.x - P1.x
        self.y = P2.y - P1.y
  
    def magnitude(self):
        return np.sqrt (self.x**2 + self.y**2)
    
    def unitnormal(self):
        """Finds the units normal vector
        """
        return -self.y/self.magnitude(),self.x/self.magnitude()
        # or return self.y/self.magnitude(),-self.x/self.magnitude()
        
    def dot(self,L2):
        """Dot product with another line
        """
        return self.x*L2.x + self.y*L2.y
    
    def angle(self,L2):
        """Angle with another line
        """
        costheta = self.dot(L2)/(self.magnitude() * L2.magnitude())
        return np.arccos(costheta)
        

def ccw(A,B,C):
	return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

def intersect(A,B,C,D):
    """
    Determines if lines connected by points A-B intersects C-D
    """
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def intersectvec(A,B,C,D):
	return op.and_( op.ne(ccwvec(A,C,D),ccwvec(B,C,D)),op.ne(ccwvec(A,B,C),ccwvec(A,B,D)) )

def ccwvec(A,B,C):
    return op.gt( (C.y-A.y)*(B.x-A.x),(B.y-A.y)*(C.x-A.x) )

def dist(x0,x1,y0,y1):
    return np.sqrt( (x0-x1)**2. + (y0-y1)**2.)

    
    
#P1a = Point(0.,0.)
#P2a = Point(0.,1.)
#P0a = Point(-3.,0.)
#L1a = Line(P1a,P0a)
#L2a = Line(P1a,P2a)
#print L1a.angle(L2a)
