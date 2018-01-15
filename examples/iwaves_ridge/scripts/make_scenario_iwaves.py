# -*- coding: utf-8 -*-
"""
Create a one cell quad grid

Created on Mon Nov 18 10:23:26 2013

@author: mrayson
"""

import os
from shutil import copyfile
import numpy as np
import matplotlib.pyplot as plt
import operator

from soda.dataio.ugrid.ugridgen import cartesian_ugrid_gen
from soda.dataio.suntans.sunboundary import modifyBCmarker, Boundary, InitialCond
from soda.dataio.suntans.sunpy import Grid

PI = np.pi
GRAV=9.81

def make_suntans(suntanspath):
    ####################################################
    # Inputs
    beta = 7e-4
    N = 0.01 # buoyancy freqncy
    H = 1000.0
    wave_period=12*3600.

    phi0 = 0.0      # Internal wave amplitude
    U0 = 0.1        # Steady flow rate
    dudz = 0e-4      # Vertical Shear

    # Size of domain
    ny = 1
    nx = 100
    nz = 40

    htopo = 50.

    #suntanspath = 'data'

    starttime = '20000101.000000'
    endtime = '20000130.000000'
    dt = 3600.

    icfile = 'IWave_IC.nc'
    bcfile = 'IWave_BC.nc'
    ####################################################
    def gaussian(x,x0,a,b,pow=2.):
	return b*np.exp(-(x-x0)**pow/(2*a**pow))


    #######
    # Compute wave parameters to get the domain size etc
    k_z = PI / H        
    omega = 2*PI/wave_period

    k_x = k_z*omega / N

    lambda_x = 2*PI / k_x

    # Grid size
    L = 2*lambda_x
    dx = L/nx

    W = dx

    #####
    # Create the grid
    if not os.path.isdir(suntanspath):
	print 'Creating new directory: %s'%suntanspath
	os.mkdir(suntanspath)
	copyfile('rundata/suntans.dat','%s/suntans.dat'%suntanspath)

    xlims = [0,L]
    ylims = [0,W]


    # Create the grid
    grd= cartesian_ugrid_gen(xlims, ylims, dx, suntanspath=suntanspath)

    # Load the grid
    grd = Grid(suntanspath)

    #grd.dv = H*np.ones_like(grd.xv)
    grd.dv = H - gaussian(grd.xv, L/2, L/16, htopo)
    grd.saveBathy('%s/depth.dat-voro'%suntanspath)

    grd.dz = grd.calcVertSpace(nz, 1.0, H)
    grd.saveVertspace('%s/vertspace.dat'%suntanspath)
    grd.setDepth(grd.dz)

    # Salinity field
    dsdz = N**2 / (GRAV * beta)
    salt = dsdz*grd.z_r

    # Create the boundary conditions


    ##########
    # Modify the boundary markers and create the boundary condition file
    ##########
    # This changes the edge types based on the polygons in the shapefile
    #modifyBCmarker(suntanspath,bcpolygonfile)

    # Modify the left and right edges and convert to type 2
    #hgrd = grd.convert2hybrid()

    grd.mark[grd.mark>0]=1 # reset all edges to type-1

    ## convert edges +/- half a grid cell from the edge
    #dx = hgrd.dg.max()
    xmin = grd.xv.min()-dx/4.0
    xmax = grd.xv.max()+dx/4.0

    grd.calc_edgecoord()

    indleft = operator.and_(grd.mark==1, grd.xe < xmin) # all boundaries
    indright = operator.and_(grd.mark==1, grd.xe > xmax) # all boundaries

    ## Free-surface boundaries
    ##grd.mark[indleft]=3
    ##grd.mark[indright]=3
    #
    ## River boundaries
    grd.mark[indleft]=2
    grd.mark[indright]=2
    ##grd.edge_id[indleft]=1
    ##grd.edge_id[indright]=2
    #
    edgefile = suntanspath+'/edges.dat'
    grd.saveEdges(edgefile)
    #print 'Updated markers written to: %s'%(edgefile)
    #grd.write2suntans(suntanspath)

    #hgrd.dv=grd.dv
    #hgrd.z_r = grd.z_r
    #grd = hgrd

    #Load the boundary object from the grid
    #   Note that this zeros all of the boundary arrays
    bnd = Boundary(suntanspath,(starttime,endtime,dt))

    bnd.setDepth(grd.dv)
    #
    ## Try out a linear interpolation 
    ##y0 = bnd.ye.min()
    ##y1 = bnd.ye.max()
    ##du = 0.0
    #       #dy = y1-y0
    ##Utst = (bnd.ye.ravel()-y0)/dy*du
    #
    ## Normal components
    #nx = hgrd.n1[grd.mark==2]
    #ny = hgrd.n2[grd.mark==2]
    #
    t = bnd.tsec-bnd.tsec[0]
    #
    ##indleft = grd.xv[bnd.cellp] < x0
    ##indright = grd.xv[bnd.cellp]>x0
    #if bnd.N3>0:
    #    indleft = bnd.xv < x0
    #    indright = bnd.xv>x0
    #
    ## Height boundary
    #for ii in range(bnd.N3):
    #    hleft = h0 * np.sin(omega*t)
    #    hright = h0 * np.sin(omega*t) # small lag
    #    #hright = 0
    #    if indleft[ii]:
    #        bnd.h[:,ii] = hleft
    #    else:
    #        bnd.h[:,ii] = hright
    #
    #if bnd.N2>0:
    #    indleft = bnd.xe < x0
    #    indright = bnd.xe>x0
    #
    #
    #
    # Velocity boundary
    for k in range(bnd.Nk):
       for ii, eptr in enumerate(bnd.edgep.tolist()):
	   # Steady shear flow
	   usteady = U0 + dudz * (H - grd.z_r[k])
	   if indleft[eptr]:
	       #bnd.boundary_u[:,k,ii] = phi0*np.cos(k_z*grd.z_r[k])*np.sin(omega*t)
	       bnd.boundary_u[:,k,ii] = phi0*np.sin(omega*t) + usteady
	       bnd.boundary_S[:,k,ii] = salt[k]# - phi0*np.sin(k_z*grd.z_r[k])*np.sin(omega*t)
	   elif indright[eptr]:
	       bnd.boundary_u[:,k,ii] = phi0*np.sin(omega*t) + usteady
	       bnd.boundary_S[:,k,ii] = salt[k]
	       #bnd.boundary_u[:,k,ii] = u*nx[ii]
	       #bnd.boundary_v[:,k,ii] = u*ny[ii]
	       #bnd.boundary_h[:,ii] = h
    #
    #
    ## River flow boundary
    #for k in range(bnd.Nk):
    #   for ii in range(bnd.Nseg):
    #
    #       Q = hmax*W*u0 * np.sin(omega*t) # Flow rate
    #       sfac = 1.
    #       if ii == 1:
    #           sfac = -1.
    #       
    #       bnd.boundary_Q[:,ii] = sfac*Q

    # type-3 
    #h = calc_fs(bnd.xv)

    #bnd.h[:] = np.ones_like(bnd.h) * h



    # Write the boundary file
    bnd.write2NC(suntanspath+'/'+bcfile)

    #########
    # Create the initial conditions file
    #########
    IC = InitialCond(suntanspath,starttime)

    #seiche=Seiche(L,W,H,A)
    #u,IC.h[:] = seiche(IC.xv,0)
    IC.h[:] = 0

    IC.T[:,0:nz//2, :] = 1.0 # set T = 1 in the upper water column
    IC.S[:,:,:] = salt[np.newaxis,:,np.newaxis]

    # Set T as a circle
    #x0 = L/2
    #y0 = W/2
    #R = W/10. # Radius
    #
    #dist = np.sqrt( (IC.xv-x0)**2 + (IC.yv-y0)**2.)
    #
    #IC.h[...,dist<=R] = Amid



    # Write the initial condition file
    IC.writeNC(suntanspath+'/'+icfile,dv=grd.dv)



    #grd.plotmesh()
    #plt.show()

if __name__=='__main__':
    import sys
    sunpath = sys.argv[1]

    make_suntans(sunpath)

