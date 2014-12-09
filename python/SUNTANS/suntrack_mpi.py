# -*- coding: utf-8 -*-
"""
Parallel implementatation of the particle tracking model

See this example:
    http://jeremybejarano.zzl.org/MPIwithPython/collectiveCom.html

Created on Mon Jan 14 14:52:43 2013

@author: mrayson
"""

from mpi4py import MPI
from suntrack import SunTrack
import numpy as np
from datetime import datetime, timedelta
import othertime

# Start the MPI communicator
comm = MPI.COMM_WORLD

size=comm.size


def runmpi(ncfile,outfile,tstart,tend,dt,dtout,x,y,z,agepoly=None,method='nearest',is3D=False):

    # Generate a list of tuples for the time info
    timevec = othertime.TimeVector(tstart,tend,dtout,timeformat ='%Y%m%d.%H%M%S')
    timevec_sec = othertime.SecondsSince(timevec)
    timeinfos=[]
    for ii in range(timevec.shape[0]-1):
        if ii == 0:
            timestart = datetime.strftime(timevec[ii],'%Y%m%d.%H%M%S')
        else:
            timestart = datetime.strftime(timevec[ii]+timedelta(seconds=dt),'%Y%m%d.%H%M%S')
        
        timeend = datetime.strftime(timevec[ii+1],'%Y%m%d.%H%M%S')
        timeinfos.append( (timestart,timeend,dt) )
    
    # Initialise the particle tracking object    
    print 'Initialising the particle tracking object on processor: %d...'%(comm.rank)
    sun = SunTrack(ncfile,interp_method='mesh',interp_meshmethod=method,is3D=is3D)
    
    # Initialise the age values
    if not agepoly==None:
        calcage=True
    else:
        calcage=False
        age=None
        agemax=None

    # On rank = 0 only
    if comm.rank == 0:
        n = int(x.shape[0])
        
        # Check if the number of processes divides eveny into the array length
        rem = np.mod(n,size)
        if rem == 0:
            #length of each process's portion of the original vector
            local_n = np.array([n/size])
        else:
            print 'Padding array with extra values...'
            nextra = size-rem
            xpad = np.zeros((nextra,))
            x = np.hstack((x,xpad))
            y = np.hstack((y,xpad))
            z = np.hstack((z,xpad))
            
            n = n+nextra    
            local_n = np.array([n/size])
        
        print 'Size of original vector = %d\nSize of split vector = %d'%(n,local_n)

	if calcage:
	    age = np.zeros_like(x)
	    agemax = np.zeros_like(x)
        
        # Initialise the output netcdf file
        sun.initParticleNC(outfile,n,age=calcage)
         
    else:
        x=None
        y=None
        z=None
        local_n = np.array([0])
	if calcage:
	    age=None
	    agemax=None
        
    comm.Barrier()
    t_start = MPI.Wtime()
    # Broadcast the particle tracking object everywhere
    #sun = comm.bcast(sun, root=0) # !! Doesn't work for this object !!
    
    # Scatter the x,y,z locations up amongst all processors
    #communicate local array size to all processes
    comm.Bcast(local_n, root=0)
    
    #initialize the local particle arrays as numpy arrays
    x_local = np.zeros(local_n)
    y_local = np.zeros(local_n)
    z_local = np.zeros(local_n)
    if calcage:
    	age_local = np.zeros(local_n)
    	agemax_local = np.zeros(local_n)

   
    #divide up vectors
    comm.Scatter(x, x_local, root=0)
    comm.Scatter(y, y_local, root=0)
    comm.Scatter(z, z_local, root=0)
    if calcage:
	comm.Scatter(age, age_local, root=0)
	comm.Scatter(agemax, agemax_local, root=0)
    else:
	age_local=None
	agemax_local=None

    ###
    # Testing plot of particle positions
    ###
    #if comm.rank ==0:
    #	import matplotlib.pyplot as plt
    #    plt.figure
    #    plt.plot(x_local,y_local,'.')
    #    plt.show()
    #comm.Barrier()

    ###
    # Write out the initial location to netcdf
    if comm.rank==0:
	sun.writeParticleNC(outfile,x,y,z,timevec_sec[0],0,age=age,agemax=agemax)
    comm.Barrier()

    ###
    # ... Call the particle tracking module on each processor
    for ii,timeinfo in enumerate(timeinfos):
        if comm.rank == 0:
            sun(x_local,y_local,z_local,timeinfo,agepoly=agepoly,age=age_local,agemax=agemax_local,verbose=True)    
        else:
            sun(x_local,y_local,z_local,timeinfo,agepoly=agepoly,age=age_local,agemax=agemax_local,verbose=False) 
            
        # Send the particles back to their main array    
        comm.Barrier()
    
        comm.Gather(sun.particles['X'],x,root=0)
        comm.Gather(sun.particles['Y'],y,root=0)
        comm.Gather(sun.particles['Z'],z,root=0)
	if calcage:
	    comm.Gather(sun.particles['age'],age,root=0)
	    comm.Gather(sun.particles['agemax'],agemax,root=0)
        
        # Write the output to a netcdf file
        if comm.rank==0:
            sun.writeParticleNC(outfile,x,y,z,sun.time_track_sec[-1],ii+1,age=age,agemax=agemax)
        
        comm.Barrier()
    
    t_diff = MPI.Wtime()-t_start ### Stop stopwatch ###
    if comm.rank==0:
        print 78*'='+'\n'+78*'='
        print 'Completed particle tracking using %d cores in %6.2f seconds.'%(comm.size,t_diff)
        print 78*'='+'\n'+78*'='


if __name__ == '__main__':
    """
    Example script
    
    run via:
	    mpirun -np 4 python suntrack_mpi.py
    """
    from suntrack import GridParticles

    ####
    # Input variables
    ncfile = '2010/GalvCoarse_20100107.nc'
    outfile = 'ParticleTest.nc'
    tstart = '20100107.140000'
    tend =  '20100107.200000'
    dt = 120.0
    dtout = 3600.0
    
    # initialse the particle locations
    dx = 250.0
    dy = 250.0
    nz = 5
    x,y,z = GridParticles(ncfile,dx,dy,nz)

    runmpi(ncfile,outfile,tstart,tend,dt,dtout,x,y,z)
