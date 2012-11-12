"""
Script for joining suntans netcdf processor-based output into one file

Created on Fri Oct 19 14:58:33 2012
@author: mrayson
"""

import sunpy
from netCDF4 import Dataset
import getopt, sys, time
import numpy as np


def main(suntanspath,basename,numprocs,nstep=-1):
    
    tic=time.clock()
    
    print '########################################################'
    print '     Started python suntans joining script...'
    print '########################################################'
    
    # Step 1) Read the main grid
    #print 'Loading suntans grid points...'
    grd = sunpy.Grid(suntanspath)
    
    # Step 2) Build a dictionary of the variable names, attributes and dimensions in the first ncfile
    proc=0
    ncfile='%s/%s.%d'%(suntanspath,basename,proc)
    
    
    variables, dims, globalatts = nc_info(ncfile)
    
    # Step 3) Set the dimension sizes based on the original grid
    dims['Nc']=grd.Nc
    dims['Np']=len(grd.xp)
    dims['nt']=None # Sets to unlimited
    
    # Step 4) Initialise the output netcdf file & write the coordinates from the original grid
    
    if nsteps == -1:
        outfile ='%s/%s'%(suntanspath,basename)
    else:
        outfile = '%s/%s_%03d.nc'%(suntanspath,basename[:-3],1)
        
    init_nc(outfile,variables,dims,globalatts)
    
    # Step 5) Write all of the variables that don't have dimension Nc (isCell = False)
    init_ncvars(ncfile,outfile,variables,grd)
    
    # Step 6) Write out the cell based variables based on their rank
    nowritevars = ['xv','yv','cells']
    
    if nsteps==-1:
        nc = Dataset(outfile,'a',format='NETCDF4')
    else:
        makenewfile=True
    
    ncin=[]
    ptr=[]
    for n in range(numprocs):
        ncfile='%s/%s.%d'%(suntanspath,basename,n)
        # Open all of the input files and get some info
        ncin.append(Dataset(ncfile,'r'))
        ncvars=ncin[n].variables
        nt = len(ncvars['time'][:])
        ptr.append(ncvars['mnptr'][:])

        
    for vv in variables:
        vname = vv['Name']
        ctr=1
        if vv['isCell'] and vname not in nowritevars:
            #print '     variable: %s, ndims = %s'%(vv['Name'],vv['ndims'])
            

            if nsteps == -1:
                # Write all time steps to the file
                nc.variables['time'][:]=ncin[0].variables['time'][:]
                if vv['ndims']==1: # Stationary variables
                    sz = (dims[vv['Dimensions'][0]])
                    outvar=np.zeros(sz)
                    for n in range(numprocs):
                        outvar[ptr[n]]=ncin[n].variables[vname][:]
                    nc.variables[vname][:]=outvar
                    
                elif vv['ndims']==2:
                    sz = (nt,dims[vv['Dimensions'][1]])
                    outvar=np.zeros(sz)
                    for n in range(numprocs):
                        outvar[:,ptr[n]]=ncin[n].variables[vname][:,:]
                    nc.variables[vname][:,:]=outvar
                    
                elif vv['ndims']==3:
                    sz = (nt,dims[vv['Dimensions'][1]],dims[vv['Dimensions'][2]])
                    outvar=np.zeros(sz)
                    for n in range(numprocs):
                        outvar[:,:,ptr[n]]=ncin[n].variables[vname][:,:,:]
                    nc.variables[vname][:,:,:]=outvar  
                    
                # Write the buffer to disk
                nc.sync()
            else:
                # Write "nsteps" time steps to the file
                tf=-1
                
                outfile = '%s/%s_%03d.nc'%(suntanspath,basename[:-3],ctr)
                nc = Dataset(outfile,'a',format='NETCDF4_CLASSIC')

                if vv['ndims']==1: # Stationary variables
                    sz = (dims[vv['Dimensions'][0]])
                    outvar=np.zeros(sz)
                    for n in range(numprocs):
                        outvar[ptr[n]]=ncin[n].variables[vname][:]
                    nc.variables[vname][:]=outvar
               
                t1=0
                for t in range(nt):
                    tf+=1
                    
                    if tf == nsteps:
                        tout = range(t1,t)
                        # Write all of the data
                        print '     Writing variable: %s, for t = (%d,%d)'%(vv['Name'],t1,t)
                        nc.variables['time'][:]=ncin[0].variables['time'][tout]
                        if vv['ndims']==2:
                            sz = (t-t1,dims[vv['Dimensions'][1]])
                            outvar=np.zeros(sz)
                            for n in range(numprocs):
                                outvar[:,ptr[n]]=ncin[n].variables[vname][tout,:]
                            nc.variables[vname][:,:]=outvar
                        elif vv['ndims']==3:
                            sz = (t-t1,dims[vv['Dimensions'][1]],dims[vv['Dimensions'][2]])
                            outvar=np.zeros(sz)
                            for n in range(numprocs):
                                outvar[:,:,ptr[n]]=ncin[n].variables[vname][tout,:,:]
                            nc.variables[vname][:,:,:]=outvar
                            
                        
                        t1=t+1
                        
                        # Create a new file
                        tf=0
                        ctr+=1
                        if makenewfile:
                            outfile = '%s/%s_%03d.nc'%(suntanspath,basename[:-3],ctr)
                            nc.close()
                            init_nc(outfile,variables,dims,globalatts)
                            init_ncvars(ncfile,outfile,variables,grd)
                            nc = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
                            
                            if vv['ndims']==1: # Stationary variables
                                sz = (dims[vv['Dimensions'][0]])
                                outvar=np.zeros(sz)
                                for n in range(numprocs):
                                    outvar[ptr[n]]=ncin[n].variables[vname][:]
                                nc.variables[vname][:]=outvar
                        else:
                            outfile = '%s/%s_%03d.nc'%(suntanspath,basename[:-3],ctr)
                            nc.close()
                            nc = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
                            
                            if vv['ndims']==1: # Stationary variables
                                sz = (dims[vv['Dimensions'][0]])
                                outvar=np.zeros(sz)
                                for n in range(numprocs):
                                    outvar[ptr[n]]=ncin[n].variables[vname][:]
                                nc.variables[vname][:]=outvar
                            
                    elif t == nt-1: # Last time step
                        tout = range(t1,t)
                        # Just write all of the data
                        print '     Writing variable: %s, for t = (%d,%d) '%\
                            (vv['Name'],t1,t)
                        nc.variables['time'][:]=ncin[0].variables['time'][tout]
                        if vv['ndims']==2:
                            sz = (t-t1,dims[vv['Dimensions'][1]])
                            outvar=np.zeros(sz)
                            for n in range(numprocs):
                                outvar[:,ptr[n]]=ncin[n].variables[vname][tout,:]
                            nc.variables[vname][:,:]=outvar
                        elif vv['ndims']==3:
                            sz = (t-t1,dims[vv['Dimensions'][1]],dims[vv['Dimensions'][2]])
                            outvar=np.zeros(sz)
                            for n in range(numprocs):
                                outvar[:,:,ptr[n]]=ncin[n].variables[vname][tout,:,:]
                            nc.variables[vname][:,:,:]=outvar
                    
                        
                # Flag to stop generating new files after the first variable    
                makenewfile=False
                nc.close()
                                 
    # Close the file for processor n 
    for n in range(numprocs):
        ncin[n].close()
    
    if nsteps==-1:
        nc.close()
     
    toc=time.clock()
    print 'Elapsed time %10.3f seconds.'%(toc-tic)
    print '########################################################'
    print '     Finished joining netcdf files'
    print '########################################################'
    
def nc_info(ncfile):
    """
    Returns the metadata of all variables, attribute and dimensions
    """
    nc = Dataset(ncfile, 'r')
    
    # Create a dictionary with the dimensions and their values
    dims = {}
    
    for dd in nc.dimensions.keys():
        dims.update({dd:nc.dimensions[dd].__len__()})
    
    # Create a list of dictionaries with the variable info
    variables=[]
    for vv in nc.variables.keys():
        # Attributes
        attdata = nc.variables[vv].__dict__
        
        dimdata = nc.variables[vv].dimensions
        hasNc = ('Nc' in dimdata)
        
        
        variables.append({'Name':vv,'dtype':nc.variables[vv].dtype.type,\
        'ndims':nc.variables[vv].ndim,'Attributes':attdata,'Dimensions':dimdata,\
        'isCell':hasNc})
    
    # Get the global attributes
    globalatts={}
    for aa in nc.ncattrs():
        globalatts.update({aa:getattr(nc,aa)})
    
    nc.close()
    
    return variables, dims, globalatts
    
def init_nc(outfile,variables,dims,globalatts):
    """
    Initialises the output netcdf file for writing
    """
    print "Generating file: %s..."%outfile    
    nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC') 
    
    # Write the global attributes
    for gg in globalatts.keys():
        nc.setncattr(gg,globalatts[gg])
            
    # Create the dimensions
    for dd in dims:
        nc.createDimension(dd,dims[dd])
    
    # Create the variables
    for vv in variables:
        tmpvar=nc.createVariable(vv['Name'],vv['dtype'],vv['Dimensions'],zlib=True,complevel=1,fill_value=99999)
    
        # Create the attributes
        for aa in vv['Attributes'].keys():
            tmpvar.setncattr(aa,vv['Attributes'][aa]) 
                   
    nc.close()    

def init_ncvars(infile,outfile,variables,grd):
    """
    Writes the grid based variables to the output file
    
    Only non-cell based variables
    """
    ncin = Dataset(infile,'r')
    nc = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
    
    #Write the coordinates in the ascii grid files
    nc.variables['xv'][:]=grd.xv
    nc.variables['yv'][:]=grd.yv
    nc.variables['cells'][:]=grd.cells
    nc.variables['yp'][:]=grd.yp
    nc.variables['xp'][:]=grd.xp
    
    nowritevars = ['time','xp','yp','xv','yv','cells']
    for vv in variables:
        vname = vv['Name']
        if not vv['isCell'] and not vv['Name'] in nowritevars:
            nc.variables[vname][:]=ncin.variables[vname][:]
    
    nc.close()
    ncin.close() 

def usage():
    print "--------------------------------------------------------------"
    print "sunjoin.py   -h                 # show this help message      "
    print "         -f suntans.nc        # SUNTANS output netcdf file  "
    print "         -p pathname          # Path to SUNTANS output folder      "
    print "         -n  N                # Number of processors"
    print "         -t  N                # Number of time steps to output (-1 all steps in one file)"
    print "\n\n Example Usage:"
    print "-----------"
    print " python sunjoin.py -f suntans.nc -p ./rundata -n 16 -t 48"
    print ""
    
    
if __name__ == '__main__':
    """ 
        Command line call to join function
    """
    nsteps = -1
    numprocs = 2
    
    try:
            opts,rest = getopt.getopt(sys.argv[1:],'hf:p:n:t:')
    except getopt.GetoptError,e:
        print e
        print "-"*80
        usage()
        exit(1)

    for opt,val in opts:
        if opt == '-h':
            usage()
            exit(1)
        elif opt == '-f':
            basename=val
        elif opt == '-p':
            suntanspath=val
        elif opt == '-t':
            nsteps=int(val)
        elif opt == '-n':
            numprocs=int(val)
            
#    nsteps = -1
#    numprocs = 2
#    suntanspath='C:/Projects/GOMGalveston/MODELLING/GalvestonSquare/rundata'
#    basename = 'suntan_output.nc'      
            
    main(suntanspath,basename,numprocs,nsteps)
