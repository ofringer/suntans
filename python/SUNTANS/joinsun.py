#!/usr/bin/python

"""
Script for joining suntans netcdf processor-based output into one file

Object oriented version created December 2013

@author: mrayson
"""

from sunpy import Grid
from netCDF4 import Dataset
import getopt, sys, time
import numpy as np
from multiprocessing import Pool
import pdb

# Globals
gridvars=['suntans_mesh','cells','face','edges','neigh','grad','nfaces','mnptr','eptr','xv','yv','xp','yp','xe','ye','normal','n1','n2','df','dg','def','Ac','dz','z_r','z_w','Nk','Nke','dv','time']	

nowritevars = ['face','neigh','grad'] # These need special attention

class JoinSuntans(Grid):
    """
    Class for joining suntans NetCDF file
    """
    def __init__(self,suntanspath,basename,numprocs,outvars=None):
        tic=time.clock()
        
        print '########################################################'
        print '     Initializing python SUNTANS NetCDF joining script...'
        print '########################################################'
        
        # Step 1) Read the main grid
        #print 'Loading suntans grid points...'
        Grid.__init__(self,suntanspath)

        self.numprocs = numprocs
        self.suntanspath=suntanspath
        self.basename=basename
    
        # Step 2) Build a dictionary of the variable names, attributes and dimensions in the first ncfile
        proc=0
        self.ncfile='%s/%s.%d'%(suntanspath,basename,proc)
        
        self.variables, self.dims, self.globalatts = nc_info(self.ncfile)

        # If outvars have been set only write those variables and the grid variables
        if not outvars==None:
            self.outvars=gridvars+outvars
            self.outvars=outvars
            newvariables = [vv for vv in self.variables if vv['Name'] in self.outvars]
            self.variables = newvariables

        # Step 3) Set the dimension sizes based on the original grid
        self.dims['Nc']=self.Nc
        self.dims['Np']=self.xp.shape[0]
        self.dims['Ne']=self.edges.shape[0]
        self.dims['time']=None # Sets to unlimited

        # Open each input file and load the cell and edge pointers that go from
        # the "localgrid" to the "maingrid"
        self.ncin=[]
        self.cptr=[] # cell pointer
        self.eptr=[] # edge pointer
        for n in range(numprocs):
            ncfile='%s/%s.%d'%(suntanspath,basename,n)
            # Open all of the input files and get some info
            self.ncin.append(Dataset(ncfile,'r'))
            ncvars= self.ncin[n].variables
            self.nt = len(ncvars['time'][:])
            self.cptr.append(ncvars['mnptr'][:])
            self.eptr.append(ncvars['eptr'][:])

    def __call__(self,nstep=-1,numthreads=4):
        """
        Call to run the joining class
        """
        tic=time.clock()
        # Work out the time steps that need to go into each file
        self.get_filetsteps(nstep) 

        # Initialize the output files
        self.nc=[]
        for outfile in self.outfiles:
            self.nc.append(self.init_outfile(outfile))
            

        # Write the non-time varying variables  
        for outfile,nc in zip(self.outfiles,self.nc):
            print 'Writing non-time varying variables to file:\n\t%s...'%outfile
            self.write_var_notime(nc)

        # Write the other variables
        for outfile,nc,t1,t2 in zip(self.outfiles,self.nc,self.t1,self.t2):
            print 'Writing time varying variables to file:\n\t%s...'%outfile
            self.write_var(nc,t1,t2)

        # Multi-core attempt...
        #   ...doesn't work.
        #numpools = min(numthreads,len(self.outfiles))
        #if numpools==1:
        #    map(self.write_var,self.nc,self.t1,self.t2)
        #else:
        #    print 'Joining files with %d threads...'%numpools
        #    p = Pool(numpools)
        #    p.map(self.write_var,(self.nc,self.t1,self.t2))

        # Close all of the open files
        #self.close_all()
        toc=time.clock()
        print 'Elapsed time %10.3f seconds.'%(toc-tic)
        print '########################################################'
        print '     Finished joining netcdf files'
        print '########################################################'


    def write_var(self,nc,t1,t2):
        """
        Write the time-varying variables
        """
        tout = range(t1,t2)
        nc.variables['time'][:]=self.ncin[0].variables['time'][tout]

        for vv in self.variables:
            vname = vv['Name']
            isCellEdge=False
            if vv['isCell']:
                isCellEdge=True
            elif vv['isEdge']:
                isCellEdge=True            

            tout = range(t1,t2)

            if vv['ndims']==2 and vv['isTime'] and isCellEdge and vname not in nowritevars:
                print '\t\t%s - steps %d to %d of %d'%(vname,t1,t2,self.nt)
                sz = (t2-t1,self.dims[vv['Dimensions'][1]])
                outvar=np.zeros(sz)
                for n in range(self.numprocs):
                    if vv['isCell']:
                        outvar[:,self.cptr[n]]=self.ncin[n].variables[vname][tout,:]
                    elif vv['isEdge']:
                        outvar[:,self.eptr[n]]=self.ncin[n].variables[vname][tout,:]
                nc.variables[vname][:,:]=outvar

            elif vv['ndims']==3 and vv['isTime'] and isCellEdge and vname not in nowritevars:
                print '\t\t%s - steps %d to %d of %d'%(vname,t1,t2,self.nt)
                sz = (t2-t1,self.dims[vv['Dimensions'][1]],self.dims[vv['Dimensions'][2]])
                outvar=np.zeros(sz)
                for n in range(self.numprocs):
                    if vv['isCell']:
                        outvar[:,:,self.cptr[n]]=self.ncin[n].variables[vname][tout,:,:]
                    elif vv['isEdge']:
                        outvar[:,:,self.eptr[n]]=self.ncin[n].variables[vname][tout,:,:]
                            
                nc.variables[vname][:,:,:]=outvar
            # Write the buffer to disk
            nc.sync()


    def write_var_notime(self,nc):
        """
        Write all of the non-time varying files to each file
        """
        for vv in self.variables:
            vname = vv['Name']
            isCellEdge=False
            if vv['isCell']:
                isCellEdge=True
            elif vv['isEdge']:
                isCellEdge=True            

            if vv['ndims']==1 and isCellEdge and vname not in nowritevars: 
                sz = (self.dims[vv['Dimensions'][0]])
                outvar=np.zeros(sz)
                for n in range(self.numprocs):
                    if vv['isCell']:
                        outvar[self.cptr[n]]=self.ncin[n].variables[vname][:]
                    elif vv['isEdge']:
                        outvar[self.eptr[n]]=self.ncin[n].variables[vname][:]
                nc.variables[vname][:]=outvar

            elif vv['ndims']==2 and not vv['isTime'] and isCellEdge and vname not in nowritevars:# Stationary 2D variables 
                sz = (self.dims[vv['Dimensions'][0]],self.dims[vv['Dimensions'][1]])
                outvar=np.zeros(sz)
                for n in range(self.numprocs):
                    if vv['isCell']:
                        outvar[self.cptr[n],:]=self.ncin[n].variables[vname][:,:]
                    elif vv['isEdge']:
                        outvar[self.eptr[n],:]=self.ncin[n].variables[vname][:,:]
                nc.variables[vname][:,:]=outvar  

            elif vname in nowritevars:
                #print '     variable: %s, ndims = %s'%(vv['Name'],vv['ndims'])
                sz = (self.dims[vv['Dimensions'][0]],self.dims[vv['Dimensions'][1]])
                outvar=np.zeros(sz)
                for n in range(self.numprocs):
                    if vname == 'grad': # edge_face_connectivity
                        outvar[self.eptr[n],:] = self.cptr[n][self.ncin[n].variables[vname][:,:]]
                    elif vname == 'face': # face_edge_connectivity
                        # Doesn't work for arbitrary polygons
                        #outvar[cptr[n],:] = eptr[n][ncin[n].variables[vname][:,:]]
                        # Need to do it cell by cell
                        for ii in range(self.cptr[n].shape[0]):
                            nf = self.nfaces[self.cptr[n][ii]]
                            outvar[self.cptr[n][ii],0:nf] = self.eptr[n][self.ncin[n].variables[vname][ii,0:nf]]
                    elif vname == 'neigh': # face_face_connectivity
                        #outvar[cptr[n],:] = cptr[n][ncin[n].variables[vname][:,:]]
                        for ii in range(self.cptr[n].shape[0]):
                            nf = self.nfaces[self.cptr[n][ii]]
                            outvar[self.cptr[n][ii],0:nf] = self.cptr[n][self.ncin[n].variables[vname][ii,0:nf]]

                nc.variables[vname][:,:]=outvar 

            # Write the buffer to disk
            nc.sync()

    def get_filetsteps(self,nstep):
        """
        Finds the time steps that need to go into each output file
        """
        if nstep==-1:
            nstep = self.nt

        nfiles = int(np.ceil(self.nt/nstep))

        self.outfiles=[]
        self.t1=[] #Start time index
        self.t2=[] #End time index
        t1=0
        for n in range(nfiles):
            self.t1.append(t1)
            self.t2.append(min(self.nt,t1+nstep))
            t1+=nstep

            #Initialise the file names
            self.outfiles.append('%s/%s_%03d.nc'%(self.suntanspath,self.basename[:-3],n+1))
        

    def close_all(self):
        """Closes all open netcdf files"""

        for nc in self.nc:
            nc.close()

        for ncin in self.ncin:
            ncin.close()

    def init_outfile(self,outfile):
        """
        Create an output netcdf files

        Returns a pointer to the output netcdf file object
        """

        # Initialise the output netcdf file & write the coordinates from the original grid
        init_nc(outfile,self.variables,self.dims,self.globalatts)
        
        # Write all of the variables that don't have dimension Nc (isCell = False)
        self.init_ncvars(self.ncfile,outfile,self.variables)
        
        # Open the current output file
        return Dataset(outfile,'a',format='NETCDF4')


    def init_ncvars(self,infile,outfile,variables):
        """
        Writes the grid based variables to the output file
        
        Only non-cell based variables
        """
        ncin = Dataset(infile,'r')
        nc = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        
        #Write the coordinates in the ascii grid files
        nc.variables['xv'][:]=self.xv
        nc.variables['yv'][:]=self.yv
        nc.variables['cells'][:]=self.cells
        nc.variables['yp'][:]=self.yp
        nc.variables['xp'][:]=self.xp
        
        # Initialise all of the other variables
        nowritevars = ['time','xp','yp','xv','yv','cells']
        for vv in variables:
            vname = vv['Name']
            if not vv['isCell'] and not vv['isEdge'] and not vv['Name'] in nowritevars:
                nc.variables[vname][:]=ncin.variables[vname][:]
        
        nc.close()
        ncin.close() 


def nc_info(ncfile):
    """
    Returns the metadata of all variables, attribute and dimensions
    """
    print 'Reading: %s...'%ncfile
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
        hasNe = ('Ne' in dimdata)
        hasTime = ('time' in dimdata)
        
        if '_FillValue' in  nc.variables[vv].ncattrs():
            isFilled = True
        else:
            isFilled = False
        #try:
        #    isFilled = nc.variables[vv]._DeflateLevel > 0
        #except:
        #    isFilled = False
        
        
        variables.append({'Name':vv,'dtype':nc.variables[vv].dtype.type,\
        'ndims':nc.variables[vv].ndim,'Attributes':attdata,'Dimensions':dimdata,\
        'isCell':hasNc,'isEdge':hasNe,'isTime':hasTime,'isFilled':isFilled})
    
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
        if vv['isFilled']:
            tmpvar=nc.createVariable(vv['Name'],vv['dtype'],vv['Dimensions'],zlib=True,complevel=2,fill_value=99999.0)
        else:
            tmpvar=nc.createVariable(vv['Name'],vv['dtype'],vv['Dimensions'])
    
        # Create the attributes
        for aa in vv['Attributes'].keys():
            tmpvar.setncattr(aa,vv['Attributes'][aa]) 
                   
    nc.close()    

def usage():
    print "--------------------------------------------------------------"
    print "sunjoin.py   -h                 # show this help message      "
    print "         -f suntans.nc        # SUNTANS output netcdf file  "
    print "         -p pathname          # Path to SUNTANS output folder      "
    print "         -n  N                # Number of processors"
    print "         -t  N                # Number of time steps to output (-1 all steps in one file)"
    print " 	    -v  'var1 var2 ...'  # List of variales to write (default: all)"
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
    outvars=None
    
    try:
        opts,rest = getopt.getopt(sys.argv[1:],'hf:p:n:t:v:')
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
	elif opt == '-v':
	     outvars=val.split(' ')

    sun = JoinSuntans(suntanspath,basename,numprocs)
    sun(nstep=nsteps)     

#	# Testing only	
#    #nsteps = 4
#    #numprocs = 64
#    #suntanspath='C:\\Projects\\GOMGalveston\\MODELLING\\GalvestonCoarse\\rundata'
#    #basename = 'GalvCoarse_ROMS_Coupling.nc'      
#            
#    #nsteps = 4
#    #numprocs = 64
#    #suntanspath='C:\\Projects\\GOMGalveston\\MODELLING\\GalvestonCoarse\\rundata'
#    #basename = 'GalvCoarse_ROMS_Coupling.nc'      
#            
#    main(suntanspath,basename,numprocs,nsteps,outvars=outvars)
#
# Testing
#nstep = -1
#numprocs = 2
#suntanspath='/home/mrayson/suntans-merge/quad-netcdf/examples_update/riverplume/data'
#basename = 'Plume.nc'      
#
#sun = JoinSuntans(suntanspath,basename,numprocs)
#sun(nstep=nstep)
