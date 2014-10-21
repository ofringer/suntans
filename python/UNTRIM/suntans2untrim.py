"""
Tool for converting from a suntans to untrim netcdf formats
"""
from sunpy import Spatial, Grid
import othertime
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
from untrim_tools import untrim_ugrid as ugrid 

import pdb

#Dictionary containing the suntans-untrim equivalent grid variables
untrim_gridvars = {\
        'xp':'Mesh2_node_x',\
        'yp':'Mesh2_node_y',\
        'xv':'Mesh2_face_x',\
        'yv':'Mesh2_face_y',\
        'xe':'Mesh2_edge_x',\
        'ye':'Mesh2_edge_y',\
        'mark':'Mesh2_edge_bc',\
        'edges':'Mesh2_edge_nodes',\
        'grad':'Mesh2_edge_faces',\
        'cells':'Mesh2_face_nodes',\
        'face':'Mesh2_face_edges',\
        'dv':'Mesh2_face_depth',\
        'de':'Mesh2_edge_depth',\
        'z_r':'Mesh2_layer_3d',\
        'time':'Mesh2_data_time',\
        'mark':'Mesh2_edge_bc',\
        'facemark':'Mesh2_face_bc'\
    }

# Dictionary containing the suntans-untrim  equivalent grid dimensions
untrim_griddims = {\
        'Np':'nMesh2_node',\
        'Ne':'nMesh2_edge',\
        'Nc':'nMesh2_face',\
        'Nkmax':'nMesh2_layer_3d',\
        'numsides':'nMaxMesh2_face_nodes',\
        'time':'nMesh2_data_time'\
    }

# Dimensions with hard-wired values
other_dims = {\
        'Three':3,\
        'Two':2,\
        'nsp':1,\
        'date_string_length':19,\
        'nMesh2_time':1\
        }

# physical variables that are directly transferable
untrim_vars = {\
        'salt':'Mesh2_salinity_3d',\
        'nu_v':'Mesh2_vertical_diffusivity',\
        'eta':'Mesh2_sea_surface_elevation'
        }

varnames = ['Mesh2_salinity_3d',\
            'Mesh2_vertical_diffusivity_3d',\
            'Mesh2_sea_surface_elevation',\
            'h_flow_avg',\
            'v_flow_avg',\
            'Mesh2_edge_wet_area',\
            'Mesh2_face_wet_area',\
            'Mesh2_edge_bottom_layer',\
            'Mesh2_edge_top_layer',\
            'Mesh2_face_bottom_layer',\
            'Mesh2_face_top_layer',\
            'Mesh2_face_water_volume',\
            ]

FILLVALUE=-9999

def suntans2untrim(ncfile,outfile,tstart,tend,grdfile=None):
    """
    Converts a suntans averages netcdf file into untrim format
    for use in particle tracking
    """
    ####
    # Step 1: Load the suntans data object
    ####
    sun = Spatial(ncfile,klayer=[-99])

    # Calculate some other variables
    sun.de = sun.get_edgevar(sun.dv,method='min')
    sun.mark[sun.mark==5]=0
    sun.mark[sun.mark==3]=2
    sun.facemark = np.zeros((sun.Nc,),dtype=np.int)

    # Update the grad variable from the ascii grid file if supplied
    if not grdfile == None:
        print 'Updating grid with ascii values...'
        grd = Grid(grdfile)
        sun.grad = grd.grad[:,::-1]

    ###
    # Step 2: Write the grid variables to a netcdf file
    ###
    nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

    # Global variable
    nc.Description = 'UnTRIM history file converted from SUNTANS output'

    # Write the dimensions
    for dd in untrim_griddims.keys():
        if dd == 'time':
            nc.createDimension(untrim_griddims[dd],0)
        elif dd =='numsides':
            nc.createDimension(untrim_griddims[dd],sun.maxfaces)
        else:
            nc.createDimension(untrim_griddims[dd],sun[dd])


    for dd in other_dims:
        nc.createDimension(dd,other_dims[dd])

    ###
    # Step 3: Initialize all of the grid variables
    ###
    def create_nc_var(name, dimensions, attdict,data=None, \
        dtype='f8',zlib=False,complevel=0,fill_value=999999.0):
            
        tmp=nc.createVariable(name, dtype, dimensions,\
            zlib=zlib,complevel=complevel,fill_value=fill_value)

        for aa in attdict.keys():
            tmp.setncattr(aa,attdict[aa])
        
        if not data==None:
            nc.variables[name][:] = data
     
    # Make sure the masked cells have a value of -1
    mask = sun['cells'].mask.copy()
    sun['cells'][mask]=FILLVALUE
    sun['face'][mask]=FILLVALUE

    for vv in untrim_gridvars.keys():
        vname = untrim_gridvars[vv]
        print 'Writing grid variable %s (%s)...'%(vname,vv)

        if vv=='time':
            continue

        # add dz_min attribute to z_r variable
        if vv == 'z_r':
            ugrid[vname]['attributes'].update({'dz_min':1e-5})
            #sun[vv][:]=sun[vv][::-1]
            sun[vv][:]=sun['z_w'][0:-1][::-1]

        # Reverse the order of grad(???)
        if vv=='grad':
            sun[vv][:]=sun[vv][:,::-1]

        ## Fix one-based indexing
        #if vv in ['cells','edges','grad']:
        #    mask = sun[vv][:]==-1
        #    tmp = sun[vv][:]+1
        #    tmp[mask]=-1
        #    #sun[vv][:]=sun[vv][:]+1
        #    create_nc_var(vname,ugrid[vname]['dimensions'],ugrid[vname]['attributes'],\
        #        data=tmp,dtype=ugrid[vname]['dtype'])

        create_nc_var(vname,ugrid[vname]['dimensions'],ugrid[vname]['attributes'],\
            data=sun[vv],dtype=ugrid[vname]['dtype'])

            
    # Initialize the two time variables
    vname=untrim_gridvars['time']
    create_nc_var(vname,ugrid[vname]['dimensions'],ugrid[vname]['attributes'],\
            dtype=ugrid[vname]['dtype'])
    vname = 'Mesh2_data_time_string'
    create_nc_var(vname,ugrid[vname]['dimensions'],ugrid[vname]['attributes'],\
            dtype=ugrid[vname]['dtype'])

    ###
    # Step 4: Initialize all of the time-varying variables (but don't write) 
    ###
    for vname  in varnames:
        print 'Creating variable %s...'%(vname)

        create_nc_var(vname,ugrid[vname]['dimensions'],ugrid[vname]['attributes'],\
            dtype=ugrid[vname]['dtype'],zlib=True,complevel=1,fill_value=999999.)

    ###
    # Step 5: Loop through all of the time steps and write the variables
    ###
    tsteps = sun.getTstep(tstart,tend)
    tdays = othertime.DaysSince(sun.time,basetime=datetime(1899,12,31))
    for ii, tt in enumerate(tsteps):
        # Convert the time to the untrim formats
        timestr = datetime.strftime(sun.time[tt],'%Y-%m-%d %H:%M:%S')

        print 'Writing data at time %s (%d of %d)...'%(timestr,tt,tsteps[-1]) 

        #Write the time variables
        nc.variables['Mesh2_data_time'][ii]=tdays[ii]
        nc.variables['Mesh2_data_time_string'][:,ii]=timestr

        # Load each variable or calculate it and convert it to the untrim format
        sun.tstep=[tt]

        ###
        # Compute a few terms first
        eta = sun.loadData(variable='eta' )  
        U = sun.loadData(variable='U_F' )  
        dzz = sun.getdzz(eta)
        dzf = sun.getdzf(eta)

        
        vname='Mesh2_sea_surface_elevation'
        #print '\tVariable: %s...'%vname
        nc.variables[vname][:,ii]=eta

        vname = 'Mesh2_salinity_3d'
        #print '\tVariable: %s...'%vname
        tmp3d = sun.loadData(variable='salt' )  
        nc.variables[vname][:,:,ii]=tmp3d.swapaxes(0,1)[:,::-1]

        vname = 'Mesh2_vertical_diffusivity_3d'
        #print '\tVariable: %s...'%vname
        tmp3d = sun.loadData(variable='nu_v' )  
        nc.variables[vname][:,:,ii]=tmp3d.swapaxes(0,1)[:,::-1]

        vname = 'h_flow_avg'
        #print '\tVariable: %s...'%vname
        nc.variables[vname][:,:,ii]=U.swapaxes(0,1)[:,::-1]


        vname = 'v_flow_avg'
        #print '\tVariable: %s...'%vname
        tmp3d = sun.loadData(variable='w' ) * sun.Ac # m^3/s
        nc.variables[vname][:,:,ii]=tmp3d.swapaxes(0,1)[:,::-1]

        # Need to calculate a few terms for the other variables

        vname = 'Mesh2_edge_wet_area'
        #print '\tVariable: %s...'%vname
        #dzf = sun.loadData(variable='dzf')
        tmp3d = dzf*sun.df
        nc.variables[vname][:,:,ii]=tmp3d.swapaxes(0,1)[:,::-1]

        vname = 'Mesh2_face_water_volume'
        #print '\tVariable: %s...'%vname
        #dzz = sun.loadData(variable='dzz')
        tmp3d = dzz*sun.Ac
        nc.variables[vname][:,:,ii]=tmp3d.swapaxes(0,1)[:,::-1]

        vname = 'Mesh2_face_wet_area'
        #print '\tVariable: %s...'%vname
        tmp3d = np.repeat(sun.Ac[np.newaxis,...],sun.Nkmax,axis=0)
        nc.variables[vname][:,:,ii]=tmp3d.swapaxes(0,1)[:,::-1]

        # UnTRIM references from bottom to top i.e.
        # k = 0 @ bed ; k = Nkmax-1 @ top

        vname = 'Mesh2_edge_bottom_layer'
        #print '\tVariable: %s...'%vname
        #tmp2d = sun.Nkmax-sun.Nke # zero based
        tmp2d = sun.Nkmax-sun.Nke+1 # one based
        nc.variables[vname][:,ii]=tmp2d

        vname = 'Mesh2_edge_top_layer'
        #print '\tVariable: %s...'%vname
        etop = sun.loadData(variable='etop')
        #tmp2d = sun.Nkmax-etop-1 # zero based
        tmp2d = sun.Nkmax-etop # one based
        nc.variables[vname][:,ii]=tmp2d

        vname = 'Mesh2_face_bottom_layer'
        #print '\tVariable: %s...'%vname
        #tmp2d = sun.Nkmax-sun.Nk + 1 # zero based
        tmp2d = sun.Nkmax-sun.Nk # one based
        nc.variables[vname][:,ii]=tmp2d

        vname = 'Mesh2_face_top_layer'
        #print '\tVariable: %s...'%vname
        ctop = sun.loadData(variable='ctop')
        #tmp2d = sun.Nkmax-ctop-1 # zero based
        tmp2d = sun.Nkmax-ctop # one based
        nc.variables[vname][:,ii]=tmp2d

    print 72*'#'
    print '\t Finished SUNTANS->UnTRIM conversion'
    print 72*'#'






    # close the file
    nc.close()


##############
# Testing
##############
                                                                            ##################
if __name__=="__main__":
    # Inputs
    ncfile = '../data/Racetrack_AVG_0000.nc'
    outfile = '../data/Racetrack_untrim.nc'
    tstart = '20020629.0000'
    tend = '20060701.1200'

    suntans2untrim(ncfile,outfile,tstart,tend)
