"""
Tool for handling UnTRIM data specific to the SUNTANS model

M.Rayson
Stanford University
March 2014
"""

from sunpy import Grid, Spatial
import numpy as np
import matplotlib.pyplot as plt

import pdb

# Dictionary containing the suntans-untrim equivalent grid variables
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
    'faces':'Mesh2_face_edges',\
    'dv':'Mesh2_face_depth',\
    'z_r':'Mesh2_layer_3d',\
    'time':'Mesh2_data_time'\
    }

# Dictionary containing the suntans-untrim equivalent grid dimensions
untrim_griddims = {\
    'Np':'nMesh2_node',\
    'Ne':'nMesh2_edge',\
    'Nc':'nMesh2_face',\
    'Nkmax':'nMesh2_layer_3d',\
    'Nk':'nMesh2_layer_3d',\
    'numsides':'nMaxMesh2_face_nodes',\
    'Two':'Two'\
    }

######
# The untrim netcdf data format
######
untrim_ugrid={}
fillval=999999.0

vname = 'Mesh2'
dimensions = ()
attributes = {
        'dimension': 2,\
		'cf_role': "mesh_topology" ,\
		'long_name': "Topology data of 2D untrim mesh" ,\
		'node_coordinatesi': "Mesh2_node_x Mesh2_node_y",\
		'edge_coordinates': "Mesh2_edge_x Mesh2_edge_y",\
		'face_coordinates': "Mesh2_face_x Mesh2_face_y",\
		'edge_node_connectivity': "Mesh2_edge_nodes",\
		'edge_face_connectivity': "Mesh2_edge_faces",\
		'face_node_connectivity': "Mesh2_face_nodes",\
		'face_edge_connectivity': "Mesh2_face_edges" \
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_node_x'
dimensions = ('nMesh2_node',)
attributes = {
		'long_name': "x-Coordinate of untrim grid node" ,\
        'units':'m'
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_node_y'
dimensions = ('nMesh2_node',)
attributes = {
		'long_name': "y-Coordinate of untrim grid node" ,\
        'units':'m'
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_face_x'
dimensions = ('nMesh2_face',)
attributes = {
		'long_name': "x-Coordinate of untrim grid face" ,\
        'units':'m'
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_face_y'
dimensions = ('nMesh2_face',)
attributes = {
		'long_name': "y-Coordinate of untrim grid face" ,\
        'units':'m'
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_edge_x'
dimensions = ('nMesh2_edge',)
attributes = {
		'long_name': "x-Coordinate of untrim polygon edge" ,\
        'units':'m'
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_edge_y'
dimensions = ('nMesh2_edge',)
attributes = {
		'long_name': "y-Coordinate of untrim polygon edge" ,\
        'units':'m'
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_edge_bc'
dimensions = ('nMesh2_edge',)
attributes = {
		'long_name': "untrim polygon edge boundary condition" ,\
        'flags':'none closed dirichlet'
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_face_bc'
dimensions = ('nMesh2_face',)
attributes = {
		'long_name': "untrim polygon boundary condition" ,\
        'flags':'none water_level'
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_edge_faces'
dimensions = ('nMesh2_edge','Two',)
attributes = {
		'long_name': "Maps every edge to its bounding faces" ,\
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_edge_nodes'
dimensions = ('nMesh2_edge','Two',)
attributes = {
		'long_name': "Maps every edge to its end nodes " ,\
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_face_nodes'
dimensions = ('nMesh2_face','nMaxMesh2_face_nodes',)
attributes = {
		'long_name': "Maps every polygon face its corner nodes" ,\
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_face_edges'
dimensions = ('nMesh2_face','nMaxMesh2_face_nodes',)
attributes = {
		'long_name': "Maps every polygon face its edges" ,\
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_edge_depth'
dimensions = ('nMesh2_edge',)
attributes = {
		'long_name': "Maximum depth of edge" ,\
        'units':'m',\
        'positive':'down',\
        'coordinates': "Mesh2_edge_x Mesh2_edge_y" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_face_depth'
dimensions = ('nMesh2_face',)
attributes = {
		'long_name': "Maximum depth of face" ,\
        'units':'m',\
        'positive':'down',\
        'coordinates': "Mesh2_face_x Mesh2_face_y" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_layer_3d'
dimensions = ('nMesh2_layer_3d',)
attributes = {
		'long_name': "elevation of layer" ,\
        'units':'m',\
        'positive':'up',\
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_data_time'
dimensions = ('nMesh2_data_time',)
attributes = {
        'units': "days since 1899-12-31 00:00:00.0" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_data_time_string'
dimensions = ('date_string_length','nMesh2_data_time')
attributes = {
        'long_name':'data time string',\
        'units': "yyyy-mm-dd HH:MM:SS" \
        }

dtype = 'c'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'h_flow_avg'
dimensions = ('nMesh2_edge','nMesh2_layer_3d','nMesh2_data_time')
attributes = {
		'long_name': "Horizontal volume flux averaged over integration interval" ,\
        'units':'m3 s-1',\
        'coordinates': "Mesh2_edge_x Mesh2_edge_y Mesh2_edge_z_3d" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'v_flow_avg'
dimensions = ('nMesh2_face','nMesh2_layer_3d','nMesh2_data_time')
attributes = {
		'long_name': "Vertical volume flux averaged over integration interval" ,\
        'units':'m3 s-1',\
        'coordinates': "Mesh2_face_x Mesh2_face_y Mesh2_edge_z_3d" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_edge_wet_area'
dimensions = ('nMesh2_edge','nMesh2_layer_3d','nMesh2_data_time')
attributes = {
		'long_name': "sea area" ,\
        'units':'m2',\
        'coordinates': "Mesh2_edge_x Mesh2_edge_y Mesh2_edge_z_3d" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_face_wet_area'
dimensions = ('nMesh2_face','nMesh2_layer_3d','nMesh2_data_time')
attributes = {
		'long_name': "sea area" ,\
        'units':'m2',\
        'coordinates': "Mesh2_face_x Mesh2_face_y Mesh2_edge_z_3d" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_edge_bottom_layer'
dimensions = ('nMesh2_edge','nMesh2_data_time')
attributes = {
		'long_name': "bottom most active layer (from either side of edge)" ,\
        'units':'',\
        'coordinates': "Mesh2_edge_x Mesh2_edge_y" \
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_edge_top_layer'
dimensions = ('nMesh2_edge','nMesh2_data_time')
attributes = {
		'long_name': "top most active layer (from either side of edge)" ,\
        'units':'',\
        'coordinates': "Mesh2_edge_x Mesh2_edge_y" \
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_face_bottom_layer'
dimensions = ('nMesh2_face','nMesh2_data_time')
attributes = {
		'long_name': "bottom most active layer (from either side of edge)" ,\
        'units':'',\
        'coordinates': "Mesh2_face_x Mesh2_face_y" \
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_face_top_layer'
dimensions = ('nMesh2_face','nMesh2_data_time')
attributes = {
		'long_name': "top most active layer (from either side of edge)" ,\
        'units':'',\
        'coordinates': "Mesh2_face_x Mesh2_face_y" \
        }

dtype = 'i4'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Mesh2_face_water_volume'
dimensions = ('nMesh2_face','nMesh2_layer_3d','nMesh2_data_time')
attributes = {
		'long_name': "Water prism volume" ,\
        'units':'m3',\
        'coordinates': "Mesh2_face_x Mesh2_face_y Mesh2_edge_z_3d" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_salinity_3d'
dimensions = ('nMesh2_face','nMesh2_layer_3d','nMesh2_data_time')
attributes = {
		'long_name': "salinity" ,\
        'units':'psu',\
        'coordinates': "Mesh2_face_x Mesh2_face_y Mesh2_edge_z_3d" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_vertical_diffusivity_3d'
dimensions = ('nMesh2_face','nMesh2_layer_3d','nMesh2_data_time')
attributes = {
		'long_name': "vertical diffusivity" ,\
        'units':'m2 s-1',\
        'coordinates': "Mesh2_face_x Mesh2_face_y Mesh2_edge_z_3d" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'Mesh2_sea_surface_elevation'
dimensions = ('nMesh2_face','nMesh2_data_time')
attributes = {
		'long_name': "sea surface elevation" ,\
        'units':'m',\
        'coordinates': "Mesh2_face_x Mesh2_face_y" \
        }

dtype = 'f8'		
untrim_ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})









class UNTRIMGrid(Grid):
    """
    UnTRIM grid class for loading from a netcdf file
    """ 
    def __init__(self,ncfile):
        self.ncfile=ncfile
        Grid.__init__(self,ncfile,gridvars=untrim_gridvars,griddims=untrim_griddims)

class UNTRIMSpatial(Spatial):
    """
    Class for handling UnTRIM netcdf hydrodynamic output
    """
    def __init__(self,ncfile,**kwargs):
        Spatial.__init__(self,ncfile,gridvars=untrim_gridvars,griddims=untrim_griddims,**kwargs)

    
    def loadDataRaw(self,variable=None):
        """ 
        Overloaded method - Untrim variables are ordered in a different order
        Untrim: [Nc,Nk,Nt] or [Nc,Nt]
        SUNTANS: [Nt,Nk,Nc] or [Nt,Nc]
        Load the specified suntans variable data as a vector
        """
        if variable==None:
            variable=self.variable
	
        if self.hasDim(variable,self.griddims['Ne']) and self.j==None:
            j=range(self.Ne)
        elif self.hasDim(variable,self.griddims['Nc']) and self.j==None:
            j=range(self.Nc)
        else:
            j = self.j
            
        nc = self.nc
        try:
            self.long_name = nc.variables[variable].long_name
        except:
            self.long_name = ''
        self.units= nc.variables[variable].units
        #        ndims = len(nc.variables[variable].dimensions)
        ndim = nc.variables[variable].ndim
        if ndim==1:
            self.data=nc.variables[variable][j]
        elif ndim==2:
            #print self.j
            data=nc.variables[variable][j,self.tstep]
            self.data = data.swapaxes(0,1)
        else:
            if self.klayer[0]==-1: # grab the seabed values
                raise Exception, 'Seabed extraction not implemented for UnTRIM'
                #klayer = np.arange(0,self.Nkmax)

                ##if type(self.tstep)==int:
                #if isinstance(self.tstep,(int,long)):
                #    data=nc.variables[variable][self.tstep,klayer,j]
                #    self.data = data[self.Nk[j],j]
                #else: # need to extract timestep by timestep for animations to save memory
                #    self.data=np.zeros((len(self.tstep),len(j)))
                #    i=-1
                #    for t in self.tstep:
                #        i+=1
                #        data=nc.variables[variable][t,klayer,j]
                #        self.data[i,:] = data[self.Nk[j],j]
            elif self.klayer[0]==-99: # Grab all layers
                klayer = np.arange(0,self.Nkmax) 
                #self.data=nc.variables[variable][self.tstep,klayer,j]
                data=nc.variables[variable][j,klayer,self.tstep]
                self.data = data.swapaxes(0,2)
            
            else:
                klayer = self.klayer
                data=nc.variables[variable][j,klayer,self.tstep]
                if data.ndim==3:
                    self.data = data.swapaxes(0,2)
                else:
                    self.data=data
        
        # Mask the data
        fillval = 999999.0
        
        self.mask = self.data==fillval
        self.data[self.mask]=0.
        self.data = self.data.squeeze()
        
        return self.data
 

