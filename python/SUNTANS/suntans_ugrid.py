"""
SUNTANS UGRID format in a python dictionary

To load the dictionary into your code type:
    >>from suntans_ugrid import ugrid

Matt Rayson
Stanford University
March 2013
"""

ugrid={}

fillval=999999.0

vname = 'suntans_mesh'
dimensions = ()
attributes = {'cf_role':'mesh_topology',\
		'long_name':'Topology data of 2D unstructured mesh',\
		'topology_dimension':'2',\
		'node_coordinates' : 'xp yp' ,\
		'face_node_connectivity':'cells' ,\
		'edge_node_connectivity': 'edges' ,\
		'face_coordinates':'xv yv' ,\
		'edge_coordinates': 'xe ye' ,\
		'face_edge_connectivity': 'face' ,\
		'edge_face_connectivity': 'grad' }
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'cells'
dimensions = ('Nc','numsides')
attributes = {'cf_role':'face_node_connectivity',\
		'long_name':'Maps every face to its corner nodes'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'face'
dimensions = ('Nc','numsides')
attributes = {'cf_role':'face_edge_connectivity',\
		'long_name':'Maps every face to its neighbouring faces'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'edges'
dimensions = ('Ne','two')
attributes = {'cf_role':'edge_node_connectivity',\
		'long_name':'Maps every edge to the two nodes it connects'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'neigh'
dimensions = ('Nc','numsides')
attributes = {'cf_role':'face_face_connectivity',\
		'long_name':'Maps every face to its neighbouring faces'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'grad'
dimensions = ('Ne','two')
attributes = {'cf_role':'edge_face_connectivity',\
		'long_name':'Maps every edge to the two faces it connects'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'xv'
dimensions = ('Nc',)
attributes = {'standard_name':'Easting',\
		'long_name':'Easting of 2D mesh face'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'yv'
dimensions = ('Nc',)
attributef = {'standard_name':'Northing',\
		'long_name':'Northing of 2D mesh face'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'xp'
dimensions = ('Np',)
attributes = {'standard_name':'Easting',\
		'long_name':'Easting of 2D mesh node'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'yp'
dimensions = ('Np',)
attributes = {'standard_name':'Northing',\
		'long_name':'Northing of 2D mesh nodes'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'xe'
dimensions = ('Ne',)
attributes = {'standard_name':'Easting',\
		'long_name':'Easting of 2D mesh edge'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'ye'
dimensions = ('Ne',)
attributes = {'standard_name':'Northing',\
		'long_name':'Northing of 2D mesh edge'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'normal'
dimensions = ('Nc','numsides')
attributes = {'long_name':'Dot product of unique normal with outward normal of each edge'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'n1'
dimensions = ('Ne',)
attributes = {'long_name':'x-component of the edge normal'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'n2'
dimensions = ('Ne',)
attributes = {'long_name':'y-component of the edge normal'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'df'
dimensions = ('Ne',)
attributes = {'long_name':'edge length',\
		'units':'m'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'dg'
dimensions = ('Ne',)
attributes = {'long_name':'Distance between faces on either side of edge',\
		'units':'m'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'def'
dimensions = ('Nc','numsides')
attributes = {'long_name':'Distance between faces and edges',\
		'units':'m'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Ac'
dimensions = ('Nc',)
attributes = {'long_name':'Horizontal area of 2D mesh',\
		'units':'m2',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'xv yv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'dz'
dimensions = ('Nk',)
attributes = {'long_name':'z layer spacing',\
		'units':'m'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'z_r'
dimensions = ('Nk',)
attributes = {'standard_name':'ocean_z_coordinate',\
		'long_name':'depth at layer mid points',\
		'units':'m',\
		'positive':'up'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'z_w'
dimensions = ('Nkw',)
attributes = {'standard_name':'ocean_z_coordinate',\
		'long_name':'depth at layer edges',\
		'units':'m',\
		'positive':'up'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Nk'
dimensions = ('Nc',)
attributes = {'long_name':'Number of layers at face'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'Nke'
dimensions = ('Ne',)
attributes = {'long_name':'Number of layers at edge'}
dtype = 'i4'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'dv'
dimensions = ('Nc',)
attributes = {'standard_name':'sea_floor_depth_below_geoid',\
		'long_name':'seafloor depth',\
		'units':'m',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'xv yv',\
		'positive':'down'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

vname = 'time'
dimensions = ('time',)
attributes = {'long_name':'time',\
		'units':'seconds since 1990-01-01 00:00:00'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False}})

### Model variables presently in use ###
vname = 'eta'
dimensions = ('time','Nc')
attributes = {'long_name':'Sea surface elevation',\
		'units':'m',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':False,'complevel':1,'fill_value':fillval}})

vname = 'uc'
dimensions = ('time','Nk','Nc')
attributes = {'long_name':'Eastward water velocity component',\
		'units':'m s-1',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_r yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'vc'
dimensions = ('time','Nk','Nc')
attributes = {'long_name':'Northward water velocity component',\
		'units':'m s-1',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_r yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'nu_v'
dimensions = ('time','Nk','Nc')
attributes = {'long_name':'Vertical eddy viscosity',\
		'units':'m2 s-1',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_r yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'salt'
dimensions = ('time','Nk','Nc')
attributes = {'long_name':'Salinity',\
		'units':'ppt',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_r yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'temp'
dimensions = ('time','Nk','Nc')
attributes = {'long_name':'Water temperature',\
		'units':'degrees C',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_r yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'rho'
dimensions = ('time','Nk','Nc')
attributes = {'long_name':'Water density',\
		'units':'kg m-3',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_r yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'w'
dimensions = ('time','Nkw','Nc')
attributes = {'long_name':'Vertical water velocity component',\
		'units':'m s-1',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_w yv xv'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

vname = 'U'
dimensions = ('time','Nk','Ne')
attributes = {'long_name':'Edge normal velocity',\
		'units':'m s-1',\
		'mesh':'suntans_mesh',\
		'location':'face',\
		'coordinates':'time z_r ye xe'}
dtype = 'f8'		
ugrid.update({vname:{'dimensions':dimensions,'attributes':attributes,'dtype':dtype,'zlib':True,'complevel':1,'fill_value':fillval}})

### Atmospheric Flux Variables ###
