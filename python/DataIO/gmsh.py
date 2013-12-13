# -*- coding: utf-8 -*-
"""
Tools for interacting with GMSH grid files

Example 1) Create a gmsh input file, generate grid (with gmsh) and convert to SUNTANS format
-------

```
from gmsh import create_gmsh_geo, gmsh2suntans
import os

gridfolder = './'
coastfile = 'coast.shp'
scalefile = 'scale.shp'
geofile = 'gmsh_test.geo'
mshfile = 'gmsh_test.msh'

# Step 1) Create the input geo file from two shapefiles
create_gmsh_geo(gridfolder+coastfile,gridfolder+scalefile,gridfolder+geofile)

# Step 2) Command-line call to gmsh from python
os.system('gmsh %s -2 -o %s'%(gridfolder+geofile,gridfolder+mshfile))

# Step 3) Convert to suntans ascii format
gmsh2suntans(gridfolder+mshfile,gridfolder)

```

Matt Rayson
Stanford University
October 2013

"""

import numpy as np
from maptools import readShpPoly, readShpPointLine
from hybridgrid import HybridGrid

import pdb

def create_gmsh_geo(shpfile,scalefile,geofile,startpoly=0,ndmin=8,scalefac=1.0,\
    r=1.08,lcmax=2000.0,sigmoid=0):
    """
    Generate a gmsh *.geo file using a boundary from a shpfile
    
    Inputs:
    ---
        - shpfile : a polygon shapefile with the domain boundaries. Make sure 
        that the first polygon is the outer ring. If not change 'startpoly' parameter.
        - scalefile : a line shapefile with field 'scale' specifying the target grid
        resolution for that region.
        - geofile: output geo ascii file.
    *Optional*
        - startpoly [default=0] - the index of the outer polygon in the shapefile.
        - ndmin [default=8] - the minimum range of the scale = ndmin*scale
        - r : expansion factor
        - lcmax [default=2000] - the maximum grid cell size.
    """
    # Load the polygon shape file 
    xy,field = readShpPoly(shpfile,FIELDNAME='FID')
    
    # Load the scale shape file
    xyscale,scale = readShpPointLine(scalefile,FIELDNAME='scale')
    
    # Load the 'embed' flag to see if the scale layer should be embedded
    try:
        xyscale,embed = readShpPointLine(scalefile,FIELDNAME='embed')
    except:
        print 'Warning - could not find "embed" field in shapefile. Setting embed = 0.'
        embed = [ss*0 for ss in scale]
        
    
    ## output geo and svg file
    fgeo = open(geofile,'w')
    
    fgeo.write("""IP = newp;
    IL = newl;
    IS = news;
    IF = newf;
    """ )
    
    ip = 0 # Point counter
    il = 0 # Line counter
    rp = 0 
    lines=[]
    for loop in xy:
        firstp = ip
        for p in range(loop.shape[0]):  
            rp = rp + 1
            fgeo.write("Point(IP + %i) = {%.16e, %.16e, %.16e}; // %i\n" % (ip, loop[p,0], loop[p,1], 0., rp -1))
            ip = ip + 1
        
        fgeo.write("BSpline(IL + %i) = {IP + %i : IP + %i, IP + %i};\n" % (il, firstp, ip - 1, firstp))
    
        il = il + 1
        lines.append(il-1)
    
    
    # Create the surface polygon
    fgeo.write('\n//%s\n//Surface Polygon Definition\n//%s\n\n'%(72*'#',72*'#'))
    
    surfstart=il
    fgeo.write("Line Loop(IL + %i) = {IL + %s};\n" % (il,startpoly) )
    il += 1
    for ll in lines:
        if ll != startpoly:
            fgeo.write("Line Loop(IL + %i) = {IL + %s};\n" % (il,ll) )
            il += 1
    surfend = il - 1
        
    fgeo.write("Plane Surface(IL + %i) = {IL + %i : IL + %i};\n"%(il,surfstart,surfend))
    fgeo.write('Physical Surface("Ocean") = {IL + %i};\n'%il)
    surface_id = il # Keep this to embed lines into
    il += 1
    
    # Create the grid scale lines and fields
    fgeo.write('\n//%s\n//Grid Scale Definition\n//%s\n\n'%(72*'#',72*'#'))
    slines=[] # Reference to scale lines
    for loop,ss,ee in zip(xyscale,scale,embed):
        firstp = ip
        
        ss *= scalefac # Applies scale factor
        for p in range(loop.shape[0]):  
            rp = rp + 1
            #fgeo.write("Point(IP + %i) = {%.16e, %.16e, %.16e}; // %i\n" % (ip, loop[p,0], loop[p,1], 0., rp -1))
            fgeo.write("Point(IP + %i) = {%.16e, %.16e, %.16e,%.16e}; // %i\n" % (ip, loop[p,0], loop[p,1], 0., rp -1,float(ss)))
            ip = ip + 1
        
        if ee:
            #fgeo.write("BSpline(IL + %i) = {IP + %i : IP + %i};\n" % (il, firstp, ip - 1))
            fgeo.write("Spline(IL + %i) = {IP + %i : IP + %i};\n" % (il, firstp, ip - 1))
            # Embed this line in the main surface so that cells are aligned
            fgeo.write("Line{IL + %i} In Surface{IL + %i};\n"%(il,surface_id))
        else:
            # Don't embed (can set as BSpline)
            fgeo.write("BSpline(IL + %i) = {IP + %i : IP + %i};\n" % (il, firstp, ip - 1))
            
        slines.append(il)
        il = il + 1
    
    ifield = 0 # Field counter    
    fids = []
    ii=0
    for ss,line in zip(scale,slines):
        fgeo.write("Field[IF + %i] = Attractor;\n"%ifield)
        fgeo.write("Field[IF + %i].EdgesList = {IL + %i};\n"%(ifield,line))
        nodesperedge = get_nnodes_from_line(xyscale[ii],2.0*float(ss))
        fgeo.write("Field[IF + %i].NNodesByEdge = %i;\n"%(ifield,nodesperedge))
        
        ifield+=1
        fgeo.write("Field[IF + %i] = Threshold;\n"%ifield)
        
        # Find the maximum distance
        Nk =  np.log(lcmax/ss)/np.log(r)
        print ss,Nk
        lmin = float(ndmin)*float(ss)
        lmax = lmin + Nk * ss
        fgeo.write("Field[IF + %i].DistMax = %6.2f;\n"%(ifield,lmax))
        fgeo.write("Field[IF + %i].DistMin = %6.2f;\n"%(ifield,lmin))
        fgeo.write("Field[IF + %i].IField = IF + %i;\n"%(ifield,ifield-1))
        fgeo.write("Field[IF + %i].LcMax = %6.2f;\n"%(ifield,lcmax))
        fgeo.write("Field[IF + %i].LcMin = %6.2f;\n"%(ifield,float(ss)))
        fgeo.write("Field[IF + %i].Sigmoid = %d;\n"%(ifield,sigmoid))
        fids.append(ifield)
        ifield+=1
        ii+=1
    
     
    fieldlist = '{'
    for ff in fids:
        fieldlist += 'IF + %d, '%ff
    fieldlist = fieldlist[:-2]+'}'
    
    fgeo.write("Field[IF + %i] = Min;\n"%ifield)
    fgeo.write("Field[IF + %i].FieldsList = %s;\n"%(ifield,fieldlist))
    
    fgeo.write("Background Field = IF + %i;\n"%ifield)
    
    # Set Quad options by default
    fgeo.write("Mesh.CharacteristicLengthMax=%.16e;//Max cell size\n"%lcmax)
    fgeo.write("Mesh.RecombineAll=1;//recombine all defined surfaces\n")
    fgeo.write("Mesh.Algorithm=8;//delquad mesher\n")
    fgeo.write("Mesh.Smoothing=1;//10 smoothing steps\n")
    fgeo.close()
    
    print 'Complete. GMSH geo file written to:\n\t %s'%geofile

def gmsh2suntans(mshfile,suntanspath):
    """
    Convert a gmsh format grid to suntans edges.dat, cells.dat, points.dat
    """
    grd = gmsh2hybridgrid(mshfile)

    grd2suntans(grd,suntanspath)
    
    
def grd2suntans(grd,suntanspath):
    ### Save cells.dat into hybrid grid format
    f = open(suntanspath+'/cells.dat','w')

    for ii in range(grd.Ncells()):
        outstr = '%d %10.6f %10.6f '%(grd.nfaces[ii],grd.xv[ii],grd.yv[ii])
        for nn in range(grd.nfaces[ii]):
            outstr += '%d '%grd.cells[ii,nn]

        for nn in range(grd.nfaces[ii]):
            outstr += '%d '%grd.neigh[ii,nn]

        outstr += '\n'
        f.write(outstr)
    
    f.close()
    
    # Save edges.dat
    f = open(suntanspath+'/edges.dat','w')

    for ee,m,gg in zip(grd.edges,grd.mark,grd.grad):
        e1=ee[0]
        e2=ee[1]
        g1=gg[0]
        g2=gg[1]
        f.write('%d %d  %d  %d  %d  0\n'%(e1,e2,m,g1,g2))

    f.close()
    
    # Save to points.dat
    f = open(suntanspath+'/points.dat','w')

    for x,y in zip(grd.xp,grd.yp):
        f.write('%10.6f %10.6f  0\n'%(x,y))

    f.close()
    print 'Completed gmsh to suntans conversion.'
    
def read_msh(mshfile):
    """
    Reads a gmsh *.msh file
    
    Returns:
        - nodes: a dictionary with nodal coordinates
        - elements: a dictionary with element dataq
        
    """
    # Code modified from here:
    #   http://users.monash.edu.au/~bburn/src/gmsh2sem.py
    
    offset = -1 # Go to python style indexing
    ifile = open(mshfile,'r')
    while 1:
        line = ifile.readline()
        if '$Nodes' in line:
            line = ifile.readline()
            Nnode = int (line)
            nodes={}
            for i in range (1, Nnode+1):
                line = ifile.readline()
                ll = line.split()
                nodes.update({int(ll[0]):{'X':float(ll[1]),'Y':float(ll[2]),'Z':float(ll[3])}})
    
                
        elif '$Elements' in line:
            line = ifile.readline()
            Nel = int (line)
            
            elements={}
            for i in range (1, Nel+1):
                line = ifile.readline()
                ll = line.split()
                eid = int(ll[0])
                etype = int(ll[1])
                ntags = int(ll[2])
                tags=[]
                for n in range(3,3+ntags):
                    tags.append(int(ll[n]))
                
                cells=[]
                nn=len(ll)
                for n in range(3+ntags,nn):
                    cells.append(int(ll[n]) + offset) 
                
                elements.update({eid:{'type':etype,'ntags':ntags,'tags':tags,'cells':cells}})
    
            break
    
    ifile.close()
    
    return nodes, elements

def gmsh2hybridgrid(mshfile):
    """
    Convert the gmsh msh file to the hybrid grid class
    """

    celltype = {2:3,3:4} # lookup table for cell type

    # Read the raw ascii file    
    nodes, elements  = read_msh(mshfile)
    
    # Get the coordinates
    x = [nodes[ii]['X'] for ii in nodes.keys()]
    y = [nodes[ii]['Y'] for ii in nodes.keys()]
    
    cells = [elements[ii]['cells'] for ii in elements.keys() if elements[ii]['type'] in celltype.keys()]
    nfaces = [len(cells[ii]) for ii in range(len(cells))]
    
    Nc = len(cells)
    cellsin = np.zeros((Nc,max(nfaces)),np.int)-999999
    for ii in range(Nc):
        cellsin[ii,0:nfaces[ii]]=cells[ii]

    return HybridGrid(x,y,cellsin,nfaces=nfaces)    


def get_nnodes_from_line( xy, dx):
    """
    Find the number of nodes along an edge to resolve it
    """
    dist = np.sqrt( (xy[1::,0]-xy[:-1,0])**2 + (xy[1::,1]-xy[:-1,1])**2)
    maxdist = np.max(dist.cumsum())
    return np.ceil(maxdist/dx).astype(int)

def write_pos_file(posfile,X,Y,Z):
    """
    Writes gridded xyz data to a background *.pos file
    """
    nj,ni = X.shape

    fpos = open(posfile,'w')

    fpos.write('View "background mesh" {\n')
    
    for j in range(nj-1):
        for i in range(ni-1):
            fpos.write("SQ(%.16e,%.16e,0,%.16e,%.16e,0,%.16e,%.16e,0,%.16e,%.16e,0){%.16e,%.16e,%.16e,%.16e};\n"%\
                (X[j,i],Y[j,i],X[j,i+1],Y[j,i+1],X[j+1,i+1],Y[j+1,i+1],X[j+1,i],Y[j+1,i],\
                Z[j,i],Z[j,i+1],Z[j+1,i+1],Z[j+1,i]))
    
    fpos.write("};\n")
    fpos.close()
    print 'Complete. pos file written to:\n\t%s'%posfile


def create_pos_file(posfile,scalefile, xlims,ylims,dx,\
    geofile=None,ndmin=5, lcmax=2000.,r=1.05, scalefac=1.0):
    """
    Generates a gmsh background scale file (*.pos)
    
    If a geofile is specified the mesh is embedded
    """
    from shapely import geometry, speedups

    if speedups.available:
        speedups.enable()
        
    X,Y = np.meshgrid(np.arange(xlims[0],xlims[1],dx),np.arange(ylims[0],ylims[1],dx))
    xy = np.vstack((X.ravel(),Y.ravel())).T
    Np = xy.shape[0]
    nj,ni=X.shape
    # Load the scalefile
    xyscale,gridscale = readShpPointLine(scalefile,FIELDNAME='scale')
    
    # Load all of the points into shapely type geometry
    
    # Distance method won't work with numpy array
    #P = geometry.asPoint(xy)
    
    P = [geometry.Point(xy[i,0],xy[i,1]) for i in range(Np)]
    
    L=[]
    for ll in xyscale:
        L.append(geometry.asLineString(ll))
     
    nlines = len(L)
    scale_all = np.zeros((nj,ni,nlines))
    for n in range(nlines):
        print 'Calculating distance from line %d...'%n
        ss = gridscale[n] * scalefac
        lmin = ndmin * ss
        
        # Find the maximum distance
        Nk =  np.log(lcmax/ss)/np.log(r)
        print ss,Nk
        lmax = lmin + Nk * ss
        
        dist = [L[n].distance(P[i]) for i in range(Np)]
        dist = np.array(dist).reshape((nj,ni))
        
        # Calculate the scale
        N = (dist-lmin)/ss
        scale = ss*r**N
        
        ind = dist<=lmin
        if ind.any():
            scale[ind] = ss
            
        ind = scale>lcmax
        if ind.any():
            scale[ind] = lcmax
        
        scale_all[:,:,n] = scale
        
    scale_min = scale_all.min(axis=-1)
    
    write_pos_file(posfile,X,Y,scale_min)  
    
    if not geofile == None:
        fgeo = open(geofile,'a')
        fgeo.write("// Merge a post-processing view containing the target mesh sizes\n")
        fgeo.write('Merge "%s";'%posfile)
        fgeo.write("// Apply the view as the current background mesh\n")
        fgeo.write("Background Mesh View[0];\n")
        fgeo.close()
    