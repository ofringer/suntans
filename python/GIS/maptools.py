# -*- coding: utf-8 -*-
"""
    Various useful tools for mapping.
	
	Mainly calls the GDAL library
	
	Uses include:
		Coordinate conversion
		Reading various raster formats
		Reading shapefiles
          Writing contours to a shapefile
		...
"""
from osgeo import osr
import ogr
import gdal
from gdalconst import * 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.nxutils as nxutils


import pdb


def ll2utm(LL,zone,CS='WGS84',north=True):
    """ Convert from lat/long coordinates to utm"""
    
    srs = osr.SpatialReference()
    if north:
        proj = "UTM %d (%s) in northern hemisphere."%(zone,CS)
    else:
        proj = "UTM %d (%s) in southern hemisphere."%(zone,CS)
        
    srs.SetProjCS( proj );
    srs.SetWellKnownGeogCS( CS );
    srs.SetUTM( zone, north );
    
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srsLatLong,srs )
    
    npt=np.size(LL,0)
    if npt > 2:
        XY=np.zeros((npt,2))
        for ii in range(0,npt):
            X,Y,z =  ct.TransformPoint(LL[ii,0],LL[ii,1])
            XY[ii,0]=X
            XY[ii,1]=Y
    else:
        XY=np.zeros((1,2))
        X,Y,z = ct.TransformPoint(LL[0],LL[1])
        XY[0,0]=X
        XY[0,1]=Y
    
    return XY
    
def utm2ll(XY,zone,CS='WGS84',north=True):
    """ Convert from utm coordinates to lat/long"""
    
    srs = osr.SpatialReference()
    if north:
        proj = "UTM %d (%s) in northern hemisphere."%(zone,CS)
    else:
        proj = "UTM %d (%s) in southern hemisphere."%(zone,CS)
        
    srs.SetProjCS( proj );
    srs.SetWellKnownGeogCS( CS );
    srs.SetUTM( zone, north );
    
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srs,srsLatLong )
    
    npt=np.size(XY,0)
    if npt > 2:
        LL=np.zeros((npt,2))
        for ii in range(0,npt):
            X,Y,z =  ct.TransformPoint(XY[ii,0],XY[ii,1])
            LL[ii,0]=X
            LL[ii,1]=Y
    else:
        LL=np.zeros((1,2))
        X,Y,z = ct.TransformPoint(XY[0],XY[1])
        LL[0,0]=X
        LL[0,1]=Y
    
    return LL
    
def readDEM(bathyfile,returnvec=False):
    """ Loads the data from a DEM file"""
    # register all of the drivers
    gdal.AllRegister()
    # open the image
    ds = gdal.Open(bathyfile, GA_ReadOnly)
    
    # Read the x and y coordinates
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    bands = ds.RasterCount
    
    geotransform = ds.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    
    x = originX + np.linspace(0,cols-1,cols)*pixelWidth
    y = originY + np.linspace(0,rows-1,rows)*pixelHeight
    X,Y=np.meshgrid(x,y)
    
    # Read the actual data
    data = ds.ReadAsArray(0,0,cols,rows)
    
    # Remove missing points
    data[data==-32767]=np.nan

    if returnvec:
        x=np.ravel(X)
        y=np.ravel(Y)
        z=np.ravel(data)
        ind = z != np.nan
        nc = np.sum(ind)
        XY = np.concatenate((np.reshape(x[ind],(nc,1)),np.reshape(y[ind],(nc,1))),axis=1)
        data = z[ind]
        return XY,data
    else:
        return X, Y, data     
    
def readShpBathy(shpfile,FIELDNAME = 'CONTOUR'):
    """ Reads a shapefile with line or point geometry and returns x,y,z
    
    See this tutorial:
        http://www.gis.usu.edu/~chrisg/python/2009/lectures/ospy_slides1.pdf
    """
    # Open the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    
    shp = driver.Open(shpfile, 0)
    
    lyr = shp.GetLayer()
    
    lyr.ResetReading()
    X=[]
    Y=[]
    Z=[]
    for feat in lyr:
        feat_defn = lyr.GetLayerDefn()
        for i in range(feat_defn.GetFieldCount()):
            field_defn = feat_defn.GetFieldDefn(i)
            if field_defn.GetName() == FIELDNAME:
                geom = feat.GetGeometryRef()
                ztmp = float(feat.GetField(i))
                if geom.GetGeometryType() == ogr.wkbPoint: # point
                    X.append(geom.getX())
                    Y.append(geom.getY())
                    Z.append(ztmp)
                elif geom.GetGeometryType() == 2:  # line
                    xyall=geom.GetPoints()
                    for xy in xyall:
                        X.append(xy[0])
                        Y.append(xy[1])
                        Z.append(ztmp)
                        
                elif geom.GetGeometryType() == 5:  # multiline
#                    pdb.set_trace()
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        xyall=geom2.GetPoints()
                        for xy in xyall:
                            X.append(xy[0])
                            Y.append(xy[1])
                            Z.append(ztmp)
#                    print geom.GetGeometryName()
#                    print geom.GetGeometryType()
    shp=None
    nc = len(X)
    XY = np.concatenate((np.reshape(np.array(X),(nc,1)),(np.reshape(np.array(Y),(nc,1)))),axis=1)
    del X
    del Y
    return XY,np.array(Z)
    
def readShpPoly(shpfile,FIELDNAME = 'marker'):
    """ Reads a shapefile with polygon geometry and returns x,y and FIELDNAME value
    
    See this tutorial:
        http://www.gis.usu.edu/~chrisg/python/2009/lectures/ospy_slides1.pdf
    """
    # Open the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    
    shp = driver.Open(shpfile, 0)
    
    lyr = shp.GetLayer()
    
    lyr.ResetReading()
    XY=[]
    field=[]
    for feat in lyr:
        feat_defn = lyr.GetLayerDefn()
        for i in range(feat_defn.GetFieldCount()):
            field_defn = feat_defn.GetFieldDefn(i)
            if field_defn.GetName() == FIELDNAME:
                geom = feat.GetGeometryRef()
                ztmp = feat.GetField(i)
               
                if geom.GetGeometryType() == ogr.wkbPolygon:  # Polygon
                
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        xyall=geom2.GetPoints()
                        
                        XY.append(np.asarray(xyall))
                        field.append(ztmp)
                        
                if geom.GetGeometryType() == ogr.wkbMultiPolygon:  # Multi Polygon
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        for jj in range(0,geom2.GetGeometryCount()):
                            geom3 = geom2.GetGeometryRef(jj)
                            xyall=geom3.GetPoints()
                            
                            XY.append(np.asarray(xyall))
                            field.append(ztmp)  
        
    shp=None
    return XY,field

def maskShpPoly(X,Y,shpfile,FIELDNAME = 'marker'):
    """
    Return a mask array of size(X/Y) 
    """
    
    # Read the polygon from the shape file
    XY,edge_id = readShpPoly(shpfile,FIELDNAME=FIELDNAME)
    
    # Reshape the input array
    sz = X.shape
    X = X.ravel()
    Y = Y.ravel()
        
    ind = nxutils.points_inside_poly(np.vstack((X,Y)).T,XY[0])
    
    mask = np.zeros(sz,dtype=np.int)
    ind = ind.reshape(sz)
    mask[ind]=1
    
    return mask,XY[0]
    
    
def readraster(infile):
    """ Loads the data from any raster-type file
        eg. *.dem, *.grd,...    
    """
    # register all of the drivers
    gdal.AllRegister()
    # open the image
    ds = gdal.Open(infile, GA_ReadOnly)
    
    # Read the x and y coordinates
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    bands = ds.RasterCount
    
    geotransform = ds.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    
    x = originX + np.linspace(0,cols-1,cols)*pixelWidth
    y = originY + np.linspace(0,rows-1,rows)*pixelHeight
    
    # Read the actual data
    data = ds.ReadAsArray(0,0,cols,rows)
    
    # Remove missing points
    data[data==-32767]=np.nan

def plotmap(shpfile,color='0.5',fieldname='FID',convert=None,zone=15):
    """
    Plots a map layer from a polygon or multipolygon shapefile
    
    Usage:
        plotmap(shpfile,color='0.5',fieldname='FID',convert=None,zone=15)
        
        convert = 'll2utm' or 'utm2ll'
        zone = utm zone number
    """
    
    from matplotlib.collections import PolyCollection
    
    # Read the shapefile
    xy,marker = readShpPoly(shpfile,FIELDNAME=fieldname)
    
    if convert=='utm2ll':
        ll=[]
        for xytmp in xy:
            ll.append(utm2ll(xytmp,zone,CS='WGS84',north=True))
        xy=ll
    elif convert=='ll2utm':
        ll=[]
        for xytmp in xy:
            ll.append(ll2utm(xttmp,zone,CS='WGS84',north=True))
        xy=ll
        
    # Add the polygons to the current axis as a series of patches
    fig = plt.gcf()
    ax = fig.gca()
    
    collection = PolyCollection(xy)
    collection.set_facecolor(color)
    ax.add_collection(collection)
    ax.axis('equal')
    
    return ax, collection

def Contour2Shp(C,outfile):
    """
    Converts a matplotlib contour object to a shapefile
    
    The shapefile has the filed Contour corresponding to the contour levels set
    in matplotlib.pyplot.contour
    
    **The projection is hardwired to WGS84 for now. This needs updating.
    """
    # Convert the contour values to a shapefile
    # see this example:
    #    http://invisibleroads.com/tutorials/gdal-shapefile-points-save.html
    
    import osgeo.ogr, osgeo.osr
    
    # Create the projection
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    
    # Create the shape file
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    shapeData = driver.CreateDataSource(outfile)
    
    # Create the layer
    layer = shapeData.CreateLayer('Contour', spatialReference, osgeo.ogr.wkbLineString)
    
    # Create a field containing the contour level
    field_def = osgeo.ogr.FieldDefn('Contour', osgeo.ogr.OFTReal)
    layer.CreateField(field_def)
    layerDefinition = layer.GetLayerDefn()
    
    
    # Loop through the contour object to get the coordinates of each layer
    ctr=0
    for coll,lev in zip(C.collections,C.levels):
        coords = coll.get_paths()
        for xyObj in coords:
            ctr+=1
            line = osgeo.ogr.Geometry(osgeo.ogr.wkbLineString)
            # Add points individually to the line
            for xy in xyObj.vertices:
                line.AddPoint_2D(xy[0],xy[1])
            
            # Update the feature with the line data
            featureIndex = ctr
            feature = osgeo.ogr.Feature(layerDefinition)
            feature.SetGeometry(line)
            feature.SetFID(featureIndex)
            feature.SetGeometryDirectly(line)
            # Set the contour level
            feature.SetField('Contour',lev)
            
            layer.CreateFeature(feature)
    
    # Close the shape file
    shapeData.Destroy()
    print 'Complete - shapefile written to:\n      %s'%outfile
 
def Polygon2GIS(xynodes,shpfile,zone,CS='WGS84',north=True):
    """
    Convert a list of polygons (Nx2 arrays) to a GIS vector format: KML or SHP.
    
    Each item in the list should contain an Nx2 array where column one is the x 
    coordinates and column two is the y coordinates of a polygon. Polygon should be closed.
    
    The output format is based on the outfile extension.
    
    The projection is assumed to be UTM. 
    """
    import osgeo.ogr, osgeo.osr
    import os
    
    # Set the projection
    srs = osgeo.osr.SpatialReference()
    if north:
        proj = "UTM %d (%s) in northern hemisphere."%(zone,CS)
    else:
        proj = "UTM %d (%s) in southern hemisphere."%(zone,CS)
        
    srs.SetProjCS( proj );
    srs.SetWellKnownGeogCS( CS );
    srs.SetUTM( zone, True );
        
    # Create the shape file
    ext=shpfile[-3:]
    if ext.lower()=='shp':
        driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    elif ext.lower()=='kml':
        driver = osgeo.ogr.GetDriverByName('KML')
    else:
        print 'Error. Unknown file extension: %s'%ext.lower()
        return
    
    
    if os.path.exists(shpfile):
        os.unlink(shpfile)
    shapeData = driver.CreateDataSource(shpfile)
    
    # Create a layer
    layer = shapeData.CreateLayer('Grid', srs, osgeo.ogr.wkbPolygon)
    layerDefinition = layer.GetLayerDefn()    
    
    # Loop through the list of nodes to get the coordinates of each polygon
    ctr=0
    for xy in xynodes:
        ctr+=1
        ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
        
        # Add points individually to the polygon
        for nodes in xy:
            ring.AddPoint(nodes[0],nodes[1])
            
        poly = osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
        poly.AddGeometry(ring)
    
        # Update the feature with the polygon data
        featureIndex = ctr
        feature = osgeo.ogr.Feature(layerDefinition)
        feature.SetGeometry(poly)
        feature.SetFID(featureIndex)
        layer.CreateFeature(feature)
        feature.Destroy()
    # Close the shape file
    shapeData.Destroy()
    print 'Complete - file written to:\n      %s'%shpfile
    
###Testing###
#LL=[-94.2,27.0]
#zone=15
#ll2utm(LL,zone,CS='WGS84',north=True)

#shpfile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/TNRIS/BathyTopoTXCoast_V0.2/Shapefiles/Bathymetry.shp'
#XY,Z = readShpBathy(shpfile)
#
#fig= plt.figure(figsize=(8,9))
##h.imshow(np.flipud(self.Z),extent=[bbox[0],bbox[1],bbox[3],bbox[2]])
#plt.plot(XY[:,0],XY[:,1],'.')
#plt.axis('equal')
#plt.show()

#shpfile= 'C:\Projects\GOMGalveston\DATA\Shoreline/GalvestonBasemap.shp'
#
#plt.figure()
#xy = plotmap(shpfile,convert='utm2ll')
#plt.show()
