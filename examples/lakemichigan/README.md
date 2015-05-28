# Lake Michigan SUNTANS tutorial

## Generating a grid

### Step 1) 

 - Unpack the gis shoreline data:

```
tar -xvf gis.tar gis
```

### Step 2)

 - Download the bathymetry data:

```
wget http://geo.glin.net/gis/shps/glin_bathymetry_lake_michigan.zip
unzip glin_bathymetry_lake_michigan.zip gis/
```

## Downloading NARR data and convert to SUNTANS format

The following code demonstrates how to download North American Regional Reanalysis (NARR) atmospheric data and convert it into the SUNTANS meteorological netcdf input format.

```python
from createSunMetFile import narr2suntans

bbox = [-96.5, -92.0, 41.5,46.0]
utmzone = 16
tstart = '20120101'
tend = '20120201'
ncfile = 'data/LkMichigan_NARR_2012.nc'

narr2suntans(ncfile,tstart,tend,bbox,utmzone)
```


