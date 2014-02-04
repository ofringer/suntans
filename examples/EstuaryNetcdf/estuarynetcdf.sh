#!/bin/sh
########################################################################
#
# Create the estuary test case netcdf input files 
#
########################################################################

pythonhome=../../../../python
maindatadir=grid
datadir=rundata

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
fi

echo Running python script...
python scripts/suntans_driver.py
