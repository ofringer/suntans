#!/usr/bin/python
"""
Converts the last time step of a suntans history netcdf file to an initial condition format
"""

from sunpy import Spatial
from sunboundary import InitialCond
from datetime import datetime
import sys

def sunhis2initial(hisfile,icfile):
    """
    Main function
    """

    print '#####\nCreating initial condition file (%s) from history file (%s)...'%(icfile,hisfile)

    # Load the history file
    sunhis = Spatial(hisfile, tstep=-1, klayer=[-99])

    # Load the initial condition object
    timeic = datetime.strftime(sunhis.time[-1],'%Y%m%d.%H%M%S')
    sunic = InitialCond(hisfile,timeic)

    # Load the last time step into each variable
    sunic.h = sunhis.loadData(variable='eta').reshape((1,sunic.h.shape))
    sunic.T = sunhis.loadData(variable='temp').reshape((1,)+sunic.T.shape)
    sunic.S = sunhis.loadData(variable='salt').reshape((1,)+sunic.S.shape)
    sunic.uc = sunhis.loadData(variable='uc').reshape((1,)+sunic.uc.shape)
    sunic.vc = sunhis.loadData(variable='vc').reshape((1,)+sunic.vc.shape)

    print sunic.h.shape
    print sunic.T.shape

    # Load the age
    if sunhis.hasVar('agec'):
        sunic.agec = sunhis.loadData(variable='agec').reshape((1,)+sunic.agec.shape)       
        sunic.agealpha = sunhis.loadData(variable='agealpha').reshape((1,)+sunic.agealpha.shape)       

    # Write to the output file
    sunic.writeNC(icfile)

def usage():
    print 'Error. Correct usage:\npython sunhis2initial.py "hisfile.nc" "icfile.nc"'

if __name__=='__main__':
    """ 
    Command line call
    """
    try:
    	hisfile = sys.argv[1]
	icfile = sys.argv[2]
    except:
    	usage()
	exit(1)

    sunhis2initial(hisfile,icfile)

