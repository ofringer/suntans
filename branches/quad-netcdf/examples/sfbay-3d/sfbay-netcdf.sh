#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

SUNTANSHOME=../../main
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot
PYTHONHOME=/home/mrayson/suntans/SourceForge/python

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
datadir=data

NUMPROCS=$1
DIM=$2

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS $SUN"
fi

if [ -z $DIM ] ; then
    DIM=2
fi
cp $maindatadir/suntans.dat-${DIM}d_netcdf $maindatadir/suntans.dat

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir

    echo Generating input files...
    python scripts/suntans_driver_${DIM}d.py

    echo Creating grid...

    $EXEC -g --datadir=$datadir
else
    cp $maindatadir/suntans.dat $datadir/.
fi


Nkmax=`cat data/suntans.dat | awk '{ print $1" "$2 }' | grep Nkmax | awk '{ print $2 }'`
NkmaxData=`wc -l data/vertspace.dat | awk '{ print $1 }'`

if [ $Nkmax -ne $NkmaxData ] ; then
    echo You specified to run with Nkmax=$Nkmax, although the grid in the $datadir directory is for Nkmax=$NkmaxData.
    if [ $NkmaxData -eq 1 ] ; then
	echo Please use "make test2d" instead.
    else
	echo Please use "make test3d" instead.
    fi
    echo Or, alternatively, use "make clobber",  followed by "make test2d" or "make test3d".
    exit 1
fi

if [ $Nkmax -eq 1 ] ; then    
    echo Running example in 2D...
else
    echo Running example in 3D with $Nkmax levels...
fi    
$EXEC -s -vv --datadir=$datadir

echo Joining output files...

python $PYTHONHOME/SUNTANS/joinsun.py -f SFBay${DIM}D.nc -p $datadir -t 24 -n $NUMPROCS
