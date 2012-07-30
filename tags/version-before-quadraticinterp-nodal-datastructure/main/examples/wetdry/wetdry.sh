#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

SUNTANSHOME=../..
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
datadir=data

NUMPROCS=$1

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS $SUN"
fi

if [ -z "$TRIANGLEHOME" ] ; then
    echo Error: This example will not run without the triangle libraries.
    echo Make sure TRIANGLEHOME is set in $SUNTANSHOME/Makefile.in
    exit 1
fi

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
    echo Creating grid...
    $EXEC -t -g --datadir=$datadir
else
    cp $maindatadir/suntans.dat $datadir/.
fi

echo "Running suntans (need to run with -g -s when wetting/drying is employed!)"
$EXEC -g -s -vv --datadir=$datadir

