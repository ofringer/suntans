#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case and then play a movie
# of the results.
#
########################################################################

SUNTANSHOME=../..
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
datadir=data

NUMPROCS=$1

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
    echo Creating grid...
    $MPIHOME/bin/mpirun -np $NUMPROCS $SUN -t -g --datadir=$datadir
else
    cp $maindatadir/suntans.dat $datadir/.
fi

echo "Running suntans (need to run with -g -s when wetting/drying is employed!)"
$MPIHOME/bin/mpirun -np $NUMPROCS $SUN -g -s -vv --datadir=$datadir
