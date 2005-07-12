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

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
    echo Creating grid...
    $MPIHOME/bin/mpirun -np $NUMPROCS $SUN -t -g --datadir=$datadir
else
    cp $maindatadir/suntans.dat $datadir/.
fi

echo Running suntans...
$MPIHOME/bin/mpirun -np $NUMPROCS $SUN -s -vv --datadir=$datadir

