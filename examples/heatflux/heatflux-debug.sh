#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

SUNTANSHOME=../../main
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
datadir=data

NUMPROCS=$1

if [ -z "$MPIHOME" ] ; then
    EXEC="gdb -command=gdb.commands $SUN"
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS xterm -e gdb -command=gdb.commands $SUN"
fi
#EXEC="gdb -command=gdb.commands $SUN"

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
else
    cp $maindatadir/suntans.dat $datadir/.
fi

echo Running suntans in debug mode...
$EXEC 

