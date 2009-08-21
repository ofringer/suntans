#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

SUNTANSHOME=../..
SUN=./sun
SUNPLOT=$SUNTANSHOME/sunplot
TAGNAME=.
JOBDIR=.
RUNNUMFILE=runnum.txt
MAINDATADIR=./rundata
DATADIR=$JOBDIR/data
GRIDDIR=$JOBDIR/grid
SUBCOMMAND=sh
VERBOSE=-vv

NUMPROCS=$1
MAXRUNS=$2

. $SUNTANSHOME/Makefile.in

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS $SUN"
fi

if [ ! -d $JOBDIR ] ; then
    mkdir $JOBDIR
fi

if [ ! -f $RUNNUMFILE ] ; then
    echo 1 > $RUNNUMFILE
fi
RUNNUM=`cat $RUNNUMFILE`
NEXTRUN=`expr $RUNNUM + 1`

echo cat $DATADIR$RUNNUM/step.dat > status.sh

if [ ! -f sun ] ; then
    cp $SUNTANSHOME/sun .
fi

if [ ! -z $VERBOSE ] ; then
    echo Job $RUNNUM of $MAXRUNS started on `date`
fi

TODIR=${DATADIR}${RUNNUM}
if [ $RUNNUM -eq 1 ] ; then
    if [ ! -d $TODIR ] ; then
	cp -r rundata $TODIR
    fi

    if [ ! -d $GRIDDIR ] ; then
	cp -r rundata $GRIDDIR
	$EXEC -g $VERBOSE --datadir=$GRIDDIR
    fi

    $EXEC -s $VERBOSE --datadir=$TODIR >& crash.txt
else
    $EXEC -s -r $VERBOSE --datadir=$TODIR >& crash.txt
fi

blowup=`grep -c -H blowing crash.txt | awk -F: '{ print $2}'`

if [ ! -z $VERBOSE ] ; then
    cat crash.txt
fi

if [ $RUNNUM -lt $MAXRUNS -a $blowup -eq 0 ] ; then
    if [ ! -d ${DATADIR}$NEXTRUN ] ; then
	cp -r rundata ${DATADIR}$NEXTRUN
    fi

    n=0;
    while(true)
    do
      if [ $n -lt $NUMPROCS ] ; then
	  cp ${DATADIR}$RUNNUM/store.dat.$n ${DATADIR}$NEXTRUN/start.dat.$n
	  n=`expr $n + 1`
      else
	  break;
      fi
    done

    if [ ! -z $VERBOSE ] ; then
	echo Job finished on `date`
    fi

    echo $NEXTRUN > $RUNNUMFILE
    $SUBCOMMAND tides-restart.sh $NUMPROCS $MAXRUNS
else
    if [ ! -z $VERBOSE ] ; then
	if [ $blowup -ne 0 ] ; then
	    echo Run $RUNNUM is blowing up.  Not resubmitting.
	else
	    echo Finished with $MAXRUNS runs!
	fi
    fi
fi

