#!/bin/bash
#
# Use this script to start job number N+1, where N is the job number of the
# last job in the file runnum.txt.
#
DATADIR=$1/data
NUMPROCS=$2
MAXRUNS=$3
RUNNUMFILE=runnum.txt

if [ -z $1 -o -z $2 ] ; then
    echo Need to specify ./restart.sh directory numprocs
    exit 1;
fi


if [ ! -f $RUNNUMFILE ] ; then
    echo File $RUNNUMFILE does not exist.
    exit 1;
else
    RUNNUM=`cat $RUNNUMFILE`
    NEXTRUN=`expr $RUNNUM + 1`
fi

if [ ! -d ${DATADIR}${RUNNUM} ] ; then
    echo Directory ${DATADIR}${RUNNUM} does not exist.
    exit 1;
fi

if [ ! -d ${DATADIR}$NEXTRUN ] ; then
    echo Copying rundata to ${DATADIR}$NEXTRUN
    cp -r rundata ${DATADIR}$NEXTRUN

    n=0;
    while(true)
    do
      if [ $n -lt $NUMPROCS ] ; then
	  echo Copying ${DATADIR}$RUNNUM/store.dat.$n to ${DATADIR}$NEXTRUN/start.dat.$n
	  cp ${DATADIR}$RUNNUM/store.dat.$n ${DATADIR}$NEXTRUN/start.dat.$n
	  n=`expr $n + 1`
      else
	  break;
      fi
    done

    echo $NEXTRUN > $RUNNUMFILE
else
    echo ${DATADIR}/data$NEXTRUN already exists.
fi

echo Submitting job $NEXTRUN...
sh tides-restart.sh $NUMPROCS $MAXRUNS


