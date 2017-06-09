#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

SUNTANSHOME=../../main
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot
PYTHONEXEC=python
SUBCOMMAND=sh
VERBOSE=-vv

. $SUNTANSHOME/Makefile.in

MAINDATADIR=rundata
DATADIR=data
makescript=make_scenario_iwaves.py
BCFILE=IWave_BC.nc
RUNFILE=restart-run.txt

SUNTANS_DT=30

MAXRUNS=2
NUMPROCS=$1
RUNNUM=$2
NEXTRUN=`expr $RUNNUM + 1`

JOBDIR=${DATADIR}$RUNNUM

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS $SUN"
fi

if [ $RUNNUM -eq 1 ] ; then
    if [ ! -d $JOBDIR ] ; then
	cp -r $MAINDATADIR $JOBDIR
	echo Creating input files...
	$PYTHONEXEC scripts/$makescript $JOBDIR
	echo Creating grid...
	$EXEC -g -vvv --datadir=$JOBDIR
    else
	cp $MAINDATADIR/suntans.dat $JOBDIR/.
    fi

    echo Running suntans...
    $EXEC -s $VERBOSE --datadir=$JOBDIR 

else
    $EXEC -s -r $VERBOSE --datadir=$JOBDIR 
fi

#######
# Copy files for the next run
#######
if [ $RUNNUM -lt $MAXRUNS ] ; then
    if [ ! -d ${DATADIR}$NEXTRUN ] ; then
	cp -r $MAINDATADIR ${DATADIR}$NEXTRUN
    fi

    # Copy the grid files
    cp $JOBDIR/vertspace.dat ${DATADIR}$NEXTRUN/
    cp ${JOBDIR}/$BCFILE ${DATADIR}$NEXTRUN/

    # Copy the processor-based files
    n=0;
    while(true)
    do
      if [ $n -lt $NUMPROCS ] ; then
	  cp $JOBDIR/cells.dat.$n ${DATADIR}$NEXTRUN/
	  cp $JOBDIR/celldata.dat.$n ${DATADIR}$NEXTRUN/
	  cp $JOBDIR/topology.dat.$n ${DATADIR}$NEXTRUN/
	  cp $JOBDIR/edges.dat.$n ${DATADIR}$NEXTRUN/
	  cp $JOBDIR/edgedata.dat.$n ${DATADIR}$NEXTRUN/
	  cp $JOBDIR/nodes.dat.$n ${DATADIR}$NEXTRUN/

	  cp $JOBDIR/store.dat.$n ${DATADIR}$NEXTRUN/start.dat.$n
	  n=`expr $n + 1`
      else
	  break;
      fi
    done

    if [ ! -z $VERBOSE ] ; then
	echo Job finished on `date`
    fi

    ###
    #Read in some variables from a text file
    ###
    ii=1;
    while read -r line
    do 
        name="$line"
        if [ $ii -eq $NEXTRUN ]; then
            jj=0;
            for word in ${name}; do
        	if [ $jj -eq 1 ]; then
        	   BASETIME=$word
        	elif [ $jj -eq 2 ]; then
        	   NUMDAYS=$word
        	fi
        	#jj=$(($jj+1))
        	jj=`expr $jj + 1`
            done
        fi
        ii=`expr $ii + 1`
    done < $RUNFILE

    NSTEPS=$(($NUMDAYS*86400/$SUNTANS_DT))

    # Update suntans.dat parameters
    echo Updating suntans.dat...
    PARAM=readinitialnc
    PARAMVALUE=0
    sedstr='s/^'$PARAM'.*/'$PARAM'\t\t\t'$PARAMVALUE'\t#/g' 
    sed -i $sedstr ${DATADIR}$NEXTRUN/suntans.dat

    PARAM=thetaramptime
    PARAMVALUE=0
    sedstr='s/^'$PARAM'.*/'$PARAM'\t\t\t'$PARAMVALUE'\t#/g' 
    sed -i $sedstr ${DATADIR}$NEXTRUN/suntans.dat

    #PARAM=ncfilectr
    #PARAMVALUE=1
    #sedstr='s/^'$PARAM'.*/'$PARAM'\t\t\t'$PARAMVALUE'\t#/g' 
    #sed -i $sedstr ${DATADIR}$NEXTRUN/suntans.dat

    ###
    # Do not need to adjust this for the restart time
    #PARAM=starttime
    #PARAMVALUE=$BASETIME
    #sedstr='s/^'$PARAM'.*/'$PARAM'\t\t\t'$PARAMVALUE'\t#/g' 
    #sed -i $sedstr ${DATADIR}$NEXTRUN/suntans.dat

    PARAM=nsteps
    PARAMVALUE=$NSTEPS
    sedstr='s/^'$PARAM'.*/'$PARAM'\t\t\t'$PARAMVALUE'\t#/g' 
    sed -i $sedstr ${DATADIR}$NEXTRUN/suntans.dat

    ###
    # Run the next scenario
    #echo $NEXTRUN > $RUNNUMFILE
    $SUBCOMMAND iwaves_ridge-restart.sh $NUMPROCS $NEXTRUN
else
    echo Finished with $MAXRUNS runs!
    #if [ ! -z $VERBOSE ] ; then
    #    if [ $blowup -ne 0 ] ; then
    #        echo Run $RUNNUM is blowing up.  Not resubmitting.
    #    else
    #        echo Finished with $MAXRUNS runs!
    #    fi
    #fi
fi
