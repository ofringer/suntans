#!/bin/sh
########################################################################
#
# Shell script to run through several different time step sizes
# to determine the time accuracy of the code as a test case.
#
########################################################################

function calc () {
   awk "BEGIN { print $* ; }"
}

SUNTANSHOME=../../main
SUN=$SUNTANSHOME/sun
maindatadir=rundata
datadir=data

. $SUNTANSHOME/Makefile.in

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np 1 $SUN"
fi

if [ -z $TRIANGLEHOME ] ; then
    echo Error: This example will not run without the triangle libraries.
    echo Make sure TRIANGLEHOME is set in $SUNTANSHOME/Makefile.in
    exit 1
fi

nmax=6
n=1
nsteps=10
tmax=1.0

cp -r $maindatadir $datadir

dtfile=$datadir/dt.txt
rm -f $dtfile

while (true)
do
  if [ $n -gt $nmax ] ; then
    break;
  fi

  np=$n
  dt=`calc "$tmax / $nsteps"`
  echo $dt $nsteps >> $dtfile
  echo "On run $n of $nmax (nsteps = $nsteps, dt = $dt)"
  
  if [ ! -d $datadir/t$np ] ; then
    mkdir $datadir/t$np
  fi

  'cp' $datadir/suntans.in $datadir/suntans.dat
  echo "dt $dt \#" >> $datadir/suntans.dat
  echo "nsteps $nsteps \#" >> $datadir/suntans.dat

  $EXEC -t -g -s -vvv --datadir=$datadir >& $datadir/t$np/run.out

  if [ `grep mergeArrays $maindatadir/suntans.in | grep -c 1` -eq 0 ] ; then
      'cp' $datadir/{fs.dat,q.dat,s.dat,u.dat,w.dat} $datadir/t$np
  else
      'cp' $datadir/fs.dat.0 $datadir/t$np/fs.dat
      'cp' $datadir/q.dat.0 $datadir/t$np/q.dat
      'cp' $datadir/s.dat.0 $datadir/t$np/s.dat
      'cp' $datadir/u.dat.0 $datadir/t$np/u.dat
      'cp' $datadir/w.dat.0 $datadir/t$np/w.dat
  fi
  
  nsteps=`calc "$nsteps * 2"`
  n=`expr $n + 1`
done

