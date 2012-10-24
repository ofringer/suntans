#!/bin/sh

datadir=data

if [ ! -f tides.m ] ; then
    echo tides.m does not exist.
    exit 1;
fi

if which matlab >& /dev/null ; then
    if [ -f $datadir/tidexy.dat.0 ] ; then
	echo Creating $datadir/tidecomponents.dat'.*'
	matlab < tides.m >& /dev/null
    else
	echo Tide files $datadir/tidexy.dat'.*' do not exist...
    fi
else
    echo Matlab does not exist on this machine.  Need to create tide files on another machine.
fi