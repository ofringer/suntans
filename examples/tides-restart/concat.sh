#!/bin/bash
#
# Short script to concatinate the output data from multiple runs
# into one file.
#
MAXRUNS=$1

if [ ! -d data ] ; then
    mkdir data
fi

cat data?/fs.dat.prof > data/fs.dat.prof
cat data?/u.dat.prof > data/u.dat.prof

if [ $MAXRUNS -gt 9 ] ; then
    cat data??/fs.dat.prof >> data/fs.dat.prof
    cat data??/u.dat.prof >> data/u.dat.prof
fi

if [ $MAXRUNS -gt 99 ] ; then
    cat data???/fs.dat.prof >> data/fs.dat.prof
    cat data???/u.dat.prof >> data/u.dat.prof
fi


