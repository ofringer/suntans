########################################################################
# 
# File: release.sh
# 
# This shell script puts together the files required to create a 
# tarball for a release
#
# $Id: release.sh,v 1.1 2004-06-17 02:56:00 fringer Exp $
# $Log: not supported by cvs2svn $
#
########################################################################
#!/bin/sh

lastversion=0.0
version=1.0
dir="suntans-$version"

files="boundaries.c initialization.c mympi.c fileio.c phys.c sunplot.c triangulate.c grid.c\
       memory.c report.c suntans.c util.c boundaries.h grid.h memory.h phys.h suntans.h util.h\
       fileio.h initialization.h mympi.h report.h triangulate.h Makefile cmaps suntans.dat"

if [ ! -f $dir.tgz ] ; then
  mkdir $dir
  cp -r $files $dir

  rm -rf $dir/cmaps/CVS
  tar czvf $dir.tgz $dir
  rm -rf $dir
else 
  echo "Tarball $dir.tgz already exists!"
  echo "Exiting..."
  exit 1;
fi
