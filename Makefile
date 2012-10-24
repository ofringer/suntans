#
# Makefile for SUNTANS
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#
# Need to define appropriate directories in Makefile.in first!
#
# Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior University. All Rights Reserved.
#
include Makefile.in

ifneq ($(MPIHOME),)
  # For the SGI Altix
  # CC = icc
  CC = $(MPIHOME)/bin/mpicc
  MPIFILE = 
  MPIDEF = 
else
  CC = gcc
  MPIFILE = no-mpi.c
  MPIDEF = -DNOMPI
endif
OPTFLAGS = -O2

XLDFLAGS=-lX11 -lm
XINC=/usr/include
XLIBDIR = /usr/lib64

ifneq ($(TRIANGLEHOME),)
  TRIANGLEINCLUDE = -I$(TRIANGLEHOME)
  TRIANGLELIB = $(TRIANGLEHOME)/triangle.o
  TRIANGLELIBDIR = -L$(TRIANGLEHOME)
  TRIANGLELD = 
  TRIANGLESRC = triangulate.c
else
  TRIANGLEINCLUDE =
  TRIANGLELIB =
  TRIANGLELIBDIR =
  TRIANGLELD =
  TRIANGLESRC = triangulate-notriangle.c
endif	

ifneq ($(PARMETISHOME),)
  PARMETISINCLUDE = -I$(PARMETISHOME)/ParMETISLib
  PARMETISLIB = $(PARMETISHOME)/libparmetis.a $(PARMETISHOME)/libmetis.a
  PARMETISLIBDIR = -L$(PARMETISHOME)
  PARMETISLD = -lparmetis -lmetis 
  PARMETISSRC = partition.c
else
  PARMETISINCLUDE =
  PARMETISLIB =
  PARMETISLIBDIR =
  PARMETISLD =
  PARMETISSRC = partition-noparmetis.c
endif

# For the Altix
#LD = $(CC) -lmpi
LD = $(CC)
LIBS = $(PARMETISLIB) $(TRIANGLELIB)
LIBDIR = $(PARMETISLIBDIR) $(TRIANGLELIBDIR)
LDFLAGS = -lm $(LIBDIR) $(LIBS)
INCLUDES = $(PARMETISINCLUDE) $(TRIANGLEINCLUDE)
DEFINES = $(MPIDEF)
CFLAGS = $(OPTFLAGS) $(INCLUDES) $(DEFINES)

EXEC = sun
PEXEC = sunplot

DEPFLAGS = -Y

SRCS = 	mympi.c grid.c report.c util.c fileio.c phys.c suntans.c initialization.c memory.c \
	turbulence.c boundaries.c check.c scalars.c tvd.c timer.c profiles.c state.c tides.c \
	sources.c diffusion.c $(TRIANGLESRC) $(PARMETISSRC) $(MPIFILE)
OBJS = $(SRCS:.c=.o)

PLOTSRCS = sunplot.c fileio.c
PLOTOBJS = $(PLOTSRCS:.c=.o)

all:	$(EXEC)

.c.o:	
	$(CC) $(CFLAGS) -c $*.c

$(EXEC): $(OBJS) 
	$(LD)  -o $@ $(OBJS) $(LDFLAGS)

$(PEXEC): $(PLOTOBJS)
	$(LD) -o sunplot $(PLOTOBJS) $(XLDFLAGS) -L$(XLIBDIR)

depend: 
	makedepend $(DEPFLAGS) -- $(SRCS) $(PLOTSRCS) &> /dev/null

clean:
	rm -f *.o

clobber:	clean
	rm -f *~ \#*\# PI* $(EXEC) $(PEXEC) $(DEPFILE)
	make -C examples clobber

# DO NOT DELETE THIS LINE - Dependencies are appended after it.

mympi.o: mympi.h suntans.h fileio.h
grid.o: grid.h suntans.h fileio.h mympi.h partition.h util.h initialization.h
grid.o: memory.h triangulate.h report.h timer.h
report.o: report.h mympi.h suntans.h fileio.h grid.h
util.o: grid.h suntans.h fileio.h mympi.h util.h
fileio.o: fileio.h defaults.h suntans.h
phys.o: suntans.h phys.h grid.h fileio.h mympi.h util.h initialization.h
phys.o: memory.h turbulence.h boundaries.h check.h scalars.h timer.h
phys.o: profiles.h state.h diffusion.h sources.h
suntans.o: suntans.h mympi.h fileio.h grid.h phys.h report.h
initialization.o: fileio.h suntans.h initialization.h
memory.o: memory.h
turbulence.o: phys.h suntans.h grid.h fileio.h mympi.h util.h turbulence.h
turbulence.o: boundaries.h scalars.h
boundaries.o: boundaries.h suntans.h phys.h grid.h fileio.h mympi.h
check.o: check.h grid.h suntans.h fileio.h mympi.h phys.h timer.h memory.h
scalars.o: scalars.h suntans.h grid.h fileio.h mympi.h phys.h util.h tvd.h
scalars.o: initialization.h
tvd.o: suntans.h phys.h grid.h fileio.h mympi.h tvd.h util.h
timer.o: mympi.h suntans.h fileio.h timer.h
profiles.o: util.h grid.h suntans.h fileio.h mympi.h memory.h phys.h
profiles.o: profiles.h
state.o: state.h grid.h suntans.h fileio.h mympi.h phys.h
tides.o: suntans.h mympi.h fileio.h grid.h tides.h memory.h
sources.o: phys.h suntans.h grid.h fileio.h mympi.h sources.h memory.h
diffusion.o: diffusion.h grid.h suntans.h fileio.h mympi.h phys.h util.h
triangulate-notriangle.o: suntans.h mympi.h fileio.h grid.h
partition-noparmetis.o: suntans.h partition.h grid.h fileio.h mympi.h
no-mpi.o: suntans.h no-mpi.h
sunplot.o: suntans.h fileio.h
fileio.o: fileio.h defaults.h suntans.h
