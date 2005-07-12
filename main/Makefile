#
# Makefile for SUNTANS
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#
# Need to define appropriate directories in Makefile.in first!
#
include Makefile.in

CC = $(MPIHOME)/bin/mpicc
OPTFLAGS = 
LIBDIR = -L$(PARMETISHOME) -L$(TRIANGLEHOME)
INCLUDES = -I$(PARMETISHOME)/ParMETISLib 

XLDFLAGS=-lX11 -lm
XINC=/usr/X11R6/include
XLIBDIR = /usr/X11R6/lib

TRIINCLUDE = -I$(TRIANGLEHOME)
TRIANGLELIB = $(TRIANGLEHOME)/triangle.o
PARMETISLIB = $(PARMETISHOME)/libparmetis.a $(PARMETISHOME)/libmetis.a

LD = $(CC) 
LDFLAGS = -lparmetis -lmetis -lm $(LIBDIR) $(PARMETISLIB) $(TRIANGLELIB)
CFLAGS = $(OPTFLAGS) $(INCLUDES) 

EXEC = sun
PEXEC = sunplot

DEPFLAGS = -Y

SRCS = mympi.c grid.c report.c util.c fileio.c phys.c suntans.c initialization.c memory.c \
	turbulence.c boundaries.c triangulate.c
OBJS = $(SRCS:.c=.o)

PLOTSRCS = sunplot.c fileio.c
PLOTOBJS = $(PLOTSRCS:.c=.o)

all:	$(EXEC)

.c.o:	
	$(CC) $(CFLAGS) -c $*.c

triangulate.o:	triangulate.c
	$(CC) -I$(TRIANGLEHOME) -I$(PARMETISHOME) -c triangulate.c

sunplot.o:	sunplot.c
	$(CC) -I$(XINC) -c $*.c

$(EXEC): $(OBJS) 
	$(LD)  -o $@ $(OBJS) $(LDFLAGS)

$(PEXEC): sunplot.o
	$(LD) -o sunplot $(PLOTOBJS) $(XLDFLAGS) -L$(XLIBDIR)

depend: 
	makedepend $(DEPFLAGS) -- $(SRCS) $(PLOTSRCS) &> /dev/null

clean:
	rm -f *.o

clobber:	clean
	rm -f *~ \#*\# PI* $(EXEC) $(PEXEC) $(DEPFILE)


# DO NOT DELETE THIS LINE - Dependencies are appended after it.

mympi.o: mympi.h suntans.h fileio.h
grid.o: grid.h suntans.h fileio.h mympi.h util.h initialization.h memory.h
grid.o: triangulate.h report.h
report.o: report.h grid.h suntans.h fileio.h mympi.h
util.o: grid.h suntans.h fileio.h mympi.h util.h
fileio.o: fileio.h
phys.o: suntans.h phys.h grid.h fileio.h mympi.h util.h initialization.h
phys.o: memory.h turbulence.h boundaries.h
suntans.o: suntans.h mympi.h fileio.h grid.h phys.h report.h
initialization.o: fileio.h suntans.h initialization.h
memory.o: memory.h
turbulence.o: phys.h suntans.h grid.h fileio.h mympi.h util.h turbulence.h
turbulence.o: boundaries.h
boundaries.o: boundaries.h suntans.h phys.h grid.h fileio.h mympi.h
triangulate.o: suntans.h mympi.h fileio.h grid.h triangulate.h
sunplot.o: suntans.h fileio.h
fileio.o: fileio.h
