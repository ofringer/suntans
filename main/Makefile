#
# Makefile for SUNTANS
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#
# $Id: Makefile,v 1.9 2004-09-12 22:12:00 fringer Exp $
# $Log: not supported by cvs2svn $
# Revision 1.8  2004/06/10 08:29:11  fringer
# Added PI* to the superclean command.
#
# Revision 1.7  2004/06/10 08:20:51  fringer
# Removed lapack and added MPIHOME to specify mpi main directory.
#
# Revision 1.6  2004/05/29 20:27:21  fringer
#  Revision before converting to CVS.
#
# Revision 1.5  2004/04/20 16:48:33  fringer
# 4/20/04
#
# Revision 1.4  2003/05/12 00:18:26  fringer
# Added the blas and lapack libraries to compile kriging.c.  These
# must be compiled from the CLAPACK libraries obtained from
# netlib.org
#
# Revision 1.3  2003/05/07 14:46:57  fringer
# Added PEXEC = suntans and updated superclean and tarball so that
# tarball creates a tarball in the directory $(BACKUP).
#
# Revision 1.2  2003/04/29 00:23:56  fringer
# Added triangulate.o.  In order to compile sunplot, need to type
# make sunplot.  The default is to compile sun.
#
# Revision 1.1  2002/11/03 00:14:13  fringer
# Initial revision
#
#
MPIHOME = /usr/local/mpich-1.2.5.2
SUNTANSHOME = /home/fringer/research/SUNTANS/suntans
PARMETISHOME = /home/fringer/research/SUNTANS/packages/ParMetis-2.0
TRIANGLEHOME = /home/fringer/research/SUNTANS/packages/triangle
BACKUP = /home/data

CC = $(MPIHOME)/bin/mpicc
OPTFLAGS = 
#LD = $(CC) -p -pg -fbounds-checking -B/usr/local/gcc-3.3.3
LD = $(CC) -p -pg 
LIBDIR =  
INCLUDES = -I$(PARMETISHOME)/ParMETISLib
CFLAGS = $(OPTFLAGS) $(INCLUDES) 
XLIB=-lX11 -lm
XINC=/usr/X11R6/include
XLIBDIR = /usr/X11R6/lib

MATHLIB = -lm
LIBSDIR = -L$(PARMETISHOME) 
LIBS = -lparmetis -lmetis $(MATHLIB)

PARMETISLIB = $(PARMETISHOME)/libparmetis.a $(PARMETISHOME)/libmetis.a

EXEC = sun
PEXEC = sunplot

OBJS = mympi.o grid.o report.o util.o fileio.o phys.o suntans.o initialization.o memory.o boundaries.o
          
all:	$(EXEC)

.c.o:	
	$(LD) $(CFLAGS) $(LIBSDIR) -c $*.c

kriging.o:	kriging.c
	$(LD) -c kriging.c $(MATHLIB)

triangulate.o:	triangulate.c
	$(LD) -c triangulate.c -I$(TRIANGLEHOME) -I$(PARMETISHOME)

$(EXEC): $(OBJS) $(PARMETISLIB) suntans.h triangulate.o 
	$(LD)  -o $@ $(OBJS) $(LIBSDIR) $(LIBS) triangulate.o  \
        $(TRIANGLEHOME)/triangle.o 

sunplot.o:	sunplot.c suntans.h $(OBJS)
		$(LD) $(CFLAGS) -c sunplot.c

$(PEXEC): sunplot.o suntans.h
	$(LD) -o sunplot sunplot.o fileio.o $(LIBSDIR) $(LIBS) $(XLIB) -I$(XINC) -L$(XLIBDIR)

tarball: 
	-tar czvf $(BACKUP)/suntans.`date +%m%d%y`.tgz -C .. main

initialization.o:	initialization.h 
grid.o:	grid.h initialization.h
mympi.o: mympi.h
phys.o:	phys.h initialization.h

clean:
	rm -f *.o ;\

superclean: 
	rm -f *~ \#*\# PI* $(EXEC) $(PEXEC)
