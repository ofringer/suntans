#
# Makefile for SUNTANS
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#
# $Id: Makefile,v 1.3 2003-05-07 14:46:57 fringer Exp $
# $Log: not supported by cvs2svn $
# Revision 1.2  2003/04/29 00:23:56  fringer
# Added triangulate.o.  In order to compile sunplot, need to type
# make sunplot.  The default is to compile sun.
#
# Revision 1.1  2002/11/03 00:14:13  fringer
# Initial revision
#
#
SUNTANSHOME = /home/fringer/research/SUNTANS/suntans
PARMETISHOME = /home/fringer/research/SUNTANS/packages/ParMetis-2.0
TRIANGLEHOME = /home/fringer/research/SUNTANS/packages/triangle
BACKUP = /home/data

CC = /usr/local/mpich-1.2.1/bin/mpicc
OPTFLAGS = -O3
COPTIONS = 
LD = $(CC)
LIBDIR =  
INCLUDES = -I$(PARMETISHOME)/ParMETISLib
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
XLIB=-lX11 -lm
XINC=/usr/X11R6/include
XLIBDIR = /usr/X11R6/lib

LIBSDIR = -L$(PARMETISHOME) 
LIBS = -lparmetis -lmetis -lm

PARMETISLIB = $(PARMETISHOME)/libparmetis.a $(PARMETISHOME)/libmetis.a

EXEC = sun
PEXEC = sunplot

OBJS = mympi.o grid.o report.o util.o fileio.o phys.o suntans.o initialization.o memory.o
          
all:	$(EXEC)

.c.o:	
	$(LD) $(CFLAGS) $(LIBSDIR) -c $*.c

triangulate.o:	triangulate.c
	$(LD) -c triangulate.c -I$(TRIANGLEHOME) -I$(PARMETISHOME)

$(EXEC): $(OBJS) $(PARMETISLIB) suntans.h triangulate.o
	$(LD)  -o $@ $(OBJS) $(LIBSDIR) $(LIBS) triangulate.o $(TRIANGLEHOME)/triangle.o

sunplot.o:	sunplot.c suntans.h $(OBJS)
		$(LD) $(CFLAGS) -c sunplot.c

$(PEXEC): sunplot.o suntans.h
	$(LD) -o sunplot sunplot.o fileio.o $(LIBSDIR) $(LIBS) $(XLIB) -I$(XINC) -L$(XLIBDIR)

tarball: 
	-tar czvf $(BACKUP)/suntans.`date +%m%d%y`.tgz ../main

initialization.o:	initialization.h
grid.o:	grid.h initialization.h
mympi.o: mympi.h
phys.o:	phys.h initialization.h

clean:
	rm -f *.o ;\

superclean: 
	rm -f *~ \#*\# $(EXEC) $(PEXEC)
