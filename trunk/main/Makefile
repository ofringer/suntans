#
# Makefile for SUNTANS
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#
# $Id: Makefile,v 1.1 2002-11-03 00:14:13 fringer Exp $
# $Log: not supported by cvs2svn $
#
SUNTANSHOME = /home/fringer/research/SUNTANS/suntans
PARMETISHOME = /home/fringer/research/SUNTANS/packages/ParMetis-2.0

CC = /usr/local/mpich-1.2.1/bin/mpicc
OPTFLAGS = -O3
COPTIONS = 
LD = $(CC)
LIBDIR =  
INCLUDES = -I$(PARMETISHOME)/ParMETISLib
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)


LIBSDIR = -L$(PARMETISHOME)
LIBS = -lparmetis -lmetis -lm
PARMETISLIB = $(PARMETISHOME)/libparmetis.a $(PARMETISHOME)/libmetis.a

EXEC = sun

OBJS = mympi.o grid.o report.o util.o fileio.o phys.o suntans.o
          
default:	$(EXEC)

.c.o:	
	$(LD) $(CFLAGS) $(LIBSDIR) -c $*.c

$(EXEC): $(OBJS) $(PARMETISLIB) main.h
	$(LD)  -o $@ $(OBJS) $(LIBSDIR) $(LIBS)

tarball: 
	-tar czvf pgrid.tgz *.c *.h Makefile

grid.o:	grid.h
mympi.o: mympi.h
phys.o:	phys.h

clean:
	rm -f *.o ;\

superclean: 
	rm -f *~ #*# $(EXEC)
