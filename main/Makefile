#
# Makefile for SUNTANS
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#
# $Id: Makefile,v 1.6 2004-05-29 20:27:21 fringer Exp $
# $Log: not supported by cvs2svn $
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
SUNTANSHOME = /home/fringer/research/SUNTANS/suntans
PARMETISHOME = /home/fringer/research/SUNTANS/packages/ParMetis-2.0
TRIANGLEHOME = /home/fringer/research/SUNTANS/packages/triangle
CLAPACKHOME = /home/fringer/research//CLAPACK
BACKUP = /home/data

CC = /usr/local/mpich-1.2.5.2/bin/mpicc
OPTFLAGS = -O3
#COPTIONS = -fsyntax-only  -pedantic  -pedantic-errors -w  -W  -Wall -Waggregate-return -Wcast-align  -Wcast-qual  -Wchar-subscripts  -Wcomment           -Wconversion  -Wno-deprecated-declarations -Wdisabled-optimization           -Wdiv-by-zero  -Werror -Wfloat-equal  -Wformat  -Wformat=2 -Wformat-nonliteral  -Wformat-security -Wimplicit  -Wimplicit-int -Wimplicit-function-declaration -Werror-implicit-function-declaration           -Wimport  -Winline -Wlong-long -Wmain  -Wmissing-braces -Wmissing-format-attribute  -Wmissing-noreturn -Wmultichar  -Wno-format-extra-args  -Wno-format-y2k -Wno-import  -Wpacked          -Wpadded -Wparentheses  -Wpointer-arith  -Wredundant-decls -Wreturn-type  -Wsequence-point  -Wshadow -Wsign-compare  -Wswitch           -Wsystem-headers -Wtrigraphs  -Wundef  -Wuninitialized -Wunknown-pragmas  -Wunreachable-code -Wunused  -Wunused-function           -Wunused-label  -Wunused-parameter -Wunused-value  -Wunused-variable  -Wwrite-strings
LD = $(CC) -p -pg
LIBDIR =  
INCLUDES = -I$(PARMETISHOME)/ParMETISLib
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
XLIB=-lX11 -lm
XINC=/usr/X11R6/include
XLIBDIR = /usr/X11R6/lib

MATHLIB = -lm
LIBSDIR = -L$(PARMETISHOME) 
LIBS = -lparmetis -lmetis $(MATHLIB)
LAPACKLIB = $(CLAPACKHOME)/lapack_LINUX.a $(CLAPACKHOME)/blas_LINUX.a
LAPACKINCLUDES = -I$(CLAPACKHOME) -I$(CLAPACKHOME)/F2CLIBS

PARMETISLIB = $(PARMETISHOME)/libparmetis.a $(PARMETISHOME)/libmetis.a

EXEC = sun
PEXEC = sunplot

OBJS = mympi.o grid.o report.o util.o fileio.o phys.o suntans.o initialization.o memory.o boundaries.o
          
all:	$(EXEC)

.c.o:	
	$(LD) $(CFLAGS) $(LIBSDIR) -c $*.c

kriging.o:	kriging.c
	$(LD) -c kriging.c $(MATHLIB) $(LAPACKINCLUDES)

triangulate.o:	triangulate.c
	$(LD) -c triangulate.c -I$(TRIANGLEHOME) -I$(PARMETISHOME)

#$(EXEC): $(OBJS) $(PARMETISLIB) suntans.h triangulate.o kriging.o
#	$(LD)  -o $@ $(OBJS) $(LIBSDIR) $(LIBS) triangulate.o kriging.o \
#        $(TRIANGLEHOME)/triangle.o $(LAPACKLIB) $(LAPACKINCLUDES)

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
	rm -f *~ \#*\# $(EXEC) $(PEXEC)
