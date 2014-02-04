SUNTANSHOME=../../main
include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
else
  CC = gcc
  MPIDEF = -DNOMPI
  MPIINC = 
endif

ifneq ($(PARMETISHOME),)
  PARMETISINC = -I$(PARMETISHOME)/ParMETISLib
endif

LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = cylinder
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 1

all:	$(SUN)

test:	$(SUN)
	sh $(EXEC).sh $(NUMPROCS)

test-tri: data-tri test

test-quad:	data-quad test

data-tri:	FORCE
	cp data-files/cells-tri.dat rundata/cells.dat
	cp data-files/edges-tri.dat rundata/edges.dat
	cp data-files/points-tri.dat rundata/points.dat
	cp data-files/suntans-tri.dat rundata/suntans.dat

data-quad: FORCE
	cp data-files/cells-quad.dat rundata/cells.dat
	cp data-files/edges-quad.dat rundata/edges.dat
	cp data-files/points-quad.dat rundata/points.dat
	cp data-files/suntans-quad.dat rundata/suntans.dat

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o state.o
	cp initialization.o boundaries.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

# Empty target makes it always happen
FORCE:

clean:
	rm -f *.o

clobber: clean data-tri
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~



