SUNTANSHOME=../../main
include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
# redefine gravity to make code rigid lid
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
else
  CC = gcc
  MPIDEF = -DNOMPI
  MPIINC = 
endif

ifneq ($(PARMETISHOME),)
  PARMETISINC = -I$(PARMETISHOME)/include/parmetis
endif

LD = $(CC) 
CFLAGS = -O2
MATHLIB = -lm

EXEC = cavity
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 1

all:	$(SUN)

test:	$(SUN)
	sh $(EXEC).sh $(NUMPROCS)

test-quad:	testXY-quad

testXY:	dataXY test

dataXY: FORCE
	cp boundariesXY.c boundaries.c 
	cp initialization-constant-dz.c initialization.c
	cp data-files/cellsXY.dat rundata/cells.dat
	cp data-files/edgesXY.dat rundata/edges.dat
	cp data-files/pointsXY.dat rundata/points.dat
	cp data-files/suntansXY.dat rundata/suntans.dat

testXZ:	dataXZ test

dataXZ: FORCE
	cp boundariesXZ.c boundaries.c
	cp initialization-constant-dz.c initialization.c
	cp data-files/cellsXZ.dat rundata/cells.dat
	cp data-files/edgesXZ.dat rundata/edges.dat
	cp data-files/pointsXZ.dat rundata/points.dat
	cp data-files/suntansXZ.dat rundata/suntans.dat

testXY-quad:	dataXY-quad test

dataXY-quad: FORCE
	cp boundariesXY.c boundaries.c
	cp initialization-constant-dz.c initialization.c
	cp data-files/cellsXY-quad.dat rundata/cells.dat
	cp data-files/edgesXY-quad.dat rundata/edges.dat
	cp data-files/pointsXY-quad.dat rundata/points.dat
	cp data-files/suntansXY-quad.dat rundata/suntans.dat

testXZ-quad:	dataXZ-quad test

dataXZ-quad: FORCE
	cp boundariesXZ.c boundaries.c
	cp initialization-chebychev-dz.c initialization.c
	cp data-files/cellsXZ-quad.dat rundata/cells.dat
	cp data-files/edgesXZ-quad.dat rundata/edges.dat
	cp data-files/pointsXZ-quad.dat rundata/points.dat
	cp data-files/suntansXZ-quad.dat rundata/suntans.dat

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o
	cp initialization.o boundaries.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

# Empty target makes it always happen
FORCE:

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data cmaps
	make dataXY
