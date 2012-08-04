SUNTANSHOME=../..
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
NUMPROCS =8

test:	testXY testXZ

testXY:	
	make -C ../../ clobber sunplot
	rm -rf dataXY
	cp boundariesXY.c boundaries.c
	cp rundata/cellsXY.dat rundata/cells.dat
	cp rundata/edgesXY.dat rundata/edges.dat
	cp rundata/pointsXY.dat rundata/points.dat
	cp rundata/suntansXY.dat rundata/suntans.dat
	make data
	sh $(EXEC)XY.sh $(NUMPROCS)

testXZ:	
	make -C ../../ clobber sunplot
	rm -rf dataXZ
	cp boundariesXZ.c boundaries.c
	cp rundata/cellsXZ.dat rundata/cells.dat
	cp rundata/edgesXZ.dat rundata/edges.dat
	cp rundata/pointsXZ.dat rundata/points.dat
	cp rundata/suntansXZ.dat rundata/suntans.dat
	make data
	sh $(EXEC)XZ.sh $(NUMPROCS)

data:	$(SUN)

plotXY:
		cp -rf $(SUNTANSHOME)/cmaps .
		$(SUNTANSHOME)/sunplot -np $(NUMPROCS) --datadir=dataXY

plotXZ:
		cp -rf $(SUNTANSHOME)/cmaps .
		$(SUNTANSHOME)/sunplot -np $(NUMPROCS) --datadir=dataXZ

.c.o:	
	$(LD) -c $(CFLAGS) $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o
	cp initialization.o boundaries.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data cmaps boundaries.c rundata/cells.dat rundata/points.dat rundata/edges.dat rundata/suntans.dat
