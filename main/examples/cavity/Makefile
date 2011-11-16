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
  PARMETISINC = -I$(PARMETISHOME)/ParMETISLib
endif

LD = $(CC) 
CFLAGS = -O2
MATHLIB = -lm

EXEC = cavity
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS =1

test:	testXY testXZ

testXY:	
	make clobber
	rm -rf dataXY
	tar zxvf cavityXY.tgz
	make data
	sh $(EXEC)XY.sh $(NUMPROCS)

testXZ:	
	make clobber 
	rm -rf dataXZ
	tar zxvf cavityXZ.tgz
	make data
	sh $(EXEC)XZ.sh $(NUMPROCS)

data:	$(SUN)

plotXY:
		cp -r $(SUNTANSHOME)/cmaps .
		$(SUNTANSHOME)/sunplot -np $(NUMPROCS) --datadir=dataXY

plotXZ:
		cp -r $(SUNTANSHOME)/cmaps .
		$(SUNTANSHOME)/sunplot -np $(NUMPROCS) --datadir=dataXZ

.c.o:	
	$(LD) -c $(CFLAGS) $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o
	cp initialization.o boundaries.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/ cmaps initialization.c boundaries.c cavityXY.sh cavityXZ.sh *.eps
