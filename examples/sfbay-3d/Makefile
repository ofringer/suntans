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

EXEC = sfbay
OBJS = 
SUN = $(SUNTANSHOME)/sun
SUNNC = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 4

all:	$(SUN)

test:	test2d

test2d:	data $(SUN)
	sh $(EXEC).sh $(NUMPROCS) 2

test3d:	data $(SUN)
	sh $(EXEC).sh $(NUMPROCS) 3

netcdf-test2d:	datanc $(SUN)
	sh $(EXEC)-netcdf.sh $(NUMPROCS) 2

netcdf-test3d:	datanc $(SUN)
	sh $(EXEC)-netcdf.sh $(NUMPROCS) 3

data:
	cp boundaries-no-nc.c boundaries.c

datanc:
	cp boundaries-nc.c boundaries.c

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o state.o
	cp initialization.o boundaries.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~ 

