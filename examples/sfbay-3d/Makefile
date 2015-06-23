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
NUMPROCS = 16

all:	$(SUN)

test:
	@echo Please use \"make test2d\" or \"make test3d\" or \"make netcdf-test2d\" or \"make netcdf-test3d\".

test2d:	data
	sh $(EXEC).sh $(NUMPROCS) 2

test3d:	data
	sh $(EXEC).sh $(NUMPROCS) 3

netcdf-test2d:	datanc
	sh $(EXEC)-netcdf.sh $(NUMPROCS) 2

netcdf-test3d:	datanc
	sh $(EXEC)-netcdf.sh $(NUMPROCS) 3


data:	$(SUN)
	cp boundaries-no-nc.c boundaries.c


datanc: $(SUN)
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

