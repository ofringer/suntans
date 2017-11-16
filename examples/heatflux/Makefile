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
ifneq ($(NETCDF4HOME),)
  NETCDFINCLUDE=-I$(NETCDF4HOME)/include
  NETCDFDEF = -DUSENETCDF
else
  ifeq ($(filter $(MAKECMDGOALS),clean clobber),)
    $(error Error: This test case will not run if NETCDF4HOME is empty in ../../main/Makefile.in)
  endif
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFSRC=
endif

LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = heatflux
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 1

all:	data

test:	data
	sh $(EXEC).sh $(NUMPROCS)

debug:	data
	sh $(EXEC)-debug.sh $(NUMPROCS)

data:	$(SUN)

.c.o:	
	$(LD) $(CFLAGS) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o state.o 
	cp initialization.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

plot:
	cp -r $(SUNTANSHOME)/cmaps .
	$(SUNTANSHOME)/sunplot -np $(NUMPROCS) --datadir=data


clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
