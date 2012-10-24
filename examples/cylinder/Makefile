SUNTANSHOME=../..
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

all:	data

test:	data
	sh $(EXEC).sh $(NUMPROCS)

data:	$(SUN)

.c.o:	
	$(LD) $(CFLAGS) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o state.o sources.o
	cp initialization.o boundaries.o state.o sources.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

plot:
	cp -r $(SUNTANSHOME)/cmaps .
	$(SUNTANSHOME)/sunplot -np $(NUMPROCS) --datadir=data


clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
