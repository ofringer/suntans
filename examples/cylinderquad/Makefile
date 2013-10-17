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

all:	data

test:   data
	cp rundata/trigridnew/cells.dat rundata/cells.dat
	cp rundata/trigridnew/edges.dat rundata/edges.dat
	cp rundata/trigridnew/points.dat rundata/points.dat
	#make data
	sh $(EXEC).sh $(NUMPROCS)

testQUAD:   data	
	cp rundata/hybridgridnew/cells.dat rundata/cells.dat
	cp rundata/hybridgridnew/edges.dat rundata/edges.dat
	cp rundata/hybridgridnew/points.dat rundata/points.dat
	#make data
	sh $(EXEC).sh $(NUMPROCS)

debug:	data
	/opt/intel/openmpi/bin/mpirun -np $(NUMPROCS) xterm -e gdb $(SUN)
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
	rm -f $(SUNTANSHOME)/*.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
