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
  PARMETISINC = -I$(PARMETISHOME)/include/parmetis
endif

LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = accuracy
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 1

all:	data

test:	data
	bash accuracy.sh
	./$(EXEC)

data:	$(EXEC) $(SUN)

plot:
	cp -fr $(SUNTANSHOME)/cmaps .
	$(SUNTANSHOME)/sunplot -np $(NUMPROCS) --datadir=data

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o state.o 
	cp initialization.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

$(EXEC):	accuracy.o
	$(LD) -o $@ $(MATHLIB) accuracy.o

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
