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

EXEC = accuracy
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
NUMPROCS = 1

all:	data

test:	data
	sh accuracy.sh
	./$(EXEC)

data:	$(EXEC) $(SUN)

.c.o:	
	$(LD) -c $(INCLUDES) $*.c

$(SUN):	initialization.o
	cp initialization.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

$(EXEC):	accuracy.o
	$(LD) -o $@ $(MATHLIB) accuracy.o

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
