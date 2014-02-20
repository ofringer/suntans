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

EXEC = tides-restart
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
JOBDIR = ./
NUMPROCS = 1
MAXRUNS = 10

all:	data

test:	data
	echo 1 > runnum.txt
	sh $(EXEC).sh $(NUMPROCS) $(MAXRUNS)

data:	$(SUN)

continue:
	sh continue.sh $(JOBDIR) $(NUMPROCS) $(MAXRUNS)

status:
	sh status.sh

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o state.o
	cp initialization.o boundaries.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)
	cp $(SUNTANSHOME)/sun .

clean:
	rm -f *.o sun

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) status.sh gmon.out data* rundata/*~ crash*.txt runnum.txt grid

