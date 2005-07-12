SUNTANSHOME=../..
include $(SUNTANSHOME)/Makefile.in

CC = $(MPIHOME)/bin/mpicc
LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = wetdry
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) -I$(MPIHOME) -I$(PARMETISHOME)/ParMETISLib
NUMPROCS = 1

all:	data

test1d:
	cp rundata/oned.dat rundata/pslg.dat
	cp initialization_1d.c initialization.c
	make test

test2d:
	cp rundata/twod.dat rundata/pslg.dat
	cp initialization_2d.c initialization.c
	make test

test:	data
	sh $(EXEC).sh $(NUMPROCS)

data:	$(SUN)

.c.o:	
	$(LD) -c $(INCLUDES) $*.c

$(SUN):	initialization.o 
	cp initialization.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
