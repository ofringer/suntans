SUNTANSHOME=../..
include $(SUNTANSHOME)/Makefile.in

CC = $(MPIHOME)/bin/mpicc
LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = gridding
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) -I$(MPIHOME) -I$(PARMETISHOME)/ParMETISLib
NUMPROCS = 1

all:	data

test:	data
	sh $(EXEC).sh $(NUMPROCS)

data:	$(SUN)

.c.o:	
	$(LD) -c $(INCLUDES) $*.c

$(SUN):	initialization.o boundaries.o
	cp initialization.o boundaries.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
