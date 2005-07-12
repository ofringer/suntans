MPIHOME = /usr/local/mpich-1.2.5.2
SUNTANSHOME = ../..
PARMETISHOME = /home/fringer/research/SUNTANS/packages/ParMetis-2.0

CC = $(MPIHOME)/bin/mpicc
LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = accuracy
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) -I$(MPIHOME) -I$(PARMETISHOME)/ParMETISLib
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
