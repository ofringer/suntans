SUNTANSHOME=../../main
include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
else
  CC = icc
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
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFSRC=
endif

LD = $(CC) 
CFLAGS =
MATHLIB = -lm

EXEC=sun
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC) $(NETCDFINCLUDE) 
DEFS = $(MPIDEF) $(NETCDFDEF)
NUMPROCS = 16
datadir=rundata

all:	data

test:	data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(EXEC) -g -vvv --datadir=$(datadir)  
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(EXEC) -s -vvv --datadir=$(datadir)  

debug:	data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) xterm -e gdb -command=gdbcommands $(EXEC)

valgrind:	data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(EXEC) -g -vv --datadir=$(datadir)
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) valgrind ./$(EXEC) -s -vvv --datadir=$(datadir)


data:	state.o 
	cp state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)
	mv $(SUNTANSHOME)/$(EXEC) ./



.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	make -C $(SUNTANSHOME)

clean:
	rm -f *.o
	rm -f $(SUNTANSHOME)/*.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
