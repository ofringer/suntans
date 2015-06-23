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

EXEC=headland
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC) $(NETCDFINCLUDE) 
DEFS = $(MPIDEF) $(NETCDFDEF)
NUMPROCS = 8 
datadir=data

all:	data

test:	data
	sh $(EXEC).sh $(NUMPROCS)

quad:	data
	python scripts/twod_grid_headland.py
	cp grids/quad/*.dat rundata/
	python scripts/suntans_driver_Headland.py

tri:	data
	python scripts/run_tom_headland.py
	cp grids/tri/*.dat rundata/
	python scripts/suntans_driver_Headland.py

hex:	data
	python scripts/run_tom_headland.py
	python scripts/convert_grid2hex.py
	cp grids/hex/*.dat rundata/
	python scripts/suntans_driver_Headland.py

debug:	data
	cp -r rundata data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) xterm -e gdb -command=gdbcommands $(SUN)

valgrind:	data
	cp -r rundata data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(SUN) -g -vv --datadir=$(datadir)
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) valgrind --tool=memcheck --leak-check=yes ./$(SUN) -s -vvv --datadir=$(datadir)


	
data:	initialization.o state.o 
	cp initialization.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	make -C $(SUNTANSHOME)

clean:
	rm -f *.o
	#rm -f $(SUNTANSHOME)/*.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) data gmon.out rundata/*
