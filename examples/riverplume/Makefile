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
  ifeq ($(filter $(MAKECMDGOALS),clean clobber),)
    $(error Error: This test case will not run if NETCDF4HOME is empty in ../../main/Makefile.in)
  endif
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFSRC=
endif

LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC=riverplume
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC) $(NETCDFINCLUDE) 
DEFS = $(MPIDEF) $(NETCDFDEF)
NUMPROCS = 2 
datadir=data
plotdir=plots/validation

all:	data

test:	data
	sh $(EXEC).sh $(NUMPROCS)

movies: 
	python scripts/buildMovies.py

quad:	data
	cp grids/quad/*.dat rundata/
	python scripts/suntans_driver_Plume.py

tri:	data
	cp grids/tri/*.dat rundata/
	python scripts/suntans_driver_Plume.py

hybrid:	data
	cp grids/hybrid/*.dat rundata/
	python scripts/suntans_driver_Plume.py

gmsh:	data
	cp grids/gmsh/*.dat rundata/
	python scripts/suntans_driver_Plume.py

hex:	data
	cp grids/hex/*.dat rundata/
	python scripts/suntans_driver_Plume.py

debug:	data
	cp -r rundata data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) xterm -e gdb -command=gdbcommands $(SUN)

valgrind:	data
	cp -r rundata data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(SUN) -g -vv --datadir=$(datadir)
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) valgrind --tool=memcheck --leak-check=yes ./$(SUN) -s -vvv --datadir=$(datadir)


join:  
	python ~/suntans/SourceForge/python/SUNTANS/joinsun.py -f Plume.nc -p $(datadir) -t -1 -n $(NUMPROCS)

	
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
	rm -rf *~ \#*\# PI* $(EXEC) sunjoin data gmon.out rundata/*~
