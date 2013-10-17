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
NUMPROCS = 4
datadir=rundata

all:	data

test:	data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(EXEC) -g -s -vv --datadir=$(datadir)  

	#python ~/suntans/SourceForge/python/SUNTANS/joinsun.py -f Estuary_SUNTANS.nc -p $(datadir) -t 168 -n $(NUMPROCS)

inputs:
	sh estuarynetcdf.sh

join:   
	python ~/suntans/SourceForge/python/SUNTANS/joinsun.py -f Estuary_SUNTANS.nc -p $(datadir) -t 168 -n $(NUMPROCS)
	
data:	initialization.o state.o 
	cp  initialization.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)
	mv $(SUNTANSHOME)/$(EXEC) ./

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	make -C $(SUNTANSHOME)

clean:
	rm -f *.o
	rm -f $(SUNTANSHOME)/*.o

clobber: clean
	#rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data 
