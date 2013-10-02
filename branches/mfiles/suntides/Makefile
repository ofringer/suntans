#
# Makefile to unpack data for use with suntides.
#
# Oliver Fringer
# Stanford University
# 18 Sep 07
#
all:	data

data:
	tar xzvf data.tgz

clean:
	rm -f *~ \#*\# data/tidal_files/tidecomponents.dat.*

clobber: clean
	rm -rf data
