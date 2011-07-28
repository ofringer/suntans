#
# Makefile to create tar archive of suntans.
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#

all:	../suntans.tgz

../suntans.tgz:	
	tar -C .. --exclude=*CVS* --exclude=*pac* --exclude=*m_map* --exclude=*big3dgrid* --exclude="*papers*" --exclude="*~*" --exclude="\#*\#" --exclude="*.o" --exclude="*.exe" -czvf ../suntans.tgz suntans

clean:
	rm -f *~ ../suntans.tgz

