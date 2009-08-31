all:	test

test:
	sh testall.sh

clean:
	rm -f *~ *.out

clobber:	clean
	make -C accuracy clobber
	make -C boundaries clobber
	make -C channel clobber
	make -C estuary clobber
	make -C iwaves clobber
	make -C lockexchange clobber
	make -C wetdry clobber
	make -C windstress clobber
	make -C cylinder clobber
	make -C tides clobber
	make -C tides-restart clobber

