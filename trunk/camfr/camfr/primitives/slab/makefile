include ../../../make.inc

OBJS = slab.o slabwall.o slabdisp.o slabmode.o slaboverlap.o

all: generalslab.o slabmatrixcache.o isoslab moslab

generalslab.o: generalslab.cpp generalslab.h ../../defs.h ../../waveguide.h \
	../../expression.h
	$(CC) $(FLAGS) -c generalslab.cpp

slabmatrixcache.o: slabmatrixcache.cpp slabmatrixcache.h ../../defs.h \
	../../waveguide.h ../../expression.h
	$(CC) $(FLAGS) -c slabmatrixcache.cpp

isoslab: FORCE
	cd isoslab ; make

moslab: FORCE
	cd moslab ; make

FORCE:

clean:
	-rm *.a *.o core *~
	cd isoslab ; make clean
	cd moslab ; make clean
