include ../../../make.inc

OBJS = slab.o slabwall.o slabdisp.o slabmode.o slaboverlap.o

all: generalslab.o isoslab moslab

generalslab.o: generalslab.cpp generalslab.h ../../defs.h ../../waveguide.h \
	../../expression.h
	$(CC) $(FLAGS) -c generalslab.cpp

isoslab: FORCE
	cd isoslab ; make

moslab: FORCE
	cd moslab ; make

FORCE:

clean:
	-rm *.a *.o core *~
	cd isoslab ; make clean
	cd moslab ; make clean
