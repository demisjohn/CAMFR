include ../../make.inc

all: planar slab circ

planar: FORCE
	cd planar ; $(MAKE)

slab: FORCE
	cd slab ; $(MAKE)

circ: FORCE
	cd circ ; $(MAKE)

FORCE:

clean:
	-rm *.a *.o *.so *.py* core *~
	cd planar ; $(MAKE) clean
	cd slab ; $(MAKE) clean
	cd circ ; $(MAKE) clean
