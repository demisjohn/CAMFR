include ../../../make.inc

all: root croot traceroot polyroot minimum quadrature

root: FORCE
	cd root ; $(MAKE)

croot: FORCE
	cd croot ; $(MAKE)

traceroot: FORCE
	cd traceroot ; $(MAKE)

polyroot: FORCE
	cd polyroot ; $(MAKE)

minimum: FORCE
	cd minimum ; $(MAKE)

quadrature: FORCE
	cd quadrature ; $(MAKE)

FORCE:

clean:
	-rm *.a *.o core *~
	cd root ; $(MAKE) clean
	cd croot ; $(MAKE) clean
	cd traceroot ; $(MAKE) clean
	cd polyroot ; $(MAKE) clean
	cd minimum ; $(MAKE) clean
	cd quadrature ; $(MAKE) clean

