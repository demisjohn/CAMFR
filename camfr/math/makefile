include ../../make.inc

all: bessel linalg calculus

bessel: FORCE
	cd bessel ; $(MAKE)

linalg: FORCE
	cd linalg ; $(MAKE)

calculus: FORCE
	cd calculus ; $(MAKE)

FORCE:

clean:
	-rm *.a *.o core *~
	cd bessel ; $(MAKE) clean
	cd linalg ; $(MAKE) clean
	cd calculus ; $(MAKE) clean
