include ../../../../make.inc

traceroot.o: traceroot.h traceroot.cpp #cvector
	$(CC) $(FLAGS) -c traceroot.cpp

cvector: FORCE
	cd ../../../util ; make

roottest.o: roottest.cpp
	$(CC) $(FLAGS) -c roottest.cpp

test: roottest.o traceroot.o
	$(CC) roottest.o traceroot.o ../../../defs.o \
	../croot/croot.a ../root/root.o ../../../util/cvector.o \
	../../bessel/bessel.a \
	$(LFLAGS) -o roottest

FORCE:
 
clean:
	-rm *.a *.o core roottest *~

