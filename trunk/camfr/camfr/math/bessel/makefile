include ../../../make.inc

bessel.a: bessel.o slatec
	$(ARCH) $(AFLAGS) bessel.a bessel.o slatec/*.o
	$(RANLIB) bessel.a

bessel.o: bessel.h bessel.cpp
	$(CC) $(FLAGS) -c bessel.cpp

slatec: FORCE
	cd slatec ; $(MAKE)

bench: bench.cpp bessel.a
	$(CC) $(FLAGS) -c bench.cpp 
	$(LINKER) bench.o bessel.a $(LFLAGS) -o bench
	time bench 

FORCE:

clean:
	-rm *.a *.o bench core *~
	cd slatec ; $(MAKE) clean
