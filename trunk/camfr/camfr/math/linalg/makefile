include ../../../make.inc

linalg.o: linalg.cpp linalg.h ../../defs.h
	$(CC) $(FLAGS) $(FSYMB) -c linalg.cpp

lintest.o: lintest.cpp linalg.h
	$(CC) $(FLAGS) -c lintest.cpp

test: lintest.o linalg.o
	$(LINKER) lintest.o linalg.o ../../defs.o $(LFLAGS) -o lintest

clean:
	-rm *.o core lintest *~

