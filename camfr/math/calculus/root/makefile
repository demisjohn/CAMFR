include ../../../../make.inc

root.o: root.h root.cpp
	$(CC) $(FLAGS) -c root.cpp

roottest.o: roottest.cpp
	$(CC) $(FLAGS) -c roottest.cpp

test: roottest.o root.o
	$(CC) roottest.o root.o ../../bessel/bessel.a \
	../../../defs.o $(LFLAGS) -o roottest

clean:
	-rm *.a *.o core roottest *~

