include ../../../../make.inc

all: patterson.o patterson_quad.o

patterson.o: patterson.h patterson.cpp
	$(CC) $(FLAGS) -c patterson.cpp

patterson_quad.o: patterson_quad.h patterson_quad.cpp
	$(CC) $(FLAGS) -c patterson_quad.cpp

test: patterson.o patterson_quad.o patterson_test.cpp
	$(CC) $(FLAGS) patterson_test.cpp patterson.o \
	patterson_quad.o ../../../defs.o -o patterson_test
	
clean:
	-rm patterson_test *.a *.o core *~
