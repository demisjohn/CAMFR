include ../../../../make.inc

all: mueller.o contour.o patterson_z_n.o allroots.o

mueller.o: mueller.h mueller.cpp
	$(CC) $(FLAGS) -c mueller.cpp

contour.o: contour.h contour.cpp
	$(CC) $(FLAGS) -c contour.cpp

patterson_z_n.o: patterson_z_n.h patterson_z_n.cpp
	$(CC) $(FLAGS) -c patterson_z_n.cpp

allroots.o: allroots.h allroots.cpp contour.h
	$(CC) $(FLAGS) -c allroots.cpp

test: mueller.o contour.o patterson_z_n.o allroots.o roottest.cpp
	$(CC) $(FLAGS) roottest.cpp mueller.o contour.o \
	patterson_z_n.o allroots.o  ../../../defs.o \
	../../../util/cvector.o ../../linalg/linalg.o ../polyroot/*.o \
	$(LFLAGS) -o roottest

clean:
	-rm *.a *.o core *~ roottest
