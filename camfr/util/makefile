include ../../make.inc

all: cvector.o stringutil.o index.o tracesorter.o

cvector.o: cvector.h cvector.cpp
	$(CC) $(FLAGS) -c cvector.cpp

stringutil.o: stringutil.h stringutil.cpp
	$(CC) $(FLAGS) -c stringutil.cpp

index.o: index.h index.cpp ../defs.h
	$(CC) $(FLAGS) -c index.cpp

tracesorter.o: tracesorter.h tracesorter.cpp ../math/linalg/linalg.h
	$(CC) $(FLAGS) -c tracesorter.cpp

wrap:

clean:
	-rm *.a *.o core *~
