include ../../../../make.inc

all: jenkins_traub.o polyroot.o

jenkins_traub.o: jenkins_traub.f
	$(F77) $(FFLAGS) -c jenkins_traub.f

polyroot.o: polyroot.h polyroot.cpp
	$(CC) $(FLAGS) -c polyroot.cpp

test: jenkins_traub.o polyroot.o polyroot_test.cpp
	$(CC) $(FLAGS) polyroot_test.cpp jenkins_traub.o \
	polyroot.o ../../../defs.o -o polyroot_test
	
clean:
	-rm polyroot_test *.a *.o core *~
