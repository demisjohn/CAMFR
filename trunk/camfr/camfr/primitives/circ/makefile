include ../../../make.inc

OBJS = circ.o circdisp.o circmode.o circoverlap.o

all: $(OBJS)
	$(MAKE) $(OBJS)

circ.o: circ.cpp circ.h circdisp.h ../../defs.h ../../mode.h ../../waveguide.h
	$(CC) -I../../math/calculus $(FLAGS) -c circ.cpp

circdisp.o: circdisp.cpp circdisp.h ../../math/bessel/bessel.h ../../defs.h
	$(CC) $(FLAGS) -c circdisp.cpp

circmode.o: circmode.cpp circmode.h circ.h ../../mode.h ../../defs.h
	$(CC) $(FLAGS) -c circmode.cpp

circoverlap.o: circoverlap.cpp circoverlap.h circ.h ../../mode.h \
	../../interface.h ../../defs.h
	$(CC) $(FLAGS) -c circoverlap.cpp

clean:
	-rm *.a *.o *.so *.py* *_wrap* core *~
