include ../../../make.inc

OBJS = sectiondisp.o section.o sectionmode.o sectionoverlap.o

all: $(OBJS)
	$(MAKE) $(OBJS)

section.o: section.cpp section.h sectiondisp.h ../../defs.h ../../mode.h \
	../../waveguide.h
	$(CC) $(FLAGS) -c section.cpp

sectiondisp.o: sectiondisp.cpp sectiondisp.h  ../../defs.h
	$(CC) $(FLAGS) -c sectiondisp.cpp

sectionmode.o: sectionmode.cpp sectionmode.h section.h ../../mode.h \
	../../defs.h
	$(CC) $(FLAGS) -c sectionmode.cpp

sectionoverlap.o: sectionoverlap.cpp sectionoverlap.h section.h ../../mode.h \
	../../interface.h ../../defs.h
	$(CC) $(FLAGS) -c sectionoverlap.cpp

clean:
	-rm *.a *.o *.so *.py* *_wrap* core *~
