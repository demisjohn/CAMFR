include ../make.inc

OBJS = defs.o material.o coord.o field.o mode.o waveguide.o \
       scatterer.o chunk.o interface.o icache.o expression.o \
       stack.o S_scheme.o T_scheme.o S_scheme_fields.o \
       T_scheme_fields.o cavity.o bloch.o	

libcamfr.a: FORCE
	$(MAKE) $(OBJS) util math primitives 
	$(ARCH) $(AFLAGS) libcamfr.a  $(OBJS) util/*.o \
	math/linalg/*.o math/calculus/root/*.o \
	math/calculus/minimum/*.o math/calculus/traceroot/*.o \
	math/calculus/polyroot/*.o math/calculus/croot/*.o \
	math/calculus/quadrature/*.o \
	primitives/planar/*.o primitives/slab/*.o \
	primitives/slab/isoslab/*.o primitives/section/*.o \
	primitives/circ/*.o math/bessel/*.o math/bessel/slatec/*.o
	$(RANLIB) libcamfr.a

camfr.so: FORCE
	$(MAKE) $(OBJS) camfr_wrap.o util math primitives 
	$(LINKER) -shared -o camfr.so $(OBJS) util/*.o \
	math/linalg/*.o math/calculus/root/*.o \
	math/calculus/minimum/*.o math/calculus/traceroot/*.o \
	math/calculus/polyroot/*.o math/calculus/croot/*.o \
	math/calculus/quadrature/*.o \
	primitives/planar/*.o primitives/slab/*.o \
	primitives/slab/isoslab/*.o \
	primitives/circ/*.o math/bessel/bessel.a $(LFLAGS)

camfr_wrap.o: camfr_wrap.cpp defs.h coord.h mode.h field.h material.h \
	waveguide.h scatterer.h expression.h stack.h cavity.h bloch.h \
	math/calculus/function.h primitives/planar/planar.h \
	primitives/circ/circ.h primitives/circ/circdisp.h \
	primitives/slab/generalslab.h primitives/slab/isoslab/slab.h \
	primitives/slab/isoslab/slabwall.h primitives/slab/isoslab/slabdisp.h
	$(CC) $(FLAGS) -c camfr_wrap.cpp	

defs.o: defs.cpp defs.h # Don't optimise calc_machine_eps.
	$(CC) -fPIC -c defs.cpp

material.o: material.cpp material.h defs.h
	$(CC) $(FLAGS) -c material.cpp

coord.o: coord.cpp coord.h
	$(CC) $(FLAGS) -c coord.cpp

field.o: field.cpp field.h
	$(CC) $(FLAGS) -c field.cpp

mode.o: mode.cpp mode.h defs.h
	$(CC) $(FLAGS) -c mode.cpp

waveguide.o: waveguide.cpp waveguide.h defs.h mode.h
	$(CC) $(FLAGS) -c waveguide.cpp

scatterer.o: scatterer.cpp scatterer.h defs.h waveguide.h material.h
	$(CC) $(FLAGS) -c scatterer.cpp

chunk.o: chunk.cpp chunk.h defs.h waveguide.h material.h
	$(CC) $(FLAGS) -c chunk.cpp

interface.o: interface.cpp interface.h waveguide.h mode.h defs.h
	$(CC) $(FLAGS) -c interface.cpp

icache.o: icache.cpp icache.h interface.h waveguide.h defs.h
	$(CC) $(FLAGS) -c icache.cpp

expression.o: expression.h expression.cpp scatterer.h waveguide.h
	$(CC) $(FLAGS) -c expression.cpp

stack.o: stack.h stack.cpp scatterer.h waveguide.h defs.h
	$(CC) $(FLAGS) -c stack.cpp

S_scheme.o: S_scheme.cpp S_scheme.h defs.h math/linalg/linalg.h
	$(CC) $(FLAGS) -c S_scheme.cpp

T_scheme.o: T_scheme.cpp T_scheme.h defs.h math/linalg/linalg.h
	$(CC) $(FLAGS) -c T_scheme.cpp

S_scheme_fields_old.o: S_scheme_fields_old.cpp S_scheme_fields_old.h defs.h \
	math/linalg/linalg.h
	$(CC) $(FLAGS) -c S_scheme_fields_old.cpp

S_scheme_fields.o: S_scheme_fields.cpp S_scheme_fields.h defs.h \
	math/linalg/linalg.h
	$(CC) $(FLAGS) -c S_scheme_fields.cpp

T_scheme_fields.o: T_scheme_fields.cpp T_scheme_fields.h defs.h \
	math/linalg/linalg.h
	$(CC) $(FLAGS) -c T_scheme_fields.cpp

cavity.o: cavity.h cavity.cpp defs.h stack.h material.h
	$(CC) $(FLAGS) -c cavity.cpp

bloch.o: bloch.h bloch.cpp defs.h stack.h
	$(CC) $(FLAGS) -c bloch.cpp

util: FORCE
	cd util ; $(MAKE)

math: FORCE
	cd math ; $(MAKE)

primitives: FORCE
	cd primitives ; $(MAKE)

FORCE:

test:
	cd testsuite ; make test

clean:
	-rm *.a *.so *.pyc *.o core *~
	cd util ; make clean
	cd math ; make clean
	cd primitives ; make clean

localclean:
	-rm *.a *.so *.pyc *.o core *~

