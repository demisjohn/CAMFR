include ../../../../make.inc

OBJS = d1mach.o   dbesy.o   djairy.o  xerhlt.o  zairy.o   zbunk.o   zs1s2.o  \
       zunik.o    d9b0mp.o  dbesy0.o  dlngam.o  xermsg.o  zasyi.o   zdiv.o   \
       zseri.o    zunk1.o   d9b1mp.o  dbesy1.o  dyairy.o  xerprn.o  zbesh.o  \
       zexp.o     zshch.o   zunk2.o   d9lgmc.o  dbsynu.o  fdump.o   xersve.o \
       zbesj.o    zkscl.o   zsqrt.o   zuoik.o   dasyjy.o  dcsevl.o  i1mach.o \
       xgetua.o   zbesy.o   zlog.o    zuchk.o   zwrsk.o   dbesj.o   dgamlm.o \
       initds.o   zabs.o    zbinu.o   zmlri.o   zunhj.o   dbesj0.o  dgamln.o \
       j4save.o   zacai.o   zbknu.o   zmlt.o    zuni1.o   dbesj1.o  dgamma.o \
       xercnt.o   zacon.o   zbuni.o   zrati.o   zuni2.o

all: $(OBJS) limits.o

# limits.c must be compiled without optimisation

limits.o: 
	$(cc) $(FSYMB) -c limits.c

.f.o: 
	$(F77) $(FFLAGS) -c $*.f

print_d1mach: d1mach.o limits.o
	$(F77) $(FFLAGS) print_d1mach.f d1mach.o limits.o -o print_d1mach

clean:
	-rm *.a *.o core print_d1mach *~
