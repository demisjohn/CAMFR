# This Python script contains all the machine dependent settings
# needed during the build process.

# Compilers to be used.

cc  = "/home/pbienst/local/bin/gcc"
cxx = "/home/pbienst/local/bin/g++"
f77 = "/home/pbienst/local/bin/g77"

# Compiler flags.
#
# Note: for the Fortran name definition you can define one of the following
#       preprocessor macros:
#
#           FORTRAN_SYMBOLS_WITHOUT_TRAILING_UNDERSCORES
#           FORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
#           FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES

base_flags = "-fPIC -ftemplate-depth-60 \
	      -DFORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE -DNDEBUG"

flags = base_flags + "-O3 --param max-inline-insns=600 -march=pentiumpro \
                       -ffast-math -funroll-loops -fstrict-aliasing -g"

fflags = base_flags

flags_noopt = base_flags

#flags = base_flags

# Include directories.

include_dirs ="/home/pbienst/blitz-20001213 \
	       /home/pbienst/boost_cvs/boost \
	       /home/pbienst/local/include/python2.2 \
	       /home/pbienst/arpack++/include"

# Temporary workaround waiting for shared lib support in scons.

camfrlib = "libcamfr.a"
dllsuffix = ".so"
dllcommand = cc + " -shared camfr/camfr_wrap.o \
	-L/home/pbienst/camfr_work/camfr \
	-L/home/pbienst/blitz-20001213/lib \
	-L/opt/intel/mkl/lib/32 \
	-L/home/pbienst/ARPACK \
	-lcamfr -lbpl -lblitz -larpack -lmkl_lapack -lmkl_p3 -lguide \
	-lm -lg2c -lstdc++ -o camfr/_camfr" + dllsuffix

# Command to strip library of excess symbols:

strip_command = "strip --strip-unneeded camfr/_camfr" + dllsuffix
strip_command = ""

# Extra files to copy to installation directory.

extra_files = []
