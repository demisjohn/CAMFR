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

base_flags = "-ftemplate-depth-60 \
	      -DFORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE -DNDEBUG "

flags_noopt = base_flags

flags = base_flags + " --param max-inline-insns=600  \
                       -ffast-math -funroll-loops -fstrict-aliasing -g"

fflags = flags

# Include directories.

include_dirs = ["/home/pbienst/blitz-20001213",
	        "/home/pbienst/boost_cvs/boost",
	        "/home/pbienst/local/include/python2.2"]

# Library directories.

library_dirs = ["/home/pbienst/blitz-20001213/lib",
                "/opt/intel/mkl/lib/32",
                "/home/pbienst/local/lib"]

# Library names.

libs = ["bpl", "blitz", "mkl_lapack", "mkl_p3", "guide", "g2c"]

# Command to strip library of excess symbols:

dllsuffix = ".so"
strip_command = "strip --strip-unneeded camfr/_camfr" + dllsuffix
strip_command = ""

# Extra files to copy into installation directory.

extra_files = []
