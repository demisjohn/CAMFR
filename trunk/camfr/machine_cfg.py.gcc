# This Python script contains all the machine dependent settings
# needed during the build process.

# Compilers to be used.

cc  = "gcc"
cxx = "g++"
f77 = "g77"

link = cxx
link_flags = ""

# Compiler flags.
#
# Note: for the Fortran name definition you can define one of the following
#       preprocessor macros:
#
#           FORTRAN_SYMBOLS_WITHOUT_TRAILING_UNDERSCORES
#           FORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
#           FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES

base_flags = "-DFORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE -DNDEBUG "

flags_noopt = base_flags

flags = base_flags + "-O3 -march=pentium4 -g -funroll-loops "

fflags = flags + "-fPIC "

# Include directories.

include_dirs = ["/usr/include/python2.4", "/usr/lib/python2.4/site-packages"]

# Library directories.

library_dirs = ["/opt/intel/mkl8/lib/32"]

# Library names.

libs = ["boost_python", "blitz", "mkl_lapack64", "mkl", "g2c"]

# Command to strip library of excess symbols:

dllsuffix = ".so"
strip_command = "strip --strip-unneeded camfr/_camfr" + dllsuffix
strip_command = ""

# Extra files to copy into installation directory.

extra_files = [("doc", ["docs/camfr.pdf"])]
