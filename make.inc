# compiler environment

cc     = gcc
CC     = g++
F77    = g77
LINKER = g++
MAKE   = make
ARCH   = ar
AFLAGS = cr
RANLIB = ranlib

# location of the Blitz library

BLITZ  = /home/pbienst/blitz-20001213

# location of the BLAS/LAPACK libraries

#LAPACK = /home/pbienst/ATLAS/lib/Linux_PIIISSE1/liblapack.a \
#	 /home/pbienst/ATLAS/lib/Linux_PIIISSE1/libcblas.a \
#	 /home/pbienst/ATLAS/lib/Linux_PIIISSE1/libf77blas.a \
#	 /home/pbienst/ATLAS/lib/Linux_PIIISSE1/libatlas.a

LAPACK = /opt/intel/mkl/lib/32/libmkl_lapack.a \
	 /opt/intel/mkl/lib/32/libmkl_p3.a \

# location of the Python header files

PYTHON = /home/pbienst/local/include/python2.1

# location of the Boost library include files

BOOST_INC = /home/pbienst/boost_1_25_0

# location of the Boost library

BOOST_LIB = /home/pbienst/boost_1_25/libs/python/src/libboost_python.a

# Arpack++ header files

ARPACK_INC = /home/pbienst/arpack++/include

# Arpack Fortran library

ARPACK_LIB = /home/pbienst/ARPACK/libarpack.a

# location of the math libraries (usually -lm, but can be replaced by
# other routines to improve performance)

MATH   = -lm

# compiler flags

FLAGS  = -fPIC -O3 --param max-inline-insns=600 -march=pentiumpro \
	-ftemplate-depth-30 -ffast-math -g\
      	-funroll-loops -fstrict-aliasing -I$(BLITZ)\
	-I$(PYTHON) -I$(BOOST_INC) -I$(ARPACK_INC)
#FLAGS  = -fPIC -I$(BLITZ) -I$(PYTHON) -ftemplate-depth-30 -g \
#	-march=pentiumpro -I$(BOOST_INC) -I$(ARPACK_INC)
FFLAGS = $(FLAGS) # Fortran flags

# linker flags

LFLAGS = $(BLITZ)/lib/libblitz.a $(BOOST_LIB) $(LAPACK) \
	$(MATH) $(ARPACK_LIB) -lg2c

# Format of Fortran symbol names. Possible values are 
#
#   -DFORTRAN_SYMBOLS_WITHOUT_TRAILING_UNDERSCORES
#   -DFORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
#   -DFORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES

FSYMB  = -DFORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
