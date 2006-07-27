from machine_cfg import *

SConsignFile() # Tell Scons not to write data all over the file system.

Default("camfr")
import os

# Construct build environments.

env = Environment(CPPPATH = include_dirs,
 		  LIBPATH = library_dirs,
		  CC   = cc,   CCFLAGS   = flags,
		  CXX  = cxx,  CXXFLAGS  = flags,
		  F77  = f77,  F77FLAGS  = fflags,
		  LINK = link, LINKFLAGS = link_flags,
		  LIBS = libs, SHLIBPREFIX = "", 
		  ENV  = os.environ)

env_noopt = env.Copy(CCFLAGS = flags_noopt, CXXFLAGS = flags_noopt)

Export("env", "env_noopt")
SConscript("camfr/SConscript")