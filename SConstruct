# tmp workaround for scons bug

import sys, os
sys.path = [os.getcwd()] + sys.path

from machine_cfg import *

Default("camfr")

# Construct build environments.

env = Environment(CPPPATH = include_dirs,
	          CC  = cc,  CCFLAGS  = flags,
		  CXX = cxx, CXXFLAGS = flags,
	          F77 = f77, F77FLAGS = fflags,
		  CAMFRLIB   = camfrlib,	
		  DLLCOMMAND = dllcommand, 
		  DLLSUFFIX  = dllsuffix,
		  ENV        = os.environ)

env_noopt = env.Copy(CCFLAGS = flags_noopt, CXXFLAGS = flags_noopt)

Export("env", "env_noopt")
SConscript("camfr/SConscript")