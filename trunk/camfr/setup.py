#! /usr/bin/env python

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext
import os

# Make sure we build libcamfr.a before running the standard build_ext.
# Also, after the make process, strip the library of debug symbols.

class camfr_build_ext(build_ext):
    def run(self):
        
        os.system('make')
        
        build_ext.run(self)

        fullname = self.get_ext_fullname(self.extensions[0].name)
        ext_filename = os.path.join(self.build_lib, \
                                    self.get_ext_filename(fullname))

#        os.system('strip ' + ext_filename)

# Define the camfr extension module.

camfr_extension = Extension(
    "camfr_work", ["./camfr/camfr_wrap.cpp"],                         \
    library_dirs=["./camfr",                                          \
                  "../boost_1_25_0/libs/python/src",                  \
                  "../blitz-20001213/lib",                            \
                  "/opt/intel/mkl/lib/32"],                           \
    include_dirs=["camfr", "../boost_1_25_0","../blitz-20001213"],    \
    libraries=["camfr", "boost_python", "blitz", "mkl_lapack",        \
               "mkl_p3", "m", "g2c", "stdc++"] )

# Set up the extension.

setup(name="camfr_work", version="1.0pre",                            \
      description="CAvity Modelling FRamework",                       \
      author="Peter Bienstman",                                       \
      author_email="Peter.Bienstman@rug.ac.be",                       \
      url="http://camfr.sourceforge.net/",                            \
      py_modules=['geometry', 'camfr_tk', 'TkPlotCanvas'],            \
      package_dir={'': 'camfr', '': 'visualisation'},                 \
      # FIXME: ugly hack
      data_files=[('lib/python2.1', ['camfr/geometry.py'])],           \
      ext_modules=[camfr_extension],                                  \
      cmdclass={'build_ext' : camfr_build_ext} )
