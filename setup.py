#! /usr/bin/env python

from distutils.core import setup
from distutils.util import byte_compile
from distutils.command.build import build
from distutils.command.install_data import install_data

from machine_cfg import *
from camfrversion import *

# Make sure we build the libraries before running the standard build.
# Also, after the build process, strip the library of debug symbols.

class camfr_build(build):
  def run(self):

    import os
        
    os.system("scons")
    os.system(strip_command)

    return build.run(self)



# Modified install_data, changing self.install_dir to the actual library dir.
# Also byte-compiles *.py files that are outside of the regular package
# hierarchy.

class camfr_install_data(install_data):
  def run(self):

    # Byte-compile Python files.

    scripts = []
          
    for i in self.data_files:
      for j in i[1]:
        if j[-2:] == "py":
          scripts.append(j)
          i[1].append(j+'c')

    byte_compile(scripts)
      
    # Change install dir to library dir.
    
    install_cmd = self.get_finalized_command('install')
    self.install_dir = getattr(install_cmd, 'install_lib')

    return install_data.run(self)



# Set up the module.

setup(name         = "camfr_work",
      version      =  camfr_version,
      description  = "CAvity Modelling FRamework",
      author       = "Peter Bienstman",
      author_email = "Peter.Bienstman@rug.ac.be",
      url          = "http://camfr.sourceforge.net",
      extra_path   = "camfr_work",
      packages     = ["examples.tutorial", "examples.other",
                      "visualisation.examples", "testsuite"],
      data_files   = [(".", ["COPYRIGHT", "camfrversion.py",
                             "camfr/__init__.py",
                             "camfr/_camfr" + dllsuffix,
                             "camfr/geometry.py",
                             "visualisation/camfr_PIL.py",
                             "visualisation/camfr_matlab.py",
                             "visualisation/camfr_tk.py",
                             "visualisation/TkPlotCanvas.py",
                             "visualisation/gifmaker.py"])] + extra_files,
      cmdclass     = {"install_data" : camfr_install_data,
                      "build"        : camfr_build},
      )
