# Imports.

from __future__ import division

import numpy as np

import os
if not os.environ.get('NO_CAMFR_GRAPHICS'):
    from pylab import *

from _camfr import *
from camfr_PIL import *     # converted numpy* to np.*
from geometry import *
from geometry3d import *
from material import *
from camfrversion import *

# Splash screen.

print
print "CAMFR", camfr_version,
print "- Copyright (C) 1998-2007 Peter Bienstman - Ghent University."
print
