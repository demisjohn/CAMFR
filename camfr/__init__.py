from __future__ import division

from _camfr import *
from camfr_tk import *
from geometry import *
from camfrversion import *

# Splash screen.

print
print "CAMFR", camfr_version,
print "- Copyright (C) 1998-2002 Peter Bienstman - Ghent University."
print

# Temporary workaround for lack of enum support in BPL V2.

Plus = 0
Min  = 1

ADR    = 0
track  = 1
series = 2

normal = 0
extra  = 1
SVD    = 2

T_T = 0
S_T = 1
S_S = 2

GEV = 0
T   = 1

lapack  = 0
arnoldi = 1

unknown = 0
TEM     = 1
TE      = 2
TM      = 3
HE      = 4
EH      = 5
TE_TM   = 6

cos_type = 0
sin_type = 1

slab_E_wall    = SlabWallMixed(1.0,  1.0)
slab_H_wall    = SlabWallMixed(1.0, -1.0)
slab_open_wall = SlabWallMixed(0.0,  1.0)

E_wall = 0
H_wall = 1
