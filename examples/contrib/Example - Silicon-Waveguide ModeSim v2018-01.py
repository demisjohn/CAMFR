# -*- coding: utf-8 -*-
"""
Python/Spyder script
Run Configuration should either:
    Execute in Current Interpreter
        OR
        
    Execute in New Interpreter 
        AND
    Interact with shell after completion (for plt.show() to work)


@author: Demis John
@email: demis@praevium.com


uses the following modules (import them to run this script outside of Spyder)
    import numpy as np               # NumPy (multidimensional arrays, linear algebra, ...)
    import scipy as sp               # SciPy (signal and image processing library)
    import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
    import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
    from pylab import *              # Matplotlib's pylab interface - for commands just like Matlab
"""
####################################################
# Module setup etc.

#from __future__ import division  # Fix nonsense division in Python2.x (where 1/2 = 0 )- unneeded in new python versions
import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import scipy as sp  # SciPy (signal and image processing library)

import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
#from pylab import *              # Matplotlib's pylab interface, to enable typing commands just like MatLab.
plt.ion()                            # Turned on Matplotlib's interactive mode

## For LaTeX usage (slows down plotting by maybe 5-10 sec per plot):
#from matplotlib import rc
#rc('text', usetex=True)        
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

####################################################
## Import required modules

import camfr as c

import sys
import os.path                # For handling folder hierarchy on the current OS
#import readline             # For pausing by waiting for a key to be pressed via raw_input("stuff")
import time                    # for getting current date & time w/ strftime() & pausing for a few seconds w/ time.pause(seconds)
#from numpy import append



## find & import the UCSB_includes functions
## Sub Folder names:
SubFolder_UCSB_includes = "UCSB_includes"

ScriptDir = sys.path[0]; print "ScriptDir is: %s" %ScriptDir                # Get directory of this script
(ParentDir , tail) = os.path.split(ScriptDir)                                # split off top-level directory from path
#SubFolder_UCSB_includes = os.path.join(ParentDir,SubFolder_UCSB_includes)    # append subdirectory
SubFolder_UCSB_includes = os.path.join(SubFolder_UCSB_includes)    # append subdirectory
if ( not os.path.isdir(SubFolder_UCSB_includes) ):                            # check if SubDirectory exists
        print "UCSB_includes Directory is NOT valid!!!"
        print "Expected path to UCSB Includes: %s" %SubFolder_UCSB_includes
sys.path.append(SubFolder_UCSB_includes)        # add subfolder to the python search path

from camfr_UCSB_v2018_03 import *         # import Demis' modified versions of Jared's CAMFR functions
#import n                                         # import cauchy models & film parameters for n vs. wavelength
#from SpiralLayoutSimulation import *             # import Demis' Spiral Layout simulation function


####################################################

print 'Running...'

## Choose your output filenames:

## Subfolder for Output files:
CurrentDate = time.strftime("%Y-%m-%d %H.%M.%S")           # Get current date and time as string in specified format
OutputSubFolder = "Testing - [" + CurrentDate + "]";     OutputSubFolder = os.path.join(ScriptDir,OutputSubFolder)




''' Simulation Options '''

## Some refractive indices ("RIx")
wl = 1.550      # wavelength in microns

nSiN = 3.4   # Silicon
nTOX = 1.446  # Thermal Oxide
nTEOS = nTOX


tCore = 0.300   # core thickness in Microns
tClad = 2.0

wCore = 1.00    # core width
wSides = 2.0    # lateral width of simulation space


#BendWG = False  #optionally include bend radius of waveguide
#BendR = 200.0   # Bend Radius in um





''' Setup some CAMFR simulation parameters '''

c.set_lambda(wl)
lam0 = abs(c.get_lambda())
k0 = 2*pi/lam0
set_polarisation(c.TM)
c.set_N(4)       # Computation time scales logarithmically with the number of eigenmodes used



''' Start simulations '''
TimeStart = time.time()    #record time we started the simulations
os.mkdir(str( OutputSubFolder ))        # Create the new folder


#[neff,MultiMode] = GetModalIndex2D(c.get_lambda(), nSiN, nTOX, tCore, wCore, tTOX, wSides)

core = c.Material(nSiN)
clad = c.Material(nTOX)
air = c.Material( 1.0 )

center = c.Slab( clad(tClad) + core(tCore) + clad(tClad) )   # Vertical, bottom-to-top
side = c.Slab( clad(2*tClad) + air(center.width() - 2*tClad)  )

all = side(wSides) + center(wCore) + side(wSides)    # Horizontal, left-to-right
        
s = c.Section(all, 16, 30)    
# Section(Waveguide, <# plane waves in 1st stage estimation>, <# 1D modes in each Slab during 2nd stage>)
# # planewaves should improve mode ordering and finding
# Num. Slab Modes should improve accuracy of final mode profiles
# used to be (all,12,100), reduced to (all,12,40) for speed

## CAMFR solver params
c.set_section_solver(c.L)
c.set_mode_correction(c.full)

c.set_left_wall(c.E_wall)
c.set_right_wall(c.E_wall)
c.set_lower_wall(c.slab_E_wall)
c.set_upper_wall(c.slab_E_wall)
c.set_right_PML(-0.04)
c.set_left_PML(-0.04)
c.set_upper_PML(-0.04)
c.set_lower_PML(-0.04)


s.calc()
#fig = s.plot(dx=0.3, dy=0.3)
fig = s.plot(field=['P','Ex','Ey'], mode=[0,1,2], dx=0.05, dy=0.05)

# Save the figure as an image.
#fig.savefig('Silicon_WG_-_Modesolver_example_v1.png')


print 'done.'



