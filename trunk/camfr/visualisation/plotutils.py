#! /usr/bin/env python

##############################################################################
#
# File:    plotutils.py
# Author:  Peter.Bienstman@rug.ac.be
# Date:    20020215
#
# Various small utilty functions for plotting and movie creation using
# plplot (http://plplot.sourcefoge.net), ImageMagick and gifsicle.
#
##############################################################################

from cmath import *
from pl import *
from Numeric import *
from MLab import *

##############################################################################
#
# Generates a shading (i.e flat) plot of a 2D matrix.
# Axes can be toggled and the picture can be written to a jpg file.
#
##############################################################################

def plot2D(z, x=arange(0,1,0.1), y=arange(0,1,0.1),
           axis="none", filename="", zrange=()):
  
  if filename != "":
      plsdev("jpeg")
      plsfnam(filename)
  else:
      plsdev("xwin")
      
  plscolbg(255,255,255) # A white background looks so much better...
  plinit()
  
  if len(zrange) == 0:
    zmin, zmax = min(min(z)), max(max(z))
  else:
    zmin, zmax = zrange
  
  if axis != "none":
    plenv(x[0],x[-1],y[0],y[-1],1,0)
  else:
    plenv(x[0],x[-1],y[0],y[-1],1,-2)

  ns = 64 # Number of shades
  sh_cmap = 1
  min_color = 1; min_width = 0; max_color = 0; max_width = 0
  clevel = zmin + (zmax - zmin) * (arrayrange(ns)+.5)/ns

  for i in range(ns):
    shade_min = zmin + (zmax - zmin) * i / ns
    shade_max = zmin + (zmax - zmin) * (i+1) / ns
    sh_color = i/(ns-1.)
    sh_width = 2
    plpsty(0)
    plshade(z, x[0],x[-1], y[0],y[-1],
	    shade_min, shade_max, sh_cmap, sh_color, sh_width,
	    min_color, min_width, max_color, max_width, 1)

  if axis != "none":
    plcol0(1)
    plbox( "bcnst", 0., 0, "bcnstv", 0., 0 )
    plcol0(2)

  plend()



##############################################################################
#
# Creates a movie starting from a complex matrix representing phasors.
#
##############################################################################

import os

def phasormovie(z, x=arange(0,1,0.1), y=arange(0,1,0.1),
                axis="none", filename=""):
  
  frames = 8

  zmax = max(max(abs(z)))
  zmin = -zmax
  
  for i in range(frames):
    plot2D(z.real,x,y,axis, "frame" + '%02d' % i + ".jpg", (zmin,zmax))
    z *= exp(2j*pi/frames)

  if filename == "":
    os.system("animate -delay 10 frame*.jpg")
  else:
    os.system("convert -loop 0 -delay 10 frame*.jpg " + filename + ".gif")

  os.system("rm frame*.jpg")



##############################################################################
#
# Plot the refractive index profile in a stack.
#
##############################################################################

from camfr import *

def plot_n(stack, r_x, r_z, axis="none", filename=""):
    
    n = zeros([len(r_z),len(r_x)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        n[i_z,i_x] = stack.n(Coord(r_x[i_x], 0, r_z[i_z])).real

    plot2D(n, r_z, r_x, axis, filename)



##############################################################################
#
# Plot the field profile in a stack.
#
##############################################################################

def plot_field(stack, r_x, r_z, component, axis="none", filename=""):
    
    f = zeros([len(r_z),len(r_x)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[i_z,i_x] = component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    plot2D(f, r_z, r_x, axis, filename)



##############################################################################
#
# Animate the field profile in a stack.
#
##############################################################################

def animate_field(stack, r_x, r_z, component, axis="none", filename=""):
    
    f = zeros([len(r_z),len(r_x)], Complex)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[i_z,i_x] = component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    phasormovie(f, r_z, r_x, axis, filename)
