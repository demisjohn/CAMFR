#! /usr/bin/env python

###########################################################################
#
# File:     geometry3d.py
# Author:   Peter.Bienstman@rug.ac.be, Pieter.Dumon@intec.UGent.be
# Date:     20040306
# Version:  1.0
#
# Copyright (C) 2004 Peter Bienstman - Ghent University
#
############################################################################

from camfr_work import *
from geometry import *

############################################################################
#
# Note: the coordinate system used here is the same as that of Sections:
# x & y for the transversal directions, z for the propagation direction
#
############################################################################

def zintersection_compare(xyobjs1,xyobjs2):
  
  if len(xyobjs1) != len(xyobjs2):
    return 1
  
  for i in range(len(xyobjs1)):
    if xyobjs2[i].compare(xyobjs1[i]):
      return 1
    
  return 0



############################################################################
#
# class Point3D
#
############################################################################

class Point3D:

    def __init__(self, x, y,z):
        self.x = x
        self.y = y
        self.z = z
    def __repr__(self):
        return "("+ repr(self.x) +","+ repr(self.y) +","+ repr(self.z) +")"



def sort_point3d(p1, p2):
  
    if p1.z < p2.z:
        return -1
    else:
        return 1



############################################################################
#
# The following are a number of geometric shapes that can be used to
# define a structure.
#
# The user can easily add support for new shapes by defining a class that
# has the 'intersection_at_z' method. This method should return a 2D shape 
# describing the intersection of the 3D shape with the XY plane at a given 
# z-coordinate
#
# The class should also have a 'mat' member with the material the shape is
# made of.
#
# Note that only convex shapes are supported (as the 2D z intersections must 
# be convex)
#
############################################################################

############################################################################
#
# class Cylinder
#
############################################################################

class Cylinder:

    def __init__(self, c, r, h, mat):
        self.c   = c
        self.r   = r
        self.h   = h
        self.mat = mat

    def intersection_at_z(self, z):
      
        # Intersection is a rectangle or a line
        
        D = self.r**2-(self.c.z-z)**2
        if D < 0:
          return 0
        if D == 0:
          return Line(Point(self.c.x,self.c.y-self.h/2.0), \
                      Point(self.c.x,self.c.y+self.h/2.0),self.mat) 
        else:
          return Rectangle(Point(self.c.x-sqrt(D),self.c.y-self.h/2.0), \
                           Point(self.c.x+sqrt(D),self.c.y+self.h/2.0), \
                           self.mat)
        
  

############################################################################
#
# class Box
#
############################################################################

class Box:

    def __init__(self, p1, p2, mat):
        p = [p1, p2]
        p.sort(sort_point)
        self.p1, self.p2 = p
        self.mat = mat

    def center(self):
        return Point3D((self.p1.x+self.p2.x)/2., \
                       (self.p1.y+self.p2.y)/2., \
                       (self.p1.z+self.p2.z)/2.)

    def width(self):
        return self.p2.x - self.p1.x

    def height(self):
        return self.p2.y - self.p1.y

    def length(self):
        return self.p2.z - self.p1.z
        
    def intersection_at_z(self, z):
      
        # Intersection is a rectangle
        
        if (z < self.p1.z) or (z > self.p2.z):
            return 0
        else:
            return Rectangle(Point(self.p1.x,self.p1.y), \
                             Point(self.p2.x,self.p2.y),self.mat)


       
############################################################################
#
# Geometry3D
#
#   A geometry consisting of a list of shapes. If two shapes are present
#   at the same point, the latest added shape takes precedence.
#   The number of section solver stage modes is specified when creating
#   the geometry
#   
#   A set of sections is constructed from a geometry as follows.
#   First, the geometry is discretised in slices along the z-axis 
#   between z0 and z1 in steps of dz. 
#   To construct a section from a slice, The slice is discretised along
#   the x-axis between x0 and x1 in steps of dx.
#   Subsequently, slices for different x-values are combined if the
#   interfaces for the materials they contain are no more than dy apart
#   in the vertical y-direction. This can sometimes result in symmetric
#   structures becoming unsymmetric, so you might want to experiment with
#   the values of dx and dy. Setting dy to zero will result in no slices
#   being combined, unless they are identical.
#
############################################################################

section_cache = []

def get_cached_section(xy):
  global section_cache
  for i in range(len(section_cache)):
    if not zintersection_compare(xy,section_cache[i][1]):
      return section_cache[i][0]
  return 0
      
class Geometry3D:

    def __init__(self, background_mat,M1,M2):
        self.background_mat = background_mat
        self.shapes = []
        self.M1 = M1
        self.M2 = M2

    def add(self, s):
        self.shapes.append(s)

    def __iadd__(self, s):
        self.shapes.append(s)
        return self
   
    def to_expression(self, x0, x1, dx, y0, y1, dy, z0, z1, dz, add_flipped=0):

        if (x0 > x1) or (y0 > y1) or (z0 > z1):
            print "Error: Invalid boundaries for to_expression."
            raise IndexError

        sections = []
        z = z0 + dz/2.0
        
        while z < z1:
          
          # Get a list of shapes intersecting this z position.
          
          xyobjs = []
          for i in range(len(self.shapes)):
            xyobj = self.shapes[i].intersection_at_z(z)
            if xyobj==0:
              continue
            xyobjs.append(xyobj)
          
          # Check if geometry is the same as at previous z position.
          
          secc = get_cached_section(xyobjs)
          if secc != 0:
            sections.append(secc)
            z += dz
            continue
           
          # Build the section at this z position.
          
          x = x0+dx/2.0
          slabs = []
          while x < x1:
            
            # Build slab for this x position.
            
            slab = [[y0,y1,self.background_mat]]
            
            for i in range(len(xyobjs)):  
              
              ys = xyobjs[i].intersection_at_x(x)
              if len(ys) == 0:
                continue  
              
              ys0, ys1 = ys
              if ys0 < y0:
                ys0 = y0
              if ys1 > y1:
                ys1 = y1
                  
              def same(y0,y1): return abs(y0-y1) < 1e-6
              
              new_slab = []
              j = 0
              while j< len(slab):
                    if (slab[j][1] < ys0) or same(slab[j][1], ys0) or \
                       (slab[j][0] > ys1) or same(slab[j][0], ys1) or \
                       (abs(ys0-ys1) < .001*dy):
                        new_slab.append(slab[j]) # No intersection.
                    else:
                        if not same(slab[j][0], ys0): # Old material pre.
                            new_slab.append([slab[j][0], ys0, slab[j][2]])

                        new_slab.append([ys0, ys1, xyobjs[i].mat])

                        while slab[j][1] < ys1:
                            j += 1

                        if not same(slab[j][1], ys1): # Old material post.
                            new_slab.append([ys1, slab[j][1], slab[j][2]])

                    j += 1

              slab = new_slab
                
            # Consolidate slab in y-direction.

            new_slab = []
            i = 0
            while i < len(slab):                
              i_end = i+1
              while (i_end < len(slab)) and (slab[i_end][2] == slab[i][2]):
                  i_end += 1
               
              new_slab.append([slab[i][0], slab[i_end-1][1], slab[i][2]])
              i = i_end
            
            # Go to next x position.
  
            slabs.append(new_slab)
            x += dx
          
          # Consolidate in x-direction.

          d = []
          new_slabs = []
        
          i = 0
          while i < len(slabs):
            
            i_end = i+1
            while     (i_end < len(slabs)) \
                  and similar(slabs[i_end], slabs[i], dy) \
                  and direction(slabs[i], slabs[i+1]) \
                   == direction(slabs[i_end-1], slabs[i_end]):
                i_end += 1

            if i_end == i+1:
                new_slabs.append(slabs[i])
            else:
                new_slabs.append(average_slabs(slabs, i, i_end))
                
            d.append((i_end - i) * dx)
                                              
            i = i_end

          d[-1] += (x1-x0) - sum(d)

          slabs = new_slabs 
          
          # Create slabs & section.
          
          e_sec = Expression()
          for i in range(len(slabs)):
            e_slab = Expression()
            d_y = []
            for j in range(len(slabs[i])):
                
                chunk_d = slabs[i][j][1] - slabs[i][j][0]                
                chunk_m = slabs[i][j][2]

                d_y.append(chunk_d)

                if j == len(slabs[i])-1:
                  chunk_d += (y1-y0) - sum(d_y)
                  
                e_slab.add(chunk_m(chunk_d))

            s = Slab(e_slab)

            slab_cache.append(s)
            e_sec.add(s(d[i]))
          
          sec = Section(e_sec,self.M1,self.M2)
          sections.append(sec)
          section_cache.append((sec,xyobjs))
          
          # Go to next z position.
          
          z += dz
          
        # Make an expression of all sections.
        
        e = Expression()
        d_z = []
        final_sections = []
        oldsection = sections[0]
        zsecs = 0
 
        for i in range(len(sections)):
          if sections[i]==oldsection:
            zsecs += dz
          else:
            e.add(oldsection(zsecs))
            d_z.append(zsecs)
            final_sections.append(oldsection)
            oldsection = sections[i]
            zsecs = dz

        zsecs = (z1-z0) - sum(d_z)
        e.add(oldsection(zsecs))
        d_z.append(zsecs)
        final_sections.append(oldsection)

        # Add flipped scatterer.

        if add_flipped:
          for i in range(len(final_sections)):
            e.add((final_sections[-i-1])(d_z[-i-1]))

        return e

        
def pretty_print(s):
    for i in range(s.N()):
        nr = s.mode(i).n_eff().real
        ni = s.mode(i).n_eff().imag

        if abs(nr) < 1e-6:
            nr = 0
            
        if abs(ni) < 1e-6:
            ni = 0

        print i, nr, ni
