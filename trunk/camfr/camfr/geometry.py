#! /usr/bin/env python

###########################################################################
#
# File:     geometry.py
# Author:   Peter.Bienstman@rug.ac.be
# Date:     20020204
# Version:  1.0
#
# Copyright (C) 2002 Peter Bienstman - Ghent University
#
############################################################################

from math import *
from camfr import *
from Numeric import *

############################################################################
#
# Note: the coordinate system used here is x for the horizontal progation
# direction and y for the vertical direction (which is different from
# 2D slabs).
#
############################################################################

############################################################################
#
# class Point
#
############################################################################

class Point:

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return "(" + repr(self.x) + "," + repr(self.y) + ")"



def sort_point(p1, p2):
    if p1.x < p2.x:
        return -1
    else:
        return 1



############################################################################
#
# class Line
#
############################################################################

class Line:
    def __init__(self, p1, p2):
        self.p1, self.p2 = p1, p2

    def intersection_at_x(self, x):
        if self.p1.x == self.p2.x:
            return
        else:
            return self.p1.y + (self.p2.y - self.p1.y) *  \
                   1./(self.p2.x - self.p1.x) * (x - self.p1.x)



############################################################################
#
# The following are a number of geometric shapes that can be used to
# define a structure.
#
# The user can easily add support for new shapes by defining a class that
# has the 'intersection_at_x' method. This method should return the two
# y-coordinates (if any) where the line x=0 intersects the shape, sorted
# in ascending order.
#
# The class should also have a 'mat' member with the material the shape is
# made of.
#
# Note that only convex shapes are supported.
#
############################################################################

############################################################################
#
# class Circle
#
############################################################################

class Circle:

    def __init__(self, c, r, mat):
        self.c   = c
        self.r   = r
        self.mat = mat

    def intersection_at_x(self, x):
        D = self.r**2 - (x - self.c.x)**2
        if D > 0:
            return [ self.c.y - sqrt(D), self.c.y + sqrt(D) ]
        else:
            return []
        


############################################################################
#
# class Rectangle
#
############################################################################

class Rectangle:

    def __init__(self, p1, p2, mat):
        p = [p1, p2]
        p.sort(sort_point)
        self.p1, self.p2 = p
        self.mat = mat

    def center(self):
        return Point((self.p1.x+self.p2.x)/2., (self.p1.y+self.p2.y)/2.)

    def width(self):
        return self.p2.x - self.p1.x

    def height(self):
        return self.p2.y - self.p1.y

    def intersection_at_x(self, x):
        if (x < self.p1.x) or (x > self.p2.x):
            return []
        else:
            r = [self.p1.y, self.p2.y]
            r.sort()
            return r



############################################################################
#
# class Square
#
############################################################################

class Square:

    def __init__(self, c, a, mat):
        self.c = c
        self.a = a
        self.mat = mat

    def intersection_at_x(self, x):
        if (x < self.c.x - self.a/2.) or (x > self.c.x + self.a/2.):
            return []
        else:
            return [self.c.y - self.a/2., self.c.y + self.a/2.]

    def to_rectangle(self):
        return Rectangle(Point(self.c.x - self.a/2., self.c.y - self.a/2),
                         Point(self.c.x + self.a/2., self.c.y + self.a/2))
    


############################################################################
#
# class Triangle
#
############################################################################
              
class Triangle:

    def __init__(self, p1, p2, p3, mat):
        p = [p1, p2, p3]
        p.sort(sort_point)
        self.p1, self.p2, self.p3 = p
        self.mat = mat

    def intersection_at_x(self, x):
        if (x < self.p1.x) or (x > self.p3.x):
            return []
        if (x < self.p2.x):
            r = [Line(self.p1, self.p2).intersection_at_x(x), \
                 Line(self.p1, self.p3).intersection_at_x(x)]
        else:
            r = [Line(self.p1, self.p3).intersection_at_x(x), \
                 Line(self.p2, self.p3).intersection_at_x(x)]
        r.sort()
        return r



############################################################################
#
# The following are some auxiliary functions operating on slabs. Here,
# slabs are represented as nested lists, e.g. for Slab(m1(1)+m2(2)):
#
#        [ [0,1,m1],
#          [1,2,m2] ]
#
# The coordinates are along the vertical y axis.
#
############################################################################

############################################################################
#
# similar
#
#  Determines if two slabs are similar, i.e. the interface positions they
#  contain are no more then 'dy' apart.
#
############################################################################

def similar(slab1, slab2, dy):

    # Same number of chunks?

    if len(slab1) != len(slab2):
        return 0

    # Same materials?
    
    for i in range(len(slab1)):
        if slab1[i][2] != slab2[i][2]:
            return 0

    # Same interface positions?

    if abs(dy) < 1e-12:
        dy = 1e-12

    for i in range(len(slab1)):
        if abs(slab1[i][0] - slab2[i][0]) >= dy:
            return 0
        if abs(slab1[i][1] - slab2[i][1]) >= dy:
            return 0

    return 1



############################################################################
#
# direction
#
#  The direction of the gradient between two slabs.
#  Assumes that the slabs only differ in their interface positions.
#
############################################################################

def direction(slab1, slab2):
    d = []
    for i in range(len(slab1)):
        
        if slab1[i][0] < slab2[i][0]:
            d.append(-1)
        else:
            d.append(+1)
            
        if slab1[i][1] < slab2[i][1]:
            d.append(-1)
        else:
            d.append(+1)
            
    return d



############################################################################
#
# average_slabs
#
#  Returns a slab that is the 'average' of the slabs in the list 'slabs'
#  with an index between i0 (inclusive) and i1 (exclusive).
#  Assumes that these slabs only differ in their interface positions.
#
############################################################################

def average_slabs(slabs, i0, i1):
    
    if i1 == i0 + 1: # Short cut.
        return slabs[i0]

    s = slabs[i0]

    for chunk in range(len(s)):
        y0 = y1 = 0.0
        
        for i in range(i0, i1):
            
            y0 += slabs[i][chunk][0]
            y1 += slabs[i][chunk][1]
            
        s[chunk][0] = y0 / len(range(i0, i1))
        s[chunk][1] = y1 / len(range(i0, i1))

    return s



############################################################################
#
# Geometry
#
#   A geometry consisting of a list of shapes. If two shapes are present
#   at the same point, the latest added shape takes precedence.
#
#   In order to construct a discrete set of slabs from a geometry, this
#   geometry is first discretised along the horizontal x-axis between
#   x0 and x1 in steps of dx.
#
#   Subsequently, slices for different x-values are combined if the
#   interfaces for the materials they contain are no more than dy apart
#   in the vertical y-direction. This can sometimes result in symmetric
#   structures becoming unsymmetric, so you might want to experiment with
#   the values of dx and dy. Setting dy to zero will result in no slices
#   being combined, unless they are identical.
#
#   The user can specify PML_bot and PML_top, which will be used as the
#   imaginary parts for the thicknesses at y=y0 and y=y1 respectively.
#
############################################################################

slab_cache = []

class Geometry:

    def __init__(self, background_mat):
        self.background_mat = background_mat
        self.shapes = []

    def add(self, s):
        self.shapes.append(s)

    def __iadd__(self, s):
        self.shapes.append(s)
        return self
   
    def to_expression(self, x0, x1, dx, y0, y1, dy, PML_bot, PML_top,
                      add_flipped=0):

        slabs = []
        
        x = x0 + dx/2.0
        while x < x1:

            # Build slab for this x position.

            slab = [[y0, y1, self.background_mat]]

            for i in range(len(self.shapes)):
                
                ys = self.shapes[i].intersection_at_x(x)
                if len(ys) == 0:
                    continue

                ys0, ys1 = ys
                
                if ys0 < y0:
                    ys0 = 0     
                if ys1 > y1:
                    ys1 = y1

                def same(x,y): return abs(x-y) < 1e-5

                new_slab = []      
                j = 0
                while j < len(slab):

                    if (slab[j][1] < ys0) or same(slab[j][1], ys0) or \
                       (slab[j][0] > ys1) or same(slab[j][0], ys1):
                        new_slab.append(slab[j]) # No intersection.
                    else:
                        if not same(slab[j][0], ys0): # Old material pre.
                            new_slab.append([slab[j][0], ys0, slab[j][2]])

                        new_slab.append([ys0, ys1, self.shapes[i].mat])

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
            
            new_slabs.append(average_slabs(slabs, i, i_end))
            d.append((i_end - i) * dx)
                                              
            i = i_end

        d[-1] += (x1-x0) - sum(d)

        slabs = new_slabs

        # Create expression.

        e = Expression()

        for i in range(len(slabs)):
            
            e_slab = Expression()
            for j in range(len(slabs[i])):
                
                chunk_d = slabs[i][j][1] - slabs[i][j][0]
                if j == 0:
                    chunk_d += PML_bot*1j
                if j == len(slabs[i]) - 1:
                    chunk_d += PML_top*1j
                
                chunk_m = slabs[i][j][2]
                
                e_slab.add(chunk_m(chunk_d))
                
            s = Slab(e_slab)
            slab_cache.append(s)
            e.add(s(d[i]))

        # Add flipped scatterer.

        if add_flipped:
            for i in range(len(slabs)-1, -1, -1):
                e.add((slab_cache[i])(d[i]))

        return e
