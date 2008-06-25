#! /usr/bin/env python

###########################################################################
#
# File:     geometry.py
# Author:   Peter.Bienstman@rug.ac.be
# Date:     20020204
# Version:  1.1
#
# Copyright (C) 2002 Peter Bienstman - Ghent University
#
############################################################################

from math import *
from camfr import *
from numpy import *

import ImageFile, Image, ImageDraw

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

    def compare(self,p):
        if self.x != p.x or self.y != p.y:
          return 1
        else:
          return 0

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
        self.c    = c
        self.r    = r
        self.mat  = mat
        self.type = "Circle"
        
    def intersection_at_x(self, x):
        D = self.r**2 - (x - self.c.x)**2
        if D > 0:
            return [ self.c.y - sqrt(D), self.c.y + sqrt(D) ]
        else:
            return []
    
    def compare(self,c):
      if c.type!=self.type:
        return 1
      if self.c.compare(c.c) or self.r != c.r or self.mat!=c.mat:
        return 1
      else:
        return 0  



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
        self.type = "Rectangle"
        
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

    def compare(self,r):
      if r.type!=self.type:
        return 1
      if self.p1.compare(r.p1) or self.p2.compare(r.p2) or self.mat != r.mat:
        return 1
      else:
        return 0



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
        self.type = "Square"
        
    def intersection_at_x(self, x):
        if (x < self.c.x - self.a/2.) or (x > self.c.x + self.a/2.):
            return []
        else:
            return [self.c.y - self.a/2., self.c.y + self.a/2.]

    def to_rectangle(self):
        return Rectangle(Point(self.c.x - self.a/2., self.c.y - self.a/2),
                         Point(self.c.x + self.a/2., self.c.y + self.a/2))
    
    def compare(self,s):
      if s.type != self.type:
        return 1
      elif self.c.compare(s.c) or self.a != s.a or self.mat != s.mat:
        return 1
      else:
        return 0



############################################################################
#
# class Triangle
#
############################################################################
              
class Triangle:

    def __init__(self, p1, p2, p3, mat):
        self.p = [p1, p2, p3]
        self.p.sort(sort_point)
        self.p1, self.p2, self.p3 = self.p
        self.mat = mat
        self.type = "Triangle"
        
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
  
    def compare(self,t):
      if t.type != self.type:
        return 1
      elif self.p1.compare(t.p1) or self.p2.compare(t.p2) or self.mat != t.mat:
        return 1
      else:
        return 0
      


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
# class Picture
#
#   If a Picture will be added to the Geometry, it can be converted 
#   in an expression. The Picture object has no intersection method 
#   and cannot be combined with other geometry objects.
#
############################################################################

class Picture(object):

    def __init__(self, name, n_min=0, n_max=2.5, materials=[]):
        self.imageName    = name
        self.camfr_mat    = [n_min, n_max]
        if len(materials):
           self.camfr_mat = materials
        self.camfr_slabs  = []
        self.mat_dict     = []
        self.expr         = Expression()

    def _free(self):
        self.camfr_slabs  = []
        self.mat_dict     = []
        self.expr         = Expression()

    ########################################################################
    #
    # img_to_expression
    #
    #  Rescales the image, makes a dictionary with all materials based
    #  on the colours, returns an expression.
    #
    ########################################################################

    def img_to_expression(self, z0, z1, dz=0, 
                          verbose=False, rescaling='ANTIALIAS'):

        # Rotate 270 degrees.                                           
        #                                       (0,0)                    
        #   x                                   .-------->x                
        #   ^                                   |                       
        #   |                                   |     270 rotated          
        #   |    camfr image        =>          |     PIL image    
        #   |                                   |               
        #   '------------------>z              \|/               
        #  (0,0)                                z                             

        im = Image.open(self.imageName).rotate(270)
        im = im.convert("L")               # Convert image to gray scale.
        xPixels, zPixels =  im.size[0], im.size[1]      
        self._free()


        # Get all dimensions.
        
        z = float(z1-z0)                   # Physical Length.
        if dz: 
            Zsteps = int(z/dz)             # Number of slabs.
        else: 
            Zsteps = zPixels               # Number of slabs: no of pixels.
            dz     = z/Zsteps

        factor = xPixels * 1.0 / zPixels                       
        x      = factor * z                # Physical height. 
        dx     = factor * dz
        Xsteps = int(factor * Zsteps)      # WARNING, may deform the img.  
     

        # Rescaling.
        
        if not(( Xsteps == im.size[0]) and (Zsteps == im.size[1])):
        
            if   (rescaling == 'NEAREST'):     
                im = rescale_nearest(Xsteps,Zsteps) 
            elif (rescaling == 'ANTIALIAS'):   
                im = rescale_antialias(im, Xsteps, Zsteps)      
            elif (rescaling == 'AVERAGE'):
                im = rescale_custom(im, x, z, Xsteps, Zsteps, 
                                         method='AVERAGE')  
            elif (rescaling == 'INV_AVERAGE'):
                im = rescale_custom(im, x, z, Xsteps, Zsteps, 
                                         method='INV_AVERAGE') 
            else:
                im = rescale_antialias(im, Xsteps, Zsteps)
            

        # Data on screen if requested.                                      
        
        if (verbose == 1): 
            print "                                                        "               
            print "IMAGE INFORMATION                                       "
            print "--------------------------------------------------------"
            print "Name          : ", self.imageName
            print "Original size : ", zPixels, xPixels , " (length, height)"
            print "New size      : ", Zsteps, Xsteps, " (lenght, height)   "
            im.save ('rescaledimage.png', 'png')     
            print "Image saved as:  rescaledimage.png                      "
            print "--------------------------------------------------------"
        
        colors = im.getcolors()            # e.g. 
                                           # colors= [(220, 0),(2330, 255)]: 
                                           # 220 * black, 2330 * white.
        self.mat_dict = makeRefractiveIndexList(self.camfr_mat, colors)  
        
        pix = im.load()                    # Matrix with pixels.
        slabs = []                         # List with unique slabs.
        stack = []                         # List refering to all slabs.
        
        
        # Collect all unique slabs.
        
        for z in range(Zsteps):            # With every new z, 
                                           # we scan a new slab.
          slab = []
          for x in range(Xsteps):          # Scan all materials of the slab.
             if len(slab):                 # Add the material to the slab.
                if (slab[-1][0] == pix[x,z]):
                   slab[-1][1] += dx
                else:
                   slab.append([pix[x,z],dx])        
             else:
               slab.append([pix[x,z],dx])
               
          try: 
              i = slabs.index(slab)        # Do we have this slab?
              if (i == stack[-1][0]):      # Same as last added?      
                  stack[-1][1] +=dz        # Make slab thicker.
              else:
                  stack.append([i,dz])     # Refer to same slab.
          except ValueError:               # This slab does not yet exist.
              slabs.append(slab)           # Add new unique slab.
                                           # Refer to new slab, length dz.   
              stack.append([len(slabs)-1,dz])  


        # Turn all unique slabs into camfr slabs
        
        for s in range(len(slabs)):        # Take slab: slabs[s]
            tmp_ex = Expression()
            for m in range(len(slabs[s])): # Take material: slabs[s][m].
                                           # Add Material(length) to expr. 
                tmp_ex.add( self.mat_dict[slabs[s][m][0]](slabs[s][m][1]) )
                
                                           # Add final slab expr to cache. 
            self.camfr_slabs.append(Slab(tmp_ex)) 
              
                                     
        # Make expression from the stack.
        
        for s in range(len(stack)):        # Take slab of stack.
                                           # Add camfr slab to expression.
            self.expr.add(self.camfr_slabs[stack[s][0]](stack[s][1]))    
       
       
        
        # Print all information.
        
        if (verbose == 1):
            print "                                                        "        
            print "MATERIAL INFORMATION                                    "
            print "--------------------------------------------------------"
            print self.mat_dict    
            print "--------------------------------------------------------"
            print "                                                        "        
            print "SLAB INFORMATION                                        "
            print "--------------------------------------------------------"
            print slabs
            print "--------------------------------------------------------"
            print "                                                        "            
            print "STACK INFORMATION                                       "
            print "--------------------------------------------------------"            
            print stack
            print "--------------------------------------------------------"
            
        return self.expr
    


############################################################################
#
# A number of rescaling algorithms.
#
#   nearest: wil not introduce new materials, but picks up the nearest
#   antialias: pil function, takes average. Same as custome TM.
#   custom: manipulate the picture at own whishes
#           camfr TM (PhC TE) : average
#           camfr TE (PhC TM) : inverse average   
#
############################################################################ 

def rescale_nearest(self, image, Xsteps, Zsteps):
    return image.resize((Xsteps,Zsteps),Image.NEAREST)    

def rescale_antialias(image, Xsteps, Zsteps):
    return image.resize((Xsteps,Zsteps),Image.ANTIALIAS )

def rescale_custom(image, x, z, Xsteps, Zsteps, method='AVARAGE'):
                       
    # Algorithm =
    # resize from fine grid x-z to X-Z 
    # Cel [X1Z1-X2Z2] covers [(xn,zm)-[(xn+j,zm+k)]
    # xndx < X1dX
    # sizes = [X1dX, (xn+1)dx, ..., (xn+j)dx]

    pix = image.load()
    dX = x / Xsteps                      # Physical lenghts of a x step.
    dZ = z / Zsteps                      # Physical lenghts of a z step.
    
    dx = x / image.size[0]
    dz = z / image.size[1]
    
    im = Image.new ("L", (Xsteps,Zsteps))
    L = list(im.getdata())
    Ln = 0
    for Z in range (Zsteps):
      for X in range (Xsteps):
      
        x1 = int(floor(X*dX/dx))         # First x-matching value.
        j = int(floor(dX/dx))            # Number of pixels in x direction.
        z1 = int(floor(Z*dZ/dz))         # First z-matching value
        k = int(floor(dZ/dz))            # Number of pixels in z direction.

        #list with x-lenghts
        LX = [(x1+1)*dx - X*dX] 
        for i in range(1,j): 
            LX.append(dx)
        LX.append((X+1)*dX - (x1+j)*dx) 

        #list with z-lenghts
        LZ = [(z1+1)*dz - Z*dZ] 
        for i in range(1,k): 
            LZ.append(dz)
        LZ.append((Z+1)*dZ - (z1+j)*dz)

        #calculate average:
        average=0
        newcolor=0
        for _z in range(k+1):
          for _x in range(j+1):
              opp = LZ[_z]*LX[_x]
              if (method=='AVERAGE'): 
                  average += opp * pix[x1+_x,z1+_z]   
              else:                        # Inverse average.
                  average += opp / (1+pix[x1+_x,z1+_z])   
        
        if (method=='AVERAGE'):
            newcolor = floor(average/(dX*dZ))  
        else:                              # Inverse average.
            newcolor = floor((dX*dZ)/average)-1               
        L[Ln] =(newcolor)
        Ln += 1
        
    im.putdata(L, 1.0, 0.0)  # List, scale, offset           
    return im

    

############################################################################
#
# Make material dictionary.
# Every colour is one material.
#
# Assumption: 0 is highest, 255 lowest material. Linear in between.
# Output example =
#       mat_dict = {0: Isotropic n=(2.8,0), 255: Isotropic n=(1,0)}
#
#
############################################################################ 

def makeRefractiveIndexList(camfr_mat, colors): 
    minn = camfr_mat[0]
    maxn = camfr_mat[-1]
    return dict([(x[1],Material(maxn-(maxn-minn)*x[1]/255)) for x in colors])



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
############################################################################

slab_cache = []

class Geometry:

    def __init__(self, background_mat = Material(1.0)):
        self.background_mat = background_mat
        self.shapes = []


    def add(self, s):
        self.shapes.append(s)

    def __iadd__(self, s):
        self.shapes.append(s)
        return self

    def to_expression(self, x0, x1, dx=0, y0=0, y1=0, dy=0, add_flipped=0,
                      calc_average=0, verbose=False, rescaling='ANTIALIAS'):

        # Test whether geometry contains a picture.
        
        picture, place = True, 0
        for i in range(len(self.shapes)):    
            if (type(self.shapes[i]) == Picture):
                picture, place = True, i

        # Geometry does contain a picture.
        
        if picture:
            
            return self.shapes[place].img_to_expression( 
                x0, x1, dx, verbose, rescaling)

        # Geometry does not contain a picture.        
        else:
            
            if (x0 > x1) or (y0 > y1):
                print "Error: Invalid boundaries for to_expression."
                raise IndexError
    
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
                        ys0 = y0     
                    if ys1 > y1:
                        ys1 = y1
    
                    def same(y0,y1): return abs(y0-y1) < 1e-6
    
                    new_slab = []      
                    j = 0
                    while j < len(slab):
    
                        if (slab[j][1] < ys0) or same(slab[j][1], ys0) or \
                           (slab[j][0] > ys1) or same(slab[j][0], ys1) or \
                           (abs(ys0-ys1) < .001*dy):
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
                    while (i_end<len(slab)) and (slab[i_end][2]==slab[i][2]):
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
    
            # Create expression.
    
            e = Expression()
            
            eps_average = inv_eps_average = 0
            
            for i in range(len(slabs)):
                
                e_slab = Expression()
                
                for j in range(len(slabs[i])):
                    
                    chunk_d = slabs[i][j][1] - slabs[i][j][0]                
                    chunk_m = slabs[i][j][2]
    
                    eps_average     +=    chunk_m.epsr() * chunk_d * d[i]
                    inv_eps_average += 1./chunk_m.epsr() * chunk_d * d[i]
                    
                    e_slab.add(chunk_m(chunk_d))
                    
                s = Slab(e_slab)
    
                if verbose == True:
                    print e_slab
                    s.calc()
                    pretty_print(s)
                    plot(s)
                    print "------------"
       
                slab_cache.append(s)
                e.add(s(d[i]))
    
            eps_average     /= (x1-x0) * (y1-y0)
            inv_eps_average /= (x1-x0) * (y1-y0)        
    
            # Add flipped scatterer.
    
            if add_flipped:
                for i in range(len(slabs)):
                    e.add((slab_cache[-i-1])(d[-i-1]))
    
            if calc_average == False:
                return e
            else:
                return e, eps_average, inv_eps_average

def pretty_print(s):
    for i in range(s.N()):
        nr = s.mode(i).n_eff().real
        ni = s.mode(i).n_eff().imag

        if abs(nr) < 1e-6:
            nr = 0
            
        if abs(ni) < 1e-6:
            ni = 0

        print i, nr, ni

