#!/usr/bin/env python

####################################################################
#
# Illustrates how to define complicated structures using Geometry.
#
####################################################################

from camfr import *

set_lambda(1)
set_N(40)

air = Material(1)
mat = Material(3)

g = Geometry(air)

g += Rectangle(Point(0.0,-0.5), Point(2.0, 0.5), mat)
g += Triangle (Point(2.0, 0.5), Point(2.0,-0.5), Point(3.0, 0.0), mat)
g += Circle   (Point(4.5, 0.0), 0.5, mat)

prop0,  prop1,  d_prop  =  0.0, 5.0, 0.50
trans0, trans1, d_trans = -2.5, 2.5, 0.01

set_lower_PML(-0.1)
set_upper_PML(-0.1)

exp = g.to_expression(prop0,  prop1,  d_prop,
                      trans0, trans1, d_trans)

s = Stack(exp)

s.calc()

print s.R12(0,0)
