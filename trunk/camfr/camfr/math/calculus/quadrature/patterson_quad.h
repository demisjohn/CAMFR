
/////////////////////////////////////////////////////////////////////////////
//
// File:     patterson_quad.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000320
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef PATTERSON_QUAD_H
#define PATTERSON_QUAD_H

#include "patterson.h"

/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad
//
//   Integrates f between a and b using Patterson quadrature formulas.
//
//   High level driver, using low level routines in patterson.h.
//
//   Succesive rules with index k are compared to see if the requested
//   relative precision 'eps' is achieved.
//
//     k=1 :   1-point rule
//     k=2 :   3-point rule
//     k=3 :   7-point rule
//     k=4 :  15-point rule
//     k=5 :  31-point rule
//     k=6 :  63-point rule
//     k=7 : 127-point rule
//     k=8 : 255-point rule
//
//   If the precision is not achieved for k=max_k, adaptive refinement of
//   the interval is done.
//
/////////////////////////////////////////////////////////////////////////////

Real patterson_quad(RealFunction& f, Real a, Real b,
                    Real eps, unsigned int max_k=8);



/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad_non_adapt
//
//   Same, but with non-adaptive refinement of the interval.
//
/////////////////////////////////////////////////////////////////////////////

Real patterson_quad_non_adapt(RealFunction& f, Real a, Real b,
                              Real eps, unsigned int max_k=8);



#endif
