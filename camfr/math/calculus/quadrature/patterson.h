
/////////////////////////////////////////////////////////////////////////////
//
// File:     patterson.h
// Author:   Peter.Bienstman@rug.ac.be
//           converted to C++ from the CACM algorithms
// Date:     20000320
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef PATTERSON_H
#define PATTERSON_H

#include "../function.h"

/////////////////////////////////////////////////////////////////////////////
//
// patterson
//
//   Integrates f between a and b using Patterson quadrature formulas.
//
//   Low-level formula. The higher level drivers are in patterson_quad.h.
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
//   If the precision is not achieved for k=max_k, *error_ptr will be true.
//
//   abs_error contains the difference between the last 2 rules invoked.
//
/////////////////////////////////////////////////////////////////////////////

Real patterson(RealFunction& f, Real a, Real b, Real eps,
               bool* error_ptr=NULL, unsigned int max_k=8,
               Real* abs_error=NULL);



#endif
