
/////////////////////////////////////////////////////////////////////////////
//
// File:     patterson_z_n.h
// Author:   Peter.Bienstman@rug.ac.be
//           converted to C++ from the CACM algorithms
// Date:     20000320
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef PATTERSON_Z_N_H
#define PATTERSON_Z_N_H

#include <vector>
#include "../function.h"

/////////////////////////////////////////////////////////////////////////////
//
// patterson_z_n
//
//   Integrates z^n / f(z) on a straight segment between a and b
//   using Patterson quadrature formulas, for n=0..M.
//
//   Low level routine, high level driver is patterson_quad_z_n.
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
//   The number of results returned can be < M+1, if the presence of
//   round-off errors is detected by a 'relaxed overflow condition' :
//   machine_eps()*abs(f.z^n) > mu.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> patterson_z_n(ComplexFunction& f,
                              const Complex& a, const Complex& b, int M,
                              Real eps, Real mu, bool* error_ptr=NULL,
                              unsigned int max_k=8,
                              vector<Complex>* abs_error=NULL);



/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad_z_n
//
//   Integrates z^n / f(z) on a straight segment between a and b
//   using Patterson quadrature formulas, for n=0..M.
//
//   If the precision is not achieved for k=max_k, adaptive refinement of
//   the interval is done.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> patterson_quad_z_n(ComplexFunction& f,
                                   const Complex& a, const Complex& b, int M,
                                   Real eps, Real mu, unsigned int max_k=8);



#endif
