
/////////////////////////////////////////////////////////////////////////////
//
// File:     minimum.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000125
// Version:  1.1
//
// Copyright (C) 1999-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef MINIMUM_H
#define MINIMUM_H

#include "../function.h"

/////////////////////////////////////////////////////////////////////////////
//
// brent_minimum
//
//   Finds minimum of real function f(x) in interval [ax,bx] with
//   precision eps.
//
/////////////////////////////////////////////////////////////////////////////

Real brent_minimum(Function1D<Real>& f, Real ax, Real bx, Real eps=1e-13);



/////////////////////////////////////////////////////////////////////////////
//
// brent_minimum
//
//   Same as above, but can handle different intervals at the same time.
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Real> brent_minimum(Function1D<Real>& f, std::vector<Real>& Ax,
                                std::vector<Real>& Bx, Real eps=1e-13);



/////////////////////////////////////////////////////////////////////////////
//
// bracket_all_minima
//
//   Finds all the intervals [Ax,Bx] containing the minima of the function f
//   in the interval [ax,bx].
//   The region [ax,bx] is divided in steps with size dx. It is checked if
//   the same number of minima is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

void bracket_all_minima(Function1D<Real>& f, Real ax, Real bx,
                        std::vector<Real>& Ax,std::vector<Real>& Bx, Real dx,
                        int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// bracket_N_minima
//
//   Finds the intervals N [Ax,Bx] containing the first N minima of the
//   function f larger than ax.
//   The region x>ax is scanned in steps with size dx. It is checked if the
//   same number of minima is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

void bracket_N_minima(Function1D<Real>& f, Real ax, int N, 
                      std::vector<Real>& Ax, std::vector<Real>& Bx, 
                      Real dx, int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// brent_all_minima
//
//   Finds all minima of real function f(x) in interval [ax,bx] with
//   precision eps.
//   The [ax,bx] is divided in steps with size dx. It is checked if the
//   same number of minima is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Real> brent_all_minima(Function1D<Real>& f, Real ax, Real bx, 
                                   Real dx, Real eps=1e-13, int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// brent_N_minima
//
//   Finds N minima of real function f(x) in interval with lower bound
//   ax and precision eps.
//   The interval is divided in steps with size dx. It is checked if the
//   same number of minima is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Real> brent_N_minima(Function1D<Real>& f, Real ax, int N, Real dx,
                                 Real eps=1e-13, int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// brent_refine_minima
//
//   Finds minima of real function f(x) in intervals around initial
//   estimates [(1-delta_x)*x[i],(1+delta_x)*x[i]] with precision eps.
//   The intervals are divided in steps with size dx. It is checked if the
//   same number of minima is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Real> brent_refine_minima
  (Function1D<Real>& f, std::vector<Real>& x, Real delta_x, Real dx,
   Real eps=1e-13, int sec_level=0);



#endif


