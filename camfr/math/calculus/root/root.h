
/////////////////////////////////////////////////////////////////////////////
//
// File:     root.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990513
// Version:  1.1
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_H
#define ROOT_H

#include <vector>
#include "../function.h"

/////////////////////////////////////////////////////////////////////////////
//
// brent_root
//
//   Finds root of real function f(x) in interval [ax,bx] with
//   precision eps.
//
/////////////////////////////////////////////////////////////////////////////

Real brent_root(Function1D<Real>& f, Real ax, Real bx, Real eps=1e-13);



/////////////////////////////////////////////////////////////////////////////
//
// brent_root
//
//  Same as above, but can handle different intervals at the same time.
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_root(Function1D<Real>& f, vector<Real>& Ax,
                        vector<Real>& Bx, Real eps=1e-13);



/////////////////////////////////////////////////////////////////////////////
//
// bracket_all_roots
//
//   Finds all the intervals [Ax,Bx] containing the roots of the function f
//   in the interval [ax,bx].
//   The region [ax,bx] is divided in steps with size dx. It is checked if
//   the same number of zeros is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

void bracket_all_roots(Function1D<Real>& f, Real ax, Real bx,
                       vector<Real>& Ax,vector<Real>& Bx, Real dx,
                       int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// bracket_N_roots
//
//   Finds the intervals N [Ax,Bx] containing the first N roots of the
//   function f larger than ax.
//   The region x>ax is scanned in steps with size dx. It is checked if the
//   same number of zeros is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

void bracket_N_roots(Function1D<Real>& f, Real ax, int N, vector<Real>& Ax,
                     vector<Real>& Bx, Real dx, int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// brent_all_roots
//
//   Finds all roots of real function f(x) in interval [ax,bx] with
//   precision eps.
//   The [ax,bx] is divided in steps with size dx. It is checked if the
//   same number of zeros is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_all_roots(Function1D<Real>& f, Real ax, Real bx, Real dx,
                             Real eps=1e-13, int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// brent_N_roots
//
//   Finds N roots of real function f(x) in interval with lower bound
//   ax and precision eps.
//   The interval is divided in steps with size dx. It is checked if the
//   same number of zeros is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_N_roots(Function1D<Real>& f, Real ax, int N, Real dx,
                           Real eps=1e-13, int sec_level=0);



/////////////////////////////////////////////////////////////////////////////
//
// brent_refine_roots
//
//   Finds roots of real function f(x) in intervals around initial
//   estimates [(1-delta_x)*x[i],(1+delta_x)*x[i]] with precision eps.
//   The intervals are divided in steps with size dx. It is checked if the
//   same number of roots is found using the larger step dx.2^sec_level.
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_refine_roots(Function1D<Real>& f, vector<Real>& x,
                                Real delta_x, Real dx,
                                Real eps=1e-13, int sec_level=0);




#endif
