
/////////////////////////////////////////////////////////////////////////////
//
// File:     fourier.h
// Author:   Peter.Bienstman@UGent.be, Peter.Debackere@intec.UGent.be
// Date:     20050208
// Version:  1.0
//
// Copyright (C) 2005-2007 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef FOURIER_H
#define FOURIER_H

#include <vector>
#include "../../linalg/linalg.h"

/////////////////////////////////////////////////////////////////////////////
//
// fourier
//
//  Returns fourier expansion of order m = -M,...,M of a piecewise 
//  continuous function f given by:
//
//    disc[0] -> disc[1]   : f[0]
//    ...
//    disc[n] -> disc[n+1] : f[n]
//
//  Basis functions are exp(j.m.2.pi/d.x). (d can be overridden).
//
//  Also has the option to extend the profile by adding mirror image at x=0.
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier(const std::vector<Complex>& f, 
                const std::vector<Complex>& disc, int M,
                const Complex* d=0, bool extend=false);



/////////////////////////////////////////////////////////////////////////////
//
// fourier_ASR
//
//  Returns fourier expansion of order m = -M,...,M of a piecewise 
//  continuous function f given by:
//
//    disc[0] -> disc[1]   : f[0]
//    ...
//    disc[n] -> disc[n+1] : f[n]
//
//  and multiplied by a parametric continouos function given by
//  
//  Basis functions are exp(j.m.2.pi/d.x). (d can be overridden).
//
//  Also has the option to extend the profile by adding mirror image at x=0. 
//  Also has the option to stretch the entire structure in order to reduce
//  the discontinouos jumps of the dielectric constant
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier_ASR(const std::vector<Complex>& f,
                    const std::vector<Complex>& disc, 
                    const std::vector<Complex>& udisc,
                    int M, bool extend, const Real eta);



////////////////////////////////////////////////////////////////////////////
//
// fourier_ASR_pseudo
//
//  This is the 1D fourier integration in the new ASR basis. The 
//  difference with the previous function is that this one takes u_disc 
//  as input, and the length of the structure, and some booleans
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier_ASR_pseudo(const std::vector<Complex>& f,
                           const std::vector<Complex>& disc, 
                           const std::vector<Complex>& udisc,
                           int M, const Complex* d, 
                           bool extend, bool begin, bool end,
                           const Real eta);



/////////////////////////////////////////////////////////////////////////////
//
// fourier_2D
//
//  Returns 2D fourier expansion of order m_x = -M,...,M and n_y = -N,...,N 
//  of a piecewise constant dielectric profile f consisting of:
//
//    disc_x[0] -> disc_x[1]   : strip_y[0]
//    ...
//    disc_x[n] -> disc_x[n+1] : strip_y[n]
//
//  strip_y[i] is given by:
//
//    disc_y[i][0] -> disc_y[i][1]   : f[i][0]
//    ...
//    disc_y[i][n] -> disc_y[i][n+1] : f[i][n]
//
//  Basis functions are exp(j.m.2.pi/d.x).exp(j.n.2.pi/d.y).
//
//  Also has the option to extend the profile by adding mirror images
//  at the x=0 and y=0 planes.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D(const std::vector<Complex>& disc_x,
                   const std::vector<std::vector<Complex> >& disc_y,
                   const std::vector<std::vector<Complex> >& f,
                   int M, int N, bool extend=false);



/////////////////////////////////////////////////////////////////////////////
//
// fourier_2D_ASR
//
//  Returns 2D fourier expansion of order m_x = -M,...,M and n_y = -N,...,N 
//  of a piecewise constant dielectric profile f.
//	This version is used in the adaptive spatial resolution formulation
//
//  The profile f consists of:
//
//    disc_x[0] -> disc_x[1]   : strip_y[0]
//    ...
//    disc_x[n] -> disc_x[n+1] : strip_y[n]
//
//  strip_y[i] is given by:
//
//    disc_y[i][0] -> disc_y[i][1]   : f[i][0]
//    ...
//    disc_y[i][n] -> disc_y[i][n+1] : f[i][n]
//
//  Basis functions are exp(j.m.2.pi/d.x).exp(j.n.2.pi/d.y).
//  
//  Also has the option to extend the profile by adding mirror images
//  at the x=0 and y=0 planes.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D_ASR(const std::vector<Complex>& disc_x,
                       const std::vector<std::vector<Complex> >& disc_y,
                       const std::vector<Complex>& disc_u,
                       const std::vector<std::vector<Complex> >& disc_v,
                       const std::vector<std::vector<Complex> >& f,
                       int M, int N, bool extend, 
                       const Real eta);



/////////////////////////////////////////////////////////////////////////////
//
//  Returns 'split' 2D fourier expansion, i.e. transform in x, invert
//  toeplitz matrix, followed by a transform in y.
//
//  Also has the option to extend the profile by adding mirror images
//  at the x=0 and y=0 planes.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D_split(const std::vector<Complex>& disc_x,
                         const std::vector<std::vector<Complex> >& disc_y,
                         const std::vector<std::vector<Complex> >& f,
                         int M, int N, bool extend=false);



/////////////////////////////////////////////////////////////////////////////
//
// fourier_2D_split_ASR
//
//  Returns 'split' 2D fourier expansion, i.e. transform in x or y, invert
//  toeplitz matrix, followed by a transform in y or x. This version 
//  can calculate both, and uses the parametric formulation for the 
//  coordinate axes.
//
//  Also has the option to extend the profile by adding mirror images
//  at the x=0 and y=0 planes.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D_split_ASR(const std::vector<Complex>& disc_x,
                             const std::vector<std::vector<Complex> >& disc_y,
                             const std::vector<Complex>& disc_u,
                             const std::vector<std::vector<Complex> >& disc_v,
                             const std::vector<std::vector<Complex> >& f,
                             int M, int N, bool extend,
                             const Real eta);



#endif
