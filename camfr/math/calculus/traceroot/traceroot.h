
/////////////////////////////////////////////////////////////////////////////
//
// File:          traceroot.h
// Author:        Peter.Bienstman@rug.ac.be
// Date:          19990907
// Version:       1.2
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef TRACEROOT_H
#define TRACEROOT_H

#include <vector>
#include <string>
#include "../function.h"

/////////////////////////////////////////////////////////////////////////////
//
// traceroot
//
//  Start from the roots estimate1 of function f with parameters params1 
//  and return the roots of f with params params2, by gradually changing 
//  the parameters.
//  Zeros that the algorithm should not converge to can be specified in
//  forbiddenzeros.
//  The initial resolution can be specified, but the sweep parameters
//  are dynamically updated for maximum performace/precision.
//  The path of zeros can be written to a file 'fname'
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> traceroot(vector<Complex>&     estimate1,
                          Function1D<Complex>& f, 
                          vector<Complex>&     params1,
                          vector<Complex>&     params2,
                          vector<Complex>&     forbiddenzeros,
                          int resolution       = 1,
                          string* fname        = NULL);

vector<Complex> traceroot(vector<Real>&        estimate1,
                          Function1D<Complex>& f, 
                          vector<Complex>&     params1,
                          vector<Complex>&     params2,
                          vector<Complex>&     forbiddenzeros,
                          int resolution       = 1,
                          string* fname        = NULL);



/////////////////////////////////////////////////////////////////////////////
//
// traceroot_chunks
//
//  Same as traceroot, but splits the list of initial zeros in separate
//  chunks with length chunk_length + 2*overlap. Each chunk overlaps
//  'overlap' elements with the previous and subsequent chunk. 
//  This can be faster, since the optimal sweep parameters for each chunk
//  can be determined independently.
//  Possible danger: if one chunk's zero converges to another chunk's zero,
//  this will go unnoticed.
//  The part of the zeros in each chunk can be written to 'fname1',
//  'fname2', ..
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> traceroot_chunks(vector<Complex>&     estimate1,
                                 Function1D<Complex>& f, 
                                 vector<Complex>&     params1,
                                 vector<Complex>&     params2,
                                 vector<Complex>&     forbiddenzeros,
                                 int resolution       = 1,
                                 int chunk_length     = 48,
                                 int overlap          = 4,
                                 string* fname        = NULL);

vector<Complex> traceroot_chunks(vector<Real>&        estimate1,
                                 Function1D<Complex>& f, 
                                 vector<Complex>&     params1,
                                 vector<Complex>&     params2,
                                 vector<Complex>&     forbiddenzeros,
                                 int resolution       = 1,
                                 int chunk_length     = 48,
                                 int overlap          = 4,
                                 string* fname        = NULL);



#endif
