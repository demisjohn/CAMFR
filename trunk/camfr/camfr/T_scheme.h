
/////////////////////////////////////////////////////////////////////////////
//
// File:     T_scheme.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000213
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef T_SCHEME_H
#define T_SCHEME_H

#include "chunk.h"

/////////////////////////////////////////////////////////////////////////////
//
// T-scheme for calculating forward and backward waves at the top end of a
// stack of chunks, given a set of forward and backward waves (fw_in,
// bw_in) at the bottom of the stack.
//
//   fw_out = fw_out_scaled * scaling_factor
//   bw_out = bw_out_scaled * scaling_factor
//
// 'chunks' should always by passed by reference, in order for the
//  parameter sweeping mechanism to work.
//
// The exit medium of each chunk should match the incidence medium of
// the next one.
//
// Currently only implemented for MonoScatterers.
//
// If in fw/bw a vector<> with a single element is provided, this is treated
// as the expansion of the incident fw/bw field. Upon exit, fw/bw will
// contain the fields at each of the interfaces.
//
/////////////////////////////////////////////////////////////////////////////

void T_scheme
  (const vector<Chunk>& chunks,vector<Complex>* fw, vector<Complex>* bw);




#endif
