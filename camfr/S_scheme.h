
/////////////////////////////////////////////////////////////////////////////
//
// File:     S_scheme.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991203
// Version:  1.1
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef S_SCHEME_H
#define S_SCHEME_H

#include <vector>
#include "scatterer.h"
#include "chunk.h"

/////////////////////////////////////////////////////////////////////////////
//
// S-scheme for calculating R and T of stack of chunks.
// This scattering matrix scheme is stable in presence of evanescent waves.
//
// Results will be stored in the Scatterer argument.
//
// 'chunks' should always by passed by reference, in order for the
//  parameter sweeping mechanism to work.
//
// The exit medium of each chunk should match the incidence medium of
// the next one.
//
// Different variants are optimised for structures with diagonal matrices
// or monomode structures.
//
//
/////////////////////////////////////////////////////////////////////////////

void S_scheme(const vector<Chunk>& chunks, DenseScatterer* result);
void S_scheme(const vector<Chunk>& chunks,  DiagScatterer* result);
void S_scheme(const vector<Chunk>& chunks,  MonoScatterer* result);



#endif
