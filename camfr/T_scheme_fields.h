
/////////////////////////////////////////////////////////////////////////////
//
// File:     T_scheme_fields.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010702
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef T_SCHEME_FIELDS_H
#define T_SCHEME_FIELDS_H

#include <vector>
#include "scatterer.h"
#include "chunk.h"
#include "field.h"

/////////////////////////////////////////////////////////////////////////////
//
// Returns field expansion in a set of chunks for a given fw/bw field
// in field[0].
//
//   element 0: field before chunk[0]
//   element 1: field after interface in chunk[0]
//   element 2: field after propagation in chunk[0]
//   ...
//
// Chunks is supposed to contain scatteres that discribe only interfaces,
// such that it is safe to invert it T21 matrix.
// Stability is improved by factoring out exponentials.
//
/////////////////////////////////////////////////////////////////////////////

void T_scheme_fields
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field);



#endif
