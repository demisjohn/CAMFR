
/////////////////////////////////////////////////////////////////////////////
//
// File:     S_scheme_fields.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000612
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef S_SCHEME_FIELDS_H
#define S_SCHEME_FIELDS_H

#include <vector>
#include "scatterer.h"
#include "chunk.h"
#include "field.h"

/////////////////////////////////////////////////////////////////////////////
//
// Returns field expansion in a set of chunks for a given entrance field
// in field[0].
//
//   element 0: field before chunk[0]
//   element 1: field after interface in chunk[0]
//   element 2: field after propagation in chunk[0]
//   ...
//
// All the submatrices are calculated using the S scheme, but the excitation
// fields can be given in a T formalism (obsolete, less accurate) or an S
// formalism.
//
// For the latter, a backward incident field form the right can also be
// included.
//
/////////////////////////////////////////////////////////////////////////////

void S_scheme_fields_T
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field);

void S_scheme_fields_S
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field,
   cVector* inc_right_bw=NULL);



#endif
