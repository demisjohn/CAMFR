
/////////////////////////////////////////////////////////////////////////////
//
// File:     S_scheme_fields_old.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000612
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef S_SCHEME_FIELDS_OLD_H
#define S_SCHEME_FIELDS_OLD_H

#include <vector>
#include "scatterer.h"
#include "chunk.h"
#include "field.h"

/////////////////////////////////////////////////////////////////////////////
//
// Returns field expansion in a set of chunks for a given excitation field
// in field[0]. No excitation is assumed from the back side.
//
//   element 0: field before chunk[0]
//   element 1: field after interface in chunk[0]
//   element 2: field after propagation in chunk[0]
//   ...
//
// All the submatrices are calculated using the S scheme, but the excitation
// fields can be given in a T formalism (faster, less accurate) or an S
// formalism (slower, but accurate).
//
/////////////////////////////////////////////////////////////////////////////

void S_scheme_fields_S_old
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field);



#endif
