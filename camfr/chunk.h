
/////////////////////////////////////////////////////////////////////////////
//
// File:     chunk.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000509
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CHUNK_H
#define CHUNK_H

#include <vector>
#include "defs.h"

/////////////////////////////////////////////////////////////////////////////
//
// STRUCT: Chunk
//
//   Low-level structure, describing a Scatterer including propagation
//   over a distance d in its exit medium.
//
/////////////////////////////////////////////////////////////////////////////

class Scatterer; // forward definition: see scatterer.h

struct Chunk
{
    Chunk(Scatterer* sc_, Complex d_=0.0) : sc(sc_), d(d_) {}

    Scatterer* sc;
    Complex    d;
};



/////////////////////////////////////////////////////////////////////////////
//
// split_chunks
//
//   Splits a vector of chunks into a left and a right part.
//   If split=after_prop, the left chunks include chunks 0 through i.
//   If split=before_prop, the interface part of chunk i is moved to the
//   left chunks and the propagation part of medium i is moved to the
//   right chunks.
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {before_prop, after_prop} Split_type;
  
void split_chunks(const vector<Chunk>& chunks,
                  vector<Chunk>* left_chunks, vector<Chunk>* right_chunks,
                  unsigned int i, Split_type split=after_prop);



#endif

