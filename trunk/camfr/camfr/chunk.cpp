
/////////////////////////////////////////////////////////////////////////////
//
// File:     chunk.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000509
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "chunk.h"
#include "scatterer.h"
#include "icache.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// split_chunks
//
/////////////////////////////////////////////////////////////////////////////
  
void split_chunks(const vector<Chunk>& chunks,
                  vector<Chunk>* left_chunks, vector<Chunk>* right_chunks,
                  unsigned int i, Split_type split)
{
  // Check values.

  if ( (!left_chunks) || (!right_chunks) || (i>=chunks.size()) )
  {
    py_error("Error: invalid arguments in split_chunks.");
    exit (-1);
  }

  // Chunks before i.

  left_chunks->clear();
  for (unsigned int j=0; j<i; j++)
    left_chunks->push_back(chunks[j]);

  // Chunk i.
  
  right_chunks->clear();
  if (split == after_prop)
    left_chunks->push_back(chunks[i]);
  else
  {
    Chunk i_noprop(chunks[i].sc, 0.0);
    left_chunks->push_back(i_noprop);

    Scatterer* transparent = interface_cache.get_interface
      (chunks[i].sc->get_ext(), chunks[i].sc->get_ext());
    
    Chunk i_prop(transparent, chunks[i].d);
    right_chunks->push_back(i_prop);
  }

  // Chunks after i.

  for (unsigned int j=i+1; j<chunks.size(); j++)
    right_chunks->push_back(chunks[j]);
}
