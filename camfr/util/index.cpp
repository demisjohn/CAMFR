
/////////////////////////////////////////////////////////////////////////////
//
// File:     index.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000204
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "index.h"

/////////////////////////////////////////////////////////////////////////////
//
// index_lookup
//
/////////////////////////////////////////////////////////////////////////////

unsigned int index_lookup(const Complex& x, Limit limit,
                          const vector<Complex>& x_table)
{
  for (unsigned int i=0; i<x_table.size(); i++)
  {
    if (abs(real(x) - real(x_table[i])) < 1e-10)
      return (limit == Min) ? i : i+1;
    
    if (real(x) < real(x_table[i]))
      return i;
  }
  
  return x_table.size();
}
