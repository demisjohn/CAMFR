
/////////////////////////////////////////////////////////////////////////////
//
// File:     index.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000204
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef INDEX_H
#define INDEX_H

#include <vector>
#include "../defs.h" // for Limit

/////////////////////////////////////////////////////////////////////////////
//
// index_lookup
//
//   Performs a lookup of an index in a piecewise constant function:
//
//                         x <= x_table[0](minus) : 0
//     x_table[0](plus) <= x <= x_table[1](minus) : 1
//     ...
//
//   Note: for complex values, only the real part is considered.
//
/////////////////////////////////////////////////////////////////////////////

unsigned int index_lookup(const Complex& x, Limit limit,
                          const vector<Complex>& x_table);



#endif
