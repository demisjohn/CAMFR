
/////////////////////////////////////////////////////////////////////////////
//
// File:     tracesorter.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     2001618
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef TRACESORTER_H
#define TRACESORTER_H

#include <vector>
#include "../math/linalg/linalg.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Tracesorter
//
//   When adding a row  [1.0 3.0 2.0], to a tracesorter and
//   subsequently a row [3.1 3.1 1.1], this last row will be reshuffled,
//   to give            [1.1 3.1 2.1].
//   The distance between each element and the corresponding one on the
//   row above will be minimised.
//   Also, no sudden turns will be allowed in each columns, except when
//   near a point in the turns_allowed list.
//   Useful when plotting e.g. bandstructures.
//
/////////////////////////////////////////////////////////////////////////////

class Tracesorter
{
  public:
  
    Tracesorter() {}

    void add_row(const cVector& row);

    void add_turning_point(const Real& r) {turning_points.push_back(r);}

    Complex operator()(int i, int j) const {return matrix[i-1](j);}

    int rows() const {return matrix.size();}
    int cols() const {return matrix[0].rows();}
    
  protected:

    bool bad_turn(int i, const Complex& c) const;
    int find_best_index(int i, const cVector& row, bool* placed,
                        bool bad_turns_allowed=false) const;

    vector<cVector> matrix;
    vector<Real> turning_points;
};



#endif
