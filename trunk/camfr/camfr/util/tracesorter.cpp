
/////////////////////////////////////////////////////////////////////////////
//
// File:     tracesorter.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     2001618
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "tracesorter.h"

/////////////////////////////////////////////////////////////////////////////
//
// Tracesorter::bad_turn
//
//  Determine whether the addition of c at index 'i' in a new row
//  breaks the trend set by the previous two rows.
//
/////////////////////////////////////////////////////////////////////////////

bool Tracesorter::bad_turn(int i, const Complex& c) const
{
  // Enough data to determine trend?

  if (matrix.size() <= 1)
    return false;

  Complex Delta01 =                c             - matrix[matrix.size()-1](i);
  Complex Delta12 =   matrix[matrix.size()-1](i) - matrix[matrix.size()-2](i);
  
  // At a turning point?

  for (unsigned j=0; j<turning_points.size(); j++)
    if (abs(real(c)-turning_points[j]) < 0.05)
      return false;
  
  // Delta's

  Complex delta01 =                c             - matrix[matrix.size()-1](i);
  Complex delta12 =   matrix[matrix.size()-1](i) - matrix[matrix.size()-2](i);

  if (    (abs(arg(delta12)-arg(delta01)     ) > pi/4.)
       && (abs(arg(delta12)-arg(delta01)-2*pi) > pi/4.) )
    return true;
  else
    return false;
}



/////////////////////////////////////////////////////////////////////////////
//
// Tracesorter::find_best_index
//
/////////////////////////////////////////////////////////////////////////////

int Tracesorter::find_best_index(int i, const cVector& row, bool* placed,
                                 bool bad_turns_allowed=false) const
{
  // Initialise.

  const int N = row.rows();
  
  const cVector& prev_row(matrix[matrix.size()-1]);

  int first_valid;
  for (first_valid=1; first_valid<=N; first_valid++)
    if (!placed[first_valid] && !bad_turn(i, row(first_valid)))
      break;

  Real best_distance = abs(prev_row(i)-row(first_valid));
  int best_index = first_valid;

  // Loop in row.

  for (int j=first_valid+1; j<=N; j++)
  {
    Real distance = abs(prev_row(i)-row(j));

    if ( (distance < best_distance) && (!placed[j]) )
      if (      bad_turns_allowed
           || (!bad_turns_allowed) && (!bad_turn(i,row(j))) )
      {
        best_distance = distance;
        best_index = j;
      }
  }

  return best_index;
}



/////////////////////////////////////////////////////////////////////////////
//
// Tracesorter::add_row
//
/////////////////////////////////////////////////////////////////////////////

void Tracesorter::add_row(const cVector& row)
{
  // First row.
  
  if (matrix.size() == 0)
  {
    matrix.push_back(row);
    return;
  }

  // Other row.

  const int N = row.rows();

  bool placed[N+1];
  for (int i=1; i<=N; i++)
    placed[i] = false;
  
  cVector ordered_row(N, fortranArray);
  
  for (int i=1; i<=N; i++)
  {
    int best_index = find_best_index(i,row,placed,false);

    if (best_index==N+1)
      best_index = find_best_index(i,row,placed,true);
    
    ordered_row(i) = row(best_index);
    placed[best_index] = true;
  }

  matrix.push_back(ordered_row);
}
