
/////////////////////////////////////////////////////////////////////////////
//
// File:     polyroot_test.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000327
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include "polyroot.h"

/////////////////////////////////////////////////////////////////////////////
//
// Main driver
//
/////////////////////////////////////////////////////////////////////////////

int main()
{
  cout << setprecision(15);

  // Polynomial with roots 1..10.
    
  vector<Complex> p1;
  p1.push_back(        1.);
  p1.push_back(      -55.);
  p1.push_back(     1320.);
  p1.push_back(   -18150.);
  p1.push_back(   157773.);
  p1.push_back(  -902055.);
  p1.push_back(  3416930.);
  p1.push_back( -8409500.);
  p1.push_back( 12753576.);
  p1.push_back(-10628640.);
  p1.push_back(  3628800.);
  
  vector<Complex> r1 = polyroot(p1);
  
  for (unsigned int i=0; i<r1.size(); i++)
    cout << i << " " << r1[i] << endl;
  
  return 0;
}
