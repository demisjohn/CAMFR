
/////////////////////////////////////////////////////////////////////////////
//
// File:     roottest.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010322
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "allroots.h"

/////////////////////////////////////////////////////////////////////////////
//
// Test function F
//
/////////////////////////////////////////////////////////////////////////////

class F : public ComplexFunction
{
  public:
    Complex operator()(const Complex& z)
      {counter++; return 1./z/z/sqrt(z)*sin(pi*z);}
};



/////////////////////////////////////////////////////////////////////////////
//
// Test function
//
/////////////////////////////////////////////////////////////////////////////

int main()
{
  F f;
  
  //allroots(f, -5-1*I, 1+1e-9+1*I);

  Real eps = 0.2;
  
  vector<Complex> roots = N_roots(f, 10, 2-eps-(1+eps)*I, 3+eps+(5+eps)*I);

  cout << "Iterations: " << f.times_called() << endl;

  cout << "All roots: " << endl;
  for (unsigned int i=0; i<roots.size(); i++)
    cout << i << " " << roots[i] << endl;
  
  return 0;
}
