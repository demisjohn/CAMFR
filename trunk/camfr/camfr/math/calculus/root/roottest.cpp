
/////////////////////////////////////////////////////////////////////////////
//
// File:     roottest.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990514
// Version:  1.0
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include "../../bessel/bessel.h"
#include "root.h"

/////////////////////////////////////////////////////////////////////////////
//
// Simple function object for roottester: Bessel function
//
/////////////////////////////////////////////////////////////////////////////

class F : public Function1D<Real>
{
  public:
    
    Real operator()(const Real& x) 
    {
      counter++;
      return J(1,x);
    }   
};



/////////////////////////////////////////////////////////////////////////////
//
// roottester
//
/////////////////////////////////////////////////////////////////////////////

int main()
{  
  F f;

  vector<Real> zeros = brent_N_roots(f, 0, 10, pi/2, 1e-15, 1);

  for (int i=0; i<zeros.size(); i++)
    cout << i << " : " << setprecision(15) << zeros[i] << endl;
  
  cout << "Iterations : " << f.times_called() << endl;
  
  return 0;
}
