
/////////////////////////////////////////////////////////////////////////////
//
// File:          roottest.cpp
// Author:        Peter.Bienstman@rug.ac.be
// Date:          19990712
// Version:       1.0
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "traceroot.h"
#include "../calculus.h"
#include "../../bessel/bessel.h"
#include "../../../util/cvector.h"

/////////////////////////////////////////////////////////////////////////////
//
// Simple function object for roottester : Bessel function + offset
//
/////////////////////////////////////////////////////////////////////////////

class F : public Function1D<Complex>
{
  public:

    F(const Complex& a) : alpha(a) {}
    
    Complex operator()(const Complex& x) 
    {
      counter++;
      return J(1,x-alpha);
    }

    vector<Complex> export_params() const
    {
      vector<Complex> param;
      param.push_back(alpha);
      return param;
    }
    
    void import_params(const vector<Complex>& param)
    {
      alpha = param[0];
    }

    Complex alpha;
    
};



/////////////////////////////////////////////////////////////////////////////
//
// roottester
//
/////////////////////////////////////////////////////////////////////////////

int main()
{

  // find roots of f1
  
  F f1(0);
  wrap_real_to_real f1_wrap(f1);
  vector<Real> zeros1 = brent_N_roots(f1_wrap, 0, 10, pi/2, 1e-15, 1);
  
  cout << "Zeros of f1:" << endl;
  for (int i=0; i<zeros1.size(); i++)
    cout << i << " : " << zeros1[i] << endl;

  // trace them to those of f2

  F f2(10*I);
  
  vector<Complex> zeros2 = traceroot(f1, zeros1, f2);
  
  cout << endl << "Zeros of f2:" << endl;
  for (int i=0; i<zeros2.size(); i++)
    cout << i << " : " << zeros2[i] << endl;
  
  return 0;
  
}
