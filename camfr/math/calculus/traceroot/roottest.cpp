
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
  // Find roots of f1.
  
  F f1(0);
  Wrap_real_to_real f1_wrap(f1);
  vector<Real> zeros1 = brent_N_roots(f1_wrap, 0, 10, pi/2, 1e-15, 1);
  
  cout << "Zeros of f1:" << endl;
  for (int i=0; i<zeros1.size(); i++)
    cout << i << " : " << zeros1[i] << endl;

  // Trace them to those of f2.

  F f2(10.0*I);

  vector<Complex> p1 = f1.get_params();
  vector<Complex> p2 = f2.get_params();
  vector<Complex> forbidden;
  
  vector<Complex> zeros2 = traceroot(zeros1, f1, p1, p2, forbidden);
  
  cout << endl << "Zeros of f2:" << endl;
  for (int i=0; i<zeros2.size(); i++)
    cout << i << " : " << zeros2[i] << endl;
  
  return 0; 
}
