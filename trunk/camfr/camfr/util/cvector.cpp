
/////////////////////////////////////////////////////////////////////////////
//
// File:     cvector.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990713
// Version:  1.02
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "cvector.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// Addition.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> operator+(const vector<Complex>& a,
                          const vector<Complex>& b)
{  
  vector<Complex> result;

  for (unsigned int i=0; i<a.size(); i++)
    result.push_back(a[i]+b[i]);

  return result; 
}

vector<Complex>& operator+=(      vector<Complex>& a,
                            const vector<Complex>& b)
{
  a = a+b;
  return a;
}



/////////////////////////////////////////////////////////////////////////////
//
// Substraction.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> operator-(const vector<Complex>& a,
                          const vector<Complex>& b)
{  
  vector<Complex> result;

  for (unsigned int i=0; i<a.size(); i++)
    result.push_back(a[i]-b[i]);

  return result; 
}

vector<Complex>& operator-=(      vector<Complex>& a,
                            const vector<Complex>& b)
{
  a = a-b;
  return a;
}



/////////////////////////////////////////////////////////////////////////////
//
// Multiplication.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> operator*(const Complex& c, const vector<Complex>& v)
{ 
  vector<Complex> result;

  for (unsigned int i=0; i<v.size(); i++)
    result.push_back(c*v[i]);

  return result; 
}

vector<Complex> operator*(const vector<Complex>& v, const Complex& c)
{
  return c*v;
}

vector<Complex>& operator*=(vector<Complex>& v, const Complex& c)
{
  v = c*v;
  return v;
}



/////////////////////////////////////////////////////////////////////////////
//
// Division.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> operator/(const vector<Complex>& v, const Complex& c)
{ 
  vector<Complex> result;

  for (unsigned int i=0; i<v.size(); i++)
    result.push_back(v[i]/c);

  return result;
}

vector<Complex>& operator/=(vector<Complex>& v, const Complex& c)
{
  v = v/c;
  return v;
}


/////////////////////////////////////////////////////////////////////////////
//
// Abs.
//
/////////////////////////////////////////////////////////////////////////////

Real abs(const vector<Complex>& v)
{
  
  Complex result = 0;

  for (unsigned int i=0; i<v.size(); i++)
    result += v[i]*v[i];

  return abs(sqrt(result));
  
}



/////////////////////////////////////////////////////////////////////////////
//
// Sign change.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> operator-(const vector<Complex>& a)
{  
  vector<Complex> result;

  for (unsigned int i=0; i<a.size(); i++)
    result.push_back(-a[i]);

  return result; 
}
