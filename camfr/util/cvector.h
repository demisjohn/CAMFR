
/////////////////////////////////////////////////////////////////////////////
//
// File:     cvector.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990713
// Version:  1.02
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CVECTOR_H
#define CVECTOR_H

#include <vector>
#include "../defs.h"

/////////////////////////////////////////////////////////////////////////////
//
// Defines a few arithmetic operations on vector<Complex>.
//
//  Note: the same functionality is also present in Blitz, but we use
//  the versions here to speed up compilation times. 
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> operator+(const vector<Complex>& a, const vector<Complex>& b);
vector<Complex>& operator+=(vector<Complex>& a, const vector<Complex>& b);

vector<Complex> operator-(const vector<Complex>& a, const vector<Complex>& b);
vector<Complex>& operator-=(vector<Complex>& a, const vector<Complex>& b);

vector<Complex> operator*(const Complex& c, const vector<Complex>& v);
vector<Complex> operator*(const vector<Complex>& v, const Complex& c);
vector<Complex>& operator*=(vector<Complex>& v, const Complex& c);

vector<Complex> operator/(const vector<Complex>& v, const Complex& c);
vector<Complex>& operator/=(vector<Complex>& v, const Complex& c);

vector<Complex> operator-(const vector<Complex>& a);
  
Real abs(const vector<Complex>& v);



#endif
