
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

std::vector<Complex> operator+(const std::vector<Complex>& a, 
                               const std::vector<Complex>& b);

std::vector<Complex>& operator+=(std::vector<Complex>& a, 
                                 const std::vector<Complex>& b);

std::vector<Complex> operator-(const std::vector<Complex>& a, 
                               const std::vector<Complex>& b);

std::vector<Complex>& operator-=(std::vector<Complex>& a, 
                                 const std::vector<Complex>& b);

std::vector<Complex> operator*(const Complex& c, 
                               const std::vector<Complex>& v);

std::vector<Complex> operator*(const std::vector<Complex>& v,const Complex& c);
std::vector<Complex>& operator*=(std::vector<Complex>& v, const Complex& c);

std::vector<Complex> operator/(const std::vector<Complex>& v,const Complex& c);
std::vector<Complex>& operator/=(std::vector<Complex>& v, const Complex& c);

std::vector<Complex> operator-(const std::vector<Complex>& a);
  
Real abs(const std::vector<Complex>& v);



#endif
