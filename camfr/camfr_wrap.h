
/////////////////////////////////////////////////////////////////////////////
//
// File:     camfr_wrap.h
// Author:   Peter.Bienstman@UGent.be
//
// Copyright (C) 1998-2006 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CAMFR_WRAP_H
#define CAMFR_WRAP_H

#include "defs.h"
#include "waveguide.h"
#include "math/calculus/function.h"

/////////////////////////////////////////////////////////////////////////////
//
// check_index
//
/////////////////////////////////////////////////////////////////////////////

inline void check_index(int i)
{
  if ( (i<0) || (i>=int(global.N)) )
  {
    PyErr_SetString(PyExc_IndexError, "index out of bounds.");
    throw std::out_of_range("index out of bounds.");
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// check_wg_index
//
/////////////////////////////////////////////////////////////////////////////

inline void check_wg_index(const Waveguide& w, int i)
{
  if ( (i<0) || (i>=int(w.N())) )
  {
    PyErr_SetString(PyExc_IndexError, "index out of bounds.");
    throw std::out_of_range("index out of bounds.");
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// The following classes are used when expanding an abritrarily shaped 
// field in slabmodes.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// PythonFunction
//
/////////////////////////////////////////////////////////////////////////////

class PythonFunction : public ComplexFunction
{
  public:

    PythonFunction(PyObject* f_): f(f_) {}

    Complex operator()(const Complex& z)
      {counter++; return boost::python::call<Complex>(f, z);}

  protected:

    PyObject* f;
};



/////////////////////////////////////////////////////////////////////////////
//
// GaussianFunction
//
/////////////////////////////////////////////////////////////////////////////

class GaussianFunction : public ComplexFunction
{
  public:

    GaussianFunction (Complex height, Complex width, Complex position)
      : h(height), w(width), p(position) {}

    Complex operator()(const Complex& x)
	{counter++; return h*exp(-(x-p)*(x-p)/(w*w*2.0));}

  protected:

    Complex h, w, p;
};



/////////////////////////////////////////////////////////////////////////////
//
// PlaneWaveFunction
//
/////////////////////////////////////////////////////////////////////////////

class PlaneWaveFunction : public ComplexFunction
{
  public:

    PlaneWaveFunction (Complex amplitude, Complex angle, Complex index)
      : am(amplitude), an(angle), n(index) {}

    Complex operator()(const Complex& x)
      {counter++; return am*exp(-2.0*I*pi*sin(an)*n*x/global.lambda);}

  protected:

    Complex am, an, n;
};



#endif


