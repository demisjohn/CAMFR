
/////////////////////////////////////////////////////////////////////////////
//
// File:     function.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000125
// Version:  1.2
//
// Copyright (C) 1998-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef FUNCTION_H
#define FUNCTION_H

#include <vector>
#include "../../defs.h"

/////////////////////////////////////////////////////////////////////////////
//
// Template for univariate function.
//
/////////////////////////////////////////////////////////////////////////////

template<class T> class Function1D
{
  public:

    Function1D<T>()                        {counter=0;}
    Function1D<T>(const std::vector<T>& p) {counter=0; set_params(p);}
    virtual ~Function1D<T>()          {}
    
    virtual T operator()(const T& t) = 0; // should contain 'counter++'
    int times_called() {return counter;}

    virtual std::vector<T> get_params() const {std::vector<T> t; return t;}
    virtual void           set_params  (const  std::vector<T>&) {}
  
  protected:

    int counter;
};



/////////////////////////////////////////////////////////////////////////////
//
// Some typedefs to hide the templates from SWIG.
//
/////////////////////////////////////////////////////////////////////////////

typedef Function1D<Real>    RealFunction;
typedef Function1D<Complex> ComplexFunction;



/////////////////////////////////////////////////////////////////////////////
//
// Wrapper routines for partially converting complex functions to real
// functions:
//
//   Wrap_real_to_real : real(f(  x))
//   Wrap_real_to_imag : imag(f(  x))
//   Wrap_imag_to_real : real(f(i.x))
//   Wrap_imag_to_imag : imag(f(i.x))
//   Wrap_real_to_abs  :  abs(f(  x))
//   Wrap_imag_to_abs  :  abs(f(i.x))
//   Wrap_real_to_arg  :  abs(f(  x))
//   Wrap_imag_to_arg  :  abs(f(i.x))
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// Wrap_real_to_real
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_real_to_real : public RealFunction
{
  public:

    Wrap_real_to_real(ComplexFunction& general_)
      : general(&general_) {}

    Real operator()(const Real& x) {return real((*general)(Complex(x,0.0)));}

  protected:

    ComplexFunction* general;
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrap_real_to_imag
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_real_to_imag : public RealFunction
{ 
  public:

    Wrap_real_to_imag(ComplexFunction& general_)
      : general(&general_) {}
   
    Real operator()(const Real& x) {return imag((*general)(Complex(x,0.0)));}

  protected:

    ComplexFunction* general;       
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrap_imag_to_real
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_imag_to_real : public RealFunction
{
  public:

    Wrap_imag_to_real(ComplexFunction& general_)
      : general(&general_) {}
   
    Real operator()(const Real& x) {return real((*general)(Complex(0.0,x)));}

  protected:

    ComplexFunction* general;      
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrap_imag_to_imag
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_imag_to_imag : public RealFunction
{
  public:

    Wrap_imag_to_imag(ComplexFunction& general_)
      : general(&general_) {}
   
    Real operator()(const Real& x) {return imag((*general)(Complex(0.0,x)));}

  protected:

    ComplexFunction* general;   
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrap_real_to_abs
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_real_to_abs : public RealFunction
{ 
  public:

    Wrap_real_to_abs(ComplexFunction& general_)
      : general(&general_) {}
   
    Real operator()(const Real& x) {return abs((*general)(Complex(x,0.0)));}

  protected:

    ComplexFunction* general;       
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrap_imag_to_abs
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_imag_to_abs : public RealFunction
{
  public:

    Wrap_imag_to_abs(ComplexFunction& general_)
      : general(&general_) {}
   
    Real operator()(const Real& x) {return abs((*general)(Complex(0.0,x)));}

  protected:

    ComplexFunction* general;      
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrap_real_to_arg
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_real_to_arg : public RealFunction
{ 
  public:

    Wrap_real_to_arg(ComplexFunction& general_)
      : general(&general_) {}
   
    Real operator()(const Real& x) 
      {return std::arg((*general)(Complex(x,0.0)));}

  protected:

    ComplexFunction* general;       
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrap_imag_to_arg
//
/////////////////////////////////////////////////////////////////////////////

class Wrap_imag_to_arg : public RealFunction
{
  public:

    Wrap_imag_to_arg(ComplexFunction& general_)
      : general(&general_) {}
   
    Real operator()(const Real& x) 
      {return std::arg((*general)(Complex(0.0,x)));}

  protected:

    ComplexFunction* general;      
};



#endif
