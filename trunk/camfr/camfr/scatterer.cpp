
/////////////////////////////////////////////////////////////////////////////
//
// File:     scatterer.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000511
// Version:  1.23
//
// Copyright (C) 1999-2000 Peter Bienstman - Ghent University
//
////////////////////////////////////////////////////////////////////////////

#include "scatterer.h"

/////////////////////////////////////////////////////////////////////////////
//
// MultiScatterer::recalc_needed
//  
/////////////////////////////////////////////////////////////////////////////

bool MultiScatterer::recalc_needed() const
{  
  const Real eps = 1e-10;

  if (abs(global.lambda - last_lambda) > eps)
    return true;

  if (abs(global.slab_ky - last_slab_ky) > eps)
    return true;

  if (!global.gain_mat)
    return false;

  if (!contains(*global.gain_mat))
    return false;

  if (abs(global.gain_mat-> n()  - last_gain_mat. n())  > eps)
    return true;
  
  if (abs(global.gain_mat->mur() - last_gain_mat.mur()) > eps)
    return true;

  return false;
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseScatterer::DenseScatterer
//  
/////////////////////////////////////////////////////////////////////////////

DenseScatterer::DenseScatterer()
  : MultiScatterer(),
    R12(fortranArray), R21(fortranArray),
    T12(fortranArray), T21(fortranArray) {}



/////////////////////////////////////////////////////////////////////////////
//
// DenseScatterer::DenseScatterer
//  
/////////////////////////////////////////////////////////////////////////////

DenseScatterer::DenseScatterer(Waveguide& inc, Waveguide& ext)
  : MultiScatterer(inc, ext),
    R12(fortranArray), R21(fortranArray),
    T12(fortranArray), T21(fortranArray) {}



/////////////////////////////////////////////////////////////////////////////
//
// DenseScatterer::allocRT
//  
/////////////////////////////////////////////////////////////////////////////

void DenseScatterer::allocRT()
{
  const int N = global.N;

  R12.resize(N,N);
  R21.resize(N,N);
  T12.resize(N,N);
  T21.resize(N,N);

  // We haven't calculated the matrices yet for any wavelength or gain.

  last_lambda   = 0.0;
  last_gain_mat = Material(0.0);
  last_slab_ky  = 0.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseScatterer::freeRT
//  
/////////////////////////////////////////////////////////////////////////////

void DenseScatterer::freeRT()
{ 
  R12.free();
  R21.free();
  T12.free();
  T21.free();

  last_lambda   = 0.0;
  last_gain_mat = Material(0.0);
  last_slab_ky  = 0.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseScatterer::copy_RT_from
//  
/////////////////////////////////////////////////////////////////////////////

void DenseScatterer::copy_RT_from(const DenseScatterer& sc)
{
  allocRT();

  R12 = sc.get_R12();
  R21 = sc.get_R21();
  T12 = sc.get_T12();
  T21 = sc.get_T21();

  // Remember wavelength and gain these matrices were calculated for.

  last_lambda = sc.last_lambda;
  if (global.gain_mat) 
    last_gain_mat = sc.last_gain_mat;
  last_slab_ky = sc.last_slab_ky;
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseScatterer::swap_RT_with
//  
/////////////////////////////////////////////////////////////////////////////

void DenseScatterer::swap_RT_with(DenseScatterer& sc)
{
  blitz::cycleArrays(R12, sc.R12);
  blitz::cycleArrays(R21, sc.R21);
  blitz::cycleArrays(T12, sc.T12);
  blitz::cycleArrays(T21, sc.T21);
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::DiagScatterer
//  
/////////////////////////////////////////////////////////////////////////////

DiagScatterer::DiagScatterer()
  : MultiScatterer(),
    R12(fortranArray), R21(fortranArray),
    T12(fortranArray), T21(fortranArray),
    R12_dense(NULL),   R21_dense(NULL),
    T12_dense(NULL),   T21_dense(NULL) {}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::DiagScatterer
//  
/////////////////////////////////////////////////////////////////////////////

DiagScatterer::DiagScatterer(Waveguide& inc, Waveguide& ext)
  : MultiScatterer(inc, ext),
    R12(fortranArray), R21(fortranArray),
    T12(fortranArray), T21(fortranArray),
    R12_dense(NULL),   R21_dense(NULL),
    T12_dense(NULL),   T21_dense(NULL)
{  
  if ( (inc != ext) && ( !inc.is_uniform() || !ext.is_uniform() ) )
  {
    py_error("Error: inc and exit media should be uniform in DiagScatterer.");
    exit (-1);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::DiagScatterer
//  
/////////////////////////////////////////////////////////////////////////////

DiagScatterer::DiagScatterer(const DiagScatterer& sc_d)
  : MultiScatterer(sc_d)
{
  copy_RT_from(sc_d);
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::operator=
//  
/////////////////////////////////////////////////////////////////////////////

DiagScatterer& DiagScatterer::operator=(const DiagScatterer& sc_d)
{
  if (this == &sc_d)
    return *this;

  freeRT();

  Scatterer::operator=(sc_d);

  copy_RT_from(sc_d);

  return *this;
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::allocRT
//  
/////////////////////////////////////////////////////////////////////////////

void DiagScatterer::allocRT()
{
  const int N = global.N;

  R12.resize(N);
  R21.resize(N);
  T12.resize(N);
  T21.resize(N);

  delete R12_dense; R12_dense = NULL;
  delete R21_dense; R21_dense = NULL;
  delete T12_dense; T12_dense = NULL;
  delete T21_dense; T21_dense = NULL;

  // We haven't calculated the matrices yet for any wavelength or gain.

  last_lambda   = 0.0;
  last_gain_mat = Material(0.0);
  last_slab_ky  = 0.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::freeRT
//  
/////////////////////////////////////////////////////////////////////////////

void DiagScatterer::freeRT()
{
  R12.free();
  R21.free();
  T12.free();
  T21.free();

  delete R12_dense; R12_dense = NULL;
  delete R21_dense; R21_dense = NULL;
  delete T12_dense; T12_dense = NULL;
  delete T21_dense; T21_dense = NULL;

  // We haven't recalculated the matrices yet for any wavelength or gain.

  last_lambda   = 0.0;
  last_gain_mat = Material(0.0);
  last_slab_ky  = 0.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::copy_RT_from
//  
/////////////////////////////////////////////////////////////////////////////

void DiagScatterer::copy_RT_from(const DiagScatterer& sc_d)
{
  allocRT();

  R12 = sc_d.get_diag_R12();
  R21 = sc_d.get_diag_R21();
  T12 = sc_d.get_diag_T12();
  T21 = sc_d.get_diag_T21();

  // Don't copy dense matrices (waste of space).

  R12_dense = NULL;
  R21_dense = NULL;
  T12_dense = NULL;
  T21_dense = NULL;

  // Remember wavelength and gain these matrices were calculated for.

  last_lambda = sc_d.last_lambda;
  if (global.gain_mat) 
    last_gain_mat = sc_d.last_gain_mat;
  last_slab_ky = sc_d.last_slab_ky;
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::swap_RT_with
//  
/////////////////////////////////////////////////////////////////////////////

void DiagScatterer::swap_RT_with(DiagScatterer& sc_d)
{
  blitz::cycleArrays(R12, sc_d.R12);
  blitz::cycleArrays(R21, sc_d.R21);
  blitz::cycleArrays(T12, sc_d.T12);
  blitz::cycleArrays(T21, sc_d.T21);
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagScatterer::convert_to_dense
//  
/////////////////////////////////////////////////////////////////////////////

void DiagScatterer::convert_to_dense() const
{
  if (!recalc_needed() && R12_dense)
    return;

  const_cast<DiagScatterer*>(this)->calcRT();
  
  const int N = global.N;

  R12_dense = new cMatrix(N,N,fortranArray); *R12_dense = 0.0;
  R21_dense = new cMatrix(N,N,fortranArray); *R21_dense = 0.0;
  T12_dense = new cMatrix(N,N,fortranArray); *T12_dense = 0.0;
  T21_dense = new cMatrix(N,N,fortranArray); *T21_dense = 0.0;
  
  for (int i=1; i<=N; i++)
  {
    (*R12_dense)(i,i) = R12(i);
    (*R21_dense)(i,i) = R21(i);
    (*T12_dense)(i,i) = T12(i);
    (*T21_dense)(i,i) = T21(i);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoScatterer::MonoScatterer
//  
/////////////////////////////////////////////////////////////////////////////

MonoScatterer::MonoScatterer(Waveguide& inc, Waveguide& ext)
  : Scatterer(inc, ext)
{
  if (    !dynamic_cast<MonoWaveguide*>(&inc)
       || !dynamic_cast<MonoWaveguide*>(&ext) )
  {
    py_error("Error: only MonoWaveguides allowed in MonoScatterer.");
    exit (-1);
  }        
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoScatterer::copy_RT_from
//  
/////////////////////////////////////////////////////////////////////////////

void MonoScatterer::copy_RT_from(const MonoScatterer& sc_m)
{
  R12 = sc_m.get_R12();
  R21 = sc_m.get_R21();
  T12 = sc_m.get_T12();
  T21 = sc_m.get_T21();
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoScatterer::swap_RT_with
//  
/////////////////////////////////////////////////////////////////////////////

void MonoScatterer::swap_RT_with(MonoScatterer& sc_m)
{
  Complex tmp;

  tmp=R12; R12=sc_m.R12; sc_m.R12=tmp;
  tmp=R21; R21=sc_m.R21; sc_m.R21=tmp;
  tmp=T12; T12=sc_m.T12; sc_m.T12=tmp;
  tmp=T21; T21=sc_m.T21; sc_m.T21=tmp;
}



/////////////////////////////////////////////////////////////////////////////
//
// TransparentScatterer::calcRT
//  
/////////////////////////////////////////////////////////////////////////////

void TransparentScatterer::calcRT()
{
  inc->find_modes();
  
  allocRT(); 

  R12 = 0.0; R21 = 0.0; 
  T12 = 1.0; T21 = 1.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// E_Wall::calcRT
//
//  Note: T is set to unity to avoid spurious warnings in field_calc.
//  
/////////////////////////////////////////////////////////////////////////////

void E_Wall::calcRT()
{
  inc->find_modes();

  allocRT(); 

  R12 = -1.0; R21 = -1.0; 
  T12 =  0.0; T21 =  1.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// H_Wall::calcRT
//  
/////////////////////////////////////////////////////////////////////////////

void H_Wall::calcRT()
{
  inc->find_modes();

  allocRT(); 

  R12 = 1.0; R21 = 1.0; 
  T12 = 0.0; T21 = 1.0;
}

