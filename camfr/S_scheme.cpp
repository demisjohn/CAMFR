
/////////////////////////////////////////////////////////////////////////////
//
// File:     S_scheme.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000612
// Version:  1.3
//
// Copyright (C) 1998-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "S_scheme.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// calc_tilde
//
//   Calculates 'tilde' matrices used in S_scheme.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_tilde(const Chunk& chunk,
                cMatrix* r12, cMatrix* r21, cMatrix* t12, cMatrix* t21)
{
  // No propagation needed?

  MultiScatterer* s = dynamic_cast<MultiScatterer*>(chunk.sc);
  
  if (abs(chunk.d) == 0)
  {
    *r12 = s->get_R12(); *r21 = s->get_R21();
    *t12 = s->get_T12(); *t21 = s->get_T21();

    return;
  }

  // Create vector with propagation factors.

  cVector prop(global.N,fortranArray);
  for (int i=1; i<=global.N; i++)
    prop(i) = exp( -I * s->get_ext()->get_mode(i)->get_kz() * chunk.d );
  
  // Transparent scatterer?

  if (dynamic_cast<TransparentScatterer*>(s))
  {
    *r12 = 0; *r21 = 0; *t12 = 0; *t21 = 0;
    
    for (int i=1; i<=global.N; i++)
      (*t12)(i,i) = (*t21)(i,i) = prop(i);
  }
  else
  {
    blitz::firstIndex i; blitz::secondIndex j;

    *r12 =           (s->get_R12());
    *r21 = prop(i) * (s->get_R21())(i,j) * prop(j);
    *t12 = prop(i) * (s->get_T12())(i,j);
    *t21 =           (s->get_T21())(i,j) * prop(j);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// S_scheme
//  
/////////////////////////////////////////////////////////////////////////////

void S_scheme(const vector<Chunk>& chunks, DenseScatterer* result)
{  
  result->allocRT();

  // Matrices to store result from previous iteration.

  DenseScatterer prev; prev.allocRT();

  const cMatrix& pR12(prev.get_R12()); const cMatrix& pR21(prev.get_R21());
  const cMatrix& pT12(prev.get_T12()); const cMatrix& pT21(prev.get_T21());

  // 'Tilde' matrices.

  const int N = global.N;

  cMatrix r12(N,N,fortranArray); cMatrix r21(N,N,fortranArray);
  cMatrix t12(N,N,fortranArray); cMatrix t21(N,N,fortranArray);
  
  // Unit and auxiliary matrices.
  
  cMatrix U1(N,N,fortranArray), M(N,N,fortranArray), tmp(N,N,fortranArray);

  U1 = 0.0;
  for (int i=1; i<=N; i++)
    U1(i,i) = Complex(1.0, 0.0);
  
  // Initialise matrices from first chunk.
  // Make sure 'result' contains separate data from other matrices,
  // since it will be overwritten.
  
  calc_tilde(chunks[0], &r12, &r21, &t12, &t21);

  result->copy_R12(r12); result->copy_R21(r21);
  result->copy_T12(t12); result->copy_T21(t21);
  
  // Loop from second chunk to last one.
  
  for (unsigned int k=1; k<chunks.size(); k++)
  {
    // Previous result = current result (fast shallow copy).

    prev.swap_RT_with(*result);
    
    // Calculate new result matrices.

    calc_tilde(chunks[k], &r12, &r21, &t12, &t21);
    
    tmp = U1 - multiply(r12, pR21);
    if (global.stability != SVD)
      M.reference(invert    (tmp));
    else
      M.reference(invert_svd(tmp));

    tmp = multiply(pT21, M, r12, pT12) + pR12;
    result->copy_R12(tmp);
    result-> set_T21(multiply(pT21, M, t21));
    
    tmp = U1 - multiply(pR21, r12);
    if(global.stability != SVD)
      M.reference(invert    (tmp));
    else
      M.reference(invert_svd(tmp));

    tmp = multiply(t12, M, pR21, t21) + r21;

    result->copy_R21(tmp);
    result-> set_T12(multiply(t12, M, pT12));
  } 
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_tilde
//
//   Calculates 'tilde' vectors used in S_scheme.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_tilde(const Chunk& chunk,
                cVector* r12, cVector* r21, cVector* t12, cVector* t21)
{
  // No propagation needed?

  DiagScatterer* s = dynamic_cast<DiagScatterer*>(chunk.sc);
  
  if (abs(chunk.d) == 0)
  {
    *r12 = s->get_diag_R12(); *r21 = s->get_diag_R21();
    *t12 = s->get_diag_T12(); *t21 = s->get_diag_T21();

    return;
  }

  // Create vector with propagation constants.

  cVector prop(global.N,fortranArray);
  for (int i=1; i<=global.N; i++)
    prop(i) = exp( -I * s->get_ext()->get_mode(i)->get_kz() * chunk.d );

  // Transparent scatterer?

  if (dynamic_cast<TransparentScatterer*>(s))
  {
    *r12 = 0; *r21 = 0; *t12 = 0; *t21 = 0;
    
    for (int i=1; i<=global.N; i++)
      (*t12)(i) = (*t21)(i) = prop(i);
  }
  else
  {
    *r12 =        s->get_diag_R12();
    *r21 = prop * s->get_diag_R21() * prop;
    *t12 = prop * s->get_diag_T12();
    *t21 =        s->get_diag_T21() * prop;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// S_scheme for diagonal Scatterers.
//  
/////////////////////////////////////////////////////////////////////////////

void S_scheme(const vector<Chunk>& chunks, DiagScatterer* result)
{ 
  result->allocRT();
  
  // Vectors to store result from previous iteration.

  DiagScatterer prev; prev.allocRT();

  const cVector& pR12(prev.get_diag_R12());
  const cVector& pR21(prev.get_diag_R21());
  const cVector& pT12(prev.get_diag_T12());
  const cVector& pT21(prev.get_diag_T21());
  
  // 'Tilde' vectors.

  const int N = global.N;
    
  cVector r12(N,fortranArray), r21(N,fortranArray);
  cVector t21(N,fortranArray), t12(N,fortranArray);
    
  // Auxiliary and temporary vectors.
  
  cVector M(N,fortranArray), tmp(N,fortranArray);
  
  // Initialise vectors from first chunk.
  // Make sure 'result' contains separate data from other matrices,
  // since it will be overwritten.
  
  calc_tilde(chunks[0], &r12, &r21, &t12, &t21);
  
  result->copy_diag_R12(r12); result->copy_diag_R21(r21);
  result->copy_diag_T12(t12); result->copy_diag_T21(t21);
  
  // Loop from second chunk to last one.
  
  for (unsigned int k=1; k<chunks.size(); k++)
  { 
    // Previous result = current result (fast shallow copy).

    prev.swap_RT_with(*result);
    
    // Calculate new result vectors.

    calc_tilde(chunks[k], &r12, &r21, &t12, &t21);
    
    M = 1.0/(1.0-r12*pR21);

    tmp = pT21 * M *  r12 * pT12  +  pR12; result->copy_diag_R12(tmp);
    tmp = pT21 * M *  t21;                 result->copy_diag_T21(tmp);
    tmp =  t12 * M * pR21 *  t21  +   r21; result->copy_diag_R21(tmp);
    tmp =  t12 * M * pT12;                 result->copy_diag_T12(tmp);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_tilde
//
//   Calculates 'tilde' variables used in S_scheme.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_tilde(const Chunk& chunk,
                Complex* r12, Complex* r21, Complex* t12, Complex* t21)
{ 
  // No propagation needed?

  MonoScatterer* s = dynamic_cast<MonoScatterer*>(chunk.sc);
  
  if (abs(chunk.d) == 0)
  {
    *r12 = s->get_R12(); *r21 = s->get_R21();
    *t12 = s->get_T12(); *t21 = s->get_T21();
    
    return;
  }

  // Propagation needed.
  
  Complex kz = s->get_ext()->get_mode(1)->get_kz();

  Complex M;
  if (imag(kz) > 1000)
  {
    py_print("Warning: S-scheme not suited for extremely high gain.");
    M = 1e6;
  }
  else
    M = exp(-I*kz*chunk.d);

  *r12 =     s->get_R12();
  *r21 = M * s->get_R21() * M;
  *t12 = M * s->get_T12();
  *t21 =     s->get_T21() * M;
}



/////////////////////////////////////////////////////////////////////////////
//
// S_scheme for Monoscatterers
//  
/////////////////////////////////////////////////////////////////////////////

void S_scheme(const vector<Chunk>& chunks, MonoScatterer* result)
{
  // Variables to store result from previous iteration.
  
  Complex pR12, pR21, pT12, pT21;
  
  // 'Tilde' variables.

  Complex r12, r21, t21, t12;
  
  // Initialise variables from first chunk.

  calc_tilde(chunks[0], &r12, &r21, &t12, &t21);
  
  result->set_R12(r12); result->set_R21(r21);
  result->set_T12(t12); result->set_T21(t21);
  
  // Loop from second chunk to last one.
  
  for (unsigned int k=1; k<chunks.size(); k++)
  {
    // Previous result = current result.

    pR12 = result->get_R12(); pR21 = result->get_R21();
    pT12 = result->get_T12(); pT21 = result->get_T21();
    
    // Calculate new result variables.

    calc_tilde(chunks[k], &r12, &r21, &t12, &t21);
    
    Complex res = 1.0-r12*pR21;
    Complex M = 1.0/res;

    if (abs(res) < 1e-12) // Regularise in case of resonance.
      M = 0.0;

    result->set_R12(pT21 * M *  r12 * pT12  +  pR12);
    result->set_T21(pT21 * M *  t21);
    result->set_R21( t12 * M * pR21 *  t21  +   r21);
    result->set_T12( t12 * M * pT12);
  }
}
