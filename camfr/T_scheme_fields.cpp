
/////////////////////////////////////////////////////////////////////////////
//
// File:     T_scheme_fields.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010702
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "T_scheme_fields.h"

/////////////////////////////////////////////////////////////////////////////
//
// cross_interface
//
//   Given forward and backward fields 'fw' and 'bw' at the beginning of
//   an interface scatterer, calculate the field at the other side of the
//   scatterer. Results are calculated in place in 'fw' and 'bw'.
//  
/////////////////////////////////////////////////////////////////////////////

void cross_interface(MultiScatterer* sc, cVector* fw, cVector* bw)
{
  // Define matrices.

  const cMatrix& R12(sc->get_R12());
  const cMatrix& R21(sc->get_R21());
  const cMatrix& T12(sc->get_T12());
  const cMatrix& T21(sc->get_T21());

  const int N = global.N;

  cVector fw0(N,fortranArray); fw0 = *fw;
  cVector bw0(N,fortranArray); bw0 = *bw;
  
  cMatrix ff(N,N,fortranArray); cMatrix fb(N,N,fortranArray);
  cMatrix bf(N,N,fortranArray); cMatrix bb(N,N,fortranArray);

  cMatrix  M(N,N,fortranArray); cMatrix tmp(N,N,fortranArray);

  // Calculate field.

  if (global.stability == extra)
    M.reference(invert_x(T21));
  else
    M.reference(invert  (T21));
  
  ff = T12 - multiply(R21,M,R12);
  fb.reference(multiply(R21,M));

  tmp = -M;
  bf.reference(multiply(tmp, R12));
  bb.reference(M);
  
  *fw = multiply(ff, fw0) + multiply(fb, bw0);
  *bw = multiply(bf, fw0) + multiply(bb, bw0);
}



/////////////////////////////////////////////////////////////////////////////
//
// cross_interface for diagonal scatterers
//  
/////////////////////////////////////////////////////////////////////////////

void cross_interface(DiagScatterer* sc, cVector* fw, cVector* bw)
{
  // Define vectors.

  const int N = global.N;

  cVector fw0(N,fortranArray); fw0 = *fw;
  cVector bw0(N,fortranArray); bw0 = *bw;

  const cVector& R12(sc->get_diag_R12());
  const cVector& R21(sc->get_diag_R21());
  const cVector& T12(sc->get_diag_T12());
  const cVector& T21(sc->get_diag_T21());

  // Calculate field.

  *fw = ( (T12*T21 - R12*R21) * fw0  +  R21 * bw0 ) / T21;
  *bw = (               -R12  * fw0  +        bw0 ) / T21;
}



/////////////////////////////////////////////////////////////////////////////
//
// cross_interface for monoscatterers.
//  
/////////////////////////////////////////////////////////////////////////////

void cross_interface(MonoScatterer* sc, cVector* fw, cVector* bw)
{
  // Define vectors.

  cVector fw0(1,fortranArray); fw0 = *fw;
  cVector bw0(1,fortranArray); bw0 = *bw;

  const Complex& R12(sc->get_R12());
  const Complex& R21(sc->get_R21());
  const Complex& T12(sc->get_T12());
  const Complex& T21(sc->get_T21());
  
  // Calculate field.

  *fw = ( (T12*T21 - R12*R21) * fw0  +  R21 * bw0 ) / T21;
  *bw = (               -R12  * fw0  +        bw0 ) / T21;
}



/////////////////////////////////////////////////////////////////////////////
//
// T_calc_fields
//
//   Calculates the fields in each chunks at two positions, before and
//   after the propagation.
//  
/////////////////////////////////////////////////////////////////////////////

template<class T>
void T_calc_fields
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field)
{
  const int N = global.N;
  
  cVector fw_scaled(N, fortranArray); fw_scaled = (*field)[0].fw;
  cVector bw_scaled(N, fortranArray); bw_scaled = (*field)[0].bw;
  
  // Loop over chunks.
  
  Complex corr = 0.0;
  for (unsigned int k=0; k<chunks.size(); k++)
  {
    // Cross interface.
    
    Waveguide* wg = chunks[k].sc->get_ext();
        
    cross_interface(dynamic_cast<T>(chunks[k].sc), &fw_scaled, &bw_scaled);

    field->push_back(exp(corr) * FieldExpansion(wg, fw_scaled, bw_scaled));

    // Find exponent to be factored out in propagation.

    Complex kz_f = wg->get_mode(1)->kz;
    for (unsigned int i=2; i<=N; i++)
    {
      const Complex kz = wg->get_mode(i)->kz;
      
      if (imag(kz) < imag(kz_f))
        kz_f = kz;
    }

    // Propagate and rescale.
    
    corr += I * kz_f;

    for (int i=1; i<=N; i++)
    {
      const Complex kz = wg->get_mode(i)->kz;
    
      fw_scaled(i) *= exp(-I * (kz+kz_f) * chunks[k].d);
      bw_scaled(i) *= exp( I * (kz-kz_f) * chunks[k].d);
    }

    field->push_back(exp(corr) * FieldExpansion(wg, fw_scaled, bw_scaled));
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// T_scheme_fields
//  
/////////////////////////////////////////////////////////////////////////////

void T_scheme_fields
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field)
{
  // Check input.

  if (field->size() != 1)
  {
    cerr << "Error: invalid source fields in T_scheme_fields." << endl;
    exit (-1);
  }
  
  // Check type of chunks and call corresponding function.

  for (unsigned int i=0; i<chunks.size(); i++)
    if (dynamic_cast<DenseScatterer*>(chunks[i].sc))
      return T_calc_fields<MultiScatterer*>(chunks, field);

  for (unsigned int i=0; i<chunks.size(); i++)
    if (dynamic_cast<DiagScatterer*>(chunks[i].sc))
      return T_calc_fields<DiagScatterer*>(chunks, field);

  return T_calc_fields<MonoScatterer*>(chunks, field);
}

