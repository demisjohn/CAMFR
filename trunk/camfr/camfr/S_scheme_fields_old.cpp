
/////////////////////////////////////////////////////////////////////////////
//
// File:     S_scheme_fields_old.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000612
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "S_scheme_fields_old.h"
#include "S_scheme.h"
#include "icache.h"

/////////////////////////////////////////////////////////////////////////////
//
// calc_field_S
//
//   Given fields incident from the outside left (1) and no field
//   incident from the right (3), calculate the field at a position (2)
//   inside the stack.
//   Add field expansion to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_S(vector<FieldExpansion>* field,
                  MultiScatterer* sc12, MultiScatterer* sc23, Waveguide* wg)
{
  // Define matrices.

  const cVector& fw1((*field)[0].fw);

  const cMatrix& R21(sc12->get_R21());
  const cMatrix& T12(sc12->get_T12());
  const cMatrix& R23(sc23->get_R12());

  const int N = global.N;

  cVector fw2(N,fortranArray), bw2(N,fortranArray);
  
  cMatrix U1(N,N,fortranArray), M(N,N,fortranArray), tmp(N,N,fortranArray);

  U1 = 0.0;
  for (int i=1; i<=N; i++)
    U1(i,i) = Complex(1.0, 0.0);

  // Calculate field.
  
  tmp = U1 - multiply(R21, R23);
  if(global.stability != SVD)
    M.reference(invert    (tmp));
  else
    M.reference(invert_svd(tmp));
  
  fw2.reference(multiply(M, T12, fw1));
  bw2.reference(multiply(R23, fw2));
    
  field->push_back(FieldExpansion(*wg, fw2, bw2));
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_field_S_end
//
//   Given fields incident from the outside left (1) and no field
//   incident from the right (3), calculate the field at a position (3).
//   Add field expansion to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_S_end(vector<FieldExpansion>* field,
                      MultiScatterer* sc12, Waveguide* wg)
{
  // Define matrices.

  const cVector& fw1((*field)[0].fw);
  const cMatrix& T12(sc12->get_T12());

  cVector fw3(global.N,fortranArray), bw3(global.N,fortranArray);

  // Calculate field.
  
  fw3.reference(multiply(T12, fw1));
  bw3 = 0;
    
  field->push_back(FieldExpansion(*wg, fw3, bw3));
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_field_S
//
//   Given fields incident from the outside left (1) and right (3),
//   calculate the field at a position (2) inside the stack.
//   Add field expansion to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_S(vector<FieldExpansion>* field,
                  DiagScatterer* sc12, DiagScatterer* sc23, Waveguide* wg)
{
  // Define vectors.

  const cVector& fw1((*field)[0].fw);
  
  const cVector& R21(sc12->get_diag_R21());
  const cVector& T12(sc12->get_diag_T12());
  const cVector& R23(sc23->get_diag_R12());
  
  cVector fw2(global.N,fortranArray), bw2(global.N,fortranArray);

  // Calculate field.
  
  fw2 = T12 * fw1 / (1.0 - R21 * R23);
  bw2 = R23 * fw2;
    
  field->push_back(FieldExpansion(*wg, fw2, bw2));
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_field_S_end
//
//   Given fields incident from the outside left (1) and no field
//   incident from the right (3), calculate the field at a position (3).
//   Add field expansion to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_S_end(vector<FieldExpansion>* field,
                      DiagScatterer* sc12, Waveguide* wg)
{
  // Define matrices.

  const cVector& fw1((*field)[0].fw);
  const cVector& T12(sc12->get_diag_T12());

  cVector fw3(global.N,fortranArray), bw3(global.N,fortranArray);

  // Calculate field.
  
  fw3 = T12 * fw1;
  bw3 = 0;
    
  field->push_back(FieldExpansion(*wg, fw3, bw3));
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_field_S
//
//   Given fields incident from the outside left (1) and right (3),
//   calculate the field at a position (2) inside the stack.
//   Add field expansion to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_S(vector<FieldExpansion>* field,
                  MonoScatterer* sc12, MonoScatterer* sc23, Waveguide* wg)
{
  // Define vector.

  const cVector& fw1((*field)[0].fw);
  
  const Complex& R21(sc12->get_R21());
  const Complex& T12(sc12->get_T12());
  const Complex& R23(sc23->get_R12());

  cVector fw2(global.N,fortranArray), bw2(global.N,fortranArray);

  // Calculate field.

  fw2 = T12 * fw1 / (1.0 - R21 * R23);
  bw2 = R23 * fw2;
    
  field->push_back(FieldExpansion(*wg, fw2, bw2));
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_field_S_end
//
//   Given fields incident from the outside left (1) and no field
//   incident from the right (3), calculate the field at a position (3).
//   Add field expansion to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_S_end(vector<FieldExpansion>* field,
                      MonoScatterer* sc12, Waveguide* wg)
{
  // Define variables.

  const cVector& fw1((*field)[0].fw);
  const Complex& T12(sc12->get_T12());

  cVector fw3(global.N,fortranArray), bw3(global.N,fortranArray);

  // Calculate field.
  
  fw3 = T12 * fw1;
  bw3 = 0;
    
  field->push_back(FieldExpansion(*wg, fw3, bw3));
}



/////////////////////////////////////////////////////////////////////////////
//
// split_and_calc_S
//
//   Templated function on Scatterer type, that splits up the chunks and
//   calculates the fields.
//  
/////////////////////////////////////////////////////////////////////////////

template<class T>
void split_and_calc_S
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field)
{
  // We want to calculate the fields in two places per chunk, once before
  // and once after the propagation. Create a new series of chunks to
  // accomodate this.

  vector<Chunk> new_chunks;
  for (unsigned int i=0; i<chunks.size(); i++)
  {
    chunks[i].sc->calcRT();
    
    Chunk i_noprop(chunks[i].sc, 0.0);
    new_chunks.push_back(i_noprop);

    Scatterer* transparent = interface_cache.get_interface
      (chunks[i].sc->get_ext(), chunks[i].sc->get_ext());

    transparent->calcRT();
    
    Chunk i_prop(transparent, chunks[i].d);
    new_chunks.push_back(i_prop);
  }

  // Set up storage.
  
  T prev12; prev12.allocRT();
  T  res12;  res12.allocRT();
  T  res23;  res23.allocRT();
  
  // Loop over chunks.

  for (unsigned int i=0; i<new_chunks.size(); i++)
  {
    vector<Chunk> left_chunks;
    if (i != 0)
      left_chunks.push_back(Chunk(&prev12, 0.0));
    left_chunks.push_back(new_chunks[i]);
    S_scheme(left_chunks, &res12);

    if (i != new_chunks.size()-1)
    {
      vector<Chunk> right_chunks;
      for (unsigned int j=i+1; j<new_chunks.size(); j++)
        right_chunks.push_back(new_chunks[j]);
      S_scheme(right_chunks, &res23);

      calc_field_S(field, &res12, &res23, new_chunks[i].sc->get_ext());
    }
    else
      calc_field_S_end(field, &res12, new_chunks[i].sc->get_ext());

    prev12.swap_RT_with(res12);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// S_scheme_fields_S_old
//  
/////////////////////////////////////////////////////////////////////////////

void S_scheme_fields_S_old
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field)
{
  // Check input.

  if (field->size() != 1)
  {
    cerr << "Error: invalid source fields in S_scheme_fields_S." << endl;
    exit (-1);
  }
  
  // Check type of chunks and call corresponding function.

  for (unsigned int i=0; i<chunks.size(); i++)
    if (dynamic_cast<DenseScatterer*>(chunks[i].sc))
      return split_and_calc_S<DenseScatterer>(chunks, field);
  
  for (unsigned int i=0; i<chunks.size(); i++)
    if (dynamic_cast<DiagScatterer*>(chunks[i].sc))
      return split_and_calc_S<DiagScatterer>(chunks, field);

  return split_and_calc_S<MonoScatterer>(chunks, field);
}
