
/////////////////////////////////////////////////////////////////////////////
//
// File:     S_scheme_fields.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000612
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "S_scheme_fields.h"
#include "S_scheme.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// calc_field_T
//
//   Given an incident field in field[0], calculate the exit field at the
//   other side of the scatterer 'sc' and add it to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_T
  (vector<FieldExpansion>* field, MultiScatterer* sc, Waveguide* wg)
{
  // Define matrices.
  
  const cVector& fw0((*field)[0].fw);
  const cVector& bw0((*field)[0].bw);

  const cMatrix& R12(sc->get_R12());
  const cMatrix& R21(sc->get_R21());
  const cMatrix& T12(sc->get_T12());
  const cMatrix& T21(sc->get_T21());

  const int N = global.N;

  cVector fw(N,fortranArray);   cVector bw(N,fortranArray);
  
  cMatrix ff(N,N,fortranArray); cMatrix fb(N,N,fortranArray);
  cMatrix bf(N,N,fortranArray); cMatrix bb(N,N,fortranArray);

  cMatrix  M(N,N,fortranArray); cMatrix tmp(N,N,fortranArray);

  // Calculate field.
  
  M.reference(invert_svd(T21));
  
  ff = T12 - multiply(R21,M,R12);
  fb.reference(multiply(R21,M));

  tmp = -M;
  bf.reference(multiply(tmp, R12));
  bb.reference(M);
  
  fw = multiply(ff, fw0) + multiply(fb, bw0);
  bw = multiply(bf, fw0) + multiply(bb, bw0);
    
  field->push_back(FieldExpansion(wg, fw, bw));
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_field_T for diagonal Scatterers
//
//   Given an incident field in field[0], calculate the exit field at the
//   other side of the scatterer 'sc' and add it to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_T
  (vector<FieldExpansion>* field, DiagScatterer* sc, Waveguide* wg)
{
  // Define vectors.
  
  const cVector& fw0((*field)[0].fw);
  const cVector& bw0((*field)[0].bw);

  const cVector& R12(sc->get_diag_R12());
  const cVector& R21(sc->get_diag_R21());
  const cVector& T12(sc->get_diag_T12());
  const cVector& T21(sc->get_diag_T21());

  cVector fw(global.N,fortranArray);
  cVector bw(global.N,fortranArray);

  // Calculate field.

  fw = ( (T12*T21 - R12*R21) * fw0  +  R21 * bw0 ) / T21;
  bw = (               -R12  * fw0  +        bw0 ) / T21;

  for (int i=1; i<=global.N; i++)
    if (abs(T21(i)) < 1e-10)
    {
      //cout << "Warning: small T21: setting transmitted field to zero."
      //     << endl;

      fw(i) = 0;
      bw(i) = 0;
    }
  
  field->push_back(FieldExpansion(wg, fw, bw));
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_field_T for Monoscatterers.
//
//   Given an incident field in field[0], calculate the exit field at the
//   other side of the scatterer 'sc' and add it to 'field'.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_field_T
  (vector<FieldExpansion>* field, MonoScatterer* sc, Waveguide* wg)
{
  // Define variables.
  
  const Complex& fw0((*field)[0].fw(1));
  const Complex& bw0((*field)[0].bw(1));

  const Complex& R12(sc->get_R12());
  const Complex& R21(sc->get_R21());
  const Complex& T12(sc->get_T12());
  const Complex& T21(sc->get_T21());

  cVector fw(1,fortranArray);
  cVector bw(1,fortranArray);

  // Calculate field.

  fw = ( (T12*T21 - R12*R21) * fw0  +  R21 * bw0 ) / T21;
  bw = (               -R12  * fw0  +        bw0 ) / T21;

  if (abs(T21) < 1e-10)
  {
    //cout << "Warning: small T21: setting transmitted field to zero."
    //     << endl;

    fw=0;
    bw=0;
  }

  field->push_back(FieldExpansion(wg, fw, bw));
}



/////////////////////////////////////////////////////////////////////////////
//
// split_and_calc_T
//
//   Templated function on Scatterer type, that splits up the chunks and
//   calculates the fields.
//  
/////////////////////////////////////////////////////////////////////////////

template<class T>
void split_and_calc_T
  (const vector<Chunk>& chunks,vector<FieldExpansion>* field)
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
  
  T   prev;   prev.allocRT();
  T result; result.allocRT();

  // Loop over chunks.

  for (unsigned int i=0; i<new_chunks.size(); i++)
  {
    vector<Chunk> i_chunks;
    if (i != 0)
      i_chunks.push_back(Chunk(&prev, 0.0));
    i_chunks.push_back(new_chunks[i]);
    S_scheme(i_chunks, &result);

    calc_field_T(field, &result, new_chunks[i].sc->get_ext());

    prev.swap_RT_with(result);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// S_scheme_fields_T
//  
/////////////////////////////////////////////////////////////////////////////

void S_scheme_fields_T
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field)
{
  // Check input.

  if (field->size() != 1)
  {
    py_error("Error: invalid source fields in S_scheme_fields_T.");
    exit (-1);
  }
  
  // Check type of chunks and call corresponding function.

  for (unsigned int i=0; i<chunks.size(); i++)
    if (dynamic_cast<DenseScatterer*>(chunks[i].sc))
      return split_and_calc_T<DenseScatterer>(chunks, field);

  for (unsigned int i=0; i<chunks.size(); i++)
    if (dynamic_cast<DiagScatterer*>(chunks[i].sc))
      return split_and_calc_T<DiagScatterer>(chunks, field);

  return split_and_calc_T<MonoScatterer>(chunks, field);
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_S_S
//  
/////////////////////////////////////////////////////////////////////////////

void calc_S_S(const vector<Chunk>& chunks,
              vector<FieldExpansion>* field,
              cVector* inc_right_bw=NULL)
{
  // Fields at left side.

  const int N = global.N;
  
  cVector fw0(N,fortranArray); fw0 = (*field)[0].fw;
  cVector bw0(N,fortranArray); bw0 = (*field)[0].bw;
  
  // Loop over chunks.

  for (unsigned int k=0; k<chunks.size(); k++)
  {
    // Calculate fw field after interface.

    MultiScatterer* sc = dynamic_cast<MultiScatterer*>(chunks[k].sc);

    const cMatrix& R12(sc->get_R12()); const cMatrix& R21(sc->get_R21());
    const cMatrix& T12(sc->get_T12()); const cMatrix& T21(sc->get_T21());

    cVector tmp(N,fortranArray);
    tmp = bw0 - multiply(R12,fw0);

    cMatrix inv_T21(N,N,fortranArray);
    if(global.stability != SVD)
      inv_T21.reference(invert    (T21));
    else
      inv_T21.reference(invert_svd(T21));

    cVector fw_int(N,fortranArray);
    fw_int = multiply(T12,fw0) + multiply(R21,inv_T21,tmp);

    // Calculate fw field after propagation.

    Waveguide* wg = chunks[k].sc->get_ext();

    cVector fw_prop(N,fortranArray);
    for (int i=1; i<=N; i++)    
      fw_prop(i) = fw_int(i) 
        * exp(-I * wg->get_mode(i)->get_kz() * chunks[k].d);

    // Calculate bw field after propagation.
    
    vector<Chunk> right_chunks;
    for (unsigned int j=k+1; j<chunks.size(); j++)
      right_chunks.push_back(chunks[j]);

    cVector bw_prop(N,fortranArray);
    if (right_chunks.size() != 0)
    {
      DenseScatterer right; right.allocRT();
      S_scheme(right_chunks, &right);

      bw_prop.reference(multiply(right.get_R12(),fw_prop));

      if (inc_right_bw)
        bw_prop = bw_prop + multiply(right.get_T21(), *inc_right_bw);
    }
    else
      if (inc_right_bw)
        bw_prop = *inc_right_bw;
      else
        bw_prop = 0.0;

    // Calculate bw field after interface.

    cVector bw_int(N,fortranArray);
    for (int i=1; i<=N; i++)
      bw_int(i) = bw_prop(i) 
        * exp(-I * wg->get_mode(i)->get_kz() * chunks[k].d);

    // Update field vector.

    field->push_back(FieldExpansion(wg, fw_int,  bw_int));
    field->push_back(FieldExpansion(wg, fw_prop, bw_prop));

    fw0 = fw_prop;
    bw0 = bw_prop;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_S_S_diag
//  
/////////////////////////////////////////////////////////////////////////////

void calc_S_S_diag(const vector<Chunk>& chunks,
                   vector<FieldExpansion>* field,
                   cVector* inc_right_bw)
{
  // Fields at left side.

  const int N = global.N;
  
  cVector fw0(N,fortranArray); fw0 = (*field)[0].fw;
  cVector bw0(N,fortranArray); bw0 = (*field)[0].bw;

  // A transparent diagonal scatterer is not automatically calculated.

  if (chunks[0].sc->get_inc() == chunks[0].sc->get_ext())
    chunks[0].sc->calcRT();

  // Loop over chunks.
  
  for (unsigned int k=0; k<chunks.size(); k++)
  {
    // Calculate fw field after interface.

    DiagScatterer* sc = dynamic_cast<DiagScatterer*>(chunks[k].sc);

    const cVector& R12(sc->get_diag_R12());
    const cVector& R21(sc->get_diag_R21());
    const cVector& T12(sc->get_diag_T12());
    const cVector& T21(sc->get_diag_T21());
    
    cVector fw_int(N,fortranArray);
    fw_int = T12*fw0 + R21*(bw0 - R12*fw0) / T21;

    // Calculate fw field after propagation.

    Waveguide* wg = chunks[k].sc->get_ext();

    cVector fw_prop(N,fortranArray);
    for (int i=1; i<=N; i++)    
      fw_prop(i) = fw_int(i) 
        * exp(-I * wg->get_mode(i)->get_kz() * chunks[k].d);

    // Calculate bw field after propagation.
    
    vector<Chunk> right_chunks;
    for (unsigned int j=k+1; j<chunks.size(); j++)
      right_chunks.push_back(chunks[j]);

    cVector bw_prop(N,fortranArray);
    if (right_chunks.size() != 0)
    {
      DiagScatterer right; right.allocRT();

      S_scheme(right_chunks, &right);
      
      bw_prop = right.get_diag_R12() * fw_prop;

      if (inc_right_bw)
        bw_prop = bw_prop + right.get_diag_T21() * (*inc_right_bw);
    }
    else
      if (inc_right_bw)
        bw_prop = *inc_right_bw;
      else
        bw_prop = 0.0;

    // Calculate bw field after interface.

    cVector bw_int(N,fortranArray);
    for (int i=1; i<=N; i++)
      bw_int(i) = bw_prop(i) 
        * exp(-I * wg->get_mode(i)->get_kz() * chunks[k].d);

    // Update field vector.

    field->push_back(FieldExpansion(wg, fw_int,  bw_int));
    field->push_back(FieldExpansion(wg, fw_prop, bw_prop));

    fw0 = fw_prop;
    bw0 = bw_prop;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_S_S_mono
//  
/////////////////////////////////////////////////////////////////////////////

void calc_S_S_mono(const vector<Chunk>& chunks,
                   vector<FieldExpansion>* field,
                   cVector* inc_right_bw)
{
  // Fields at left side

  cVector fw0(1,fortranArray); fw0 = (*field)[0].fw;
  cVector bw0(1,fortranArray); bw0 = (*field)[0].bw;

  // A transparent diagonal scatterer is not automatically calculated.

  if (chunks[0].sc->get_inc() == chunks[0].sc->get_ext())
    chunks[0].sc->calcRT();

  // Loop over chunks.

  for (unsigned int k=0; k<chunks.size(); k++)
  {
    // Calculate fw field after interface.

    MonoScatterer* sc = dynamic_cast<MonoScatterer*>(chunks[k].sc);

    const Complex& R12(sc->get_R12());
    const Complex& R21(sc->get_R21());
    const Complex& T12(sc->get_T12());
    const Complex& T21(sc->get_T21());
    
    cVector fw_int(1,fortranArray);
    fw_int = T12*fw0 + R21*(bw0 - R12*fw0) / T21;

    // Calculate fw field after propagation.

    Waveguide* wg = chunks[k].sc->get_ext();

    cVector fw_prop(1,fortranArray);
    fw_prop(1) = fw_int(1) * exp(-I * wg->get_mode(1)->get_kz() * chunks[k].d);

    // Calculate bw field after propagation.
    
    vector<Chunk> right_chunks;
    for (unsigned int j=k+1; j<chunks.size(); j++)
    {
      right_chunks.push_back(chunks[j]);
      chunks[j].sc->calcRT(); // Transparent scatterer isn't auto-calculateded.
    }

    cVector bw_prop(1,fortranArray);
    if (right_chunks.size() != 0)
    {
      MonoScatterer right; right.allocRT();
      S_scheme(right_chunks, &right);
      
      bw_prop = right.get_R12() * fw_prop;

      if (inc_right_bw)
        bw_prop = bw_prop + right.get_T21() * (*inc_right_bw);
    }
    else
      if (inc_right_bw)
        bw_prop = *inc_right_bw;
      else
        bw_prop = 0.0;

    // Calculate bw field after interface.

    cVector bw_int(1,fortranArray);
    bw_int(1) = bw_prop(1) * exp(-I * wg->get_mode(1)->get_kz() * chunks[k].d);

    // Update field vector.

    field->push_back(FieldExpansion(wg, fw_int,  bw_int));
    field->push_back(FieldExpansion(wg, fw_prop, bw_prop));

    fw0 = fw_prop;
    bw0 = bw_prop;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// S_scheme_fields_S
//  
/////////////////////////////////////////////////////////////////////////////

void S_scheme_fields_S
  (const vector<Chunk>& chunks, vector<FieldExpansion>* field,
   cVector* inc_right_bw)
{
  // Check input.

  if (field->size() != 1)
  {
    py_error("Error: invalid source fields in S_scheme_fields_S.");
    exit (-1);
  }

  // Check type of chunks and call corresponding function.

  for (unsigned int i=0; i<chunks.size(); i++) 
    if (    ! dynamic_cast<DiagScatterer*>(chunks[i].sc) 
         && ! chunks[i].sc->is_mono())
      return calc_S_S(chunks, field, inc_right_bw);
  
  for (unsigned int i=0; i<chunks.size(); i++)
    if (! chunks[i].sc->is_mono())
      return calc_S_S_diag(chunks, field, inc_right_bw);

  return calc_S_S_mono(chunks, field, inc_right_bw);
}
