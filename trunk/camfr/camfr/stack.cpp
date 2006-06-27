
/////////////////////////////////////////////////////////////////////////////
//
// File:     stack.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000121
// Version:  2.1
//
// Copyright (C) 1998-2000 Peter Bienstman - Ghent University
//
////////////////////////////////////////////////////////////////////////////

#include "stack.h"
#include "util/index.h"
#include "S_scheme.h"
#include "S_scheme_fields.h"
#include "T_scheme_fields.h"
#include "bloch.h"
#include "primitives/blochsection/blochsection.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::StackImpl
//  
/////////////////////////////////////////////////////////////////////////////

StackImpl::StackImpl
(const vector<Chunk>& chunks_, unsigned int no_of_periods_)
  : chunks(chunks_), no_of_periods(no_of_periods_)
{
  // Check validity of chunks.

  for (unsigned int i=0; i<chunks.size()-1; i++)
  {
    Waveguide* wg1 = chunks[i]  .sc->get_ext();
    Waveguide* wg2 = chunks[i+1].sc->get_inc();
    
    if ( wg1 && wg2 && (*wg1 != *wg2) )
    {
      py_error("Error: intermediate incidence and exit media don't match.");
      return;
    }
  }
  
  Waveguide* inc = chunks[0].sc->get_inc();
  Waveguide* ext = chunks[chunks.size()-1].sc->get_ext();

  if ( (no_of_periods > 1) && (*inc != *ext) )
  {
    py_error("Error: inc and exit media don't match for periodic extension.");
    return;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::StackImpl
//
//   Create chunks from an expression where the exit medium of each term
//   matches the incidence medium of the next.
//
//   Each waveguide is joined in a chunk with its preceding scatterer.
//  
/////////////////////////////////////////////////////////////////////////////

StackImpl::StackImpl(const Expression& e_, unsigned int no_of_periods_)
  : no_of_periods(no_of_periods_)
{
  Expression e = e_; // Create working copy.

  // Check validity of expression.

  for (unsigned int i=0; i<e.get_size(); i++)
  {    
    if (    (e.get_term(i)->get_type() == MATERIAL)
         || (e.get_term(i)->get_type() == MAT_EXPRESSION) )
    {
      py_error("Error: unexpected material term in stack expression.");
      return;
    }
  } 
  
  for (unsigned int i=0; i<e.get_size()-1; i++)
  {    
    if (e.get_term(i)->get_ext() != e.get_term(i+1)->get_inc())
    {
      py_error("Error: intermediate incidence and exit media don't match.");
      return;
    }
  }  

  if ( (no_of_periods > 1) && (e.get_inc() != e.get_ext()) )
  {
    py_error("Error: inc and exit media don't match for periodic extension.");
    return;
  }

  // Corner case of expression consisting of single waveguide.

  if ((e.get_size() == 1) && (e.get_term(0)->get_type() == WAVEGUIDE))
  {
    Waveguide* wg = e.get_term(0)->get_wg();
    Scatterer* sc = interface_cache.get_interface(wg, wg);
    e.insert_term_front(Term(*sc));
  }
  
  // Create chunks.
  
  for (unsigned int i=0; i<e.get_size(); i++)
  {
    Term* t1 = e.get_term(i);
    Term* t2 = e.get_term(i+1);

    // No waveguide.
    
    if (t1->get_type() != WAVEGUIDE)
      if ( (i+1 < e.get_size()) && (t2->get_type() == WAVEGUIDE) )
      {
        Complex d = t2->get_d();

        // Join subsequent same waveguides.

        Term* t;
        unsigned int j = i+2;
        while ( (j < e.get_size()) && ((t=e.get_term(j++)) != NULL)
                 && ( t->get_type() == WAVEGUIDE)
                 && (*(t->get_wg()) == *(t2->get_wg())) )
        {
          d += t->get_d();
          i++;
        }
        
        chunks.push_back(Chunk(t1->get_sc(), d));
        i++;
      }
      else
        chunks.push_back(Chunk(t1->get_sc(), 0.0));
    
    // Waveguide.

    if (t1->get_type() == WAVEGUIDE)
      if (    (i>0)
           || ( (i==0) && abs(t1->get_d()) ))
      {
        py_error("Error: unexpected waveguide.");
        return;
      }
  }  
}



/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::stack_get_total_thickness
//
/////////////////////////////////////////////////////////////////////////////

Complex StackImpl::stack_get_total_thickness() const
{
  Complex d = 0.0;
  
  for (unsigned int i=0; i<chunks.size(); i++)
    d += chunks[i].sc->get_total_thickness() + chunks[i].d;

  return Real(no_of_periods) * d;
}



/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::stack_get_materials
//  
/////////////////////////////////////////////////////////////////////////////

vector<Material*> StackImpl::stack_get_materials() const
{
  vector<Material*> materials;

  for (unsigned int i=0; i<chunks.size(); i++)
  {
    vector<Material*> ch_materials = chunks[i].sc->get_materials();
    for (unsigned int j=0; j<ch_materials.size(); j++)
      materials.push_back(ch_materials[j]);
  }

  return materials;
}



/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::stack_contains
//  
/////////////////////////////////////////////////////////////////////////////

bool StackImpl::stack_contains(const Material& m) const
{
  for (unsigned int i=0; i<chunks.size(); i++)
    if (chunks[i].sc->contains(m))
      return true;

  return false;
}



/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::stack_no_gain_present
//  
/////////////////////////////////////////////////////////////////////////////

bool StackImpl::stack_no_gain_present() const
{
  for (unsigned int i=0; i<chunks.size(); i++)
  {
    if (!chunks[i].sc->no_gain_present())
      return false;
        
    if (imag(chunks[i].d) > 1e-3)
      return false;
  }

  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::stack_all_layers_uniform
//  
/////////////////////////////////////////////////////////////////////////////

bool StackImpl::stack_all_layers_uniform() const
{
  for (unsigned int i=0; i<chunks.size(); i++)
    if (!chunks[i].sc->all_layers_uniform())
      return false;

  return true; 
}



/////////////////////////////////////////////////////////////////////////////
//
// StackImpl::get_thicknesses
//  
/////////////////////////////////////////////////////////////////////////////

vector<Complex*> StackImpl::get_thicknesses() const
{
  vector<Complex*> thicknesses;

  for (unsigned int i=0; i<chunks.size(); i++)
    thicknesses.push_back(const_cast<Complex*>(&(chunks[i].d)));
  
  return thicknesses;
}



/////////////////////////////////////////////////////////////////////////////
//
// stack_calcRT
//
//   Template for code reuse in calcRT in the different Stacks with
//   T = DenseStack, DiagStack or MonoStack.
//  
/////////////////////////////////////////////////////////////////////////////

template <class T>
void stack_calcRT(T* stack)
{ 
  if (stack->no_of_periods == 1)
  {
    S_scheme(stack->chunks, stack);
    return;
  }

  // Calculate basic period.

  T period(stack->chunks);
  period.calcRT();

  // Do periodic extension.

  // The straightforward way: O(n)

  // vector<Chunk> extension;
  // for (unsigned int i=1; i<=stack->no_of_periods; i++)
  //   extension.push_back(period);
  //
  // S_scheme(extension, stack);

  // The smart way: O(log(n))
  // Based on the observation that
  //   x^(2n)   = (x^2)^n
  //   x^(2n+1) =  x.x^2n

  // Set up temporary storage.

  T tmp_stack;
  tmp_stack.copy_RT_from(period);

  // We will manipulate pointers instead of copying matrices.

  T* result(&tmp_stack);
  T* storage(stack);
  T* swap;

  // Determine bit mask for MSB in no_of_periods.

  unsigned int mask = 1;
  for (unsigned int i = stack->no_of_periods; i != 1; i >>= 1)
    mask <<= 1;

  // Loop from from the bit to the right of MSB to the LSB.

  while (mask != 1)
  {
    mask >>= 1;

    // Square result calculated so far.

    vector<Chunk> extension;
    extension.push_back(result);
    extension.push_back(result);

    S_scheme(extension, storage);

    swap = result; result = storage; storage = swap;

    if (stack->no_of_periods & mask) // bit is 1
    {
      // Add orginal period to result calculated so far.

      extension.clear();
      extension.push_back(result);
      extension.push_back(&period);

      S_scheme(extension, storage);

      swap = result; result = storage; storage = swap;
    }
  }

  // Place final results in the right place.

  if (result != stack)
    stack->copy_RT_from(*result);
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseStack::DenseStack
//  
/////////////////////////////////////////////////////////////////////////////

DenseStack::DenseStack
 (const vector<Chunk>& chunks, unsigned int no_of_periods)
  : StackImpl(chunks, no_of_periods)
{ 
  inc = chunks[0].sc->get_inc();
  ext = chunks[chunks.size()-1].sc->get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseStack::DenseStack
//  
/////////////////////////////////////////////////////////////////////////////

DenseStack::DenseStack(const Expression& e, unsigned int no_of_periods)
  : StackImpl(e, no_of_periods)
{ 
  inc = chunks[0].sc->get_inc();
  ext = chunks[chunks.size()-1].sc->get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseStack::calcRT
//  
/////////////////////////////////////////////////////////////////////////////

void DenseStack::calcRT()
{
  if (!recalc_needed())
    return;

  allocRT();

  for (unsigned int i=0; i<chunks.size(); i++)
    chunks[i].sc->calcRT();
  
  stack_calcRT<DenseStack>(this);

  // Remember wavelength and gain these matrices were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat_n = global.gain_mat->n();
  last_slab_ky = global.slab_ky;
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseStack::freeRT
//  
/////////////////////////////////////////////////////////////////////////////

void DenseStack::freeRT()
{
  DenseScatterer::freeRT();
  
  for (unsigned int i=0; i<chunks.size(); i++)
    chunks[i].sc->freeRT();
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagStack::DiagStack
//  
/////////////////////////////////////////////////////////////////////////////

DiagStack::DiagStack
 (const vector<Chunk>& chunks, unsigned int no_of_periods)
  : StackImpl(chunks, no_of_periods)
{ 
  inc = chunks[0].sc->get_inc();
  ext = chunks[chunks.size()-1].sc->get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagStack::DiagStack
//  
/////////////////////////////////////////////////////////////////////////////

DiagStack::DiagStack(const Expression& e, unsigned int no_of_periods)
  : StackImpl(e, no_of_periods)
{  
  inc = chunks[0].sc->get_inc();
  ext = chunks[chunks.size()-1].sc->get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagStack::calcRT
//  
/////////////////////////////////////////////////////////////////////////////

void DiagStack::calcRT()
{ 
  if (!recalc_needed())
    return;

  allocRT();

  for (unsigned int i=0; i<chunks.size(); i++)
    chunks[i].sc->calcRT();

  stack_calcRT<DiagStack>(this);
  
  // Remember wavelength and gain these matrices were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat_n = global.gain_mat->n();
  last_slab_ky = global.slab_ky;
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagStack::freeRT
//  
/////////////////////////////////////////////////////////////////////////////

void DiagStack::freeRT()
{
  DiagScatterer::freeRT();
  
  for (unsigned int i=0; i<chunks.size(); i++)
    chunks[i].sc->freeRT();
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoStack::MonoStack
//  
/////////////////////////////////////////////////////////////////////////////

MonoStack::MonoStack
 (const vector<Chunk>& chunks, unsigned int no_of_periods)
  : StackImpl(chunks, no_of_periods)
{
  inc = chunks[0].sc->get_inc();
  ext = chunks[chunks.size()-1].sc->get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoStack::MonoStack
//  
/////////////////////////////////////////////////////////////////////////////

MonoStack::MonoStack(const Expression& e, unsigned int no_of_periods)
  : StackImpl(e, no_of_periods)
{   
  inc = chunks[0].sc->get_inc();
  ext = chunks[chunks.size()-1].sc->get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoStack::calcRT
//  
/////////////////////////////////////////////////////////////////////////////

void MonoStack::calcRT()
{
  int old_N = global.N;
  global.N = 1;

  if (!recalc_needed())
    return;

  for (unsigned int i=0; i<chunks.size(); i++)
    chunks[i].sc->calcRT();
  
  stack_calcRT<MonoStack>(this);

  global.N = old_N;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::Stack
//  
/////////////////////////////////////////////////////////////////////////////

Stack::Stack(const Expression& e, unsigned int no_of_periods_)
  : expression(e), no_of_periods(no_of_periods_), 
    inc_field(fortranArray), inc_field_bw(fortranArray)
{
  sc = create_sc(expression, no_of_periods);
  flat_sc = create_sc(expression.flatten());
  calc_interface_positions();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::Stack
//  
/////////////////////////////////////////////////////////////////////////////

Stack::Stack(const Term& t)
  : expression(Expression(t)), no_of_periods(1), 
    inc_field(fortranArray), inc_field_bw(fortranArray)
{
  sc = create_sc(expression, no_of_periods);
  flat_sc = create_sc(expression.flatten());
  calc_interface_positions();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::Stack
//  
/////////////////////////////////////////////////////////////////////////////

Stack::Stack(const Stack& s)
  : expression(s.expression), no_of_periods(s.no_of_periods),
    interface_positions(s.interface_positions),
    interface_field(s.interface_field),
    inc_field(fortranArray), inc_field_bw(fortranArray)
{
  inc_field.resize(s.inc_field.shape());
  inc_field = s.inc_field;

  inc_field_bw.resize(s.inc_field_bw.shape());
  inc_field_bw = s.inc_field_bw;
  
  sc = create_sc(expression, no_of_periods);
  flat_sc = create_sc(expression.flatten());
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::operator=
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::operator=(const Expression& e)
{
  delete sc;
  delete flat_sc;

  expression    = e;
  no_of_periods = 1;
  
  sc = create_sc(expression, no_of_periods);
  flat_sc = create_sc(expression.flatten());
  calc_interface_positions();

  inc_field = 0;
  inc_field_bw = 0;

  interface_field.clear();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::operator=
//  
/////////////////////////////////////////////////////////////////////////////

Stack& Stack::operator=(const Stack& s)
{
  if (this == &s)
    return *this;
  
  delete sc;
  delete flat_sc;
  
  expression    = s.expression;
  no_of_periods = s.no_of_periods;
  
  sc = create_sc(expression, no_of_periods);
  flat_sc = create_sc(expression.flatten());
  calc_interface_positions();
  
  inc_field = s.inc_field;
  inc_field_bw = s.inc_field_bw;
  interface_field = s.interface_field;
  
  return *this;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::single_material
//
/////////////////////////////////////////////////////////////////////////////

bool Stack::single_material() const
{
  vector<Material*> materials = get_materials();

  for (unsigned int i=0; i<materials.size(); i++)
    for (unsigned int j=i+1; j<materials.size(); j++)
      if ( (*materials[i]) != (*materials[j]) )
        return false;

  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::calcRT
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::calcRT()
{
  if (sc)
    sc->calcRT();
  else
    py_error("No scatterer defined.");
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::allocRT
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::allocRT()
{
  if (sc)
    sc->allocRT();
  else
    py_error("No scatterer defined.");

}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::freeRT
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::freeRT()
{
  MultiScatterer* sc_m = dynamic_cast<MultiScatterer*>(sc);
  
  if (sc_m)
    sc_m->freeRT();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::set_inc_field
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::set_inc_field(const cVector& inc_field_, cVector* inc_field_bw_)
{
  inc_field.resize(inc_field_.shape());
  inc_field = inc_field_;

  inc_field_bw.resize(inc_field.shape());
  if (inc_field_bw_)
  {
    inc_field_bw.resize(inc_field_bw_->shape());
    inc_field_bw = *inc_field_bw_;
    bw_inc = true;
  }
  else
  {
    inc_field_bw = 0.0;
    bw_inc = false;
  }

  interface_field.clear();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_refl_field
//  
/////////////////////////////////////////////////////////////////////////////

cVector Stack::get_refl_field()
{
  if (interface_field.size())
    return interface_field[0].bw;

  calcRT();

  cVector refl_field(inc_field.rows(), fortranArray);
  refl_field.reference(multiply(as_multi()->get_R12(), inc_field));

  if (bw_inc)
    refl_field += multiply(as_multi()->get_T21(), inc_field_bw);

  FieldExpansion inc(get_inc(), inc_field, refl_field);
  interface_field.push_back(inc);
  
  return refl_field;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_trans_field
//  
/////////////////////////////////////////////////////////////////////////////

cVector Stack::get_trans_field()
{  
  if (interface_field.size() > 1)
    return interface_field.back().fw;

  calcRT();

  cVector trans_field(inc_field.rows(), fortranArray);
  trans_field.reference(multiply(as_multi()->get_T12(), inc_field));

  if (bw_inc)
    trans_field += multiply(as_multi()->get_R21(), inc_field_bw);
  
  return trans_field;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::inc_field_expansion
//  
/////////////////////////////////////////////////////////////////////////////

FieldExpansion Stack::inc_field_expansion()
{
  if (inc_field.rows() == 0)
  {
    py_error("Error: incident field not set.");
    exit(-1);
  }

  if (interface_field.size())
    return interface_field[0];

  calcRT();

  cVector refl_field(inc_field.rows(), fortranArray);
  refl_field.reference(multiply(as_multi()->get_R12(), inc_field));

  if (bw_inc)
    refl_field += multiply(as_multi()->get_T21(), inc_field_bw);

  FieldExpansion inc_field_exp(get_inc(), inc_field, refl_field);
  interface_field.push_back(inc_field_exp);

  return inc_field_exp;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::ext_field_expansion
//  
/////////////////////////////////////////////////////////////////////////////

FieldExpansion Stack::ext_field_expansion()
{
  if (inc_field.rows() == 0)
  {
    py_error("Error: incident field not set.");
    exit(-1);
  }

  if (interface_field.size() > 1)
    return interface_field.back();

  calcRT();

  cVector trans_field(inc_field.rows(), fortranArray);
  trans_field.reference(multiply(as_multi()->get_T12(), inc_field));

  if (bw_inc)
    trans_field += multiply(as_multi()->get_R21(), inc_field_bw);

  return FieldExpansion(get_ext(), trans_field, inc_field_bw);
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::set_interface_field
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::set_interface_field(const vector<FieldExpansion>& field)
{
  interface_field.clear();
  for (unsigned int i=0; i<field.size(); i++)
    interface_field.push_back(field[i]);

  inc_field = interface_field[0].fw;

  if (bw_inc) // FIXME: not entirely general.
    inc_field_bw = interface_field.back().bw;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_interface_field
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::get_interface_field(vector<FieldExpansion>* field)
{
  field->clear();

  if (interface_field.size() <= 1)
    calc_interface_fields();

  for (unsigned int i=0; i<interface_field.size(); i++)
    field->push_back(interface_field[i]);
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::field
//  
/////////////////////////////////////////////////////////////////////////////

Field Stack::field(const Coord& coord)
{
  // If needed, calculate field at each interface.
  
  if (interface_field.size() <= 1)
    calc_interface_fields();

  // Calculate field expansion.

  unsigned int index =
    index_lookup(coord.z, coord.z_limit, interface_positions);

  const vector<Chunk>* chunks
    = dynamic_cast<StackImpl*>(flat_sc)->get_chunks();
  
  if (index == chunks->size())
    index--;
  
  Waveguide* wg = (*chunks)[index].sc->get_ext();

  // fw and bw matrices for blochstacks have only global.N components.

  cVector fw(global.N,fortranArray);
  cVector bw(global.N,fortranArray);

  fw_bw_field(coord, &fw, &bw);

  // Calculate total field. Note that the z-dependence has already been
  // taken care of.

  FieldExpansion field_expansion(wg, fw, bw);

  const Coord c(coord.c1,       coord.c2,       0,
                coord.c1_limit, coord.c2_limit, coord.z_limit);

  return field_expansion.field(c);
}



/////////////////////////////////////////////////////////////////////////////
//
// fw_bw_field
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::fw_bw_field(const Coord& coord, cVector* fw, cVector* bw)
{ 
  // Check if we are calculating the field in an internal point or not.

  const Complex last_z = interface_positions.back();

  bool coord_in_inc_medium = (real(coord.z) < 0)
    || ((abs(coord.z) < 1e-10) && (coord.z_limit == Min));
  
  bool coord_in_ext_medium = (real(coord.z) > real(last_z))
    || ((abs(coord.z - last_z) < 1e-10) && (coord.z_limit == Plus));

  if (coord_in_inc_medium || coord_in_ext_medium)
  {
    // If we already have all the interface fields, we can use those.
  
    if (interface_field.size() > 1)
    {
      if (coord_in_inc_medium)
      {    
        FieldExpansion f(interface_field[0].propagate(coord.z));

        *fw = f.fw; *bw = f.bw;

        return;
      }

      if (coord_in_ext_medium)
      {    
        FieldExpansion f(interface_field.back().
                         propagate(coord.z-get_total_thickness()));
    
        *fw = f.fw; *bw = f.bw;
    
        return;
      }
    }

    // Else, no need to calculate all the interface fields.

    if (inc_field.rows() == 0)
    {
      py_error("Error: incident field not set.");
      return;
    }

    calcRT();

    if (coord_in_inc_medium)
    { 
      cVector left_bw(inc_field.rows(), fortranArray);

      if (as_multi())
      {
        left_bw.reference(multiply(as_multi()->get_R12(), inc_field));

        if (bw_inc)
          left_bw += multiply(as_multi()->get_T21(), inc_field_bw);
      }
      else
      {
        left_bw(1) = R12(1,1) * inc_field(1);

        if (bw_inc)
          left_bw(1) += T21(1,1) * inc_field_bw(1);
      }

      FieldExpansion inc(get_inc(), inc_field, left_bw);
      FieldExpansion f(inc.propagate(coord.z));
      *fw = f.fw; *bw = f.bw;

      return;
    }

    if (coord_in_ext_medium)
    {
      cVector right_fw(inc_field.rows(), fortranArray);

      if (as_multi())
      {
        right_fw.reference(multiply(as_multi()->get_T12(), inc_field));
        
        if (bw_inc)
          right_fw += multiply(as_multi()->get_R21(), inc_field_bw);
      }
      else
      {
        right_fw(1) = T12(1,1) * inc_field(1);

        if (bw_inc)
          right_fw(1) += R21(1,1) * inc_field_bw(1);
      }

      FieldExpansion ext(get_ext(), right_fw, inc_field_bw);
      FieldExpansion f(ext.propagate(coord.z-get_total_thickness()));
      *fw = f.fw; *bw = f.bw;

      return;
    }
  }

  // If needed, calculate field at each interface.
  
  if (interface_field.size() <= 1)
    calc_interface_fields();

  // Calculate index in chunk vector (first chunk is zero).
  
  const unsigned int index =
    index_lookup(coord.z, coord.z_limit, interface_positions);

  const vector<Chunk>* chunks
    = dynamic_cast<StackImpl*>(flat_sc)->get_chunks();

  Waveguide* wg = (*chunks)[index].sc->get_ext();

  // Calculate (positive) distance from previous discontinuity.
  
  const Complex d_prev =
    (index==0) ? coord.z : coord.z - interface_positions[index-1];
  
  // Calculate (positive) distance to next discontinuity.
  
  const Complex d_next = interface_positions[index] - coord.z;

  // Calculate indices in interface_field.

  const unsigned int prev_index = 2*index+1;
  const unsigned int next_index = 2*index+2;

  // Calculate field expansion at desired z position.

  for (int i=1; i<=wg->N(); i++)
  {
    const Complex kz = wg->get_mode(i)->get_kz();
    
    if (imag(kz) < 0) // Propagate forwards.
    {
      (*fw)(i) = interface_field[prev_index].fw(i) * exp(-I * kz * d_prev);
      (*bw)(i) = interface_field[next_index].bw(i) * exp(-I * kz * d_next);
    }
    else // Propagate backwards.
    {
      (*fw)(i) = interface_field[next_index].fw(i) * exp( I * kz * d_next);
      (*bw)(i) = interface_field[prev_index].bw(i) * exp( I * kz * d_prev);
    }
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// prop_int
//
//  Calculates 'propagation integral' Int(exp(-I*k*x),x=0...d).
//  
/////////////////////////////////////////////////////////////////////////////

inline Complex prop_int(const Complex& k, const Complex& d)
{
  return (abs(k) < 1e-10) ? d : (exp(-I*k*d) - 1.0) / (-I*k);
}



/////////////////////////////////////////////////////////////////////////////
//
// safe_mult
//
//  Safe multiplication of number with propagation integral.
//  
/////////////////////////////////////////////////////////////////////////////

Complex safe_mult(const Complex& a, const Complex& prop_int)
{
  if (    (abs(a) < global.unstable_exp_threshold)
       && (abs(prop_int) > 1) )
      return 0;

  return a*prop_int;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::lateral_S_flux
//
//  Returns lateral flux of Poynting vector at a given lateral coordinate c1.
//  S_k gives the flux in each of the chunks.
//  
/////////////////////////////////////////////////////////////////////////////

Complex Stack::lateral_S_flux(const Complex& c1, vector<Complex>* S_k)
{  
  // If needed, calculate field at each interface.

  if (interface_field.size() <= 1)
    calc_interface_fields();
  
  // Loop over all chunks.
  
  const vector<Chunk>* chunks
    = dynamic_cast<StackImpl*>(flat_sc)->get_chunks();
  
  Complex S_total=0;
  for (unsigned int k=0; k<chunks->size(); k++)
  {
    const Waveguide* wg = (*chunks)[k].sc->get_ext();
    const Complex     d = (*chunks)[k].d;
    
    // Calculate field values of each forward mode.
    // Note that this vector offsets at 0, while mode numbers offset at 1.

    vector<Field> f;
    for (unsigned int i=1; i<=global.N; i++)
      f.push_back(wg->get_mode(i)->field(Coord(c1,0,0)));

    // Calculate lateral power flux.

    Complex S=0;
    for (unsigned int i=1; i<=global.N; i++)
      for (unsigned int j=1; j<=global.N; j++)
      { 
        const Complex A_i = interface_field[2*k+1].fw(i);
        const Complex A_j = interface_field[2*k+1].fw(j);

        const Complex B_i = interface_field[2*k+1].bw(i);
        const Complex B_j = interface_field[2*k+1].bw(j);

        const Complex k_i = wg->get_mode(i)->get_kz();
        const Complex k_j = wg->get_mode(j)->get_kz();        

        const Complex fw_fw =     A_i*f[i-1].E2  *   conj(A_j*f[j-1].Hz)
                              -   A_i*f[i-1].Ez  *   conj(A_j*f[j-1].H2);

        const Complex fw_bw =     A_i*f[i-1].E2  *   conj(B_j*f[j-1].Hz)
                              -   A_i*f[i-1].Ez  * (-conj(B_j*f[j-1].H2));

        const Complex bw_fw =     B_i*f[i-1].E2  *   conj(A_j*f[j-1].Hz)
                              - (-B_i*f[i-1].Ez) *   conj(A_j*f[j-1].H2);

        const Complex bw_bw =     B_i*f[i-1].E2  *   conj(B_j*f[j-1].Hz)
                              - (-B_i*f[i-1].Ez) * (-conj(B_j*f[j-1].H2));
        
        S +=   safe_mult( fw_fw, prop_int( k_i + conj(-k_j), d) )
             + safe_mult( fw_bw, prop_int( k_i - conj(-k_j), d) )
             + safe_mult( bw_fw, prop_int(-k_i + conj(-k_j), d) )
             + safe_mult( bw_bw, prop_int(-k_i - conj(-k_j), d) );
      }
    
    // Apply geometry dependent correction factor and update results.

    S *= wg->lateral_S_corr(Coord(c1,0,0));
    
    S_total += S;
    if (S_k)
      S_k->push_back(S); 
  }

  return S_total;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::material_at
//  
/////////////////////////////////////////////////////////////////////////////

Material* Stack::material_at(const Coord& coord) const
{  
  if ( (abs(coord.z) < 1e-10) && (coord.z_limit == Min))
    return get_inc()->material_at(coord);

  unsigned int index =
    index_lookup(coord.z, coord.z_limit, interface_positions);

  if (index == interface_positions.size())
    return get_ext()->material_at(coord);

  const vector<Chunk>* chunks
    = dynamic_cast<StackImpl*>(flat_sc)->get_chunks();
  
  return (*chunks)[index].sc->get_ext()->material_at(coord);
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::R12
//  
/////////////////////////////////////////////////////////////////////////////

const Complex Stack::R12(int i, int j) const
{
  MultiScatterer* multi = as_multi();

  if (multi)
  {
    Complex R12ij = (multi->get_R12())(i,j);

    BlochStack* bs_inc = dynamic_cast<BlochStack*>(multi->get_inc());

    if (bs_inc) // Normalise.
    { 
      BlochMode* fm = dynamic_cast<BlochMode*>(bs_inc->get_fw_mode(j));
      BlochMode* bm = dynamic_cast<BlochMode*>(bs_inc->get_bw_mode(i));

      Real w = bs_inc->c1_size().real();

      Complex fflux = fm->S_flux(0,w,1e-10);
      Complex bflux = bm->S_flux(0,w,1e-10);

      R12ij *= sqrt(bflux/fflux);    
    }

    return R12ij;
  }
  else
    return as_mono()->get_R12();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::R21
//  
/////////////////////////////////////////////////////////////////////////////

const Complex Stack::R21(int i, int j) const
{
  MultiScatterer* multi = as_multi();

  if (multi)
  {
    Complex R21ij = (multi->get_R21())(i,j);

    BlochStack* bs_ext = dynamic_cast<BlochStack*>(multi->get_ext());
    
    if (bs_ext) // Normalise.
    { 
      BlochMode* fm = dynamic_cast<BlochMode*>(bs_ext->get_fw_mode(i));
      BlochMode* bm = dynamic_cast<BlochMode*>(bs_ext->get_bw_mode(j));
      
      Real w = bs_ext->c1_size().real();

      Complex fflux = fm->S_flux(0,w,1e-10);
      Complex bflux = bm->S_flux(0,w,1e-10);

      R21ij *= sqrt(fflux/bflux);
    }

    return R21ij;
  } 
  else
    return as_mono()->get_R21();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::T12
//  
/////////////////////////////////////////////////////////////////////////////

const Complex Stack::T12(int i, int j) const
{
  MultiScatterer* multi = as_multi();

  if (multi)
  {
    Complex T12ij = (multi->get_T12())(i,j);

    BlochStack* bs_inc = dynamic_cast<BlochStack*>(multi->get_inc());
    BlochStack* bs_ext = dynamic_cast<BlochStack*>(multi->get_ext());

    Complex fflux1, fflux2;

    if (bs_inc)
    {
      Real w = bs_inc->c1_size().real();
      BlochMode* fm  = dynamic_cast<BlochMode*>(bs_inc->get_fw_mode(j));
      fflux1 = fm->S_flux(0,w,1e-10);

      T12ij /= sqrt(fflux1);
    }

    if (bs_ext)
    {
      if (bs_ext==bs_inc && i==j)
        fflux2 = fflux1;
      else
      {
        Real w = bs_ext->c1_size().real();
        BlochMode* fm = dynamic_cast<BlochMode*>(bs_ext->get_fw_mode(i));
        fflux2 = fm->S_flux(0,w,1e-10);
       }

       T12ij *= sqrt(fflux2);
    }
    return T12ij;
  }
  else
    return as_mono()->get_T12();
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::T21
//  
/////////////////////////////////////////////////////////////////////////////

const Complex Stack::T21(int i, int j) const
{
  MultiScatterer* multi = as_multi();

  if (multi)
  {
    Complex T21ij = (multi->get_T21())(i,j);

    BlochStack* bs_inc = dynamic_cast<BlochStack*>(multi->get_inc());
    BlochStack* bs_ext = dynamic_cast<BlochStack*>(multi->get_ext());

    Complex bflux1, bflux2;

    if (bs_inc)
    {
      Real w = bs_inc->c1_size().real();
      BlochMode* bm = dynamic_cast<BlochMode*>(bs_inc->get_bw_mode(i));

      bflux2 = bm->S_flux(0,w,1e-10);

      T21ij *= sqrt(bflux2);
    }

    if (bs_ext)
    {
      if (bs_ext==bs_inc && i==j)
        bflux1 = bflux2;
      else
      {
        Real w = bs_ext->c1_size().real();
        BlochMode* bm = dynamic_cast<BlochMode*>(bs_ext->get_bw_mode(j));
        bflux1 = bm->S_flux(0,w,1e-10);
       }
      
       T21ij /= sqrt(bflux1);
    }
    return T21ij;
  }
  else
    return as_mono()->get_T21();
}



/////////////////////////////////////////////////////////////////////////////
//
// In the case of a stack starting and/or ending with a semi-infinite 
// BlochStack, the following methods get the normalised R/T elements, 
// as opposed to DenseScatterer->get_Rxx/Txx.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_R12
//  
/////////////////////////////////////////////////////////////////////////////

const cMatrix Stack::get_R12() const
{
  BlochStack* bs_inc = dynamic_cast<BlochStack*>(as_multi()->get_inc());

  if (!bs_inc)
    return as_multi()->get_R12();
  
  cMatrix R12p(global.N,global.N,fortranArray);

  for (int i=1; i<=global.N; i++)
    for (int j=1; j<=global.N; j++)
      (R12p)(i,j) = R12(i,j);

  return R12p;	
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_R21
//  
/////////////////////////////////////////////////////////////////////////////

const cMatrix Stack::get_R21() const
{
  BlochStack* bs_ext = dynamic_cast<BlochStack*>(as_multi()->get_ext());

  if (!bs_ext)
    return as_multi()->get_R21();

  cMatrix R21p(global.N,global.N,fortranArray);

  for (int i=1; i<=global.N; i++)
    for (int j=1; j<=global.N; j++)
      (R21p)(i,j) = R21(i,j);

   return R21p;	
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_T12
//  
/////////////////////////////////////////////////////////////////////////////

const cMatrix Stack::get_T12() const
{
  BlochStack* bs_inc = dynamic_cast<BlochStack*>(as_multi()->get_inc());
  BlochStack* bs_ext = dynamic_cast<BlochStack*>(as_multi()->get_ext());

  if (!bs_inc && !bs_ext)
    return as_multi()->get_T12();  

  cMatrix T12p(global.N,global.N,fortranArray);

  for (int i=1; i<=global.N; i++)
    for (int j=1; j<=global.N; j++)
      (T12p)(i,j) = T12(i,j);

   return T12p;	
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_T21
//  
/////////////////////////////////////////////////////////////////////////////

const cMatrix Stack::get_T21() const
{
  BlochStack* bs_inc = dynamic_cast<BlochStack*>(as_multi()->get_inc());
  BlochStack* bs_ext = dynamic_cast<BlochStack*>(as_multi()->get_ext());

  if (!bs_inc && !bs_ext)
    return as_multi()->get_T21();  

  cMatrix T21p(global.N,global.N,fortranArray);

  for (int i=1; i<=global.N; i++)
    for (int j=1; j<=global.N; j++)
      (T21p)(i,j) = T21(i,j);

  return T21p;	
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_R12_power
//
//  The main reasons these functions are faster than their Python 
//  equivalents is because they avoid calls to convert_to_dense() which 
//  are expensive when always_recalculate is true.
//  
/////////////////////////////////////////////////////////////////////////////

const cMatrix Stack::get_R12_power() const
{
  MultiScatterer* multi = as_multi();

  if (! dynamic_cast<BlochSection*>(multi->get_inc()))
  {
    py_error("Power functions only implemented for BlochSections.");
    exit(-1);
  }

  cMatrix R(global.N,global.N,fortranArray);

  R = blitz::pow2(blitz::abs(multi->get_R12()));

  for (int i=1; i<=global.N; i++)
    for (int j=1; j<=global.N; j++)
    {
      Complex P_in  = multi->get_inc()->get_mode(j)->field(Coord(0,0,0)).Sz();
      Complex P_out = multi->get_inc()->get_mode(i)->field(Coord(0,0,0)).Sz();

      if (abs(real(P_in)) > 1e-10)
        R(i,j) *= (P_out) / (P_in);
      else
        R(i,j) = 0.0;
    }

  return R;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::get_T12_power
//  
/////////////////////////////////////////////////////////////////////////////

const cMatrix Stack::get_T12_power() const
{
  MultiScatterer* multi = as_multi();

  if (! dynamic_cast<BlochSection*>(multi->get_inc()))
  {
    py_error("Power functions only implemented for BlochSections.");
    exit(-1);
  }

  cMatrix T(global.N,global.N,fortranArray);

  T = blitz::pow2(blitz::abs(multi->get_T12()));

  for (int i=1; i<=global.N; i++)
    for (int j=1; j<=global.N; j++)
    {
      Complex P_in  = multi->get_inc()->get_mode(j)->field(Coord(0,0,0)).Sz();
      Complex P_out = multi->get_ext()->get_mode(i)->field(Coord(0,0,0)).Sz();

      if (abs(real(P_in)) > 1e-10)
        T(i,j) *= (P_out) / (P_in);
      else
        T(i,j) = 0.0;
    }

  return T;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::create_sc
//  
/////////////////////////////////////////////////////////////////////////////

Scatterer* Stack::create_sc(const Expression& e, unsigned int no_of_periods)
{  
  Scatterer* sc;
  
  if (e.all_layers_uniform())
  { 
    if (e.is_mono())
      sc = new MonoStack(e, no_of_periods);    
    else
      sc = new DiagStack(e, no_of_periods);
  }
  else
    sc = new DenseStack(e, no_of_periods);

  return sc;
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::calc_interface_positions
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::calc_interface_positions()
{
  interface_positions.clear();
  
  vector<Complex*> thicknesses
    = dynamic_cast<const StackImpl*>(flat_sc)->get_thicknesses();
  
  Complex current_z = 0.0; 
  for (unsigned int i=0; i<thicknesses.size(); i++)
  {
    current_z += *(thicknesses[i]);
    interface_positions.push_back(current_z);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Stack::calc_interface_fields
//  
/////////////////////////////////////////////////////////////////////////////

void Stack::calc_interface_fields()
{ 
  // Build incident field if necessary.

  if (interface_field.size() == 0)
  {
    if (inc_field.rows() == 0)
    {
      py_error("Error: incident field not set.");
      return;
    }

    calcRT();

    cVector left_bw(inc_field.rows(), fortranArray);

    if (as_multi())
    {
      left_bw.reference(multiply(as_multi()->get_R12(), inc_field));

      if (bw_inc)
        left_bw += multiply(as_multi()->get_T21(), inc_field_bw);
    }
    else
    {
      left_bw(1) = R12(1,1) * inc_field(1);

      if (bw_inc)
        left_bw(1) += T21(1,1) * inc_field_bw(1);
    }

    FieldExpansion inc(get_inc(), inc_field, left_bw);
    interface_field.push_back(inc);
  }

  if (interface_field.size() > 1)
    return;
  
  // Calculate field at each interface.

  const vector<Chunk>* chunks
    = dynamic_cast<StackImpl*>(flat_sc)->get_chunks();

  return S_scheme_fields_S(*chunks, &interface_field, &inc_field_bw);
}
