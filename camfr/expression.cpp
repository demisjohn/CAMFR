
/////////////////////////////////////////////////////////////////////////////
//
// File:     expression.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000511
// Version:  1.2
//
// Copyright (C) 1999-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "expression.h"
#include "stack.h"
#include "interface.h"
#include "primitives/slab/slabmatrixcache.h"

/////////////////////////////////////////////////////////////////////////////
//
// Static variables.
//
/////////////////////////////////////////////////////////////////////////////

TmpStorage<Stack>      Expression::tmp_stacks;
TmpStorage<Expression> Expression::tmp_exprs;



/////////////////////////////////////////////////////////////////////////////
//
// free_tmps
//
/////////////////////////////////////////////////////////////////////////////

void free_tmps()
{
  Expression::tmp_stacks.clear();
  Expression::tmp_exprs.clear();
  interface_cache.clear();
  slabmatrix_cache.clear();
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::Expression
//
/////////////////////////////////////////////////////////////////////////////

Expression::Expression(const Term& t_)
{
  Term* t = new Term(t_);
  terms.push_back(t);
  
  transparent_dummy = NULL;
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::Expression
//
/////////////////////////////////////////////////////////////////////////////

Expression::Expression(const Expression& e)
{ 
  for (unsigned int i=0; i<e.terms.size(); i++)
  {
    Term* t = new Term(*(e.terms[i]));
    terms.push_back(t);
  }

  transparent_dummy = e.transparent_dummy;
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::~Expression
//
/////////////////////////////////////////////////////////////////////////////

Expression::~Expression()
{
  for (unsigned int i=0; i<terms.size(); i++)
    delete terms[i];
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::add_term
//
/////////////////////////////////////////////////////////////////////////////

void Expression::add_term(const Term& t_) const
{
  Term* t = new Term(t_);
  terms.push_back(t);
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::insert_term_front
//
/////////////////////////////////////////////////////////////////////////////

void Expression::insert_term_front(const Term& t_) const
{
  Term* t = new Term(t_);
  terms.insert(terms.begin(), t);
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::operator+=
//
/////////////////////////////////////////////////////////////////////////////

void Expression::operator+=(const Term& t)
{
  if ( (t.get_type() == MATERIAL) || (t.get_type() == MAT_EXPRESSION) )
  {
    add_term(t);
    return;
  }

  // For waveguides, insert interfaces as needed.
  
  if (terms.size() == 0)
  {
    // Special case of propagation in incidence medium.
    
    if ( (t.get_type() == WAVEGUIDE) && (t.get_d() != 0.0) )
    {
      Scatterer* sc = interface_cache.get_interface(t.get_wg(), t.get_wg());
      add_term(Term(*sc));
      set_transparent_dummy(sc);
    }
  }
  else
  {
    // Interface needed?

    Waveguide* wg1 = get_ext();
    Waveguide* wg2 = t.get_inc();

    if (wg1 != wg2)
    {
      Scatterer* sc = interface_cache.get_interface(wg1, wg2);
      add_term(Term(*sc));
    }
  }
  
  add_term(t);
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::operator=
//
/////////////////////////////////////////////////////////////////////////////

const Expression& Expression::operator=(const Expression& e)
{
  if (this == &e)
    return e;

  for (unsigned int i=0; i<terms.size(); i++)
    delete terms[i];

  terms.clear();

  for (unsigned int i=0; i<e.terms.size(); i++)
  {
    Term* t = new Term(*(e.terms[i]));
    terms.push_back(t);
  }

  transparent_dummy = e.transparent_dummy;
  
  return *this;
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::flatten
//
//   Returns a single level view of the expression, i.e. without
//   subexpressions and without explicit periodicity.
//
/////////////////////////////////////////////////////////////////////////////

Expression Expression::flatten() const
{  
  // Pass 1: eliminate recursion and periodicity.
  
  Expression flat;
  
  for (unsigned int i=0; i<terms.size(); i++)
  {
    if (    (terms[i]->get_type() == STACK_EXPRESSION)
         || (terms[i]->get_type() == MAT_EXPRESSION)   )
    {
      Expression flat_subexpression
        = terms[i]->get_expression()->flatten();

      for (unsigned int j=0; j<terms[i]->get_no_of_periods(); j++)
        for (unsigned int k=0; k<flat_subexpression.get_size(); k++)
          flat.add_term(*(flat_subexpression.get_term(k)));
    }
    else
      flat.add_term(*(terms[i]));
  }

  // In case of Material expression, don't optimise any further.

  if (    (terms[0]->get_type() == MATERIAL)
       || (terms[0]->get_type() == MAT_EXPRESSION) )
    return flat;

  // Pass 2: optimise wg(d1) + interface(wg,x) + interface(x,wg) + wg(d2)
  // to wg(d1+d2).

  Expression flat_optim;
  for (unsigned int i=0; i<flat.get_size(); i++)
  {
    Term* t = flat.get_term(i);
    if ((t->get_type() == WAVEGUIDE) && (i+3<flat.get_size()))
    { 
      Term* t1 = flat.get_term(i+1);
      Term* t2 = flat.get_term(i+2);          
      Term* t3 = flat.get_term(i+3);
      
      if (t1->is_interface() && t2->is_interface()
          && (t3->get_type()   == WAVEGUIDE)
          && (*( t->get_wg())  == *(t1->get_inc()))
          && (*(t1->get_ext()) == *(t2->get_inc()))
          && (*(t2->get_ext()) == *(t3->get_wg()))
          && (*( t->get_wg())  == *(t3->get_wg())))
      {
        Term optim_term( (*(t->get_wg()))(t->get_d() + t3->get_d()) );
        flat_optim.add_term(optim_term);
        
        i += 3;
      }
      else
        flat_optim.add_term(*t);
    }
    else
      flat_optim.add_term(*t);
  }

  // Corner case if optimisation exposed a waveguide at the front of
  // the expression, without a scatterer.

  Term* t0 = flat_optim.get_term(0);
  if (t0->get_type() == WAVEGUIDE)
  { 
    Scatterer* sc = interface_cache.get_interface(t0->get_wg(), t0->get_wg());
    flat_optim.insert_term_front(Term(*sc));
  }
  
  return flat_optim;  
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::all_layers_uniform
//
/////////////////////////////////////////////////////////////////////////////

bool Expression::all_layers_uniform() const
{
  for (unsigned int i=0; i<terms.size(); i++)
    if (terms[i]->all_layers_uniform() == false)
      return false;
  
  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::no_gain_present
//
/////////////////////////////////////////////////////////////////////////////

bool Expression::no_gain_present() const
{
  for (unsigned int i=0; i<terms.size(); i++)
    if (terms[i]->no_gain_present() == false)
      return false;
  
  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::is_mono
//
/////////////////////////////////////////////////////////////////////////////

bool Expression::is_mono() const
{
  for (unsigned int i=0; i<terms.size(); i++)
    if (terms[i]->is_mono() == false)
      return false;
  
  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::get_inc
//
/////////////////////////////////////////////////////////////////////////////

Waveguide* Expression::get_inc() const
{
  return terms[0]->get_inc();
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::get_ext
//
/////////////////////////////////////////////////////////////////////////////

Waveguide* Expression::get_ext() const
{
  return terms[terms.size()-1]->get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// Expression::repr
//
/////////////////////////////////////////////////////////////////////////////

string Expression::repr() const
{
  ostringstream s;
  
  for (unsigned int i=0; i<get_size(); i++)
  {
    s << get_term(i)->repr();
    
    if (i != get_size()-1)
      s << endl;
  }
  
  return s.str();
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::Term(Material_length)
//
/////////////////////////////////////////////////////////////////////////////

Term::Term(const Material_length& m_l)
  : type(MATERIAL), mat(m_l.mat), wg(NULL), d(m_l.d),
    sc(NULL), ex(NULL), N(1), st(NULL) {}



/////////////////////////////////////////////////////////////////////////////
//
// Term::Term(Waveguide_length)
//
/////////////////////////////////////////////////////////////////////////////

Term::Term(const Waveguide_length& wg_l)
  : type(WAVEGUIDE), mat(NULL), wg(wg_l.wg), d(wg_l.d),
    sc(NULL), ex(NULL), N(1), st(NULL) {}



/////////////////////////////////////////////////////////////////////////////
//
// Term::Term(Scatterer)
//
/////////////////////////////////////////////////////////////////////////////

Term::Term(Scatterer& s)
  : type(SCATTERER), mat(NULL), wg(NULL), d(0.0),
    sc(&s), ex(NULL), N(1), st(NULL) {}



/////////////////////////////////////////////////////////////////////////////
//
// Term::Term(Stack)
//
/////////////////////////////////////////////////////////////////////////////

Term::Term(Stack& st_)
  : type(SCATTERER), mat(NULL), wg(NULL), d(0.0), sc(st_.get_sc()),
    N(st_.get_no_of_periods()), st(&st_)
{
  ex = new Expression(st_.get_expression());
  Expression::tmp_exprs.store(ex);
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::Term(Expression)
//
/////////////////////////////////////////////////////////////////////////////

Term::Term(const Expression& e, unsigned int N_=1)
  : mat(NULL), wg(NULL), d(0.0), N(N_)
{
  ex = new Expression(e);
  Expression::tmp_exprs.store(ex);
    
  if (    (e.get_term(0)->get_type() == MATERIAL)
       || (e.get_term(0)->get_type() == MAT_EXPRESSION) )
  {
    type = MAT_EXPRESSION;

    sc = NULL;
    st = NULL;
  }
  else
  {
    type = STACK_EXPRESSION;
    
    st = new Stack(e, N);
    sc = st->get_sc();
    Expression::tmp_stacks.store(st);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::is_interface
//
/////////////////////////////////////////////////////////////////////////////

bool Term::is_interface() const
{
  if (type == SCATTERER)
    if (    dynamic_cast<DenseInterface*>(sc)
         || dynamic_cast< DiagInterface*>(sc)
         || dynamic_cast< MonoInterface*>(sc) )
      return true;
  
  return false;
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::all_layers_uniform
//
/////////////////////////////////////////////////////////////////////////////

bool Term::all_layers_uniform() const
{
  switch (type)
  {
    case (MATERIAL):
    case (MAT_EXPRESSION):
      return true;
      
    case (WAVEGUIDE):
      return wg->is_uniform();

    case (SCATTERER):
      return sc->all_layers_uniform();

    case (STACK_EXPRESSION):
      return st->all_layers_uniform();

  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::no_gain_present
//
/////////////////////////////////////////////////////////////////////////////

bool Term::no_gain_present() const
{
  switch (type)
  {
    case (MATERIAL):
      return mat->no_gain_present();
      
    case (WAVEGUIDE):
      return ( (!wg->no_gain_present()) || (imag(d) > 1e-3) )
        ? false : true;

    case (SCATTERER):
      return sc->no_gain_present();
 
    case (MAT_EXPRESSION):
      return get_expression()->no_gain_present();

    case (STACK_EXPRESSION):
      return st->no_gain_present();
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::is_mono
//
/////////////////////////////////////////////////////////////////////////////

bool Term::is_mono() const
{
  switch (type)
  {
    case (MATERIAL):
    case (MAT_EXPRESSION):
      return true;
      
    case (WAVEGUIDE):
      return dynamic_cast<MonoWaveguide*>(wg) ? true : false;

    case (SCATTERER):
      return sc->is_mono();
      
    case (STACK_EXPRESSION):
      return st->is_mono();
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::get_inc
//
/////////////////////////////////////////////////////////////////////////////

Waveguide* Term::get_inc() const
{
  switch (type)
  {
    case (MATERIAL):
    case (MAT_EXPRESSION):
      return NULL;
    
    case (WAVEGUIDE):
      return wg;

    case (SCATTERER):
      return sc->get_inc();

    case (STACK_EXPRESSION):
      return st->get_inc();
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::get_ext
//
/////////////////////////////////////////////////////////////////////////////

Waveguide* Term::get_ext() const
{
  switch (type)
  {
    case (MATERIAL):
    case (MAT_EXPRESSION):
      return NULL;
      
    case (WAVEGUIDE):
      return wg;

    case (SCATTERER):
      return sc->get_ext();

    case (STACK_EXPRESSION):
      return st->get_ext();
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Term::repr
//
/////////////////////////////////////////////////////////////////////////////

string Term::repr() const
{
  ostringstream s;
  
  switch (type)
  {
    case (MATERIAL):
      s << "Mat " << d << " " << mat->repr();
      break;
      
    case (WAVEGUIDE):
      s << "WG " << d << " " << wg << ":" << wg->get_core()->n();
      break;

    case (SCATTERER):
      s << "SC "
        << sc->get_inc() << ":" << sc->get_inc()->get_core()->n()
        << " " << sc << " "
        << sc->get_ext() << ":" << sc->get_ext()->get_core()->n();     
      break;
      
    case (STACK_EXPRESSION):
    case (MAT_EXPRESSION):
      s << "EX " << N << " * " << endl << "(" << endl << *ex << ")";
  }

  return s.str();
}



/////////////////////////////////////////////////////////////////////////////
//
// operators+(Term, Term)
//
//   Starts the expression parsing.
//
/////////////////////////////////////////////////////////////////////////////

const Expression operator+(const Term& L, const Term& R)
{ 
  Expression e;

  // Special case of propagation in incidence medium.
  
  if ( (L.get_type() == WAVEGUIDE) && (L.get_d() != 0.0) )
  {
    Scatterer* sc = interface_cache.get_interface(L.get_wg(), L.get_wg());
    e.add_term(Term(*sc));
    e.set_transparent_dummy(sc);
  }

  e.add_term(L);

  // Interface needed between L and R?
  
  Waveguide* wg1 = L.get_ext();
  Waveguide* wg2 = R.get_inc();
  
  if (wg1 != wg2)
  {
    Scatterer* sc = interface_cache.get_interface(wg1, wg2);
    e.add_term(Term(*sc));
  }
  
  e.add_term(R);
  
  return e;
}



/////////////////////////////////////////////////////////////////////////////
//
// operators+(Expression, Term)
//
//   Continue expression parsing, adding interfaces along the way if needed.
//
/////////////////////////////////////////////////////////////////////////////

const Expression& operator+(const Expression& L, const Term& R)
{ 
  Waveguide* wg1 = L.get_ext();
  Waveguide* wg2 = R.get_inc();
  
  if (wg1 != wg2)
  {
    Scatterer* sc = interface_cache.get_interface(wg1, wg2);
    L.add_term(Term(*sc));
  }

  L.add_term(R);
  
  return L;
}



/////////////////////////////////////////////////////////////////////////////
//
// operator*(N, Expression)
//
/////////////////////////////////////////////////////////////////////////////

const Term operator*(unsigned int N, const Expression& e)
{
  if (N == 0)
  {
    cerr << "Error: zero is an invalid number of periods." << endl;
    exit (-1);
  }

  Expression new_e;

  // Add extra interface add the beginning (uses less stack iterations
  // than when doing it at the end).

  if (    (e.get_term(0)->get_type() != MATERIAL)
       && (e.get_term(0)->get_type() != MAT_EXPRESSION) )
  {
    Scatterer* sc = interface_cache.get_interface(e.get_ext(), e.get_inc());
    new_e.add_term(Term(*sc));
  }
  
  // Add other terms, but remove transparent dummy scatterer in beginning.

  for (unsigned int i=0; i<e.get_size(); i++)
  {
    const Term* t = e.get_term(i);

    if (!(    (t->get_type() == SCATTERER)
           && (t->get_sc()   == e.get_transparent_dummy()) ))
      new_e.add_term(*t);
  }

  return Term(new_e, N);
}



/////////////////////////////////////////////////////////////////////////////
//
// material_expression_to_table
//
/////////////////////////////////////////////////////////////////////////////

void material_expression_to_table(const Expression& e, 
                                  vector<Complex>* eps, 
                                  vector<Complex>* mu,
                                  vector<Complex>* d)
{
  Expression ex = e.flatten();
  
  for (unsigned int i=0; i<ex.get_size(); i++)
  {
    Material* m = dynamic_cast<Material*>(ex.get_term(i)->get_mat());

    if (!m)
    {
      cerr << "Error: expression contains non-material term." << endl;
      exit (-1);
    }

    Complex thickness = ex.get_term(i)->get_d();

    // Combine two succesive terms containing the same material.

    if ( (i+1 < ex.get_size()) && (m == ex.get_term(i+1)->get_mat()) )
    {
      thickness += ex.get_term(i+1)->get_d();
      i++;
    }

    if (abs(thickness))
    {
      eps->push_back(m->eps());
       mu->push_back(m->mu());
        d->push_back(thickness);
    }
  }
}
