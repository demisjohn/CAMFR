
/////////////////////////////////////////////////////////////////////////////
//
// File:     expression.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000124
// Version:  1.1
//
// Copyright (C) 1999-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef EXPRESSION_H
#define EXPRESSION_H

#include <string>
#include <vector>
#include "defs.h"
#include "waveguide.h"
#include "scatterer.h"
#include "icache.h"
#include "util/storage.h"

/////////////////////////////////////////////////////////////////////////////
//
// These following classes implement expressions describing a sequence of
// structures using the grammar:
//
//      expr -> expr + term
//      term -> waveguide | scatterer | N*(expr) | (expr)*N | (expr)
//
// This grammar is implemented by:
//
//      expr -> expr + term : operator+ in expression.h
//      term -> waveguide   : constructor in class Term
//      term -> scatterer   : constructor in class Term
//      term -> N*(expr)    : operator* in expression.h
//      term -> (expr)*N    : operator* in expression.h
//      term -> (expr)      : constructor in class Term
//
// Note: gaas(10) + alas(20) has gaas as inc medium and alas as exit,
// while 10*(gaas(10) + alas(20)) has alas as inc medium and alas as exit. 
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Expression
//
//   Holds a list of pointers to Terms. Pointers need to be used, since
//   Term is forward declared. This implies that for every Term, a copy
//   needs to be constructed to avoid losing data if Term is a temporary
//   object.
//
/////////////////////////////////////////////////////////////////////////////

class Term; // Forward declaration
class Stack; // Forward declaration

class Expression
{
  public:

    Expression() : transparent_dummy(NULL) {}
    Expression(const Term& t);
    Expression(const Expression& e);
    ~Expression();

    // Two low level routines, don't check if interfaces need to be added.
    
    void add_term(const Term& t) const;
    
    void insert_term_front(const Term& t) const;

    // Same as add term, but adds an interface if needed.
    
    void operator+=(const Term& t);

    const Expression& operator=(const Expression& e);
    
    Term* get_term(int i) const
      {return terms[i];}
    
    unsigned int get_size() const
      {return terms.size();}
    
    Scatterer* get_transparent_dummy() const
      {return transparent_dummy;}
    
    void set_transparent_dummy(Scatterer* sc)
      {transparent_dummy=sc;}

    Expression      flatten() const;
    bool all_layers_uniform() const;
    bool    no_gain_present() const;
    bool            is_mono() const;
    Waveguide*      get_inc() const;
    Waveguide*      get_ext() const;
    std::string        repr() const;
    
    static TmpStorage<Stack>      tmp_stacks;
    static TmpStorage<Expression> tmp_exprs;
      
  protected:
    
    mutable std::vector<Term*> terms;
        
    // 'mutable' allows 'add_term' to be declared 'const', which
    // in turn allows operator+ to take Expressions passed as
    // 'const Expression&' instead of 'Expression'.
    // Note: 'Expression&' is not an option since temporaries only
    // bind to const references.

    Scatterer* transparent_dummy;
};

void free_tmps();



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Term
//
//   Note: we use a type field instead of a class hierachy to avoid
//   having to use two succesive type conversions.
//
//   Currently bad design, using the same class to handle both
//   material expressions and waveguide/stack expressions.
//  
/////////////////////////////////////////////////////////////////////////////

typedef enum
  {MATERIAL, WAVEGUIDE, SCATTERER, STACK_EXPRESSION, MAT_EXPRESSION} Termtype;
  
class Term
{
  public:

    Term() {}
    Term(const Material_length& m_l);
    Term(const Waveguide_length& wg_l);
    Term(Scatterer& s);
    Term(Stack& st);
    Term(const Expression& e, unsigned int N=1);
    ~Term() {}

    bool       is_interface() const;
    bool all_layers_uniform() const;
    bool    no_gain_present() const;
    bool            is_mono() const;
    Waveguide*      get_inc() const;
    Waveguide*      get_ext() const;
    std::string        repr() const;

    Termtype              get_type() const {return type;}
    BaseMaterial*          get_mat() const {return mat;}
    Waveguide*              get_wg() const {return wg;}
    Complex                  get_d() const {return d;}
    Scatterer*              get_sc() const {return sc;}
    Expression*     get_expression() const {return ex;}
    unsigned int get_no_of_periods() const {return N;}
    Stack*                  get_st() const {return st;}

  protected:
    
    Termtype type;

    BaseMaterial* mat;// only valid if Termtype == MATERIAL
    Waveguide*    wg; // only valid if Termtype == WAVEGUIDE
    Complex       d;  // only valid if Termtype == MATERIAL | WAVEGUIDE
    Scatterer*    sc; // only valid if Termtype == SCATTERER | STACK_EXPRESSION
    Expression*   ex; // only valid if Termtype == {STACK|MAT}_EXPRESSION
    unsigned int  N;  // only valid if Termtype == {STACK|MAT}_EXPRESSION
    Stack*        st; // only valid if Termtype == STACK_EXPRESSION

    friend class Stack;
};



/////////////////////////////////////////////////////////////////////////////
//
// operators+ implementing expr -> expr + term
//
/////////////////////////////////////////////////////////////////////////////

const Expression operator+(const Term& L,       const Term& R);
const Expression operator+(const Expression& L, const Term& R);



/////////////////////////////////////////////////////////////////////////////
//
// operators* implementing term -> N * (expr)
//                         term -> (expr) * N
//
//   The Term versions short-cut the need for two succesive type conversions.
//
/////////////////////////////////////////////////////////////////////////////

const Term operator*(unsigned int N, const Expression& e);

inline const Term operator*(unsigned int N, const Term& t)
    {return N*Expression(t);}

inline const Term operator*(const Expression& e, unsigned int N)
    {return N*e;}

inline const Term operator*(const Term& t, unsigned int N)
    {return N*Expression(t);}



/////////////////////////////////////////////////////////////////////////////
//
// operator<<
//
/////////////////////////////////////////////////////////////////////////////

inline std::ostream& operator<<(std::ostream& s, const Expression& e)
  {return s << e.repr() << std::endl;}



/////////////////////////////////////////////////////////////////////////////
//
// Convert material expressions to tables with eps, mu and d's.
//
/////////////////////////////////////////////////////////////////////////////

void material_expression_to_table(const Expression& e, 
                                  std::vector<Complex>* eps, 
                                  std::vector<Complex>* mu,
                                  std::vector<Complex>* d);



#endif
