
/////////////////////////////////////////////////////////////////////////////
//
// File:     stack.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000121
// Version:  2.1
//
// Copyright (C) 1998-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef STACK_H
#define STACK_H

#include <vector>
#include "field.h"
#include "chunk.h"
#include "S_scheme.h"
#include "expression.h"

/////////////////////////////////////////////////////////////////////////////
//
// Note: the coordinate axes are chosen such that the interfaces between
// the different layers are perpendicular the z-axis.
// The modes propagate in the z-direction.
//  
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: StackImpl
//
//   Class containing the implementation shared by the different types of
//   Stacks (DenseStack, DiagStack, MonoStack).
//  
/////////////////////////////////////////////////////////////////////////////

class StackImpl
{
  public:

    StackImpl() {}
  
    StackImpl(const std::vector<Chunk>& chunks, unsigned int no_of_periods=1);
    StackImpl(const Expression& e,              unsigned int no_of_periods=1);
    virtual ~StackImpl() {}

    Complex stack_get_total_thickness()          const;
    std::vector<Material*> stack_get_materials() const;
    bool stack_contains(const Material& m)       const;
    bool stack_no_gain_present()                 const;
    bool stack_all_layers_uniform()              const;
    
    std::vector<Complex*> get_thicknesses() const;

    const std::vector<Chunk>* get_chunks() const {return &chunks;}
    
    template <class T> friend void stack_calcRT(T* stack);
    
  protected:

    std::vector<Chunk> chunks;
    unsigned int  no_of_periods;

  private:
    
    void create_from(const Expression& e);
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: DenseStack
//
//   A stack of chunks, that can be periodically extended.
//   Assumes incidence and exit media of chunks are matched.
//  
/////////////////////////////////////////////////////////////////////////////

class DenseStack : public DenseScatterer, public StackImpl
{
  public:

    DenseStack() {}
    DenseStack(const std::vector<Chunk>& chunks, unsigned int no_of_periods=1);
    DenseStack(const Expression& e,              unsigned int no_of_periods=1);
    ~DenseStack() {}

    Complex get_total_thickness() const
      {return stack_get_total_thickness();}
    
    std::vector<Material*> get_materials() const 
      {return stack_get_materials();}
    
    bool contains(const Material& m) const
      {return stack_contains(m);}
    
    bool no_gain_present() const
      {return stack_no_gain_present();}

    bool all_layers_uniform() const
      {return stack_all_layers_uniform();}
    
    void calcRT();
    void freeRT(); // Also frees chunk memory.
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: DiagStack
//
//   DiagScatterer version of DenseStack.
//  
/////////////////////////////////////////////////////////////////////////////

class DiagStack : public DiagScatterer,  public StackImpl
{
  public:

    DiagStack() {}
    DiagStack(const std::vector<Chunk>& chunks, unsigned int no_of_periods=1); 
    DiagStack(const Expression& e,              unsigned int no_of_periods=1);
    ~DiagStack() {}

    Complex get_total_thickness() const
      {return stack_get_total_thickness();}
    
    std::vector<Material*> get_materials() const 
      {return stack_get_materials();}
    
    bool contains(const Material& m) const
      {return stack_contains(m);}
    
    bool no_gain_present() const
      {return stack_no_gain_present();}

    bool all_layers_uniform() const
      {return stack_all_layers_uniform();}
    
    void calcRT();
    void freeRT(); // Also frees chunk memory.
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: MonoStack
//
//   MonoScatterer version of DenseStack.
//  
/////////////////////////////////////////////////////////////////////////////

class MonoStack : public MonoScatterer, public StackImpl
{
  public:
  
    MonoStack() {}
    MonoStack(const std::vector<Chunk>& chunks, unsigned int no_of_periods=1); 
    MonoStack(const Expression& e,              unsigned int no_of_periods=1);
    ~MonoStack() {}

    Complex get_total_thickness() const
      {return stack_get_total_thickness();}
    
    std::vector<Material*> get_materials() const 
      {return stack_get_materials();}
    
    bool contains(const Material& m) const
      {return stack_contains(m);}
    
    bool no_gain_present() const
      {return stack_no_gain_present();}
    
    void calcRT();
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Stack
//
//   Encapsulates different kinds of stacks through a pointer.
//   Can be initialised by an expression.
//   
/////////////////////////////////////////////////////////////////////////////

class Stack
{
  public:

    Stack() : sc(NULL), flat_sc(NULL), 
              inc_field(fortranArray), inc_field_bw(fortranArray) {}

    Stack(const Expression& e, unsigned int no_of_periods=1);
    Stack(const Term& t);
    Stack(const Stack& s);
    ~Stack() {delete sc; delete flat_sc;}

    Waveguide* get_inc() const {return sc->get_inc();}
    Waveguide* get_ext() const {return sc->get_ext();}    
    
    void   operator=(const Expression& e);
    void   operator=(const Term& t) {operator=(Expression(t));}
    Stack& operator=(const Stack& s);
    
    bool            is_mono() const {return expression.is_mono();}
    bool all_layers_uniform() const {return expression.all_layers_uniform();}
    bool    single_material() const;

    Complex get_total_thickness() const {return sc->get_total_thickness();}
    
    std::vector<Material*> get_materials() const {return sc->get_materials();}

    bool contains(const Material& m) const {return sc->contains(m);}
    bool no_gain_present()           const {return sc->no_gain_present();}
    bool recalc_needed()             const {return sc->recalc_needed();}
      
    void calcRT();
    void allocRT();
    void freeRT();

    void set_inc_field(const cVector& inc_field, cVector* inc_field_bw=NULL);
    cVector get_inc_field() {return inc_field;}
    cVector get_refl_field();
    cVector get_trans_field();

    FieldExpansion inc_field_expansion();
    FieldExpansion ext_field_expansion();
    
    Field field(const Coord& coord);
    void fw_bw_field(const Coord& coord, cVector* fw, cVector* bw);

    void set_interface_field(const std::vector<FieldExpansion>& field);
    void get_interface_field(      std::vector<FieldExpansion>* field);

    Complex lateral_S_flux(const Complex& c1, std::vector<Complex>* S_k=NULL);

    Complex eps_at(const Coord& coord) const;
    Complex  mu_at(const Coord& coord) const;
    Complex   n_at(const Coord& coord) const
      {return sqrt(eps_at(coord)/eps0);}  

    // Easier interface to get R and T matrix elements.

    const Complex R12(int i, int j) const;
    const Complex R21(int i, int j) const;
    const Complex T12(int i, int j) const;
    const Complex T21(int i, int j) const;
    
    // Low level functions.
    
    Scatterer*   get_sc()            const {return sc;}
    Scatterer*   get_flat_sc()       const {return flat_sc;}
    Expression   get_expression()    const {return expression;}
    unsigned int get_no_of_periods() const {return no_of_periods;}

    MultiScatterer* as_multi() const 
      {return dynamic_cast<MultiScatterer*>(sc);}

    MonoScatterer*   as_mono() const 
      {return dynamic_cast< MonoScatterer*>(sc);}
    
  protected:
    
    Scatterer* sc;

    Expression expression;      // These are kept to allow rebuilding sc
    unsigned int no_of_periods; // in copy ctor and operator=.

    // A flat view of the stack, for calculating field profiles.

    Scatterer* flat_sc;
    std::vector<Complex> interface_positions;

    cVector inc_field, inc_field_bw;
    bool bw_inc;
    
    std::vector<FieldExpansion> interface_field;

  private:

    Scatterer* create_sc(const Expression& e, unsigned int no_of_periods=1);

    void calc_interface_positions();

    void calc_interface_fields();
};



#endif
