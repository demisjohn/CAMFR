
/////////////////////////////////////////////////////////////////////////////
//
// File:     scatterer.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000511
// Version:  1.3
//
// Copyright (C) 1998-2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SCATTERER_H
#define SCATTERER_H

#include "math/linalg/linalg.h"
#include "waveguide.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Scatterer
//
//   Abstract black box scatterer from an incidence waveguide to an exit
//   waveguide.
//  
/////////////////////////////////////////////////////////////////////////////

class Scatterer
{
  public:

    Scatterer() : inc(NULL), ext(NULL) {}
    Scatterer(Waveguide& inc_, Waveguide& ext_) : inc(&inc_), ext(&ext_) {}
    virtual ~Scatterer() {}

    Waveguide* get_inc() const {return inc;}
    Waveguide* get_ext() const {return ext;}

    virtual std::vector<Material*> get_materials() 
      const {std::vector<Material*> m; return m;}

    virtual Complex get_total_thickness()     const {return 0.0;}
    virtual bool contains(const Material& m)  const {return false;}
    virtual bool no_gain_present()            const {return true;}
    virtual bool recalc_needed()              const {return true;}
    
    virtual bool all_layers_uniform() const {return false;}
    virtual bool is_mono() const {return false;}
    
    virtual void  calcRT() {}
    virtual void allocRT() {}
    virtual void  freeRT() {}
     
  protected:
    
    Waveguide* inc;
    Waveguide* ext;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: MultiScatterer
//
//   Scatterer where multiple modes are considered, therefore having
//   (dense or diagonal) R and T matrices.
//
//   Note that it is not yet specified how the storage should be handled,
//   only that a cMatrix can be returned.
//  
/////////////////////////////////////////////////////////////////////////////

class MultiScatterer : public Scatterer
{
  public:

    MultiScatterer()
      : last_lambda(0.0), last_gain_mat(Material(0.0)) {}
    MultiScatterer(Waveguide& inc, Waveguide& ext)
      : Scatterer(inc, ext), last_lambda(0.0), last_gain_mat(Material(0.0)) {}

    bool recalc_needed() const;

    virtual const cMatrix& get_R12() const = 0;
    virtual const cMatrix& get_R21() const = 0;
    virtual const cMatrix& get_T12() const = 0;
    virtual const cMatrix& get_T21() const = 0;
     
  protected:

    // The wavelength and gain the matrices were last calculated for,
    // are used to determine if recalculation is needed.

    Complex  last_lambda;
    Material last_gain_mat;
    Complex  last_slab_ky;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: DenseScatterer
//
//   MultiScatterer with dense R and T matrices.
//  
/////////////////////////////////////////////////////////////////////////////

class DenseScatterer : public MultiScatterer
{
  public:

    DenseScatterer();
    DenseScatterer(Waveguide& inc, Waveguide& ext);
    ~DenseScatterer() {freeRT();}

    void allocRT();
    void  freeRT();

    const cMatrix& get_R12() const {return R12;}
    const cMatrix& get_R21() const {return R21;}
    const cMatrix& get_T12() const {return T12;}
    const cMatrix& get_T21() const {return T21;}

    void set_R12(const cMatrix& M) {R12.reference(M);}
    void set_R21(const cMatrix& M) {R21.reference(M);}
    void set_T12(const cMatrix& M) {T12.reference(M);}
    void set_T21(const cMatrix& M) {T21.reference(M);}

    void copy_R12(const cMatrix& M) {R12=M;}
    void copy_R21(const cMatrix& M) {R21=M;}
    void copy_T12(const cMatrix& M) {T12=M;}
    void copy_T21(const cMatrix& M) {T21=M;}
    
    void copy_RT_from(const DenseScatterer& sc);
    void swap_RT_with(      DenseScatterer& sc);
    
  protected:

    cMatrix R12, R21, T12, T21;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: DiagScatterer
//
//   MultiScatterer with diagonal R and T matrices.
//   Implies uniform incidence and exit media.
//   Internal storage is a cVector, but a cMatrix can be returned for
//   interoperability with DenseScatterers
//  
/////////////////////////////////////////////////////////////////////////////

class DiagScatterer : public MultiScatterer
{
  public:

    DiagScatterer();
    DiagScatterer(Waveguide& inc, Waveguide& ext);
    DiagScatterer(const DiagScatterer& sc_d);
    DiagScatterer& operator=(const DiagScatterer& sc_d);
    ~DiagScatterer() {freeRT();}

    void allocRT();
    void  freeRT();

    const cMatrix& get_R12() const {convert_to_dense(); return *R12_dense;}
    const cMatrix& get_R21() const {convert_to_dense(); return *R21_dense;}
    const cMatrix& get_T12() const {convert_to_dense(); return *T12_dense;}
    const cMatrix& get_T21() const {convert_to_dense(); return *T21_dense;}

    const cVector& get_diag_R12() const {return R12;}
    const cVector& get_diag_R21() const {return R21;}
    const cVector& get_diag_T12() const {return T12;}
    const cVector& get_diag_T21() const {return T21;}

    void set_diag_R12(const cVector& V) {R12.reference(V);}
    void set_diag_R21(const cVector& V) {R21.reference(V);}
    void set_diag_T12(const cVector& V) {T12.reference(V);}
    void set_diag_T21(const cVector& V) {T21.reference(V);}

    void copy_diag_R12(const cVector& V) {R12=V;}
    void copy_diag_R21(const cVector& V) {R21=V;}
    void copy_diag_T12(const cVector& V) {T12=V;}
    void copy_diag_T21(const cVector& V) {T21=V;}

    void copy_RT_from(const DiagScatterer& sc_d);
    void swap_RT_with(      DiagScatterer& sc_d);
       
  protected:

    cVector R12, R21, T12, T21;

    void convert_to_dense() const;

    mutable cMatrix* R12_dense;
    mutable cMatrix* R21_dense;
    mutable cMatrix* T12_dense;
    mutable cMatrix* T21_dense;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: MonoScatterer
//
//   Scatterer where only a single mode is considered, therefore having
//   'single' element R and T matrices.
//   Implies uniform infinite media, where all modes are decoupled
//   and a single mode is chosen from a continuous set.
//  
/////////////////////////////////////////////////////////////////////////////

class MonoScatterer : public Scatterer
{
  public:

    MonoScatterer() : Scatterer() {}
    MonoScatterer(Waveguide& inc, Waveguide& ext);
    ~MonoScatterer() {}
    
    // We always recalculate MonoScatterers, since the overhead of
    // checking for a needed recalc outweighs the possible gains.

    bool recalc_needed()      const {return true;}
    bool all_layers_uniform() const {return true;} 
    bool is_mono()            const {return true;}

    const Complex& get_R12() const {return R12;}
    const Complex& get_R21() const {return R21;}
    const Complex& get_T12() const {return T12;}
    const Complex& get_T21() const {return T21;}

    void set_R12(const Complex& C) {R12=C;}
    void set_R21(const Complex& C) {R21=C;}
    void set_T12(const Complex& C) {T12=C;}
    void set_T21(const Complex& C) {T21=C;}

    void copy_RT_from(const MonoScatterer& sc_m);
    void swap_RT_with(      MonoScatterer& sc_m);
       
  protected:
    
    Complex R12, R21, T12, T21;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: FlippedScatterer
//
//   A flipped version of an existing MultiScatterer, reusing its
//   calculated matrices.
//  
/////////////////////////////////////////////////////////////////////////////

class FlippedScatterer : public MultiScatterer
{
  public:
  
    FlippedScatterer(MultiScatterer& sc_)
      : MultiScatterer(*sc_.get_ext(), *sc_.get_inc()), sc(&sc_) {}

    Complex get_total_thickness() const {return sc->get_total_thickness();}
    std::vector<Material*> get_materials() const {return sc->get_materials();}
    bool contains(const Material& m) const {return sc->contains(m);}

    void allocRT() {return sc->allocRT();}
    void  freeRT() {return sc->freeRT();}
    void  calcRT() {return sc->calcRT();}

    const cMatrix& get_R12() const {return sc->get_R21();}  
    const cMatrix& get_R21() const {return sc->get_R12();} 
    const cMatrix& get_T12() const {return sc->get_T21();} 
    const cMatrix& get_T21() const {return sc->get_T12();} 
    
  protected:

    MultiScatterer* sc;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: TransparentScatterer
//
//   Transparent DiagScatterer, i.e. artificial interface between two
//   waveguides that are the same.
//  
/////////////////////////////////////////////////////////////////////////////

class TransparentScatterer : public DiagScatterer
{
  public:

    TransparentScatterer(Waveguide& w) : DiagScatterer(w, w) {}

    void calcRT();

    Complex get_total_thickness()          const {return 0.0;}
    std::vector<Material*> get_materials() const {return inc->get_materials();}
    bool contains(const Material& m)       const {return inc->contains(m);}
    bool recalc_needed()                   const {return false;}
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: E_Wall
//
//   Perfectly conducting electric wall.
//  
/////////////////////////////////////////////////////////////////////////////

class E_Wall : public DiagScatterer
{
  public:

    E_Wall(Waveguide& w) : DiagScatterer(w, w) {}
      
    void calcRT();

    Complex get_total_thickness()          const {return 0.0;}   
    std::vector<Material*> get_materials() const {return inc->get_materials();}
    bool contains(const Material& m)       const {return inc->contains(m);}
    bool recalc_needed()                   const {return false;}
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: H_Wall
//
//   Perfectly conducting magnetic wall.
//  
/////////////////////////////////////////////////////////////////////////////

class H_Wall : public DiagScatterer
{
  public:

    H_Wall(Waveguide& w) : DiagScatterer(w, w) {}

    void calcRT();

    Complex get_total_thickness()          const {return 0.0;}
    std::vector<Material*> get_materials() const {return inc->get_materials();}
    bool contains(const Material& m)       const {return inc->contains(m);}
    bool recalc_needed()                   const {return false;}
};


    
#endif
