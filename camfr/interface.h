
/////////////////////////////////////////////////////////////////////////////
//
// File:     interface.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991105
// Version:  1.0
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef INTERFACE_H
#define INTERFACE_H

#include "scatterer.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: DenseInterface
//
//   Abstract interface between an incidence and an exit waveguide.
//   Uses generalised mode matching to compute the R and T matrices.
//   Assumes mode are suitably normalised.
//  
/////////////////////////////////////////////////////////////////////////////

class DenseInterface : public DenseScatterer
{
  public:

    DenseInterface(Waveguide& inc, Waveguide& ext)
      : DenseScatterer(inc, ext) {}

    Complex get_total_thickness() const {return 0.0;}

    std::vector<Material*> get_materials() const;
 
    bool contains(const Material& m) const
      {return (inc->contains(m) || ext->contains(m));}

    bool no_gain_present() const
      {return (inc->no_gain_present() && ext->no_gain_present());}

    bool all_layers_uniform() const
      {return (inc->is_uniform() && ext->is_uniform());} 

    void calcRT();

  protected:

    void calcRT_fast();
    void calcRT_safe();
    void calcRT_non_orth_fast();
    void calcRT_non_orth_safe();
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: DiagInterface
//
//   Same, but with uniform incidence and exit media
//   and diagonal R and T matrices.
//
//   The chosen sign convention is such that r_te = r_tm for theta=0.
//  
/////////////////////////////////////////////////////////////////////////////

class DiagInterface : public DiagScatterer
{
  public:

    DiagInterface(Waveguide& inc, Waveguide& ext)
      : DiagScatterer(inc, ext) {}

    Complex get_total_thickness() const {return 0.0;}

    std::vector<Material*> get_materials() const;

    bool contains(const Material& m) const
      {return (inc->contains(m) || ext->contains(m));}

    bool no_gain_present() const
      {return (inc->no_gain_present() && ext->no_gain_present());}

    bool all_layers_uniform() const
      {return (inc->is_uniform() && ext->is_uniform());}  
    
    void calcRT();
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: MonoInterface
//
//   Same as DiagScatterer but with one-element 'matrices'.
//
//   The chosen sign convention is such that r_te = r_tm for theta=0.
//  
/////////////////////////////////////////////////////////////////////////////

class MonoInterface : public MonoScatterer
{
  public:

    MonoInterface(Waveguide& inc, Waveguide& ext)
      : MonoScatterer(inc, ext) {}

    Complex get_total_thickness() const {return 0.0;}

    std::vector<Material*> get_materials() const;

    bool contains(const Material& m) const
      {return (inc->contains(m) || ext->contains(m));}

    bool no_gain_present() const
      {return (inc->no_gain_present() && ext->no_gain_present());}
    
    void calcRT();
};



#endif
    

  
