
/////////////////////////////////////////////////////////////////////////////
//
// File:     waveguide.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991116
// Version:  1.1
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef WAVEGUIDE_H
#define WAVEGUIDE_H

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include "math/linalg/linalg.h"
#include "mode.h"
#include "material.h"
#include "icache.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Waveguide_length
//
//   A piece of waveguide with a certain length
//
/////////////////////////////////////////////////////////////////////////////

class Waveguide; // advance declaration

struct Waveguide_length
{
    Waveguide_length(Waveguide& wg_, const Complex& d_)
      : wg(&wg_), d(d_) {}

    Waveguide_length(Waveguide* wg_, const Complex& d_)
      : wg(wg_),  d(d_) {}

    Waveguide*    wg;
    const Complex d;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Waveguide
//
//   An abstract waveguide containing N modes.
//
//   Note that it is not yet specified how the storage should be handled,
//   only that a Mode* can be returned.
//
/////////////////////////////////////////////////////////////////////////////

class Waveguide
{
  public:

    Waveguide(bool uniform_=false, Material* core_=NULL)
      : uniform(uniform_), core(core_) {}
    
    virtual ~Waveguide() {}

    bool    is_uniform() const {return uniform;}
    Material* get_core() const {return core;}

    virtual Complex c1_size() const = 0;
    
    virtual Complex eps_at(const Coord& coord) const = 0;
    virtual Complex  mu_at(const Coord& coord) const = 0;

    Complex n_at(const Coord& coord) const 
      {return sqrt(eps_at(coord)/eps0);}  
    
    virtual bool  operator==(const Waveguide& w)   const = 0;
    virtual std::vector<Material*> get_materials() const = 0;
    virtual bool  contains(const Material& m)      const = 0;
    virtual bool  no_gain_present()                const = 0;
    virtual bool  recalc_needed()                  const = 0;
    virtual int   N()                              const = 0;
    virtual Mode* get_mode(int i)                  const = 0;
    virtual Mode* get_fw_mode(int i)               const {return 0;};
    virtual Mode* get_bw_mode(int i)               const {return 0;};
    virtual void  find_modes()                           = 0;
    
    const Waveguide_length operator()(const Complex& d=0.0) const;

    virtual bool operator!=(const Waveguide& w) const
      {return !(*this == w);}

    virtual Complex lateral_S_corr(const Coord& c) const
      {return 1.0;}

    virtual std::string repr() const = 0;
    
  protected:

    bool      uniform;
    Material* core;
};

inline std::ostream& operator<<(std::ostream& s, const Waveguide& wg)
  {return s << wg.repr() << std::endl;}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: MultiWaveguide
//
//   A waveguide described by a number of modes, indexed from 1 to N.
//   No range checking is done for performance reasons.
//
//   Pointers to Modes are stored, because Mode is polymorphic.
//
//   Also provides an interface for calculating the overlap matrices
//   O_I_II (O_II_I) with elements (i,j) containing the overlapintegral
//   Int(Ei x Hj) between mode i in medium I (II) and mode j in medium II (I),
//   where medium I is the current object, and medium II is another waveguide.
// 
//   Assumes that all modes are normalised, so that O_I_I and O_II_II are
//   unit matrices.
//
/////////////////////////////////////////////////////////////////////////////

class MultiWaveguide : public Waveguide
{
  public:

    MultiWaveguide(bool uniform=false, Material* core=NULL)
      : Waveguide(uniform, core),
        last_lambda(0.0), last_gain_mat(Material(0.0)) {}

    MultiWaveguide(const MultiWaveguide&);

    ~MultiWaveguide();

    bool recalc_needed() const;

    int N() const
      {return modeset.size();}

    Mode* get_mode(int i) const
      {return modeset[i-1];}

    const FieldExpansion field_from_source
      (const Coord& pos, const Coord& orientation);

    virtual Real S_flux(const FieldExpansion& f,
                        Real c1_start, Real c1_stop,
                        Real precision = 1e-10) const {return 0.0;};

    void add_mode(Mode& m)
      {modeset.push_back(&m);}

    void sort_modes()
      {std::sort(modeset.begin(), modeset.end(), modesorter());}

    void sort_modes_bloch()
      {std::sort(modeset.begin(), modeset.end(), modesorter_bloch());}
      
    void sort_modes_BDM()
      {std::sort(modeset.begin(), modeset.end(), modesorter_BDM());}
      
    virtual void calc_overlap_matrices
      (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL) {};

    void truncate_N_modes(int N=global.N);

    std::string repr() const;
    
  protected:

    std::vector<Mode*> modeset;
    
    // The wavelength and gain the modes were last calculated for,
    // are used to determine if recalculation is needed.

    Complex  last_lambda;
    Material last_gain_mat;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: MonoWaveguide
//
//   Waveguide where only a single mode is considered from a decoupled
//   continuous set (e.g. in infinite planar stratified media).
//
/////////////////////////////////////////////////////////////////////////////

class MonoWaveguide : public Waveguide
{
  public:

    MonoWaveguide(Material* core=NULL)
      : Waveguide(true, core) {}
    
    MonoWaveguide(const MonoWaveguide& w)
      : Waveguide(true, w.core) {mode = new Mode(*w.mode);}
    
    ~MonoWaveguide() {delete mode;}

    std::vector<Material*> get_materials() const;
    
    bool contains(const Material& m) const {return core == &m;}

    bool no_gain_present() const {return core->no_gain_present();}
    
    // We always recalculate MonoWaveguides, since the overhead of
    // checking for a needed recalc outweighs the possible gains.

    bool  recalc_needed() const {return true;}
    int   N()             const {return 1;}
    Mode* get_mode(int i) const {return mode;}

    std::string repr() const {return mode->repr();}
    
  protected:

    Mode* mode;
};



#endif
