
/////////////////////////////////////////////////////////////////////////////
//
// File:     material.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990503
// Version:  1.02
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_H
#define MATERIAL_H

#include <sstream>
#include <iostream>
#include "defs.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Material_length
//
//   A piece of material with a certain length.
//
/////////////////////////////////////////////////////////////////////////////

class BaseMaterial; // Forward declaration.

struct Material_length
{
    Material_length(BaseMaterial& mat_, const Complex& d_=0.0)
      : mat(&mat_), d(d_) {}

    Material_length(BaseMaterial* mat_, const Complex& d_=0.0)
      : mat(mat_), d(d_) {}
    
    BaseMaterial* mat;
    const Complex d;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BaseMaterial
//
/////////////////////////////////////////////////////////////////////////////

class BaseMaterial
{
  public:

    const Material_length operator() (const Complex& d=0.0) const;

    virtual const Complex epsr(int) const = 0; // Diagonal tensor elements.
    virtual const Complex  mur(int) const = 0;    

    virtual bool no_gain_present() const = 0;
    virtual std::string repr()     const = 0;
};

inline std::ostream& operator<<(std::ostream& s, const BaseMaterial& m)
  {return s << m.repr() << std::endl;}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Material
//
//   The purpose of this class is to allow for easy incorporation of
//   future, more complicated material models.
//
//   For example, to add dispersive materials, you would want to subclass
//   this class and override its functions.
//
//   Note: access the current wavelength through the global variable
//   'global.lambda' (see defs.h)
//  
/////////////////////////////////////////////////////////////////////////////

class Material : public BaseMaterial
{
  public:
    
    Material(const Complex& n) : i_n(n), i_etar(i_n) {}
    Material(const Complex& n, const Complex& etar) : i_n(n), i_etar(etar) {}
    
    const Complex   epsr() const {return i_n * i_etar;}
    const Complex    mur() const {return i_n / i_etar;}
    const Complex    eps() const {return epsr() * eps0;}
    const Complex     mu() const {return mur() * mu0;}
	  
    const Complex      n() const {return i_n;}
    const Complex eps_mu() const {return epsr() * eps0 * mur() * mu0;}
    const Complex   etar() const {return i_etar;}
    const Complex    eta() const {return i_etar * sqrt(eps0 / mu0);}
    
    const Complex epsr(int) const {return epsr();}
    const Complex  mur(int) const {return mur();}     
    
    Real gain() const {return 4*imag(i_n)*pi/(real(global.lambda)*1e-4);}
    
    void set_epsr_mur(const Complex& epsr, const Complex& mur);
    void set_epsr    (const Complex& epsr) {set_epsr_mur(epsr,   mur());}
    void set_mur     (const Complex& mur)  {set_epsr_mur(epsr(), mur);}

    void set_n(Complex n)        {i_n = n;}
    void set_n_imag(Real n_imag) {i_n = Complex(real(i_n), n_imag);}
    void set_etar(Complex etar)  {i_etar = etar;}
    
    bool no_gain_present() const {return (imag(i_n) < 1e-12);}
    
    bool operator==(const Material& m) const
      {return (abs(i_n - m.i_n) < 1e-12) && (abs(i_etar - m.i_etar) < 1e-12);}
    
    bool operator!=(const Material& m) const
      {return !(*this == m);};
    
    std::string repr() const
      {std::ostringstream s; s << "Isotropic n=" << i_n; return s.str();};

  protected:
    
    Complex i_n;
    Complex i_etar;
};

inline std::ostream& operator<<(std::ostream& s, const Material& m)
  {return s << m.repr() << std::endl;}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BiaxialMaterial
//
//   Material with a diagonal eps and mu tensor.
//
//   TODO: flesh out this class and improve inheritance hierarchy.
//  
/////////////////////////////////////////////////////////////////////////////

class BiaxialMaterial : public Material
{
  public:

     BiaxialMaterial(const Complex& epsr_1,
                     const Complex& epsr_2,
                     const Complex& epsr_3,
                     const Complex& mur_1,
                     const Complex& mur_2,
                     const Complex& mur_3) : Material(sqrt(epsr_1))
       {i_epsr[0] = epsr_1; i_epsr[1] = epsr_2; i_epsr[2] = epsr_3;
         i_mur[0] =  mur_1;  i_mur[1] =  mur_2;  i_mur[2] =  mur_3;}

     const Complex epsr(int i) const {return i_epsr[i-1];}
     const Complex  mur(int i) const {return  i_mur[i-1];}

     bool no_gain_present() const {return true;}

     bool operator==(const Material& m) const
       {return false;}

     bool operator!=(const Material& m) const
       {return !(*this == m);}

     std::string repr() const
       {std::ostringstream s; s << "Biaxial material"; return s.str();}
     
  protected:

     Complex i_epsr[3];
     Complex i_mur[3];
};

inline std::ostream& operator<<(std::ostream& s, const BiaxialMaterial& m)
  {return s << m.repr() << std::endl;}



/////////////////////////////////////////////////////////////////////////////
//
// signedsqrt(kz2, material)
//   Extra layer of indirection, for experimenting with sign of sqrt
//   in relation to material loss or gain.
//
//   Will be removed in the future.
//
/////////////////////////////////////////////////////////////////////////////

inline Complex signedsqrt(const Complex& kz2, const Material& material)
   {return sqrt(kz2);}

inline Complex signedsqrt(const Complex& kz2, const Material* material)
   {return signedsqrt(kz2, *material);}

inline Complex signedsqrt(const Complex& kz2)
{
  Complex kz = sqrt(kz2);

  if (real(kz) < 0) 
      kz = -kz;

  if (abs(real(kz)) < 1e-12)
    if (imag(kz) > 0)
      kz = -kz;

  return kz;
}



#endif
