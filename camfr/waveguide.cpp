
/////////////////////////////////////////////////////////////////////////////
//
// File:     waveguide.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991105
// Version:  1.0
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "waveguide.h"

/////////////////////////////////////////////////////////////////////////////
//
// Waveguide::operator()
//
/////////////////////////////////////////////////////////////////////////////

const Waveguide_length Waveguide::operator()(const Complex& d=0.0) const
{
  if (real(d) < 0)
    cout << "Warning: negative real length of waveguide." << endl;
  
  return Waveguide_length(const_cast<Waveguide*>(this),d);
}



/////////////////////////////////////////////////////////////////////////////
//
// MultiWaveguide::MultiWaveguide(const MultiWaveguide&)
//
/////////////////////////////////////////////////////////////////////////////

MultiWaveguide::MultiWaveguide(const MultiWaveguide& w)
  : Waveguide(w.uniform, w.core), 
    last_lambda(w.last_lambda), last_gain_mat(w.last_gain_mat)
{
  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];

  modeset.clear();
  
  for (unsigned int i=0; i<w.modeset.size(); i++)
  {
    Mode* m = new Mode(*(w.modeset[i]));
    modeset.push_back(m);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// MultiWaveguide::~MultiWaveguide
//
/////////////////////////////////////////////////////////////////////////////

MultiWaveguide::~MultiWaveguide()
{
  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
}



/////////////////////////////////////////////////////////////////////////////
//
// MultiWaveguide::recalc_needed
//  
/////////////////////////////////////////////////////////////////////////////

bool MultiWaveguide::recalc_needed() const
{
  const Real eps = 1e-10;

  if (abs(global.lambda - last_lambda) > eps)
    return true;

  if (!global.gain_mat)
    return false;

  if (!contains(*global.gain_mat))
    return false;

  if (abs(global.gain_mat-> n()  - last_gain_mat. n())  > eps)
    return true;
  
  if (abs(global.gain_mat->mur() - last_gain_mat.mur()) > eps)
    return true;

  return false;
}



/////////////////////////////////////////////////////////////////////////////
//
// MultiWaveguide::field_from_source
//  
/////////////////////////////////////////////////////////////////////////////

const FieldExpansion MultiWaveguide::field_from_source
   (const Coord& pos, const Coord& orientation)
{
  find_modes();

  cVector fw(global.N,fortranArray), bw(global.N,fortranArray);

  for (unsigned int i=1; i<=global.N; i++)
  {
    Field field = get_mode(i)->field(pos);
    
    fw(i) = -0.5 * (  field.E1 * orientation.c1
                    + field.E2 * orientation.c2
                    - field.Ez * orientation.z  ); // - : bw eigenmode

    bw(i) = -0.5 * (  field.E1 * orientation.c1
                    + field.E2 * orientation.c2
                    + field.Ez * orientation.z  );
  }
  
  return FieldExpansion(*this, fw, bw);
}


/////////////////////////////////////////////////////////////////////////////
//
// MultiWaveguide::truncate_N_modes
//  
/////////////////////////////////////////////////////////////////////////////

void MultiWaveguide::truncate_N_modes(int N=global.N)
{
  if (modeset.size() <= N)
    return;

  for (unsigned int i=N; i<modeset.size(); i++)
    delete modeset[i];

  modeset.erase(modeset.begin()+N, modeset.end());
}



/////////////////////////////////////////////////////////////////////////////
//
// MultiWaveguide::repr
//
/////////////////////////////////////////////////////////////////////////////

string MultiWaveguide::repr() const
{
  ostringstream s;
  
  for (int i=1; i<=N(); i++)
  {
    s << i-1 << " " << get_mode(i)->n_eff(); // i-1: reflect Python offset.
    if (i != N())
      s << endl;
  }
  
  return s.str();
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoWaveguide::get_materials
//
/////////////////////////////////////////////////////////////////////////////

vector<Material*> MonoWaveguide::get_materials() const
{
  vector<Material*> materials;
  materials.push_back(core);
  return materials;
}
