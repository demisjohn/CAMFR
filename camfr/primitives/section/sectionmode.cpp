
/////////////////////////////////////////////////////////////////////////////
//
// File:     sectionmode.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020307
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "sectionmode.h"
#include "sectionoverlap.h"
#include "../slab/generalslab.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

/////////////////////////////////////////////////////////////////////////////
//
// Section2D_Mode::Section2D_Mode
//
/////////////////////////////////////////////////////////////////////////////

Section2D_Mode::Section2D_Mode
  (Polarisation pol, const Complex& kz, Section2D* geom) 
    : SectionMode(pol, kz, geom) 
{
  // Initialise.

  int old_N = global.N;
  int M = geom->M;
  global.N = geom->M;

  Complex old_beta = global_slab.beta;
  global_slab.beta = kz;

  // Do eigenvalue calculation with eigenvectors.

  cMatrix Q(M,M,fortranArray);
  if (! geom->symmetric)
    Q.reference(multiply(geom-> left.as_multi()->get_R12(), 
                         geom->right.as_multi()->get_R12()));
  else
    Q.reference(multiply(geom-> left.as_multi()->get_R12(), 
                         geom-> left.as_multi()->get_R12()));
  
  cVector e(M,fortranArray);
  cMatrix E(M,M,fortranArray);
  if (global.stability == normal)
    e.reference(eigenvalues  (Q, &E));
  else
    e.reference(eigenvalues_x(Q, &E));

  // Find eigenvalue closest to 1.

  int index = 1;
  for (int i=2; i<=M; i++)
    if (abs(e(i) - 1.0) < abs(e(index) - 1.0))
      index = i;
  
  cVector f(M,fortranArray);
  for (int i=1; i<=M; i++)
    f(i) = E(i, index);

  // Set fields.

  geom->right.set_inc_field(f);
  geom->left .set_inc_field(geom->right.get_refl_field());
  
  geom->left .get_interface_field( &left_interface_field);
  geom->right.get_interface_field(&right_interface_field);

  global.N = old_N;
  global_slab.beta = old_beta;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D_Mode::field
//
/////////////////////////////////////////////////////////////////////////////

Field Section2D_Mode::field(const Coord& coord) const 
{
  // Initialise.

  Section2D* section = dynamic_cast<Section2D*>(geom);

  int old_N = global.N;
  global.N = section->M;

  Complex old_beta = global_slab.beta;
  global_slab.beta = kz;

  section->left .set_interface_field( left_interface_field);
  section->right.set_interface_field(right_interface_field); 

  // Convert from Section to Stack coordinates.

  Coord stack_coord(coord);

  stack_coord.c1 = coord.c2;
  stack_coord.c1_limit = coord.c2_limit;

  stack_coord.c2 = coord.z;
  stack_coord.c2_limit = coord.z_limit;
  
  // Determine correct substack.

  Field stack_f;
  Complex s;
  if (real(coord.c1) < real(section->left.get_total_thickness()))
  {    
    stack_coord.z = section->left.get_total_thickness() - coord.c1;
    stack_coord.z_limit = (coord.c1_limit == Plus) ? Min : Plus;

    stack_f = section->left.field(stack_coord);
    s = -1.0;
  }
  else
  {  
    stack_coord.z = coord.c1 - section->left.get_total_thickness();
    stack_coord.z_limit = coord.c1_limit;

    stack_f = section->right.field(stack_coord);
    s = 1.0;
  }

  global.N = old_N;
  global_slab.beta = old_beta;

  // Transform fields from Stack to Section coordinates.

  Field f;

  f.E1 = stack_f.Ez * s;  f.H1 = stack_f.Hz;
  f.E2 = stack_f.E1;      f.H2 = stack_f.H1 * s;
  f.Ez = stack_f.E2;      f.Hz = stack_f.H2 * s;

  return f;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void Section2D_Mode::normalise() 
{
  Complex norm = sqrt(overlap(this, this));
  if (abs(norm) < 1e-10)
  {
    py_print("Warning: section mode close to cutoff.");
    norm = 1.0;
  }

  for(unsigned int i=0; i<left_interface_field.size(); i++)
  {
    left_interface_field[i].fw /= norm;
    left_interface_field[i].bw /= norm;
  }
  
  for(unsigned int i=0; i<right_interface_field.size(); i++)
  {
    right_interface_field[i].fw /= norm;
    right_interface_field[i].bw /= norm;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D_Mode::Section1D_Mode
//
/////////////////////////////////////////////////////////////////////////////

Section1D_Mode::Section1D_Mode
  (Polarisation pol, const Complex& kz, Section1D* geom) 
    : SectionMode(pol, kz, geom) 
{
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D_Mode::field
//
/////////////////////////////////////////////////////////////////////////////

Field Section1D_Mode::field(const Coord& coord) const 
{
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void Section1D_Mode::normalise() 
{
  cout << "Warning: Section1D_Mode::normalise not yet implemented." << endl;
}
