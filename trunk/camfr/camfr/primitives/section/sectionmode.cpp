
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
#include "../slab/isoslab/slabmode.h"

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
  int M = geom->M2;
  global.N = geom->M2;

  int old_orthogonal = global.orthogonal;
  global.orthogonal = false;

  Complex old_beta = global.slab_ky;
  global.slab_ky = kz;

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
  global.slab_ky = old_beta;
  global.orthogonal = old_orthogonal;
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
  global.N = section->M2;  

  int old_orthogonal = global.orthogonal;
  global.orthogonal = false;

  Complex old_beta = global.slab_ky;
  global.slab_ky = kz;

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
    stack_coord.z = section->left.get_total_thickness() - coord.c1
      + I*global_section.left_PML;

    stack_coord.z_limit = (coord.c1_limit == Plus) ? Min : Plus;

    if (abs(coord.c1) < 1e-8) // TMP: kludge
      stack_coord.z -= 1e-8;

    stack_f = section->left.field(stack_coord);
    s = -1.0;
  }
  else
  {  
    stack_coord.z = coord.c1 - section->left.get_total_thickness()
      + I*global_section.left_PML;

    stack_coord.z_limit = coord.c1_limit;

    if (abs(coord.c1 - section->right.get_total_thickness()) < 1e-8)
      stack_coord.z -= 1e-8; // TMP: kludge

    stack_f = section->right.field(stack_coord);
    s = 1.0;
  }

  global.N = old_N;
  global.slab_ky = old_beta;
  global.orthogonal = old_orthogonal;

  // Transform fields from Stack to Section coordinates.

  Field f;

  f.E1 = stack_f.Ez * s;  f.H1 = stack_f.Hz;
  f.E2 = stack_f.E1;      f.H2 = stack_f.H1 * s;
  f.Ez = stack_f.E2;      f.Hz = stack_f.H2 * s;

  return f;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D_Mode::get_fw_bw
//
/////////////////////////////////////////////////////////////////////////////

void Section2D_Mode::get_fw_bw(const Complex& c, Limit c_limit,
                               cVector* fw, cVector* bw)
{
  // Initialise.

  Section2D* section = dynamic_cast<Section2D*>(geom);

  int old_N = global.N;
  global.N = section->M2;

  int old_orthogonal = global.orthogonal;
  global.orthogonal = false;

  Complex old_beta = global.slab_ky;
  global.slab_ky = kz;

  section->left .set_interface_field( left_interface_field);
  section->right.set_interface_field(right_interface_field);
  
  // Determine correct substack.

  Coord stack_coord;
  if (real(c) < real(section->left.get_total_thickness()))
  {    
    stack_coord.z = section->left.get_total_thickness() - c;
    stack_coord.z_limit = (c_limit == Plus) ? Min : Plus;

    if (abs(c) < 1e-8) // TMP: kludge
      stack_coord.z -= 1e-8;

    section->left.fw_bw_field(stack_coord,bw,fw);
  }
  else
  {  
    stack_coord.z = c - section->left.get_total_thickness();
    stack_coord.z_limit = c_limit;

    if (abs(c - section->right.get_total_thickness()) < 1e-8)
      stack_coord.z -= 1e-8; // TMP: kludge

    section->right.fw_bw_field(stack_coord,fw,bw);
  }

  global.N = old_N;
  global.slab_ky = old_beta;  
  global.orthogonal = old_orthogonal;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void Section2D_Mode::normalise() 
{
  vector<Complex> disc(geom->get_disc());
  disc.insert(disc.begin(), 0.0);
  
  Complex norm = 0.0;

  for (unsigned int k=0; k<disc.size()-1; k++)
    norm += overlap_slice(this, this, disc[k], disc[k+1]);
  
  norm = sqrt(norm);
  
  if (abs(norm) < 1e-10)
  {
    py_print("Warning: section mode close to cutoff.");
    norm = 1.0;
  }

  for (unsigned int i=0; i<left_interface_field.size(); i++)
  {
    left_interface_field[i].fw /= norm;
    left_interface_field[i].bw /= norm;
  }
  
  for (unsigned int i=0; i<right_interface_field.size(); i++)
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
  (Polarisation pol, const Complex& kz, SlabMode* m_, Section1D* geom) 
    : SectionMode(pol, kz, geom), m(m_)
{
  fw0 = 1.0;
  bw0 = (global_section.leftwall == E_wall) ? -1.0 : 1.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D_Mode::field
//
/////////////////////////////////////////////////////////////////////////////

Field Section1D_Mode::field(const Coord& coord) const 
{
  Complex old_beta = global.slab_ky;
  global.slab_ky = kz;

  Complex fw, bw;
  get_fw_bw(coord.c1, &fw, &bw);

  // Convert from Section to Slab coordinates.

  Coord slab_coord(coord);

  slab_coord.c1 = coord.c2;
  slab_coord.c1_limit = coord.c2_limit;

  slab_coord.c2 = 0.0;
  slab_coord.c2_limit = coord.z_limit;

  slab_coord.z = 0.0;
  slab_coord.z_limit = coord.c1_limit;

  Field f_slab = m->field(slab_coord);

  // Convert from Slab to Section coordinates.

  Field f;

  f.E1 = fw * f_slab.Ez  -  bw * f_slab.Ez;
  f.E2 = fw * f_slab.E1  +  bw * f_slab.E1;
  f.Ez = fw * f_slab.E2  +  bw * f_slab.E2;

  f.H1 = fw * f_slab.Hz  +  bw * f_slab.Hz;
  f.H2 = fw * f_slab.H1  -  bw * f_slab.H1;
  f.Hz = fw * f_slab.H2  -  bw * f_slab.H2;

  global.slab_ky = old_beta;

  return f;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D_Mode::get_fw_bw
//
/////////////////////////////////////////////////////////////////////////////

void Section1D_Mode::get_fw_bw(const Complex& c, 
                               Complex* fw, Complex* bw) const
{
  Complex old_beta = global.slab_ky;
  global.slab_ky = kz;

  Complex kx = m->get_kz();

  *fw = fw0 * exp(-I*kx*c);
  *bw = bw0 * exp(+I*kx*c);

  global.slab_ky = old_beta;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D_Mode::get_fw_bw
//
/////////////////////////////////////////////////////////////////////////////

void Section1D_Mode::get_fw_bw(const Complex& c, Limit c_limit, 
                               cVector* fw, cVector* bw)
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
