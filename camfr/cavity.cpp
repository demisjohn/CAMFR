
/////////////////////////////////////////////////////////////////////////////
//
// File:     cavity.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000104
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <vector>
#include "cavity.h"
#include "math/calculus/minimum/minimum.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// Sigma_lambda::operator()
//  
/////////////////////////////////////////////////////////////////////////////

Real Sigma_lambda::operator()(const Real& lambda)
{
  counter++;
  
  global.lambda = lambda;
  
  return cavity->calc_sigma();
}



/////////////////////////////////////////////////////////////////////////////
//
// Sigma_n_imag::operator()
//  
/////////////////////////////////////////////////////////////////////////////

Real Sigma_n_imag::operator()(const Real& n_imag)
{
  counter++;
  
  global.gain_mat->set_n_imag(n_imag);
  
  return cavity->calc_sigma();
}



/////////////////////////////////////////////////////////////////////////////
//
// Cavity::Cavity
//  
/////////////////////////////////////////////////////////////////////////////

Cavity::Cavity(Stack& top_, Stack& bot_)
  : top(&top_), bot(&bot_),
    sigma_lambda(this), sigma_n_imag(this)
{
  if (!(    dynamic_cast<MultiScatterer*>(top->get_sc())
         && dynamic_cast<MultiScatterer*>(bot->get_sc()) ))
  {
    py_error("Error: only MultiScatterers allowed in Cavity.");
    exit (-1);
  }
  
  if (top->get_inc() != bot->get_inc())
  {
    py_error(
     "Error: top and bottom half of cavity have different incidence media.");
    exit (-1);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Cavity::find_modes_in_region
//  
/////////////////////////////////////////////////////////////////////////////

void Cavity::find_modes_in_region
  (Real lambda_start, Real lambda_stop, Real delta_lambda,
   unsigned int number, Real n_imag_start, Real n_imag_stop,
   unsigned int passes)
{ 
  // Sweep wavelength and look for possible minima.

  vector<Real> Ax, Bx;
  bracket_all_minima
    (sigma_lambda, lambda_start, lambda_stop, Ax, Bx, delta_lambda);

  for (unsigned int i=0; i<Ax.size(); i++)
  {
    std::ostringstream s;
    s << "Minimum " << i << " between " << Ax[i] << " and " << Bx[i] << ".";
    py_print(s.str());
  }
  
  // Do a fine search for the minima.

  unsigned int no_of_modes = Ax.size();
  if (number && (number < no_of_modes))
    no_of_modes = number;

  if (Ax.size() == 0)
  {
    py_print("No minimum found in this region.");
    return;
  }

  int modes_found = 0;
  for (int i=Ax.size()-1; i>=0; i--)
  {
    std::ostringstream s;
    s << "Looking for mode between " << Ax[i] << " and " << Bx[i] << ".";
    py_print(s.str());
    
    find_mode(Ax[i], Bx[i],n_imag_start, n_imag_stop, passes);

    modes_found++;
    if (modes_found == no_of_modes)
      break;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Cavity::find_mode
//  
/////////////////////////////////////////////////////////////////////////////

void Cavity::find_mode(Real lambda_start, Real lambda_stop,
                       Real n_imag_start, Real n_imag_stop,
                       unsigned int passes)
{  
  if (passes < 1)
  {
    py_print("Invalid number of passes. Setting it to 1.");
    passes = 1;
  }

  Real lambda_prec = 1e-5; // Gives 0.01 nm precision. 
  Real n_imag_prec = 5e-5; // Gives 1 cm-1 precision.
  
  for (unsigned int i=1; i<=passes; i++)
  {
    // Sweep wavelength.

    Real lambda = brent_minimum
      (sigma_lambda, lambda_start, lambda_stop, lambda_prec);
    global.lambda = lambda;

    // Sweep gain.

    if (!global.gain_mat)
    {
      py_error("Error: no gain material defined for cavity.");
      exit (-1);
    }
    
    Real n_imag = brent_minimum
      (sigma_n_imag, n_imag_start, n_imag_stop, n_imag_prec);
    global.gain_mat->set_n_imag(n_imag);

    // Print diagnostics.

    std::ostringstream s;

    s << "Done pass " << i 
      << ": lambda " << global.lambda
      << ", gain "   << global.gain_mat->gain();

    py_print(s.str());
    
    // Refine intervals and precisions for next pass.

    Real delta_lambda = lambda_stop - lambda_start;
    lambda_start = global.lambda - delta_lambda/10.0;
    lambda_stop  = global.lambda + delta_lambda/10.0;

    Real delta_n_imag = n_imag_stop - n_imag_start;
    n_imag_start = n_imag - delta_n_imag/10.0;
    n_imag_stop  = n_imag + delta_n_imag/10.0;

    lambda_prec /= 10.0;
    n_imag_prec /= 10.0;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Cavity::calc_sigma
//  
/////////////////////////////////////////////////////////////////////////////

Real Cavity::calc_sigma(int* dominant_mode)
{
  int N = global.N;
  
  // Calculate matrices.

  top->calcRT();
  bot->calcRT();
  
  // Do SVD of cavity matrix Q = R_top.R_bot - U1
  
  cMatrix Q(N, N, fortranArray);
  Q.reference(multiply(top->as_multi()->get_R12(),
                       bot->as_multi()->get_R12()));
  for (int i=1; i<=N; i++)
    Q(i,i) -= 1.0;

  cMatrix Vh(N, N, fortranArray);
  rVector sigma(N, fortranArray);
  sigma.reference(svd(Q, &Vh));

  // Find dominant mode.

  Real dominant_power = 0.0;

  int mode;
  if (!dominant_mode)
     dominant_mode = &mode;
  *dominant_mode = 0;

  cVector bot_field(N, fortranArray);
  cVector top_field(N, fortranArray);

  for (int i=1; i<=N; i++)
  {
    bot_field(i) = conj(Vh(N,i));
    if (abs(Vh(N,i)) > dominant_power)
    {
      dominant_power = abs(Vh(N,i));
      *dominant_mode = i;
    }
  }

  // Set incident fields.

  bot->set_inc_field(bot_field);
  
  top_field = bot->get_refl_field();
  top->set_inc_field(top_field);
  
  // Print diagnostics.

  std::ostringstream s;

  s << "@ " << global.lambda
    << " " << imag(global.gain_mat->n())
    << " " << global.gain_mat->gain()
    << " " << sigma(N)
    << " " << *dominant_mode;
  
  py_print(s.str());
  
  current_sigma = sigma(N);
  
  return sigma(N);
}



/////////////////////////////////////////////////////////////////////////////
//
// Cavity::set_source
//  
/////////////////////////////////////////////////////////////////////////////

void Cavity::set_source(const FieldExpansion& f)
{
  int N = global.N;
  
  // Calculate matrices.

  top->calcRT(); const cMatrix& Ru(top->as_multi()->get_R12());
  bot->calcRT(); const cMatrix& Rd(bot->as_multi()->get_R12());

  cMatrix U1(N,N,fortranArray), tmp(N,N,fortranArray);
  U1 = 0.0;
  for (int i=1; i<=N; i++)
    U1(i,i) = Complex(1.0, 0.0);

  cMatrix inv_Rd_Ru(N,N,fortranArray);
  tmp = U1 - multiply(Rd, Ru);
  if (global.stability != SVD)
    inv_Rd_Ru.reference(invert    (tmp));
  else
    inv_Rd_Ru.reference(invert_svd(tmp));

  cMatrix inv_Ru_Rd(N,N,fortranArray);
  tmp = U1 - multiply(Ru, Rd);
  if (global.stability != SVD)
    inv_Ru_Rd.reference(invert    (tmp));
  else
    inv_Ru_Rd.reference(invert_svd(tmp));

  // Calculate total field incident on top mirror.

  cVector U_equivalent_source(N,fortranArray);
  U_equivalent_source = f.fw + multiply(Rd,f.bw);
  
  cVector U_top(N,fortranArray);
  U_top = multiply(inv_Rd_Ru, U_equivalent_source);
  
  top->set_inc_field(U_top);

  // Calculate total field incident on bottom mirror.

  cVector D_equivalent_source(N,fortranArray);
  D_equivalent_source = f.bw + multiply(Ru,f.fw);

  cVector D_bot(N,fortranArray);
  D_bot = multiply(inv_Ru_Rd, D_equivalent_source);
  
  bot->set_inc_field(D_bot);
}



/////////////////////////////////////////////////////////////////////////////
//
// Cavity::field
//  
/////////////////////////////////////////////////////////////////////////////

Field Cavity::field(const Coord& c)
{
  if (real(c.z) >= 0)
    return top->field(c);

  // Make necessary sign changes for bottom field (fw modes <-> bw modes).

  Complex z = -c.z;
  Limit z_limit = (c.z_limit == Plus) ? Min : Plus;
  Field f = bot->field(Coord(c.c1,c.c2,z,c.c1_limit,c.c2_limit,z_limit));

  f.Ez *= -1.0;
  f.H1 *= -1.0;
  f.H2 *= -1.0;

  return f;
}


