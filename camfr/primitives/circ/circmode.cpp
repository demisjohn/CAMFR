
/////////////////////////////////////////////////////////////////////////////
//
// File:     circmode.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     11991104
// Version:  1.2
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "circmode.h"
#include "circdisp.h"
#include "circoverlap.h"

/////////////////////////////////////////////////////////////////////////////
//
// Circ_M_Mode::Circ_M_Mode
//
/////////////////////////////////////////////////////////////////////////////

Circ_M_Mode::Circ_M_Mode(Polarisation pol, const Complex& kz,
                         const Circ_M* geom_)
  : Mode(pol, kz, -kz), geom(geom_)
{ 
  // Initialise kr for each ring.
    
  Complex ki_2;
  for (unsigned int i=0; i<geom->M; i++)
  {
    ki_2 = pow(2*pi/global.lambda,2) * geom->material[i]->epsr()
                                     * geom->material[i]->mur(); 
    kr.push_back(signedsqrt(ki_2 - kz*kz, geom->material[i]));
  } 
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void Circ_M_Mode::normalise()
{
  Complex power = overlap(this, this);

  if (abs(power) < 1e-10)
  {
    cout << "Warning: mode close to cutoff." << endl;
    power = 1;
  }

  A /= sqrt(power);
  B /= sqrt(power);
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M_Mode::field_ring
//
/////////////////////////////////////////////////////////////////////////////

Field Circ_M_Mode::field_ring(int i, const Coord& coord,
                              Complex* dEzdr=0, Complex* dHzdr=0,
                              bool ang_dep=true) const
{
  cerr << "Field profiles for structures with M rings "
       << "not yet implemented." << endl;
  exit (-1);
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2_Mode::Circ_2_Mode
//
/////////////////////////////////////////////////////////////////////////////

Circ_2_Mode::Circ_2_Mode(Polarisation pol_,   const Complex& kz_,
                         const Complex& kr1,  const Complex& kr2,
                         const Circ_M* geom_)
  : Circ_M_Mode(pol_, kz_, geom_) 
{  
  // Set kr vector.
  
  kr.clear();
  kr.push_back(kr1);
  kr.push_back(kr2);
  
  // Set constants.
  
  const Real eps      = 1e-10;
  const Real lambda   = global.lambda;
  const Real k0       = 2*pi/lambda;

  const Real order    = global_circ.order;
  
  const Complex r     = geom->radius[0];
  const Complex R     = geom->radius[1];

  const Material co = *(geom->material[0]);
  const Material cl = *(geom->material[1]);
  
  const Complex epsr1 = co.epsr();
  const Complex epsr2 = cl.epsr();
  const Complex  mur1 = co.mur();
  const Complex  mur2 = cl.mur();
  
  // Set scaling factors and check values.

  scaling_co = true;
  scaling_cl = true;

  if (abs(kr2) < 0.05)
    cout << "Warning: mode very close to cutoff. "
         << "Possible loss of precision." << endl;

  // Calculate J(z).exp(-abs(z_imag))  (scaling factors cancel here).

  Complex J1r, dJ1r = dJ(order, kr1*r, &J1r, NULL, scaling_co);
  Complex J2R, dJ2R = dJ(order, kr2*R, &J2R, NULL, scaling_cl);

  if (abs(J1r) < 1e-6)
    cout << "Warning: mode is zero of J1r. "
         << "Possible loss of precision for field profiles." << endl;

  // Calculate c_TE and c_TM.

  Complex a_TE, b_TE, a_TM, b_TM;
  Circ_2_closed disp(r, R, co, cl, lambda, int(order), guided, geom->hankel);
  disp.ab_TE_TM(&a_TE, &b_TE, &a_TM, &b_TM, kr2);
  
  const Complex c_TE = a_TE / b_TE;  
  const Complex c_TM = a_TM / b_TM;

  const Complex order_0_TE =    mur1 / kr1 * dJ1r / J1r * b_TE
                             -  mur2 / kr2              * a_TE;

  const Complex order_0_TM =   epsr1 / kr1 * dJ1r / J1r * b_TM
                             - epsr2 / kr2              * a_TM;
  
  // Calculate field expansion coefficients A and B.

  if ( (order==0) && (abs(order_0_TE) < abs(order_0_TM)) ) // TE
  {
    A = 0;
    B = J1r;
  }
  else // TM, EH or HE
  {
    const Complex K = pow(k0/kr1/kr2, 2) * (epsr2*mur2 - epsr1*mur1);

    const Real sign = (global_circ.fieldtype == cos_type) ? -1 : 1;
    
    if ( abs(order_0_TE) > abs(order_0_TM) )
    {
      A = J1r;
      B = A * kz * sign * order / r / k0 / c / mu0 / order_0_TE * b_TE * K;
    }
    else
    {
      B = J1r;
      A = B * kz * sign * order / r / k0 / c / eps0 / order_0_TM * b_TM * K;
    }
  }

  // Determine polarisation.

  // Pathological case of homogeneous medium, but necessary to be general.

  if (abs(geom->material[0]->n() - geom->material[1]->n()) < eps)
    pol = (abs(dJ2R) < abs(J2R)) ? TE : TM;
  else 
  {
    if (order == 0)
      pol = (abs(order_0_TE) < abs(order_0_TM)) ? TE : TM;
    else
    {
      const Complex C =   2.0 / kr1 *     epsr1 * mur1 * dJ1r / J1r
                        - 1.0 / kr2 * (   epsr1 * mur2 * c_TE
                                        + epsr2 * mur1 * c_TM);

      pol = (real(C) > 0) ? EH : HE;
    }
  }
  
  // Check boundary conditions.

  //cout << setprecision(15);
  //cout << "---" << endl;
  //cout << "Mode with neff " << kz/2/pi*lambda << endl;
  //cout << "pol " << Pol_string[pol] << endl;
  //cout << "kr1 : " << kr1 << endl;
  //cout << "kr2 : " << kr2 << endl;
  //cout << "A   : " << A << endl;
  //cout << "B   : " << B << endl;

  const Real theta = (order != 0) ? pi/4./order : 0.0;
        
  Field Rfield  = field_at(Coord(R,theta,0, Min));
  Field rpfield = field_at(Coord(r,theta,0, Plus));
  Field rmfield = field_at(Coord(r,theta,0, Min));
  
  Real error;
  Real worst = 0;  
  
  if ( (error = abs(  Rfield.Ez -           0.0         )) > worst )
    worst = error;
  if ( (error = abs(  Rfield.E2 -           0.0         )) > worst )
    worst = error;
  if ( (error = abs((rpfield.Ez - rmfield.Ez)/rmfield.Ez)) > worst )
    worst = error;
  if ( (error = abs((rpfield.Hz - rmfield.Hz)/rmfield.Hz)) > worst )
    worst = error;
  if ( (error = abs((rpfield.E2 - rmfield.E2)/rmfield.E2)) > worst )
    worst = error;
  if ( (error = abs((rpfield.H2 - rmfield.H2)/rmfield.H2)) > worst )
    worst = error;

  //cout << "Worst error : " << worst << endl;

  if (worst > 1e-8)
    cout << "Warning: error " << worst << " higher than 1e-8 "
         << "in boundary conditions check: kz = " << kz_ << endl; 
}


/////////////////////////////////////////////////////////////////////////////
//
// calc_phi_dependence
//
//  Auxiliary function to calculate phi-dependence of E fields (f)
//  and H fields (g).
//
/////////////////////////////////////////////////////////////////////////////

void calc_phi_dependence(int order, const Complex& phi, bool ang_dep,
                         Complex* f, Complex* df, Complex* g, Complex* dg)
{
  Complex cs, d_cs, sn, d_sn;
  Real n = order;

  if (n != 0)
  {
      cs = ang_dep ?    cos(n*phi) : 1.0;
    d_cs = ang_dep ? -n*sin(n*phi) : 1.0;

      sn = ang_dep ?    sin(n*phi) : 1.0;
    d_sn = ang_dep ?  n*cos(n*phi) : 1.0;
  }
  else
  {
      cs = 1.0;
    d_cs = 0.0;

      sn = 1.0;
    d_sn = 0.0;
  }
  
  if (global_circ.fieldtype == cos_type)
  {
    *f = cs; *df = d_cs; *g = sn; *dg = d_sn;
  }
  else
  {
    *f = sn; *df = d_sn; *g = cs; *dg = d_cs;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2_Mode::field_core
//
/////////////////////////////////////////////////////////////////////////////

Field Circ_2_Mode::field_core(const Coord& coord, Complex* dEzdr_p,
                              Complex* dHzdr_p, bool ang_dep=true) const
{
  // Set constants.

  const Real eps2   = 1e-10;

  const Real order  = global_circ.order;
  
  const Complex k0  = 2*pi/global.lambda;
  const Complex kr1 = kr[0];
  
  const Complex mu  = geom->material[0]->mu();
  const Complex eps = geom->material[0]->eps();

  const Complex rho = coord.c1;
  const Complex r   = geom->radius[0];

  // Calculate bessel functions and scaling factors.
  
  Complex J1r, J1rho, dJ1rho;
  
  J1r     = J(order, kr1*r,                 scaling_co); 
  dJ1rho = dJ(order, kr1*rho, &J1rho, NULL, scaling_co);

  Complex s = 1; 
  if (scaling_co)
    s = exp(abs(imag(kr1))*(rho-r)); // Assumes rho real in core.
  
  // Calculate phi dependence. Note: z-dependence is handled in 'field.cpp'.

  Complex f,df,g,dg;
  calc_phi_dependence(int(order), coord.c2, ang_dep, &f, &df, &g, &dg);
  
  // Resulting field.

  Field field;

  field.Ez             = A *        J1rho / J1r * s *  f;
  const Complex dEz_dr = A * kr1 * dJ1rho / J1r * s *  f;
  const Complex dEz_df = A *        J1rho / J1r * s * df;
  
  field.Hz             = B *        J1rho / J1r * s *  g;
  const Complex dHz_dr = B * kr1 * dJ1rho / J1r * s *  g;
  const Complex dHz_df = B *        J1rho / J1r * s * dg;

  Complex dEz_df_div_rho, dHz_df_div_rho;
  
  if (abs(rho) > eps2)
  {
    dEz_df_div_rho = dEz_df / rho;
    dHz_df_div_rho = dHz_df / rho;
  }
  else // rho = 0
  {
    if (order == 1)
    {
      dEz_df_div_rho = A * (kr1/2.0) / J1r * s * df;
      dHz_df_div_rho = B * (kr1/2.0) / J1r * s * dg;
    }
    else
    {
      dEz_df_div_rho = 0.0;
      dHz_df_div_rho = 0.0;
    } 
  }

  const Complex factor = I / kr1 / kr1;
    
  field.E1 = factor * ( - k0 * c * mu * dHz_df_div_rho
                        -          kz * dEz_dr );
  
  field.H1 = factor * ( + k0 * c * eps * dEz_df_div_rho
                        -           kz * dHz_dr );

  field.E2 = factor * ( + k0 * c * mu * dHz_dr
                        -          kz * dEz_df_div_rho );

  field.H2 = factor * ( - k0 * c * eps * dEz_dr
                        -           kz * dHz_df_div_rho );

  if (dEzdr_p)
    *dEzdr_p = dEz_dr;

  if (dHzdr_p)
    *dHzdr_p = dHz_dr;
  
  return field;
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2_Mode::field_cladding
//
/////////////////////////////////////////////////////////////////////////////

Field Circ_2_Mode::field_cladding(const Coord& coord, Complex* dEzdr_p,
                                  Complex* dHzdr_p, bool ang_dep=true) const
{
  // Set constants.
  
  const Real order  = global_circ.order;
  
  const Complex k0  = 2*pi/global.lambda;
  const Complex kr2 = kr[1];
  
  const Complex mu  = geom->material[1]->mu();
  const Complex eps = geom->material[1]->eps();

  const Complex rho = coord.c1;
  const Complex r   = geom->radius[0];
  const Complex R   = geom->radius[1];

  // Calculate bessel functions.
  
  Complex J2r,           H2r;
  Complex J2R,   dJ2R,   H2R,   dH2R;
  Complex J2rho, dJ2rho, H2rho, dH2rho;
  
  J2r    =  J(order, kr2*r,                 scaling_cl);
  dJ2R   = dJ(order, kr2*R,   &J2R,   NULL, scaling_cl);
  dJ2rho = dJ(order, kr2*rho, &J2rho, NULL, scaling_cl);

  if (geom->hankel == kind_1)
  {
    H2r    =  H1(order, kr2*r,                 scaling_cl);
    dH2R   = dH1(order, kr2*R,   &H2R,   NULL, scaling_cl);
    dH2rho = dH1(order, kr2*rho, &H2rho, NULL, scaling_cl);
  }
  else
  {
    H2r    =  H2(order, kr2*r,                 scaling_cl);
    dH2R   = dH2(order, kr2*R,   &H2R,   NULL, scaling_cl);
    dH2rho = dH2(order, kr2*rho, &H2rho, NULL, scaling_cl);
  }

  // Calculate scaling factors.

  Complex s_r   = 0.0;
  Complex s_rho = 0.0;
  Complex s_R   = 0.0;
  Complex s     = 1.0;

  Complex dummy;
  
  if (scaling_cl)
  {
    s_r   = scale(&J2r,   &dummy,  kr2*r,   geom->hankel);
    s_rho = scale(&J2rho, &dJ2rho, kr2*rho, geom->hankel);
    s_R   = scale(&J2R,   &dJ2R,   kr2*R,   geom->hankel);

    s = (geom->hankel == kind_1) ? exp(-I*kr2*(r-rho))
                                 : exp(+I*kr2*(r-rho));
  }

  const Complex s1 = exp(s_rho - s_R);
  const Complex s2 = exp(s_r   - s_R);

  // Calculate phi dependence. Note: z-dependence is handled in 'field.cpp'.

  Complex f,df,g,dg;
  calc_phi_dependence(int(order), coord.c2, ang_dep, &f, &df, &g, &dg);
  
  // Resulting field.

  Field field;
  
  // Ez + derivatives

  field.Ez       = A       * (J2R *  H2rho  -  H2R * J2rho  * s1)
                           / (J2R *  H2r    -  H2R * J2r    * s2) * s *  f;
  
  Complex dEz_dr = A * kr2 * (J2R * dH2rho  -  H2R * dJ2rho * s1)
                           / (J2R *  H2r    -  H2R *  J2r   * s2) * s *  f;
  
  Complex dEz_df = A       * (J2R *  H2rho  -  H2R * J2rho  * s1)
                           / (J2R *  H2r    -  H2R * J2r    * s2) * s * df;
  
  // Hz + derivatives

  field.Hz       = B       * (dJ2R *  H2rho - dH2R * J2rho  * s1)
                           / (dJ2R *  H2r   - dH2R * J2r    * s2) * s *  g;
  
  Complex dHz_dr = B * kr2 * (dJ2R * dH2rho - dH2R * dJ2rho * s1)
                           / (dJ2R *  H2r   - dH2R * J2r    * s2) * s *  g;
  
  Complex dHz_df = B       * (dJ2R *  H2rho - dH2R * J2rho  * s1)
                           / (dJ2R *  H2r   - dH2R * J2r    * s2) * s * dg;
  
  // Other field components.

  const Complex factor = I / kr2 / kr2;
  
  field.E1 = factor * ( - k0 * c * mu / rho * dHz_df
                        -                kz * dEz_dr );
  
  field.H1 = factor * ( + k0 * c * eps / rho * dEz_df
                        -                 kz * dHz_dr );

  field.E2 = factor * ( + k0 * c * mu  * dHz_dr
                        -     kz / rho * dEz_df );

  field.H2 = factor * ( - k0 * c * eps * dEz_dr
                        -     kz / rho * dHz_df );
  
  if (dEzdr_p)
    *dEzdr_p = dEz_dr;

  if (dHzdr_p)
    *dHzdr_p = dHz_dr;
  
  return field; 
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_1_Mode::Circ_1_Mode
//
/////////////////////////////////////////////////////////////////////////////

Circ_1_Mode::Circ_1_Mode(Polarisation pol_,  const Complex& kz_,
                         const Complex& kr_, const Circ_M*  geom_)
  : Circ_M_Mode(pol_, kz_, geom_) 
{
  // Set kr vector.

  kr.clear();
  kr.push_back(kr_);

  // Set constants and check values.
  
  const Real order = global_circ.order;
  const Complex R  = geom->radius[0];

  if (abs(kr[0]) < 0.05)
    cout << "Warning: mode very close to cutoff. "
         << "Possible loss of precision." << endl;
  
  // Determine polarisation and set field expansion coefficients A and B.

  Complex JR, dJR;
  dJR = dJ(order, kr[0]*R, &JR);
  
  if (abs(dJR) < abs(JR))
  {
    pol = TE;
    A   = 0;
    B   = 1;
  }
  else
  {
    pol = TM;
    A   = 1;
    B   = 0;
  }
    
  // Check boundary conditions.

  //cout << "---" << endl;
  //cout << "Mode with neff " << kz/2/pi*global.lambda << endl;
  //cout << "pol " << Pol_string[pol] << endl;
  //cout << "kr  : " << kr[0] << endl;
  
  const Real theta = (order != 0) ? pi/4./order : 0.0;
  
  Field Rfield = field_at(Coord(R,theta,0, Min));
  
  Real error;
  Real worst = 0;

  if ( (error = abs(Rfield.Ez - 0.0)) > worst ) worst = error;
  if ( (error = abs(Rfield.E2 - 0.0)) > worst ) worst = error;
  
  //cout << "Worst error : " << worst << endl;

  if (worst > 1e-8)
    cout << "Warning: error " << worst << " higher than 1e-8 "
         << "in boundary conditions check." << endl; 
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_1_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void Circ_1_Mode::normalise()
{
  // Faster equivalent of Complex power = overlap(this, this);
  
  // Set constants.
  
  const Real order  = global_circ.order;
  
  const Complex mu  = geom->material[0]->mu();
  const Complex eps = geom->material[0]->eps();
  
  const Real omega  = 2*pi/global.lambda * c;
  
  const Complex R   = geom->radius[0];

  Complex JR, dJR;
  dJR = dJ(order, kr[0]*R, &JR);

  // Calculate power.
  
  Complex power;

  if (pol == TE)
    power = -mu  * pi*kz*omega/2.0/kr[0]/kr[0] *  JR*JR  * B*B
                 * (R*R - order*order/kr[0]/kr[0]);
  else
    power = -eps * pi*kz*omega/2.0/kr[0]/kr[0] * dJR*dJR * A*A
                 *  R*R;

  if (order==0)
    power *= 2;
    
  // Normalise.

  if (abs(power) < 1e-10)
  {
    cout << "Warning: mode close to cutoff." << endl;
    power = 1;
  }
  
  A /= sqrt(power);
  B /= sqrt(power);
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_1_Mode::field_ring
//
/////////////////////////////////////////////////////////////////////////////

Field Circ_1_Mode::field_ring(int notused, const Coord& coord,
                              Complex* dEzdr_p=0, Complex* dHzdr_p=0,
                              bool ang_dep=true) const
{
  // Set constants.

  const Real eps2   = 1e-10;
  const Real order  = global_circ.order;
  
  const Complex k0  = 2*pi/global.lambda;
  
  const Complex mu  = geom->material[0]->mu();
  const Complex eps = geom->material[0]->eps();

  const Complex rho = coord.c1;
  
  Complex Jrho, dJrho;
  dJrho = dJ(order, kr[0]*rho, &Jrho);
  
  // Calculate phi dependence. Note: z-dependence is handled in 'field.cpp'.
  
  Complex f,df,g,dg;
  calc_phi_dependence(int(order), coord.c2, ang_dep, &f, &df, &g, &dg);
    
  // Resulting field.

  Field field;

  field.Ez             = A *          Jrho *  f;
  const Complex dEz_dr = A * kr[0] * dJrho *  f; 
  const Complex dEz_df = A *          Jrho * df;
  
  field.Hz             = B *          Jrho *  g;
  const Complex dHz_dr = B * kr[0] * dJrho *  g;
  const Complex dHz_df = B *          Jrho * dg;

  Complex dEz_df_div_rho, dHz_df_div_rho;
  
  if (abs(rho) > eps2)
  {
    dEz_df_div_rho = dEz_df / rho;
    dHz_df_div_rho = dHz_df / rho;
  }
  else // rho = 0
  {
    if (order == 1)
    {
      dEz_df_div_rho = A * (kr[0]/2.0) * df;
      dHz_df_div_rho = B * (kr[0]/2.0) * dg;
    }
    else
    {
      dEz_df_div_rho = 0.0;
      dHz_df_div_rho = 0.0;
    } 
  }
    
  const Complex factor = I / kr[0] / kr[0];
  
  field.E1 = factor * ( - k0 * c * mu * dHz_df_div_rho
                        -          kz * dEz_dr );
  
  field.H1 = factor * ( + k0 * c * eps * dEz_df_div_rho
                        -           kz * dHz_dr );

  field.E2 = factor * ( + k0 * c * mu * dHz_dr
                        -          kz * dEz_df_div_rho );

  field.H2 = factor * ( - k0 * c * eps * dEz_dr
                        -           kz * dHz_df_div_rho );

  if (dEzdr_p)
    *dEzdr_p = dEz_dr;

  if (dHzdr_p)
    *dHzdr_p = dHz_dr;
  
  return field;
}
