
/////////////////////////////////////////////////////////////////////////////
//
// File:          traceroot.cpp
// Author:        Peter.Bienstman@rug.ac.be
// Date:          20000221
// Version:       1.23
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "traceroot.h"
#include "../croot/mueller.h"

using std::vector;
using std::string;
using std::arg;

#include "../../../util/vectorutil.h"
#include "../../../util/cvector.h"

#ifdef _WIN32
#include <float.h>
#define ISNAN _isnan
#else
#define ISNAN isnan
#endif

/////////////////////////////////////////////////////////////////////////////
//
// traceroot
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> traceroot(vector<Complex>&     estimate1,
                          Function1D<Complex>& f, 
                          vector<Complex>&     params1,
                          vector<Complex>&     params2,
                          vector<Complex>&     forbiddenzeros,
                          int resolution,
                          string* fname)
{ 
  // Declare all variables before first goto-label.

  std::ofstream* s = fname ? new std::ofstream(fname->c_str()) : NULL;

  Real eps_coarse = global.eps_trace_coarse;
  Real eps_fine   = 10*machine_eps();
  Real eps_copies = (eps_coarse < 1e-6) ? 1e-6 : eps_coarse;
  
  int initial_resolution = resolution;
  int max_fine_iters     = 20;

  vector<Complex> params, step;
  vector<Complex> zeros, zeros_bis; // Two recent zero estimates.
  vector<Complex> zeros_try;        // New trial estimate.
  vector<Complex> zeros_jump;       // The last jump made by the zeros.
  vector<Complex> *oldest, *newest; // Points to either zeros or zeros_bis.

  vector<bool> valid(estimate1.size(),true);
  int fine_iters, trouble_index;
  bool trouble;

  if (estimate1.size() == 0)
    return estimate1;

  vector<Complex> params_orig = f.get_params();

  // Check if sweep needed.
  
  step = (params2 - params1) / resolution;
  
  if (abs(step) < 1e-13)
    return estimate1;

  // Write initial zero locations to file.

  restart:
  
  if (s)
    for (unsigned int i=0; i<estimate1.size();i++)
      *s << i << " "
         << real(estimate1[i]) << " "
         << imag(estimate1[i]) << " "
         << valid[i] << std::endl;

  // Initialise zeros_jump.

  zeros_jump.clear();
  if (estimate1.size() == 1)
    zeros_jump.push_back( 0.1 * estimate1[0] );
  else
  {
    for (unsigned int i=0; i<estimate1.size()-1; i++)
      zeros_jump.push_back( 0.1 * (estimate1[i] - estimate1[i+1]) );
    zeros_jump.push_back(zeros_jump[zeros_jump.size()-1]);
  }
  
  // Create a small perturbation of estimate1 as initial second estimate.
  
  zeros_bis.clear();
  for (unsigned int i=0; i<estimate1.size(); i++)
    if (abs(real(estimate1[i])) > abs(imag(estimate1[i])))
      zeros_bis.push_back(estimate1[i] + 10*eps_coarse*I);
    else
      zeros_bis.push_back(estimate1[i] + 10*eps_coarse);

  zeros      = estimate1;
  params     = params1;
  oldest     = &zeros_bis;
  newest     = &zeros;
  fine_iters = 0;

  do // Gradually sweep parameters from params1 to params2.
  {
    
    // See if we can try to increase the step.

    if (    (fine_iters++ == max_fine_iters)
         && (abs(5*step) < abs(params2-params)) )
    {
      step       *= 5;
      fine_iters  = 0;
      //cout << "Increasing step to " << abs(step) << endl;
    }
    
    params += step;

    // Remedy stagnation (i.e two estimates too close together).
      
    for (unsigned int i=0; i<zeros.size(); i++)
    {
      Real eps_stag = 2*eps_coarse;
      while (abs(zeros[i]-zeros_bis[i]) < 2*eps_coarse)
      {
        if (abs(real(zeros_bis[i])) > abs(imag(zeros_bis[i])))
          zeros_bis[i] += eps_stag*I;
        else
          zeros_bis[i] += eps_stag;   
        eps_stag *= 2;
      }
    }
    
    // Try to solve a new dispersion relation.

    undo:

    f.set_params(params);
    trouble       = false;
    trouble_index = -1;
    
    //cout << "Traceroot progress " << params1.size() << ":" 
    //     << abs(params-params1)/abs(params2-params1)*100 << "% " << endl;
    
    zeros_try.clear();
    for (unsigned int i=0; !trouble && i<zeros.size(); i++)
      if (valid[i])
      {        
        zeros_try.push_back(mueller(f,zeros[i],zeros_bis[i],eps_coarse,
                                    0,100,&trouble));  
        if (trouble)
        {
          //cout << "Mueller didn't converge for " << zeros[i] << endl;
          trouble_index = i;
          break;
        }
        
        //cout << setprecision(6);
        //cout << i << " " << zeros[i] << " " << zeros_bis[i]
        //<< "->" << zeros_try[i] << " " << zeros_jump[i] << endl;
      }
      else
        zeros_try.push_back(zeros[i]); // dummy fill

    // Check if mueller yielded NaN.

    for (unsigned int i=0; !trouble && i<zeros_try.size(); i++)
    
    if ( valid[i] &&
         ( ISNAN(real(zeros_try[i])) || ISNAN(imag(zeros_try[i])) ) )
        {
          //cout << "NaN detected in zero " << i << endl; 
          trouble       = true;
          trouble_index = i;
          break;
        }
    
    // Check if we lost a zero (i.e. two solutions are the same).
      
    for (unsigned int i=0; !trouble && i<zeros_try.size(); i++)
      for (unsigned int j=i+1; j<zeros_try.size(); j++)       
        if (    valid[i] && valid[j]
            && (abs(zeros_try[i] - zeros_try[j]) < eps_copies) )
        {
          //cout << "Zeros " <<i<< " and " <<j<< " are the same." << endl;
          trouble       = true;
          trouble_index = i;
          break;
        }

    // Check if we made a jump that is too large.
    
    for (unsigned int i=0; !trouble && i<zeros_try.size(); i++)
      if ( valid[i] && ( abs(zeros_jump[i]) > 1e-6 )
                    && ( abs(imag(zeros_try[i])) > 1e-4) )
        if (abs(zeros_try[i] - (*newest)[i]) > 2500 * abs(zeros_jump[i]))
        {
          //cout << "Too large jump for zero " << i << "." << endl;     
          trouble       = true;
          trouble_index = i;
          break;
        }

    // Check if we made a medium sized jump with too sharp a bend.

    for (unsigned int i=0; !trouble && i<zeros_try.size(); i++)         
      if (      valid[i]
           && ( abs(zeros_jump[i]) > 1e-6 )
           && ( abs(imag(zeros_try[i])) > 1e-4)     
           && ( abs(zeros_try[i] - (*newest)[i]) > 10 * abs(zeros_jump[i]) ) )
        if (    (abs(   arg(zeros_try[i] - (*newest)[i])
                      - arg(zeros_jump[i]))                  > pi )
             && (abs(abs(   arg(zeros_try[i] - (*newest)[i])
                          - arg(zeros_jump[i]))   - 2*pi   ) > pi ) )
        {
          //cout << "Too sharp a bend in jump of zero " << i << "." <<endl; 
          trouble       = true;
          trouble_index = i;
          break;
        }

    // Check if we converged to a forbidden zero.

    for (unsigned int i=0; !trouble && i<zeros_try.size(); i++)
      for (unsigned int j=0; j<forbiddenzeros.size(); j++)
        if (     valid[i]
             && (abs(zeros_try[i] - forbiddenzeros[j]) < eps_copies) )
        {
          //cout << "Zero " << i << " jumped to forbidden zero " << j <<endl;
          trouble       = true;
          trouble_index = i;
          break;
        }
    
    // In the presence of a convergence problem, take increasingly
    // drastic measures depending on the severity.
    
    if (trouble)
    {
      if (abs(step) > 1e-10) // Decrease step size and undo last move.
      {
        params     -= step;
        step       /= 10.0;
        params     += step;
        fine_iters  = 0;

        //cout << "Decreasing step to " << abs(step) << endl;

        if (s)
          *s << "-1 "
             << real(zeros_try[trouble_index]) << " "
             << imag(zeros_try[trouble_index]) << " "
             << valid[trouble_index] << std::endl;
        
        goto undo;
      }
      else if (resolution < pow(2.0,10)) // Refine resolution and restart.
      {
        resolution *= 2;
        step = (params2 - params1) / resolution;
        if (eps_coarse/2.0 > 10*machine_eps())
          eps_coarse /= 2.0;
        fine_iters  = 0;
        
        std::ostringstream st;
        st << "Restarting: resolution now " << resolution;
        py_print(st.str());

        // Clear file.

        if (s)
        {
          s->close();
          s->open(fname->c_str());
        }
        
        goto restart;
      }
      else // Give up on this zero and carry on.
      {
        std::ostringstream st;
        st << "Error: inrecoverable problem during sweep. "
          << "Invalidating zero " << zeros_try[trouble_index];
        py_error(st.str());

        valid[trouble_index] = false;
        step = (params2 - params) / initial_resolution;

        if (s)
          *s << "-2 "
             << real(zeros_try[trouble_index]) << " "
             << imag(zeros_try[trouble_index]) << " "
             << valid[trouble_index] << std::endl;
      }
    } // end if (trouble)
    
    // Update jumps.

    for (unsigned int i=0; i<zeros_jump.size(); i++)
      zeros_jump[i] = (*newest)[i] - (*oldest)[i];

    // Overwrite oldest with updated estimates.
    
    *oldest = zeros_try;

    // Swapping pointers is faster than swapping vectors.

    newest = oldest;
    oldest = (oldest == &zeros) ? &zeros_bis : oldest = &zeros;
    
    if (s)
      for (unsigned int i=0; i<zeros_try.size();i++)
        *s << i << " "
           << real(zeros_try[i]) << " "
           << imag(zeros_try[i]) << " "
           << valid[i] << std::endl;
    
  }
  while (abs(params-params2) >= abs(step/2));
    
  // Do a final fine search for zeros.

  f.set_params(params2);

  vector<Complex> zeros2;
  for (unsigned int i=0; i<zeros.size(); i++)
    if (valid[i])
    {
      if (eps_coarse <= eps_fine) // No need for fine search.
        zeros2.push_back((*newest)[i]);
      else
      {
        Real eps_stag = 2*eps_fine;
        while (abs(zeros[i]-zeros_bis[i]) < 2*eps_fine)
        {
          if (abs(real(zeros_bis[i])) > abs(imag(zeros_bis[i])))
            (*oldest)[i] += eps_stag*I;
          else
            (*oldest)[i] += eps_stag;
          eps_stag *= 2;
        }

        trouble = false;
        Complex final
          = mueller(f,zeros[i],zeros_bis[i],eps_fine,0,100,&trouble);

        if (trouble)
        {
          std::ostringstream st;
          st << "Mueller didn't converge on final search for " << final;
          py_print(st.str());
        }

        if (s)
          *s << i << " "
             << real(final) << " "
             << imag(final) << " "
             << valid[i] << std::endl;

        zeros2.push_back(final);
      }
    }

  delete s;

  f.set_params(params_orig); // Restore f. 
  
  return zeros2;
}



/////////////////////////////////////////////////////////////////////////////
//
// traceroot
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> traceroot(vector<Real>&        estimate1,
                          Function1D<Complex>& f, 
                          vector<Complex>&     params1,
                          vector<Complex>&     params2,
                          vector<Complex>&     forbiddenzeros,
                          int resolution,
                          string* fname)
{
  vector<Complex> c_estimate1;
  for (unsigned int i=0; i<estimate1.size(); i++)
    c_estimate1.push_back(Complex(estimate1[i], 0.0));

  return traceroot(c_estimate1, f, params1, params2, 
                   forbiddenzeros, resolution, fname); 
}



/////////////////////////////////////////////////////////////////////////////
//
// traceroot_chunks
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> traceroot_chunks(vector<Complex>&     estimate1,
                                 Function1D<Complex>& f, 
                                 vector<Complex>&     params1,
                                 vector<Complex>&     params2,
                                 vector<Complex>&     forbiddenzeros,
                                 int resolution,
                                 int chunk_length,
                                 int overlap,
                                 string* fname)
{ 
  vector<Complex> zeros2;
  int chunknumber=0;
  vector<Complex>::iterator start = estimate1.begin();

  if (estimate1.size() == 0)
    return estimate1;

  do
  {
    // Create chunk.
    
    chunknumber++;
    
    vector<Complex> chunk;    
    vector<Complex>::iterator stop;
    if (start+chunk_length+2*overlap > estimate1.end())
      stop = estimate1.end();
    else
      stop = start+chunk_length+2*overlap;
    
    chunk.insert(chunk.begin(),start,stop);

    string* chunkfname;
    if (fname)
    {
      std::ostringstream s;
      s << *fname << chunknumber;
      chunkfname = new string(s.str());
    }
    else
      chunkfname = NULL;
      
    // Find zeros of chunk.

    vector<Complex> zeros_chunk = traceroot
      (chunk, f, params1, params2, forbiddenzeros, resolution, chunkfname);

    delete chunkfname;
    
    zeros2.insert(zeros2.end(), zeros_chunk.begin(), zeros_chunk.end());

    start += chunk_length+overlap;
  }
  while (start+overlap < estimate1.end());

  const Real eps_copies = 1e-6;
  remove_copies(&zeros2, eps_copies);
  
  return zeros2;
}



/////////////////////////////////////////////////////////////////////////////
//
// traceroot_chunks
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> traceroot_chunks(vector<Real>&        estimate1,
                                 Function1D<Complex>& f, 
                                 vector<Complex>&     params1,
                                 vector<Complex>&     params2,
                                 vector<Complex>&     forbiddenzeros,
                                 int resolution,
                                 int chunk_length,
                                 int overlap,
                                 string* fname)
{
  vector<Complex> c_estimate1;
  for (unsigned int i=0; i<estimate1.size(); i++)
    c_estimate1.push_back(Complex(estimate1[i], 0.0));

  return traceroot_chunks(c_estimate1, f, params1, params2, forbiddenzeros,
                          resolution, chunk_length, overlap, fname); 
}
