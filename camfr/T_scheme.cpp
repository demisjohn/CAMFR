
/////////////////////////////////////////////////////////////////////////////
//
// File:     T_scheme.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000213
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "T_scheme.h"
#include "scatterer.h"

/////////////////////////////////////////////////////////////////////////////
//
// T_scheme for monoscatterers.
//  
/////////////////////////////////////////////////////////////////////////////

void T_scheme
  (const vector<Chunk>& chunks,vector<Complex>* fw, vector<Complex>* bw)
{
  // Initialise.

  if ( (fw->size() != 1) || (bw->size() != 1) )
  {
    cerr << "Error: incorrect setting of input fields." << endl;
    exit (-1);
  }

  Complex fw_chunk_begin_scaled = (*fw)[0];
  Complex bw_chunk_begin_scaled = (*bw)[0];

  Complex scaling_factor = 1;
  
  // Loop through chunks and relate fields at the end of each chunk to
  // those at the beginning of the chunk.

  Complex fw_chunk_end_scaled;
  Complex bw_chunk_end_scaled;

  for (unsigned int k=0; k<chunks.size(); k++)
  {
    MonoScatterer* s = dynamic_cast<MonoScatterer*>(chunks[k].sc);
    if (!s)
    {
      cerr << "Error: T-scheme only valid for MonoScatterers." << endl;
      exit (-1);
    }

    const Complex& R12(s->get_R12()); const Complex& R21(s->get_R21());
    const Complex& T12(s->get_T12()); const Complex& T21(s->get_T21());

    // Cross the interface.
    
    fw_chunk_end_scaled = (T12 * T21 - R12 * R21) * fw_chunk_begin_scaled +
                                             R21  * bw_chunk_begin_scaled;

    bw_chunk_end_scaled = (-R12) * fw_chunk_begin_scaled +
                           (1.0) * bw_chunk_begin_scaled;

    if (abs(T21) < 1e-10)
    {
      cout << "Warning: small T21: " << T21 << endl;
      cout << "Possible loss of precision." << endl;
    }
    
    fw_chunk_end_scaled /= T21;
    bw_chunk_end_scaled /= T21;

    fw->push_back(fw_chunk_end_scaled * scaling_factor);
    bw->push_back(bw_chunk_end_scaled * scaling_factor);

    // Propagate in medium and scale.
    
    Complex I_kz_d = I * s->get_ext()->get_mode(1)->get_kz() * chunks[k].d;
    
    if (real(I_kz_d) > 0)
    {
      fw_chunk_end_scaled *= exp(-2.0*I_kz_d);
      scaling_factor      *= exp(     I_kz_d);
    }
    else
    {
      bw_chunk_end_scaled *= exp(+2.0*I_kz_d);
      scaling_factor      *= exp(    -I_kz_d);
    }

    fw->push_back(fw_chunk_end_scaled * scaling_factor);
    bw->push_back(bw_chunk_end_scaled * scaling_factor);
    
    // Update values for next iteration.
    
    fw_chunk_begin_scaled = fw_chunk_end_scaled;
    bw_chunk_begin_scaled = bw_chunk_end_scaled;
  }
}
