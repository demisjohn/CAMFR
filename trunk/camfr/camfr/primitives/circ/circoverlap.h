
/////////////////////////////////////////////////////////////////////////////
//
// File:     circoverlap.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990519
// Version:  1.0
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CIRCOVERLAP_H
#define CIRCOVERLAP_H

#include "../../math/linalg/linalg.h"
#include <vector>

/////////////////////////////////////////////////////////////////////////////
//
// CircCache
//  
/////////////////////////////////////////////////////////////////////////////

struct CircCache
{
    CircCache(int N, int k)
      :  E_l(2,N,k,fortranArray),  H_l(2,N,k,fortranArray),
         E_u(2,N,k,fortranArray),  H_u(2,N,k,fortranArray),
         
        dE_l(2,N,k,fortranArray), dH_l(2,N,k,fortranArray),
        dE_u(2,N,k,fortranArray), dH_u(2,N,k,fortranArray)  {}
    
    cHyperM  E_l, H_l, E_u, H_u, dE_l, dH_l, dE_u, dH_u;
};



/////////////////////////////////////////////////////////////////////////////
//
// Calculates overlapintegral Int(E_I x H_II) between mode_I and mode_II.
//
// Optionally, the field profiles can be retrieved from a cache, indexed
// by i,j, I_index, II_index.
// i (j) is the mode number of mode_I (mode_II) in medium_I (medium_II).
// I_index is either 1 or 2, indicating which part of the cache to use.
// II_index is complementary to I_index.
// The discontinuities are stored in the vector disc.
//  
/////////////////////////////////////////////////////////////////////////////

class Circ_M_Mode; // forward declation -- see circmode.h

Complex overlap(const Circ_M_Mode* mode_I,
                const Circ_M_Mode* mode_II,
                const CircCache* cache=NULL,
                const std::vector<Complex>* disc=NULL,
                int i=0, int j=0, int I_index=0, int II_index=0);



#endif



