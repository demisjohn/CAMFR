
/////////////////////////////////////////////////////////////////////////////
//
// File:     sectionoverlap.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020225
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SECTIONOVERLAP_H
#define SECTIONOVERLAP_H

#include "sectionmode.h"

/////////////////////////////////////////////////////////////////////////////
//
// SectionCache
//  
/////////////////////////////////////////////////////////////////////////////

struct SectionCache
{
    SectionCache(int N, int k)
      : fw_l(2,N,k,fortranArray),  bw_l(2,N,k,fortranArray),
        fw_u(2,N,k,fortranArray),  bw_u(2,N,k,fortranArray) {}
    
    cHyperM fw_l, bw_l, fw_u, bw_u;
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
// The list of discontinuities is in the vector disc.
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap(const SectionMode* mode_I,
                const SectionMode* mode_II,
                const SectionCache* cache=NULL,
                const vector<Complex>* disc=NULL,
                int i=0, int j=0,  int I_index=0, int II_index=0) {};



#endif



