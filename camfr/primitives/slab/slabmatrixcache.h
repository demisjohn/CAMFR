
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabmatrixcache.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20011211
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SLABMATRIXCACHE_H
#define SLABMATRIXCACHE_H

#include "../../math/linalg/linalg.h"
#include "../../util/storage.h"

/////////////////////////////////////////////////////////////////////////////
//
// STRUCT: OverlapMatrices
//
//   Off-axis overlap matrices for slabs.
//
/////////////////////////////////////////////////////////////////////////////

class SlabImpl;  // forward declation - see generalslab.h
class SlabCache; // forward declation - see generalslab.h

struct OverlapMatrices
{
  OverlapMatrices(const SlabImpl* medium_I, const SlabImpl* medium_II, 
                  const SlabCache* cache, const vector<Complex>* disc);

  int n; // Auxiliary variable to simplify comstructor.  
    
  cMatrix TE_TE_I_II, TM_TM_I_II, TE_TE_II_I,  TM_TM_II_I;
  cMatrix Ex_Hz_I_II, Ez_Hx_I_II, Ex_Hz_II_I,  Ez_Hx_II_I;
  cMatrix Ex_Hz_I_I,  Ez_Hx_I_I,  Ex_Hz_II_II, Ez_Hx_II_II;
};


/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabMatrixCache
//
//   Cache for off-axis overlap matrices for slabs.
//
/////////////////////////////////////////////////////////////////////////////

class SlabMatrixCache
{
  public:

    SlabMatrixCache() {}
    ~SlabMatrixCache();
      
    OverlapMatrices* get_matrices(SlabImpl* wg1, SlabImpl* wg2,
       const SlabCache* cache, const vector<Complex>* disc);

    void deregister(SlabImpl* wg);
    
    void clear();

  protected:

    Cache<std::pair<SlabImpl*, SlabImpl*>, OverlapMatrices*> cache;
};



/////////////////////////////////////////////////////////////////////////////
//
// Global cache.
//
/////////////////////////////////////////////////////////////////////////////

extern SlabMatrixCache slabmatrix_cache;


#endif



