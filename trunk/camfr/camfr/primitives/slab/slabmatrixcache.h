
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
                  const SlabCache* cache, const std::vector<Complex>* disc,
                  bool calc_all = true);

  int n; // Auxiliary variable to simplify comstructor.  
    
  cHyperM TE_TE,       TM_TM;
  cHyperM Ex_Hz_cross, Ez_Hx_cross;
  cHyperM Ex_Hz_self,  Ez_Hx_self;
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
       const SlabCache* cache, const std::vector<Complex>* disc);

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



