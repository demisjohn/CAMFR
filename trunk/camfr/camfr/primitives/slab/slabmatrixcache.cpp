
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabmatrixcache.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20011211
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "isoslab/slaboverlap.h"
#include "isoslab/slabmode.h"
#include "slabmatrixcache.h"
#include "generalslab.h"

using std::vector;
using std::pair;

/////////////////////////////////////////////////////////////////////////////
//
// OverlapMatrices::OverlapMatrices
//
/////////////////////////////////////////////////////////////////////////////

OverlapMatrices::OverlapMatrices
  (const SlabImpl* medium_I, const SlabImpl* medium_II,
   const SlabCache* cache, const vector<Complex>* disc, bool calc_all)
    : n(medium_I->N()/2), 
      TE_TE      (2,n,n,fortranArray), TM_TM      (2,n,n,fortranArray),
      Ex_Hz_cross(2,n,n,fortranArray), Ez_Hx_cross(2,n,n,fortranArray),
      Ex_Hz_self (2,n,n,fortranArray), Ez_Hx_self (2,n,n,fortranArray)
{
  Complex Ex_Hz_I_II_ij, Ez_Hx_I_II_ij, Ex_Hz_II_I_ij,  Ez_Hx_II_I_ij;
  Complex Ex_Hz_I_I_ij,  Ez_Hx_I_I_ij,  Ex_Hz_II_II_ij, Ez_Hx_II_II_ij;

  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
    {
      TE_TE(1,i,j) = overlap
        (dynamic_cast<const SlabMode*>(medium_I ->get_mode(i)),
         dynamic_cast<const SlabMode*>(medium_II->get_mode(j)),
         cache, disc, i, j, 1, 2);

      TE_TE(2,i,j) = overlap
        (dynamic_cast<const SlabMode*>(medium_II->get_mode(i)),
         dynamic_cast<const SlabMode*>(medium_I ->get_mode(j)),
         cache, disc, i, j, 2, 1);

      TM_TM(1,i,j) = overlap
        (dynamic_cast<const SlabMode*>(medium_I ->get_mode(n+i)),
         dynamic_cast<const SlabMode*>(medium_II->get_mode(n+j)),
         cache, disc, n+i, n+j, 1, 2);

      TM_TM(2,i,j) = overlap
        (dynamic_cast<const SlabMode*>(medium_II->get_mode(n+i)),
         dynamic_cast<const SlabMode*>(medium_I ->get_mode(n+j)),
         cache, disc, n+i, n+j, 2, 1);

      overlap_TM_TE(dynamic_cast<const SlabMode*>(medium_I ->get_mode(n+i)),
                    dynamic_cast<const SlabMode*>(medium_II->get_mode(  j)),
                    &Ex_Hz_I_II_ij, &Ez_Hx_I_II_ij,
                    cache, disc, n+i, j, 1, 2);

      overlap_TM_TE(dynamic_cast<const SlabMode*>(medium_II->get_mode(n+i)),
                    dynamic_cast<const SlabMode*>(medium_I ->get_mode(  j)),
                    &Ex_Hz_II_I_ij, &Ez_Hx_II_I_ij,
                    cache, disc, n+i, j, 2, 1);      

      Ex_Hz_cross(1,i,j) = Ex_Hz_I_II_ij;  Ez_Hx_cross(1,i,j) = Ez_Hx_I_II_ij;
      Ex_Hz_cross(2,i,j) = Ex_Hz_II_I_ij;  Ez_Hx_cross(2,i,j) = Ez_Hx_II_I_ij;

      if (calc_all)
      {
        overlap_TM_TE(dynamic_cast<const SlabMode*>(medium_I ->get_mode(n+i)),
                      dynamic_cast<const SlabMode*>(medium_I ->get_mode(  j)),
                      &Ex_Hz_I_I_ij, &Ez_Hx_I_I_ij,
                      cache, disc, n+i, j, 1, 1);

        overlap_TM_TE(dynamic_cast<const SlabMode*>(medium_II->get_mode(n+i)),
                      dynamic_cast<const SlabMode*>(medium_II->get_mode(  j)),
                      &Ex_Hz_II_II_ij, &Ez_Hx_II_II_ij,
                      cache, disc, n+i, j, 2, 2);

        Ex_Hz_self(1,i,j)=Ex_Hz_I_I_ij;   Ez_Hx_self(1,i,j)=Ez_Hx_I_I_ij;
        Ex_Hz_self(2,i,j)=Ex_Hz_II_II_ij; Ez_Hx_self(2,i,j)=Ez_Hx_II_II_ij;
      }
    }

  if (global.davy == true)
  {
    std::cout << "Medium I core n "  << medium_I ->get_core()->n() 
              << std::endl;
    std::cout << "Medium II core n " << medium_II->get_core()->n() 
              << std::endl;

    blitz::Range r(1,n);

    write_to_textfile("TE_TE_1.txt", TE_TE(1,r,r));
    write_to_textfile("TE_TE_2.txt", TE_TE(2,r,r));
    write_to_textfile("TM_TM_1.txt", TM_TM(1,r,r));
    write_to_textfile("TM_TM_2.txt", TM_TM(2,r,r));
    write_to_textfile("Ex_Hz_cross_1.txt", Ex_Hz_cross(1,r,r));
    write_to_textfile("Ex_Hz_cross_2.txt", Ex_Hz_cross(2,r,r));
    write_to_textfile("Ez_Hx_cross_1.txt", Ez_Hx_cross(1,r,r));
    write_to_textfile("Ez_Hx_cross_2.txt", Ez_Hx_cross(2,r,r));
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabMatrixCache::get_matrices
//
/////////////////////////////////////////////////////////////////////////////

OverlapMatrices* SlabMatrixCache::get_matrices(SlabImpl* wg1, SlabImpl* wg2,
  const SlabCache* slabcache, const vector<Complex>* disc)
{
  // Already in cache?

  OverlapMatrices* m;
  
  bool found = cache.lookup(pair<SlabImpl*, SlabImpl*>(wg1, wg2), &m);

  if (found)
    return m;

  if (!slabcache || !disc)
    py_error("SlabMatrixCache: no field cache provided.");
  
  // Calculate matrices and cache them.

  m = new OverlapMatrices(wg1, wg2, slabcache, disc);
  
  cache.store(pair<SlabImpl*, SlabImpl*>(wg1, wg2), m);

  return m;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabMatrixCache::~SlabMatrixCache
//
/////////////////////////////////////////////////////////////////////////////

SlabMatrixCache::~SlabMatrixCache()
{  
  for (Cache<pair<SlabImpl*, SlabImpl*>, OverlapMatrices*>::iter
         i=cache.begin(); i!=cache.end(); ++i)
    delete i->second;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabMatrixCache::deregister
//
/////////////////////////////////////////////////////////////////////////////

void SlabMatrixCache::deregister(SlabImpl* wg)
{
  for (Cache<pair<SlabImpl*, SlabImpl*>, OverlapMatrices*>::iter
         i=cache.begin(); i!=cache.end(); ++i)
    if ( (i->first.first == wg) || (i->first.second == wg) )
      cache.erase(i->first);
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabMatrixCache::clear
//
/////////////////////////////////////////////////////////////////////////////

void SlabMatrixCache::clear()
{ 
  for (Cache<pair<SlabImpl*, SlabImpl*>, OverlapMatrices*>::iter
         i=cache.begin(); i!=cache.end(); ++i)
    delete i->second;

  cache.clear();
}



/////////////////////////////////////////////////////////////////////////////
//
// Global overlap matrix cache
//
/////////////////////////////////////////////////////////////////////////////

SlabMatrixCache slabmatrix_cache;

