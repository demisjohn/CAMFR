
/////////////////////////////////////////////////////////////////////////////
//
// File:     icache.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000608
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "icache.h"
#include "waveguide.h"
#include "scatterer.h"
#include "interface.h"
#include "bloch.h"

/////////////////////////////////////////////////////////////////////////////
//
// InterfaceCache::get_interface
//
/////////////////////////////////////////////////////////////////////////////

Scatterer* InterfaceCache::get_interface(Waveguide* wg1, Waveguide* wg2)
{
  // Interface already in cache?

  Scatterer* sc;
  
  bool found = cache.lookup(std::pair<Waveguide*, Waveguide*>(wg1, wg2), &sc);

  if (found)
  {
    
    if (    (global.always_dense == true)
         && (wg1->is_uniform() && wg2->is_uniform())
         && (wg1 != wg2) )
    {
      deregister(wg1, wg2);
      deregister(wg2, wg1);
    }
    else
      return sc;
  }
  
  // Transparent DiagScatterer?

  if ( (wg1 == wg2) && (!dynamic_cast<MonoWaveguide*>(wg1)) )
  {
    sc = new TransparentScatterer(*wg1);
    cache.store(std::pair<Waveguide*, Waveguide*>(wg1, wg2), sc);
    return sc;
  }
 
  // General case.

  if (wg1->is_uniform() && wg2->is_uniform() && (global.always_dense == false))
  {
    if (    dynamic_cast<MonoWaveguide*>(wg1)
         && dynamic_cast<MonoWaveguide*>(wg2) )
      sc = new MonoInterface(*wg1, *wg2);
    else 
      sc = new DiagInterface(*wg1, *wg2);
  }
  else
    sc = new DenseInterface(*wg1, *wg2);
  
  cache.store(std::pair<Waveguide*, Waveguide*>(wg1, wg2), sc);

  if (wg1 != wg2)
  {
    Scatterer* sc_flip;
    if (    wg1->is_uniform() && wg2->is_uniform() 
        && (global.always_dense == false))
    {
      if (    dynamic_cast<MonoWaveguide*>(wg1)
           && dynamic_cast<MonoWaveguide*>(wg2) )
        sc_flip = new MonoInterface(*wg2, *wg1);
      else
        sc_flip = new DiagInterface(*wg2, *wg1);
    }
    else
    {
      // Don't make a flipped scatterer for interfaces with a BlochStack
      
      if ( dynamic_cast<BlochStack*>(wg1) || dynamic_cast<BlochStack*>(wg2) )
	sc_flip = new DenseInterface(*wg2,*wg1);
      else
        sc_flip = new FlippedScatterer(*dynamic_cast<MultiScatterer*>(sc));
    }

    cache.store(std::pair<Waveguide*, Waveguide*>(wg2, wg1), sc_flip);
  }

  return sc;
}



/////////////////////////////////////////////////////////////////////////////
//
// InterfaceCache::~InterfaceCache
//
/////////////////////////////////////////////////////////////////////////////

InterfaceCache::~InterfaceCache()
{  
  for (Cache<std::pair<Waveguide*, Waveguide*>, Scatterer*>::iter
         i=cache.begin(); i!=cache.end(); ++i)
    delete i->second;
}



/////////////////////////////////////////////////////////////////////////////
//
// InterfaceCache::deregister
//
/////////////////////////////////////////////////////////////////////////////

void InterfaceCache::deregister(Waveguide* wg)
{
  std::vector<std::pair<Waveguide*,Waveguide*> > to_wipe;

  for (Cache<std::pair<Waveguide*, Waveguide*>, Scatterer*>::iter
         i=cache.begin(); i!=cache.end(); ++i)
    if ( (i->first.first == wg) || (i->first.second == wg) )
    {
      delete i->second;
      to_wipe.push_back(i->first);
    }

  for (int i=0; i<to_wipe.size(); i++)
    cache.erase(to_wipe[i]);
}



/////////////////////////////////////////////////////////////////////////////
//
// InterfaceCache::deregister
//
/////////////////////////////////////////////////////////////////////////////

void InterfaceCache::deregister(Waveguide* wg1, Waveguide* wg2)
{  
  std::vector<std::pair<Waveguide*,Waveguide*> > to_wipe;

  for (Cache<std::pair<Waveguide*, Waveguide*>, Scatterer*>::iter
         i=cache.begin(); i!=cache.end(); ++i)
    if ( (i->first.first == wg1) && (i->first.second == wg2) )
    {
      delete i->second;      
      to_wipe.push_back(i->first);
    }  

  for (int i=0; i<to_wipe.size(); i++)
    cache.erase(to_wipe[i]);
}



/////////////////////////////////////////////////////////////////////////////
//
// InterfaceCache::clear
//
/////////////////////////////////////////////////////////////////////////////

void InterfaceCache::clear()
{ 
  for (Cache<std::pair<Waveguide*, Waveguide*>, Scatterer*>::iter
         i=cache.begin(); i!=cache.end(); ++i)
    delete i->second;

  cache.clear();
}



/////////////////////////////////////////////////////////////////////////////
//
// Global interface cache
//
/////////////////////////////////////////////////////////////////////////////

InterfaceCache interface_cache;




