
/////////////////////////////////////////////////////////////////////////////
//
// File:     icache.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000608
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ICACHE_H
#define ICACHE_H

#include "util/storage.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: InterfaceCache
//
//   Cache for previously created Interfaces.
//
/////////////////////////////////////////////////////////////////////////////

class Waveguide; // forward declaration - see waveguide.h
class Scatterer; // forward declaration - see scatterer.h

class InterfaceCache
{
  public:

    InterfaceCache() {}
    ~InterfaceCache();
      
    Scatterer* get_interface(Waveguide* wg1, Waveguide* wg2);

    void deregister(Waveguide* wg);
    
    void clear();

  protected:

    Cache<std::pair<Waveguide*, Waveguide*>, Scatterer*> cache;
};



/////////////////////////////////////////////////////////////////////////////
//
// Global interface cache
//
/////////////////////////////////////////////////////////////////////////////

extern InterfaceCache interface_cache;


#endif



