
/////////////////////////////////////////////////////////////////////////////
//
// File:     storage.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010223
// Version:  1.1
//
// Copyright (C) 1999,2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef STORAGE_H
#define STORAGE_H

#include <map>
#include <vector>

/////////////////////////////////////////////////////////////////////////////
//
// A general lookup cache that can grow arbitrarily large.
//
/////////////////////////////////////////////////////////////////////////////

template <class Key, class Value>
class Cache
{
  public:
  
    Cache() {}
    
    ~Cache() {}
    
    typedef typename std::map<Key, Value>::const_iterator const_iter;
    typedef typename std::map<Key, Value>::      iterator       iter;
    
    bool lookup(const Key& key, Value* value) const
      {
        const_iter i = storage.find(key);

        if (i == storage.end())
          return false;
        else
        {
          *value = i->second;
          return true;
        }
      }
    
    void store(const Key& key, const Value& value)
      {
        std::pair<Key, Value> new_pair(key, value);

        std::pair<const_iter, bool> result = storage.insert(new_pair);

        if (result.second == false) // Entry was already present.
        {
          storage.erase(key);
          storage.insert(new_pair);
        }
      }

    void erase(const Key& key) {storage.erase(key);}

    void clear() {storage.clear();}

    const_iter begin() const {return storage.begin();}
    const_iter   end() const {return storage.end();}

          iter begin()       {return storage.begin();}
          iter   end()       {return storage.end();}    
    
  protected:
    
    std::map<Key, Value> storage;
};



/////////////////////////////////////////////////////////////////////////////
//
// Temporary storage for pointers to objects. The objects itself are
// destroyed when 'clear' is called or when TmpStorage goes out of scope.
//
/////////////////////////////////////////////////////////////////////////////

template <class T>
class TmpStorage
{
  public:

    TmpStorage() {}
    
    ~TmpStorage() {clear();}
    
    void store(T* t) {tmps.push_back(t);}

    void clear()
      {
        for (unsigned int i=0; i<tmps.size(); i++)
          delete tmps[i];

        tmps.clear();
      }

  protected:

    std::vector<T*> tmps;
};



#endif
