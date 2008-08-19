
/////////////////////////////////////////////////////////////////////////////
//
// File:     vectorutil.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990625
// Version:  1.0
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef VECTORUTIL_H
#define VECTORUTIL_H

#include "../defs.h"
  
/////////////////////////////////////////////////////////////////////////////
//
// remove_elems
//
//   Remove all elements from vector where abs(v[i]-t) < eps.
//
/////////////////////////////////////////////////////////////////////////////

template <class T>
void remove_elems(vector<T>* v, const T& t, Real eps=0)
{
  vector<T>    v_bis(*v);
  vector<bool> keep(v_bis.size(), true);
  
  for (unsigned int i=0; i<v_bis.size(); i++)
  {
    if (abs(v_bis[i] - t) <= eps)
      keep[i] = false;   
  }

  v->clear();
  for (unsigned int i=0; i<v_bis.size(); i++)
    if (keep[i])
      v->push_back(v_bis[i]);
}



/////////////////////////////////////////////////////////////////////////////
//
// remove_copies
//
//   Remove double elements from vector if abs(v[i]-v[j]) < eps.
//
//   Note: we do the funky cast to Complex because we also use this
//   routine to compare pointers. Oh my...
//
/////////////////////////////////////////////////////////////////////////////

template <class T>
void remove_copies(vector<T>* v, Real eps=0)
{
  vector<T>    v_bis(*v);
  vector<bool> unique(v_bis.size(), true);
  
  for (unsigned int i=0; i<v_bis.size(); i++)
  {
    if (unique[i] == false)
      continue;
    
    for (unsigned int j=i+1; j<v_bis.size(); j++)
      if (abs(Complex(v_bis[i] - v_bis[j])) <= eps)
        unique[j] = false;   
  }

  v->clear();
  for (unsigned int i=0; i<v_bis.size(); i++)
    if (unique[i])
      v->push_back(v_bis[i]);  
}



/////////////////////////////////////////////////////////////////////////////
//
// remove_opposites
//
//   Remove opposites elements from vector if abs(v[i]+v[j]) < eps.
//
/////////////////////////////////////////////////////////////////////////////

template <class T>
void remove_opposites(vector<T>* v, Real eps=0)
{
  vector<T>    v_bis(*v);
  vector<bool> unique(v_bis.size(), true);
  
  for (unsigned int i=0; i<v_bis.size(); i++)
  {
    if (unique[i] == false)
      continue;
    
    for (unsigned int j=i+1; j<v_bis.size(); j++)
      if (abs(v_bis[i] + v_bis[j]) <= eps)
        unique[j] = false;   
  }

  v->clear();
  for (unsigned int i=0; i<v_bis.size(); i++)
    if (unique[i])
      v->push_back(v_bis[i]);  
}



/////////////////////////////////////////////////////////////////////////////
//
// is_present
//
//   Check if element is present with tolarance eps.
//
/////////////////////////////////////////////////////////////////////////////

template <class T>
bool is_present(const vector<T>& v, const T& z, Real eps=0)
{
  for (unsigned int i=0; i<v.size(); i++)
    if (abs(v[i]-z) <= eps)
      return true;
  
  return false;
}


#endif
