
/////////////////////////////////////////////////////////////////////////////
//
// File:     infstack.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020206
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef INFSTACK_H
#define INFSTACK_H

#include "scatterer.h"
#include "bloch.h"

/////////////////////////////////////////////////////////////////////////////
//
// InfStack
//
//   Semi-infinite repetition of the same period.
//
/////////////////////////////////////////////////////////////////////////////

class InfStack : public DenseScatterer
{
  public:

    InfStack(const Expression& e, const Complex& W_=0);
  
    void calcRT();
    
  protected:

    BlochStack s;
    Complex W; // TMP
};



#endif
